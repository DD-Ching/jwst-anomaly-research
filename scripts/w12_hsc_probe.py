"""W2 deflector test at HST resolution from the Hubble Source Catalog (HSC v3) (D-060).

For every galaxy-scale lensed-quasar or radio-selected system of the merged lens catalogues
(``w12_lenscats.load_systems``), fetch the HSC v3 summary sources within ``radius`` (MAST
catalogs API, one request per system, thread pool) and classify: >= 2 point-like sources
(concentration index CI < ``ci_point``; sources in fewer than ``min_images``
HSC images are dropped as likely artifacts) inside ``image_radius`` are the quasar images; an
extended source (CI >= ``ci_ext``) near their centroid and not on an image is the deflector.
The test is validated on systems whose lens galaxy has a published redshift (lenscat
``lens_z_known``); its efficiency there bounds what a "none" means. Thresholds are ASSUMPTIONs
(``Params``); outputs are ``derived``. Nothing here is evidence of exotic physics.

    python scripts/w12_hsc_probe.py --out outputs/w12_hsc_probe
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table, vstack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_lenscats  # noqa: E402

from jwst_anomaly import schema  # noqa: E402

HSC = "https://catalogs.mast.stsci.edu/api/v0.1/hsc/v3/summary/magaper2.csv"


@dataclass(frozen=True)
class Params:
    radius: float = 4.0  # arcsec, HSC cone around the catalogue position
    image_radius: float = 3.0  # arcsec, quasar images within this of the position
    ci_point: float = 1.3  # HSC concentration index below this: point-like (ASSUMPTION)
    ci_ext: float = 1.5  # at or above this: extended (ASSUMPTION)
    lens_frac: float = 0.6  # deflector within lens_frac * max image separation of the centroid
    lens_min: float = 0.5  # arcsec, minimum deflector search radius
    image_exclusion: float = 0.2  # arcsec, an extended source this close to an image is the image
    min_images: int = (
        2  # HSC NumImages >= this; single-image detections are mostly artifacts (MAST)
    )
    workers: int = 8


def classify(ra: float, dec: float, src: Table, p: Params) -> tuple[str, int, float]:
    """(status, number of point images, max image separation in arcsec) for one system.

    status: "deflector", "none" (>= 2 point images, no extended source near their centroid)
    or "undecided" (fewer than two point images)."""
    if len(src) and "NumImages" in src.colnames:
        src = src[np.asarray(src["NumImages"]) >= p.min_images]
    if len(src) == 0:
        return "undecided", 0, np.nan
    x = (np.asarray(src["MatchRA"], float) - ra) * np.cos(np.radians(dec)) * 3600
    y = (np.asarray(src["MatchDec"], float) - dec) * 3600
    ci = np.asarray(src["CI"], float)
    pt = (ci < p.ci_point) & (np.hypot(x, y) < p.image_radius)
    if pt.sum() < 2:
        return "undecided", int(pt.sum()), np.nan
    px, py = x[pt], y[pt]
    sep = float(np.max(np.hypot(px[:, None] - px, py[:, None] - py)))
    d_cen = np.hypot(x - px.mean(), y - py.mean())
    d_img = np.min(np.hypot(x[:, None] - px, y[:, None] - py), axis=1)
    near = d_cen < max(p.lens_frac * sep, p.lens_min)
    ext = (ci >= p.ci_ext) & near & (d_img > p.image_exclusion)
    return ("deflector" if ext.any() else "none"), int(pt.sum()), sep


def fetch(ra: float, dec: float, p: Params) -> Table:
    err = None
    for _ in range(3):
        try:
            r = requests.get(
                HSC, params=dict(ra=ra, dec=dec, radius=p.radius / 3600, pagesize=500), timeout=60
            )
            r.raise_for_status()
            text = r.text.strip()
            return Table.read(text, format="ascii.csv") if "\n" in text else Table()
        except (requests.RequestException, ValueError) as e:
            err = e
            time.sleep(2)
    raise RuntimeError(f"HSC query failed at {ra:.5f} {dec:+.5f}: {err}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    p, lp = Params(), w12_lenscats.Params()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    s, _ = w12_lenscats.load_systems(lp)
    gal = ~np.asarray(s["any_group_scale"], bool)
    gal &= ~(np.asarray(s["theta_e"], float) > lp.theta_e_max)
    g = s[gal]
    g["selection"] = [w12_lenscats.selection_class(r) for r in g]
    q = g[np.isin(g["selection"], ("quasar", "radio"))]
    with ThreadPoolExecutor(p.workers) as ex:
        srcs = list(ex.map(lambda r: fetch(float(r["ra"]), float(r["dec"]), p), q))
    res = Table(
        {c: q[c] for c in ("system_id", "name", "ra", "dec", "selection", "lens_z_known", "z_lens")}
    )
    cls = [classify(float(r["ra"]), float(r["dec"]), t, p) for r, t in zip(q, srcs, strict=True)]
    res["n_hsc"] = [len(t) for t in srcs]
    res["status"] = [c[0] for c in cls]
    res["n_point"] = [c[1] for c in cls]
    res["max_sep"] = [c[2] for c in cls]
    res.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"merged lens catalogues (D-056); HSC v3 summary magaper2 ({HSC})",
        params=asdict(p),
    )
    res.write(out / "systems.ecsv", overwrite=True)
    hits = [t for t in srcs if len(t)]
    if hits:
        vstack(hits, metadata_conflicts="silent").write(out / "hsc_sources.ecsv", overwrite=True)

    def counts(m):
        vals, n = np.unique(np.asarray(res["status"])[m], return_counts=True)
        return {str(k): int(v) for k, v in zip(vals, n, strict=True)}

    known = np.asarray(res["lens_z_known"], bool)
    cov = np.asarray(res["n_hsc"]) > 0
    k = counts(known & cov)
    dec_known = k.get("deflector", 0) + k.get("none", 0)
    summary = {
        "systems": len(res),
        "with_hsc_sources": int(cov.sum()),
        "lens_z_known": k,
        "no_lens_z": counts(~known & cov),
        "efficiency_on_known_lenses": k.get("deflector", 0) / dec_known if dec_known else None,
        "none_names": [str(n) for n in res["name"][np.asarray(res["status"]) == "none"]],
        "params": asdict(p),
        "wall_s": round(time.time() - t0, 1),
        "provenance": schema.Provenance.DERIVED.value,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "params"}, indent=1))


if __name__ == "__main__":
    main()
