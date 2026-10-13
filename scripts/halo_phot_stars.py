"""Calibrate ``halo_phot``'s axis-ratio gate on Gaia DR3 stars, and list every catalogued image's
axis ratio and concentration on the same mosaics.

    python scripts/halo_phot_stars.py --from results/halo_phot/elgordo/halo_phot.json \\
        --out results/halo_phot/stars/elgordo.json

Reads the model and ``_i2d`` URIs of an earlier ``halo_phot`` run. One VizieR cone (I/355/gaiadr3,
4' around the mean image position); stars are Gaia sources with proper motion or parallax > 3 sigma
(ASSUMPTION: that removes compact galaxies). A star is kept in a band when its stamp core (0.3")
has no masked pixel (saturation), and it is counted towards the gate when S/N >= ``MIN_SNR`` and
concentration >= ``MIN_CONCENTRATION``, as for images. ``MIN_AXIS_RATIO`` is the 5th percentile
of the pooled counted stars per band (printed by ``--pool``). Provenance: Gaia rows observed, the
rest derived.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import halo_phot as hp  # noqa: E402
import lens_consistency as lc  # noqa: E402

GAIA_RADIUS_ARCMIN = 4.0
PM_SIGMA = 3.0


def gaia_stars(ra: float, dec: float):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    cols = ["RA_ICRS", "DE_ICRS", "Gmag", "pmRA", "pmDE", "e_pmRA", "e_pmDE", "Plx", "e_Plx"]
    g = Vizier(columns=cols, row_limit=-1).query_region(
        SkyCoord(ra, dec, unit="deg"), radius=GAIA_RADIUS_ARCMIN * u.arcmin, catalog="I/355/gaiadr3"
    )[0]

    def f(name, fill):
        return np.asarray(g[name].filled(fill), float)

    pm = np.hypot(f("pmRA", 0), f("pmDE", 0)) / np.hypot(f("e_pmRA", 99), f("e_pmDE", 99))
    return g[(pm > PM_SIGMA) | (f("Plx", 0) / f("e_Plx", 99) > PM_SIGMA)]


def shape(stamp_tuple, recentre=hp.RECENTRE_ARCSEC) -> dict:
    sci, err, xx, yy = stamp_tuple
    filt = hp.halo_filtered(sci, hp._pix(xx, yy))
    flux, e, conc = hp.measure(sci, err, xx, yy, filtered=filt, recentre=recentre)
    q = hp.axis_ratio(filt, xx, yy)
    return {"q": q, "conc": conc, "snr": flux / e if e > 0 else float("nan")}


def calibrate(run_json: Path) -> dict:
    d = json.loads(run_json.read_text())
    model, files, par = lc.load_model(d["model"])
    images, _ = lc.image_list(d["model"], files, par)
    lc.apply_frame_offset(d["model"], model, images)
    stars = gaia_stars(float(np.mean(images["ra"])), float(np.mean(images["dec"])))
    out = {"model": d["model"], "bands": d["bands"], "n_gaia": len(stars)}
    out |= {"stars": {}, "images": {}}
    for band, uri in d["bands"].items():
        stamp, close = lc.image_stamper(uri, 3.0)
        try:
            rows = []
            for r in stars:
                st = stamp(float(r["RA_ICRS"]), float(r["DE_ICRS"]))
                if not np.isfinite(st[0][np.hypot(st[2], st[3]) <= 0.3]).all():
                    continue  # off the mosaic, or a saturated / masked core
                rows.append({"G": float(r["Gmag"]), **shape(st)})
            out["stars"][band] = rows
            out["images"][band] = {
                str(r["image_id"]): shape(stamp(float(r["ra"]), float(r["dec"]))) for r in images
            }
        finally:
            close()
    return out


def pool(paths: list[Path]) -> dict:
    """5th / 16th / 50th percentiles of the counted stars' axis ratio per band, pooled."""
    q: dict[str, list[float]] = {}
    for p in paths:
        for band, rows in json.loads(p.read_text())["stars"].items():
            q.setdefault(band, []).extend(
                r["q"]
                for r in rows
                if (r["snr"] or 0) >= hp.MIN_SNR and (r["conc"] or 0) >= hp.MIN_CONCENTRATION
            )
    return {
        b: {"n": len(v), "pct_5_16_50": [round(float(x), 3) for x in np.percentile(v, [5, 16, 50])]}
        for b, v in q.items()
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--from", dest="src", type=Path, help="an earlier halo_phot.json")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--pool", type=Path, nargs="+", help="pool these outputs and print the gate")
    args = ap.parse_args(argv)
    if args.pool:
        print(json.dumps(pool(args.pool)))
        return 0
    res = calibrate(args.src)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    clean = json.loads(json.dumps(res, default=float).replace("NaN", "null"))
    args.out.write_text(json.dumps(clean, indent=0) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
