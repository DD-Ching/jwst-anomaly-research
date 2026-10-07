"""Combine per-band transient candidates and exclude positions near bright Gaia sources (D-027).

Takes the ``transients.ecsv`` of ``scripts/transient_search.py`` for several bands of one epoch
pair and keeps positions where the same ``kind`` is found in at least ``--min-bands`` bands within
``--radius-arcsec``. A real transient changes in several bands; a deblending split, a cosmic ray
or an asteroid (the SW filters of one epoch are taken minutes apart) usually shows up in one.

Positions near bright Gaia DR3 sources are then excluded: there, diffraction spikes that rotate
with the position angle and PSF-wing differences dominate a fixed aperture (SMACS/VENUS, D-027).
The exclusion radius grows with brightness, ``r = r0 * 10**(0.2 * (g0 - G))`` clipped to
``[r_min, r_max]`` (ASSUMPTION; defaults 1.5" at G = 20, 1.5–12"). Gaia sources include compact
galaxies, so the mask is conservative. Every output is ``derived``; the Gaia table is ``observed``.

    python scripts/transient_combine.py --out outputs/transients_sunrise \\
        --transients F090W outputs/transients_sunrise/f090w/transients.ecsv \\
        --transients F277W outputs/transients_sunrise/f277w/transients.ecsv \\
        --gaia-centre 24.355 -8.457 --reference 24.346822 -8.464511 earendel
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack


def coincident(tables: dict[str, Table], radius_arcsec: float = 0.3, min_bands: int = 2) -> Table:
    """One row per position whose ``kind`` recurs in >= ``min_bands`` bands within the radius.

    The row is the first band's (in ``tables`` order); ``bands`` lists every band with that kind
    there and ``n_bands`` counts them."""
    parts = []
    for band, t in tables.items():
        t = t.copy()
        t["band"] = band
        parts.append(t)
    if not parts or sum(len(t) for t in parts) == 0:
        return Table(
            names=("ra", "dec", "kind", "band", "bands", "n_bands"),
            dtype=(float, float, str, str, str, int),
        )
    allc = vstack(parts, metadata_conflicts="silent")
    c = SkyCoord(allc["ra"], allc["dec"], unit="deg")
    i, j, _, _ = c.search_around_sky(c, radius_arcsec * u.arcsec)
    kinds = np.asarray(allc["kind"]).astype(str)
    bands_at: dict[int, set[str]] = {}
    for a, b in zip(i, j, strict=True):
        if kinds[a] == kinds[b]:
            bands_at.setdefault(int(a), set()).add(str(allc["band"][b]))
    order = list(tables)
    keep, used = [], np.zeros(len(allc), bool)
    for k in range(len(allc)):
        bs = bands_at.get(k, set())
        if used[k] or len(bs) < min_bands:
            continue
        same = np.flatnonzero((c.separation(c[k]).arcsec <= radius_arcsec) & (kinds == kinds[k]))
        used[same] = True
        keep.append((k, ",".join(sorted(bs, key=order.index)), len(bs)))
    out = allc[[k for k, _, _ in keep]]
    out["bands"] = [b for _, b, _ in keep]
    out["n_bands"] = np.array([n for _, _, n in keep], int)
    out.meta = {"provenance": "derived", "radius_arcsec": radius_arcsec, "min_bands": min_bands}
    return out


def exclusion_radius(
    gmag, r0: float = 1.5, g0: float = 20.0, r_min: float = 1.5, r_max: float = 12.0
):
    """Exclusion radius (arcsec) around a Gaia source of G magnitude ``gmag`` (ASSUMPTION)."""
    g = np.where(np.isfinite(gmag), gmag, g0)
    return np.clip(r0 * 10 ** (0.2 * (g0 - g)), r_min, r_max)


def near_bright(ra, dec, gaia: Table, **kw) -> np.ndarray:
    """Index into ``gaia`` of the bright source each position falls near, or -1."""
    out = np.full(len(ra), -1, int)
    if len(ra) == 0 or len(gaia) == 0:
        return out
    pos = SkyCoord(ra, dec, unit="deg")
    g = SkyCoord(gaia["ra"], gaia["dec"], unit="deg")
    r = exclusion_radius(np.asarray(gaia["gmag"], float), **kw)
    sep = pos[:, None].separation(g[None, :]).arcsec
    inside = sep <= r[None, :]
    hit = inside.any(axis=1)
    out[hit] = np.argmax(inside[hit], axis=1)
    return out


def fetch_gaia(ra: float, dec: float, radius_arcmin: float) -> Table:
    """Gaia DR3 (VizieR I/355/gaiadr3) sources in one cone; columns ra, dec, gmag, source_id."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag"], row_limit=-1)
    res = v.query_region(
        SkyCoord(ra, dec, unit="deg"), radius=radius_arcmin * u.arcmin, catalog="I/355/gaiadr3"
    )
    t = res[0] if len(res) else Table(names=("Source", "RA_ICRS", "DE_ICRS", "Gmag"))
    out = Table(
        {
            "source_id": np.asarray(t["Source"], np.int64),
            "ra": np.asarray(t["RA_ICRS"], float),
            "dec": np.asarray(t["DE_ICRS"], float),
            "gmag": np.ma.filled(np.ma.asarray(t["Gmag"], float), np.nan),
        }
    )
    out.meta = {
        "provenance": "observed",
        "source": "VizieR I/355/gaiadr3 (Gaia DR3, epoch J2016.0)",
    }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--transients", nargs=2, action="append", required=True, metavar=("BAND", "PATH")
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--radius-arcsec", type=float, default=0.3)
    ap.add_argument("--min-bands", type=int, default=2)
    ap.add_argument(
        "--gaia-centre", nargs=2, type=float, metavar=("RA", "DEC"), help="query Gaia DR3 here"
    )
    ap.add_argument("--gaia-radius-arcmin", type=float, default=4.0)
    ap.add_argument(
        "--reference", nargs=3, action="append", default=[], metavar=("RA", "DEC", "NAME")
    )
    args = ap.parse_args(argv)

    tables = {band: Table.read(p) for band, p in args.transients}
    cand = coincident(tables, args.radius_arcsec, args.min_bands)
    cand["near_gaia"] = np.full(len(cand), -1, np.int64)
    if args.gaia_centre:
        gaia = fetch_gaia(*args.gaia_centre, args.gaia_radius_arcmin)
        args.out.mkdir(parents=True, exist_ok=True)
        gaia.write(args.out / "gaia.ecsv", overwrite=True)
        idx = near_bright(np.asarray(cand["ra"]), np.asarray(cand["dec"]), gaia)
        cand["near_gaia"][idx >= 0] = gaia["source_id"][idx[idx >= 0]]
    excluded = cand[cand["near_gaia"] >= 0]
    targets = cand[cand["near_gaia"] < 0]
    for ra, dec, name in args.reference:  # known positions measured alongside, e.g. a lensed star
        targets.add_row(
            {"ra": float(ra), "dec": float(dec), "kind": f"reference_{name}", "n_bands": 0}
        )
    targets.meta.update(provenance="derived", n_excluded_near_gaia=len(excluded))
    args.out.mkdir(parents=True, exist_ok=True)
    targets.write(args.out / "targets.ecsv", overwrite=True)
    excluded.write(args.out / "excluded.ecsv", overwrite=True)
    kinds = {
        k: int((np.asarray(cand["kind"]).astype(str) == k).sum())
        for k in ("variable", "appeared", "disappeared")
    }
    print(
        json.dumps(
            {
                "coincident": len(cand),
                **kinds,
                "excluded_near_gaia": len(excluded),
                "targets": len(targets),
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
