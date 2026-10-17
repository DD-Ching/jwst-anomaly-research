"""ASTRODEEP Frontier Fields photo-z (Di Criscienzo+2017, J/A+A/607/A30) as an eazy-style zout.

MACS0717 has no DJA photo-z (v7 index, 2026-10-17), so the W1 shear screen takes the HST-based
ASTRODEEP-FF redshifts: positions from table ``m0717cla``, ``zbest`` / ``e_zbest``
(semi-interquartile range, SIQR) from ``m0717clz``, joined on ``ID``.
The zout columns ``load_field`` reads are

- ``z_phot`` = ``zbest`` (a spectroscopic z where ``f_zbest`` = 1, as published);
- ``z160``, ``z840`` = ``zbest`` -/+ 1.483 SIQR (Gaussian 16/84 % points; ASSUMPTION), floored at 0;
- NaN for rows with ``relFlag`` != 1 (the catalogue's robustness flag; ASSUMPTION).

    python scripts/astrodeep_zout.py --field m0717 \\
        --out <data root>/cache/external/astrodeep/m0717clu_zout.fits
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, join

SIQR_TO_SIGMA = 1.483  # 1 / 0.6745: Gaussian sigma per semi-interquartile range (ASSUMPTION)
CATALOG = "J/A+A/607/A30"


def fetch(field: str) -> tuple[Table, Table]:
    from astroquery.vizier import Vizier

    viz = Vizier(columns=["**"], row_limit=-1)
    pos = viz.get_catalogs(f"{CATALOG}/{field}cla")[0]
    z = viz.get_catalogs(f"{CATALOG}/{field}clz")[0]
    return pos, z


def to_zout(pos: Table, z: Table) -> Table:
    """Join positions and redshifts on ``ID`` and convert to zout columns (``derived``)."""
    t = join(pos["ID", "RAJ2000", "DEJ2000"], z["ID", "zbest", "e_zbest", "relFlag"], keys="ID")
    zb = np.ma.filled(np.ma.asarray(t["zbest"], float), np.nan)
    siqr = np.ma.filled(np.ma.asarray(t["e_zbest"], float), np.nan)
    rel = np.ma.filled(np.ma.asarray(t["relFlag"], float), 0.0)
    good = np.isfinite(zb) & (zb > 0) & np.isfinite(siqr) & (rel == 1)
    half = SIQR_TO_SIGMA * siqr
    out = Table(
        {
            "id": np.asarray(t["ID"]),
            "ra": np.asarray(t["RAJ2000"], float),
            "dec": np.asarray(t["DEJ2000"], float),
            "z_phot": np.where(good, zb, np.nan),
            "z160": np.where(good, np.maximum(zb - half, 0.0), np.nan),
            "z840": np.where(good, zb + half, np.nan),
        }
    )
    out.meta.update(provenance="derived", source=CATALOG, SIQRSIG=SIQR_TO_SIGMA)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--field", default="m0717", choices=("m0717", "m1149"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    pos, z = fetch(args.field)
    out = to_zout(pos, z)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.write(args.out, overwrite=True)
    sha = hashlib.sha256(
        b"".join(np.ascontiguousarray(out[c]).tobytes() for c in out.colnames)
    ).hexdigest()
    n_ok = int(np.isfinite(out["z_phot"]).sum())
    print(
        f"{args.field}: {len(pos)} positions, {len(z)} redshifts, {len(out)} joined, {n_ok} usable"
    )
    print(f"data_sha256 (column bytes) {sha}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
