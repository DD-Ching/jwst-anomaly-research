"""ASTRODEEP Frontier Fields photo-z (Di Criscienzo+2017, J/A+A/607/A30) as an eazy-style zout.

MACS0717 has no DJA photo-z (v7 index, 2026-10-17), so the W1 shear screen takes the HST-based
ASTRODEEP-FF redshifts: positions from table ``m0717cla``, ``zbest`` / ``e_zbest``
(semi-interquartile range, SIQR) from ``m0717clz``, joined on ``ID``.
The zout columns ``load_field`` reads are

- ``z_phot`` = ``zbest`` (a spectroscopic z where ``f_zbest`` = 1, as published);
- ``z160``, ``z840`` = ``zbest`` -/+ 1.483 SIQR (Gaussian 16/84 % points; ASSUMPTION); z160 is
  floored at ``Z_FLOOR`` (> 0, so ``attach_photoz`` keeps it and the row counts as possibly
  foreground) but never exceeds ``zbest``;
- rows with ``relFlag`` != 1 (the catalogue's photo-z robustness flag; ASSUMPTION) and no
  spectroscopic z, or with no zbest, are dropped, so the nearest-neighbour match cannot land on
  them and hide a usable row (it may then reach a neighbour within 0.3 arcsec instead).

    python scripts/astrodeep_zout.py --field m0717 \\
        --out <data root>/cache/external/astrodeep/m0717clu_zout.fits
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, join

SIQR_TO_SIGMA = 1.483  # 1 / 0.6745: Gaussian sigma per semi-interquartile range (ASSUMPTION)
Z_FLOOR = 0.01  # lower bound of z160 (ASSUMPTION; attach_photoz turns values <= 0 into NaN)
CATALOG = "J/A+A/607/A30"
ZOUT_COLUMNS = ("id", "ra", "dec", "z_phot", "z160", "z840")


def fetch(field: str) -> tuple[Table, Table]:
    from astroquery.vizier import Vizier

    viz = Vizier(columns=["**"], row_limit=-1)
    pos = viz.get_catalogs(f"{CATALOG}/{field}cla")[0]
    z = viz.get_catalogs(f"{CATALOG}/{field}clz")[0]
    return pos, z


def to_zout(pos: Table, z: Table) -> Table:
    """Join positions and redshifts on ``ID`` and convert to zout columns (``derived``)."""
    cols = ["ID", "zbest", "e_zbest", "relFlag", "f_zbest"]
    t = join(pos["ID", "RAJ2000", "DEJ2000"], z[cols], keys="ID")
    zb = np.ma.filled(np.ma.asarray(t["zbest"], float), np.nan)
    siqr = np.ma.filled(np.ma.asarray(t["e_zbest"], float), np.nan)
    rel = np.ma.filled(np.ma.asarray(t["relFlag"], float), 0.0)
    spec = np.ma.filled(np.ma.asarray(t["f_zbest"], float), 0.0) == 1
    siqr = np.where(spec, 0.0, siqr)  # a spectroscopic z has no range
    good = np.isfinite(zb) & (zb > 0) & np.isfinite(siqr) & ((rel == 1) | spec)
    t, zb, half = t[good], zb[good], SIQR_TO_SIGMA * siqr[good]
    out = Table(
        {
            "id": np.asarray(t["ID"]),
            "ra": np.asarray(t["RAJ2000"], float),
            "dec": np.asarray(t["DEJ2000"], float),
            "z_phot": zb,
            "z160": np.minimum(np.maximum(zb - half, Z_FLOOR), zb),  # never above z_phot
            "z840": zb + half,
        }
    )
    out.meta.update(provenance="derived", source=CATALOG, SIQRSIG=SIQR_TO_SIGMA)
    return out


def content_sha256(t: Table) -> str:
    """sha256 of the zout columns' values (little-endian, rows in file order), not the FITS bytes.

    A new astropy or FITS header layout leaves it unchanged; a changed value, row or column fails.
    """
    h = hashlib.sha256()
    for name in ZOUT_COLUMNS:
        dtype = "<i8" if name == "id" else "<f8"
        h.update(name.encode() + b"\0")
        h.update(np.ascontiguousarray(np.asarray(t[name]), dtype=dtype).tobytes())
    return h.hexdigest()


class HashMismatch(RuntimeError):
    """The fetched table's ``content_sha256`` differs from the pinned value; nothing was written."""


def write(field: str, path: Path, expect: str | None = None) -> str:
    """Fetch, convert and write ``path``; returns ``content_sha256``.

    With ``expect``, a table whose hash differs raises ``HashMismatch`` before ``path`` is touched,
    so only a verified file reaches the final name (written via a ``.part`` file).
    """
    pos, z = fetch(field)
    out = to_zout(pos, z)
    got = content_sha256(out)
    if expect is not None and got != expect:
        raise HashMismatch(f"content sha256 {got} != {expect}")
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    try:
        out.write(part, format="fits", overwrite=True)
        os.replace(part, path)
    finally:
        part.unlink(missing_ok=True)
    print(f"{field}: {len(pos)} positions, {len(z)} redshifts, {len(out)} usable rows")
    return got


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--field", default="m0717", choices=("m0717",))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    print(f"content sha256 {write(args.field, args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
