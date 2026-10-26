"""Post hoc checks of the MACS J1931 6882 o038 → 9493 o001 grid hits (not in the frozen plan).

1. Per hit group (first member): contributing exposures (``CON``) in each epoch and band, and the
   epoch-2 ``_i2d`` peak fraction (functions of ``m1931_grid_posthoc.py``).
2. Distance to the nearest bright star: grid rows whose epoch-1 F444W flux is above the 99.5th
   percentile of all measured rows (ASSUMPTION; a few saturated-core stars in this field).
3. The roll angle of both epochs (``PA_V3`` of the F444W ``SCI`` header): diffraction spikes rotate
   by the difference.
4. Known case: SN Eos 101.1 / 101.2, nearest grid row and its significances.

All numbers are ``derived``.

    python scripts/m1931_9493_posthoc.py --data data/m1931_9493 \\
        --out results/transients/macs1931_grid_9493/posthoc.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m1931_grid_posthoc import BANDS, _local_pf, _num, _pix, known_case, n_exposures

E1 = "jw06882-o038_t041_nircam_clear-"
E2 = "jw09493-o001_t001_nircam_clear-"
BRIGHT_PERCENTILE = 99.5


def nearest_bright(
    ra: float, dec: float, bra: np.ndarray, bdec: np.ndarray
) -> tuple[int, float] | None:
    """Index and distance (arcsec) of the nearest bright row to (ra, dec); None if there is none."""
    if len(bra) == 0:
        return None
    d = np.hypot((np.asarray(bra) - ra) * math.cos(math.radians(dec)), np.asarray(bdec) - dec)
    k = int(np.argmin(d))
    return k, float(d[k] * 3600)


def bright_rows(grid: Table, band: str = "F444W", pct: float = BRIGHT_PERCENTILE) -> Table:
    """Grid rows whose epoch-1 flux in ``band`` is finite and above the ``pct`` percentile."""
    f = np.asarray(grid[f"{band}_flux1"], dtype=float)
    ok = np.isfinite(f)
    if not ok.any():
        return grid[:0]
    return grid[ok & (f > np.percentile(f[ok], pct))]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    hits = Table.read(a.data / "out" / "hits.ecsv")
    grid = Table.read(a.data / "out" / "grid.ecsv")
    bright = bright_rows(grid)
    images = {}
    for b in BANDS:
        for e, pre in (("e1", E1), ("e2", E2)):
            h = fits.open(a.data / f"{pre}{b}_i2d.fits", memmap=True)
            images[(b, e)] = (h["SCI"].data, h["CON"].data, WCS(h["SCI"].header), h)
    groups = []
    for g in sorted({int(v) for v in hits["group"]}):
        members = hits[hits["group"] == g]
        r = members[0]
        row = {
            "group": g,
            "pos_id": str(r["pos_id"]),
            "members": [str(p) for p in members["pos_id"]],
        }
        nb = nearest_bright(float(r["ra"]), float(r["dec"]), bright["ra"], bright["dec"])
        if nb is not None:
            row["nearest_bright_arcsec"] = round(nb[1], 2)
            row["nearest_bright_f444w_ujy"] = round(float(bright["F444W_flux1"][nb[0]]), 1)
        for b in BANDS:
            for e in ("e1", "e2"):
                sci, con, wcs, _ = images[(b, e)]
                x, y = _pix(wcs, float(r["ra"]), float(r["dec"]))
                row[f"{b.upper()}_n_exposures_{e}"] = n_exposures(con, x, y)
                if e == "e2":
                    row[f"{b.upper()}_peak_fraction_e2"] = _local_pf(sci, x, y)
        groups.append(row)
    pa = {e: images[("f444w", e)][3]["SCI"].header.get("PA_V3") for e in ("e1", "e2")}
    sizes = {
        p.name: p.stat().st_size
        for e in (E1, E2)
        for b in BANDS
        for p in [a.data / f"{e}{b}_i2d.fits"]
    }
    out = {
        "label": "derived, post hoc (not in PLAN.md)",
        "mosaic_bytes": sizes,
        "grid_ecsv_bytes": (a.data / "out" / "grid.ecsv").stat().st_size,
        "bright_rule": (
            f"epoch-1 F444W flux above the {BRIGHT_PERCENTILE} percentile of measured rows"
        ),
        "bright_min_flux_ujy": _num(np.min(bright["F444W_flux1"])) if len(bright) else None,
        "n_bright": len(bright),
        "pa_v3_deg": {k: _num(v) if v is not None else None for k, v in pa.items()},
        "known_case": known_case(grid),
        "groups": groups,
    }
    for _, _, _, h in images.values():
        h.close()
    a.out.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
