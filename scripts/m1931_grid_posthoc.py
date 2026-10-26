"""Post hoc sharpness and coverage check of the MACS J1931 grid hits (not in the frozen plan).

For every hit row of ``transient_grid.py`` and each LW band, reads the epoch-2 ``_i2d`` at the
hit and reports (a) the peak fraction: the brightest pixel over the background-subtracted 5 × 5
pixel sum (a PSF-sampled point source sits near the reference sources' value, a single hot pixel
or cosmic ray near 1), and (b) the number of exposures contributing there (bits set in the
``CON`` context image). Reference: the same peak fraction for up to 300 catalogued grid
positions in the band's brightest 20 % of epoch-2 fluxes. All numbers are ``derived``.

    python scripts/m1931_grid_posthoc.py --data data/m1931 \\
        --out results/transients/macs1931_grid/sharp.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

E2 = "jw12299-o001_t001_nircam_clear-"
BANDS = ("f277w", "f356w", "f444w")


def peak_fraction(sci: np.ndarray, x: int, y: int) -> float:
    """Brightest pixel over the background-subtracted 5 × 5 sum at (x, y); NaN off the image."""
    stamp = sci[y - 2 : y + 3, x - 2 : x + 3].astype(float)
    ring = sci[max(y - 15, 0) : y + 16, max(x - 15, 0) : x + 16]
    if stamp.shape != (5, 5) or not np.isfinite(stamp).all():
        return float("nan")
    stamp = stamp - np.nanmedian(ring)
    total = stamp.sum()
    return float(stamp.max() / total) if total > 0 else float("nan")


def n_exposures(con: np.ndarray, x: int, y: int) -> int:
    """Exposures contributing at (x, y): bits set across all context planes."""
    planes = con[:, y, x] if con.ndim == 3 else np.atleast_1d(con[y, x])
    return int(sum(bin(int(v) & 0xFFFFFFFF).count("1") for v in planes))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    hits = Table.read(a.data / "out" / "hits.ecsv")
    grid = Table.read(a.data / "out" / "grid.ecsv")
    out: dict = {"label": "derived, post hoc (not in PLAN.md)"}
    for b in BANDS:
        with fits.open(a.data / f"{E2}{b}_i2d.fits", memmap=True) as h:
            sci, con, wcs = h["SCI"].data, h["CON"].data, WCS(h["SCI"].header)

            def pix(ra: float, dec: float, wcs: WCS = wcs) -> tuple[int, int]:
                x, y = wcs.world_to_pixel_values(ra, dec)
                return int(round(float(x))), int(round(float(y)))

            rows = []
            for r in hits:
                x, y = pix(r["ra"], r["dec"])
                rows.append(
                    {
                        "group": int(r["group"]),
                        "pos_id": str(r["pos_id"]),
                        "peak_fraction": round(peak_fraction(sci, x, y), 3),
                        "n_exposures_e2": n_exposures(con, x, y),
                    }
                )
            f2 = np.asarray(grid[f"{b.upper()}_flux2"], dtype=float)
            ok = np.isfinite(f2) & (f2 > np.nanpercentile(f2, 80))
            ref = np.array([peak_fraction(sci, *pix(r["ra"], r["dec"])) for r in grid[ok][:300]])
            ref = ref[np.isfinite(ref)]
        out[b.upper()] = {
            "hits": rows,
            "reference": {
                "n": int(ref.size),
                "p50": round(float(np.median(ref)), 3),
                "p99": round(float(np.percentile(ref, 99)), 3),
                "max": round(float(ref.max()), 3),
            },
        }
    a.out.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
