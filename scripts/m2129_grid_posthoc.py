"""Post hoc checks of the MACS J2129 two-epoch grid hits (not in the frozen plan).

1. Coverage: for every hit group and band, the number of contributing exposures (bits set in
   ``CON``) and the relative weight (``WHT`` over the band's median positive weight) per epoch.
2. Nearest bright star: distance and position angle to the nearest grid row whose epoch-2 F444W
   flux is above the 99.5th percentile of measured rows (ASSUMPTION, as in the MACS J0308 post hoc).
3. Grid coverage: per band, the number of measured grid rows by epoch-1 exposure count (a
   diagnostic for gate G1, which counts injections at measured rows).
4. ``classify`` reads each group from 1–2: star wings, single-exposure coverage, and the
   fixed-offset star image of the MACS J0308 grid (8.33″ at 120.1° measured at ``PA_APER`` 75.12°),
   rotated here by this pair's epoch-1 ``PA_APER`` − 75.12°. Thresholds are post hoc ASSUMPTIONs.

All numbers are ``derived``.

    python scripts/m2129_grid_posthoc.py --data data/m2129 --grid <out>/grid.ecsv \\
        --hits <out>/hits.ecsv --out results/transients/macs2129_grid/posthoc.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m0308_grid_posthoc import GHOST_PA, GHOST_SEP, STAR_NEAR, STAR_PCT, n_exposures

E1 = "jw04903-o004_t004_nircam_clear-"
E2 = "jw06882-o040_t043_nircam_clear-"
BANDS = ("f277w", "f444w")
GHOST_PA_REF = 75.12  # epoch-1 PA_APER of the MACS J0308 pair, where GHOST_PA was measured


def first_rows(hits: Table) -> list:
    """One row per group: its first member."""
    seen: set[int] = set()
    out = []
    for r in hits:
        if int(r["group"]) not in seen:
            seen.add(int(r["group"]))
            out.append(r)
    return out


def coverage(data: Path, hits: Table) -> dict:
    out: dict = {}
    rows = first_rows(hits)
    for b in BANDS:
        for tag, stem in (("e1", E1), ("e2", E2)):
            with fits.open(data / f"{stem}{b}_i2d.fits", memmap=True) as h:
                w = WCS(h["SCI"].header)
                con, wht = h["CON"].data, h["WHT"].data
                sub = wht[:: max(wht.shape[0] // 32, 1)]
                ref = float(np.median(sub[sub > 0]))
                for r in rows:
                    x, y = w.world_to_pixel_values(float(r["ra"]), float(r["dec"]))
                    xi, yi = int(round(float(x))), int(round(float(y)))
                    inside = 0 <= yi < wht.shape[0] and 0 <= xi < wht.shape[1]
                    rec = out.setdefault(str(int(r["group"])), {})
                    rec[f"{b}_{tag}_n_exp"] = n_exposures(con, yi, xi) if inside else 0
                    rec[f"{b}_{tag}_wht_rel"] = (
                        round(float(wht[yi, xi]) / ref, 3) if inside else 0.0
                    )
    return out


def grid_coverage(data: Path, grid: Table) -> dict:
    out: dict = {}
    for b in BANDS:
        ok = np.isfinite(np.asarray(grid[f"{b.upper()}_flux1"], float))
        with fits.open(data / f"{E1}{b}_i2d.fits", memmap=True) as h:
            w, con = WCS(h["SCI"].header), h["CON"].data
            x, y = w.world_to_pixel_values(np.asarray(grid["ra"])[ok], np.asarray(grid["dec"])[ok])
            n = []
            for xi, yi in zip(x, y, strict=True):
                xi, yi = int(round(float(xi))), int(round(float(yi)))
                inside = 0 <= yi < con.shape[1] and 0 <= xi < con.shape[2]
                n.append(n_exposures(con, yi, xi) if inside else 0)
        vals, counts = np.unique(n, return_counts=True)
        out[b] = {str(int(v)): int(c) for v, c in zip(vals, counts, strict=True)}
    return out


def bright_stars(grid: Table, hits: Table) -> dict:
    f = np.asarray(grid["F444W_flux2"], float)
    cut = np.nanpercentile(f[np.isfinite(f)], STAR_PCT)
    stars = grid[np.isfinite(f) & (f >= cut)]
    sc = SkyCoord(stars["ra"], stars["dec"], unit="deg")
    out = {"flux_cut_uJy": round(float(cut), 2), "n_stars": len(stars), "groups": {}}
    if len(stars) == 0:
        return out
    for r in first_rows(hits):
        c = SkyCoord(r["ra"], r["dec"], unit="deg")
        d = c.separation(sc).arcsec
        d[d <= 0.05] = np.inf  # the hit itself when it is bright
        i = int(np.argmin(d))
        out["groups"][str(int(r["group"]))] = {
            "sep_arcsec": round(float(d[i]), 2),
            "pa_deg": round(float(sc[i].position_angle(c).deg), 1),
            "star_ra": round(float(stars["ra"][i]), 6),
            "star_dec": round(float(stars["dec"][i]), 6),
            "star_f444w_e2_uJy": round(float(stars["F444W_flux2"][i]), 1),
        }
    return out


def classify(cov: dict, star: dict, pa_e1: float) -> dict:
    """Every matching label per group; the fixed-offset vector is rotated to this epoch-1 roll."""
    ghost_pa = (GHOST_PA + pa_e1 - GHOST_PA_REF) % 360.0
    out = {}
    for g, c in cov.items():
        s = star["groups"].get(g)
        n_exp = min(v for k, v in c.items() if k.endswith("_n_exp"))
        labels = []
        if n_exp == 0:
            labels.append("off_footprint")
        if s and s["sep_arcsec"] < STAR_NEAR:
            labels.append("star_wings")
        if s and abs(s["sep_arcsec"] - GHOST_SEP) < 0.2:
            dpa = abs((s["pa_deg"] - ghost_pa + 180.0) % 360.0 - 180.0)
            if dpa < 1.0:
                labels.append("fixed_offset_star_image")
        if n_exp == 1:
            labels.append("single_exposure")
        out[g] = labels or ["unexplained"]
    return out


def pa_aper(data: Path) -> dict:
    return {
        tag: round(float(fits.getheader(data / f"{stem}f277w_i2d.fits", "SCI")["PA_APER"]), 2)
        for tag, stem in (("e1", E1), ("e2", E2))
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--grid", type=Path, required=True)
    ap.add_argument("--hits", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    hits, grid = Table.read(a.hits), Table.read(a.grid)
    cov, star, pa = coverage(a.data, hits), bright_stars(grid, hits), pa_aper(a.data)
    res = {
        "provenance": "derived",
        "coverage": cov,
        "bright_star": star,
        "pa_aper": pa,
        "ghost_pa_rotated_deg": round((GHOST_PA + pa["e1"] - GHOST_PA_REF) % 360.0, 2),
        "classification": classify(cov, star, pa["e1"]),
        "grid_e1_n_exp": grid_coverage(a.data, grid),
        "mosaic_bytes": {p.name: p.stat().st_size for p in sorted(a.data.glob("*_i2d.fits"))},
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res["classification"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
