"""Post hoc checks of the Abell 1689 two-epoch grid (results/transients/a1689_grid; not planned).

1. Hit groups: coverage (``CON`` exposure count, relative ``WHT``) per band and epoch, the nearest
   bright star (epoch-2 F356W flux above the 99.5th percentile of measured rows) and the
   ``m2129_grid_posthoc.classify`` labels (star wings, single exposure, off footprint, the rotated
   fixed-offset star image of #576). Thresholds are the post hoc ASSUMPTIONs of #576 / #577.
2. G1 split as in #578 (``grid_g1_coverage_posthoc.split_recovery``) and each band's floor
   significance at median errors.
3. Recovery of the frac-1.0 injections by their expected calibrated significance: the injected
   flux over the error ``transient_forced.compare`` assigns (both ERR terms, epoch 1 on the
   epoch-2 zero point, the 0.05 mag systematic floor on the brighter epoch), divided by the band's
   noise scale; per injection the smaller of the two bands. This tests whether the misses are
   injections the 5σ threshold cannot reach, rather than a tool failure.

All numbers are ``derived``.

    python scripts/a1689_grid_posthoc.py --data data/a1689_grid --run results/transients/a1689_grid
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transient_grid as tg  # noqa: E402
from grid_g1_coverage_posthoc import floor_significance, split_recovery  # noqa: E402
from m0308_grid_posthoc import STAR_PCT  # noqa: E402
from m2129_grid_posthoc import GHOST_PA, GHOST_PA_REF, classify, first_rows  # noqa: E402

E1 = "jw01840-o019_t009_nircam_clear-"
E2 = "jw05782-o001_t001_nircam_clear-"
BANDS = ("F277W", "F356W")
STAR_BAND = "F356W"
SYS_FLOOR = 0.05  # transient_grid default (mag), used by the frozen run
Z_EDGES = (0.0, 3.0, 4.0, 5.0, 6.0, 8.0, np.inf)  # post hoc bins of expected significance


def expected_sigma(f1z, e1z, f2, e2, inj, scale: float, sys_floor: float = SYS_FLOOR):
    """Calibrated significance an injected flux ``inj`` reaches with no noise draw."""
    floor = sys_floor * np.log(10) / 2.5 * np.fmax(np.abs(f1z), np.abs(f2 + inj))
    return inj / np.sqrt(e1z**2 + e2**2 + floor**2) / scale


def recovery_by_bins(z: np.ndarray, rec: np.ndarray, edges=Z_EDGES) -> list[dict]:
    out = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        m = (z >= lo) & (z < hi)
        out.append(
            {
                "z_lo": lo,
                "z_hi": None if np.isinf(hi) else hi,
                "n": int(m.sum()),
                "recovered_frac": float(rec[m].mean()) if m.any() else None,
            }
        )
    return out


def coverage(data: Path, rows: list) -> dict:
    out: dict = {}
    for b in BANDS:
        for tag, stem in (("e1", E1), ("e2", E2)):
            with fits.open(data / f"{stem}{b.lower()}_i2d.fits", memmap=True) as h:
                w, wht = WCS(h["SCI"].header), h["WHT"].data
                sub = wht[:: max(wht.shape[0] // 32, 1)]
                ref = float(np.median(sub[sub > 0]))
                ra = np.array([float(r["ra"]) for r in rows])
                dec = np.array([float(r["dec"]) for r in rows])
                n_exp = tg.n_exposures(str(data / f"{stem}{b.lower()}_i2d.fits"), ra, dec)
                for r, n in zip(rows, n_exp, strict=True):
                    x, y = w.world_to_pixel_values(float(r["ra"]), float(r["dec"]))
                    xi, yi = int(round(float(x))), int(round(float(y)))
                    inside = 0 <= yi < wht.shape[0] and 0 <= xi < wht.shape[1]
                    rec = out.setdefault(str(int(r["group"])), {})
                    rec[f"{b.lower()}_{tag}_n_exp"] = int(n) if inside else 0
                    rec[f"{b.lower()}_{tag}_wht_rel"] = (
                        round(float(wht[yi, xi]) / ref, 3) if inside else 0.0
                    )
    return out


def bright_stars(grid: Table, rows: list) -> dict:
    f = np.asarray(grid[f"{STAR_BAND}_flux2"], float)
    cut = np.nanpercentile(f[np.isfinite(f)], STAR_PCT)
    stars = grid[np.isfinite(f) & (f >= cut)]
    sc = SkyCoord(stars["ra"], stars["dec"], unit="deg")
    out = {"band": STAR_BAND, "flux_cut": round(float(cut), 2), "n_stars": len(stars), "groups": {}}
    if len(stars) == 0:
        return out
    for r in rows:
        c = SkyCoord(r["ra"], r["dec"], unit="deg")
        d = c.separation(sc).arcsec
        d[d <= 0.05] = np.inf  # the hit itself when it is bright
        i = int(np.argmin(d))
        out["groups"][str(int(r["group"]))] = {
            "sep_arcsec": round(float(d[i]), 2),
            "pa_deg": round(float(sc[i].position_angle(c).deg), 1),
        }
    return out


def injection_sigma(data: Path, run: Path) -> dict:
    summary = json.loads((run / "summary.json").read_text())
    grid = Table.read(run / "grid.ecsv")
    inj_all = Table.read(run / "injections.ecsv")
    inj = inj_all[np.isclose(inj_all["frac"], 1.0)]
    ra, dec = np.asarray(grid["ra"]), np.asarray(grid["dec"])
    row = {p: n for n, p in enumerate(grid["pos_id"])}
    z, bands, match = [], {}, {}
    for b in BANDS:
        p1, p2 = (data / f"{s}{b.lower()}_i2d.fits" for s in (E1, E2))
        with tg._open(str(p1)) as i1, tg._open(str(p2)) as i2:
            m = tg.measure_band(i1, i2, ra, dec, np.zeros(len(ra), bool), 0.15, None)
        scale = summary["bands"][b]["noise_scale"]
        k = 10 ** (-0.4 * summary["bands"][b]["zero_point_mag"])
        i = [row[p] for p in inj["pos_id"]]  # injections sit at grid rows
        f = np.asarray(inj[f"{b}_inj_flux"], float)
        # the run measured at the refined grid positions (rc None), so this must reproduce them
        match[b] = bool(np.allclose(m["f1"], grid[f"{b}_flux1"], equal_nan=True))
        z.append(expected_sigma(m["f1"][i] * k, m["e1"][i] * k, m["f2"][i], m["e2"][i], f, scale))
        ok = np.isfinite(m["e1"] + m["e2"])
        floor = float(np.min(inj_all[f"{b}_inj_flux"]))
        bands[b] = {
            "floor_flux": floor,
            "median_e1_zp": float(np.median(m["e1"][ok] * k)),
            "median_e2": float(np.median(m["e2"][ok])),
            "floor_sigma_at_median": floor_significance(m["e1"][ok] * k, m["e2"][ok], scale, floor),
        }
    zmin = np.minimum(*z)  # NaN in either band stays NaN and falls in no bin
    rec = np.asarray(inj["recovered"], bool)
    nom = np.asarray(inj["nominal"], bool)
    ge5 = zmin >= 5
    return {
        "bands": bands,
        "remeasure_matches_grid_flux1": match,
        "frac_1.0_n_z_min_nan": int(np.isnan(zmin).sum()),
        "frac_1.0_z_min_bins_all": recovery_by_bins(zmin, rec),
        "frac_1.0_z_min_bins_nominal": recovery_by_bins(zmin[nom], rec[nom]),
        "frac_1.0_n_z_min_ge_5": int(ge5.sum()),
        "frac_1.0_recovered_where_z_min_ge_5": float(rec[ge5].mean()) if ge5.any() else None,
        "split_frac_1.0": split_recovery(inj_all, 1.0),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--run", type=Path, required=True)
    a = ap.parse_args(argv)
    hits, grid = Table.read(a.run / "hits.ecsv"), Table.read(a.run / "grid.ecsv")
    rows = first_rows(hits)
    cov, star = coverage(a.data, rows), bright_stars(grid, rows)
    pa = {
        tag: round(float(fits.getheader(a.data / f"{s}f277w_i2d.fits", "SCI")["PA_APER"]), 2)
        for tag, s in (("e1", E1), ("e2", E2))
    }
    res = {
        "provenance": "derived",
        "post_hoc": True,
        "coverage": cov,
        "bright_star": star,
        "pa_aper": pa,
        "ghost_pa_rotated_deg": round((GHOST_PA + pa["e1"] - GHOST_PA_REF) % 360.0, 2),
        "classification": classify(cov, star, pa["e1"]),
        "injection": injection_sigma(a.data, a.run),
        "mosaic_bytes": {p.name: p.stat().st_size for p in sorted(a.data.glob("*_i2d.fits"))},
    }
    (a.run / "posthoc.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("classification", "injection")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
