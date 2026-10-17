"""Does PSF anisotropy measured from stars set the W1 shear-screen B-mode floor? (Abell 2744)

Plan and frozen rules: ``results/w1_psf/PLAN.md``. Stars are catalogue rows with a Gaia DR3 match
(G < ``GAIA_GMAX``, within ``MATCH_ARCSEC``), S/N > ``STAR_MIN_SNR`` and mag >= ``STAR_MIN_MAG``.
Their observed-moment distortion χ★ gives a PSF anisotropy; the screen's isotropic Gaussian
deconvolution is replaced by an anisotropic one (``exotic_screens.intrinsic_ellipticity``
``psf_chi``), constant or, when rule R1 allows a map, inverse-distance weighted over the nearest
``MAP_K`` stars. Rule R2 decides whether the correction goes on to injections. All outputs are
``derived``; thresholds are ASSUMPTIONs.

    python scripts/w1_psf_floor.py --stars gaia.ecsv --work <dir> \
        --out results/w1_psf/psf_abell2744.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject_shear as ish  # noqa: E402

es = ish.es
GAIA_GMAX = 21.0  # ASSUMPTION
MATCH_ARCSEC = 0.3  # ASSUMPTION
STAR_MIN_SNR = 50.0  # ASSUMPTION
STAR_MIN_MAG = 19.5  # ASSUMPTION: brighter stars' isophotal moments are spike-dominated
R1_MAX_MEDIAN_SEP = 10.0  # arcsec (the aperture radius), ASSUMPTION
R1_MIN_STARS = 10  # ASSUMPTION
R2_MAX_BMAX = 3.79  # 10 % below Abell 2744's 4.21, ASSUMPTION
MAP_K = 5  # ASSUMPTION
FLOOR_CENTRE = (88.0, -42.0)  # model-frame arcsec, results/w1_floor


def moments_chi(a, b, pa_deg) -> np.ndarray:
    """Distortion χ = (a² - b²)/(a² + b²) exp(2i PA) of observed moments (PA east of north), the
    convention of ``exotic_screens.intrinsic_ellipticity``."""
    a2, b2 = np.asarray(a, float) ** 2, np.asarray(b, float) ** 2
    return (a2 - b2) / (a2 + b2) * np.exp(2j * np.deg2rad(np.asarray(pa_deg, float)))


def r1_map_feasible(nearest_sep_arcsec, n_stars: int) -> bool:
    """R1: a star map resolves the aperture scale only with >= R1_MIN_STARS stars and a median
    used-source-to-nearest-star distance <= R1_MAX_MEDIAN_SEP (both must hold)."""
    sep = np.asarray(nearest_sep_arcsec, float)
    return bool(n_stars >= R1_MIN_STARS and np.median(sep) <= R1_MAX_MEDIAN_SEP)


def r2_adopt(b_max: float) -> bool:
    """R2: carry the correction to injections only if the field's max |S_×| is <= R2_MAX_BMAX."""
    return bool(np.isfinite(b_max) and b_max <= R2_MAX_BMAX)


R3_MASSES = ("2e+12", "8e+12", "2e+13")


def r3_adopt(new_limits: dict, old_limits: dict, masses=R3_MASSES) -> bool:
    """R3: adopt only if the recombined 95 % limit is strictly lower at every mass in ``masses``
    (``limits_95`` dicts of ``inject_radial.combine``); a missing or infinite limit fails."""
    for m in masses:
        a = new_limits.get(m, {}).get("upper_limit_deg2")
        b = old_limits.get(m, {}).get("upper_limit_deg2")
        if a is None or b is None or not (np.isfinite(a) and a < b):
            return False
    return True


def idw_map(xs, ys, chi_s, x, y, k: int = MAP_K) -> np.ndarray:
    """Inverse-distance mean of the ``k`` nearest stars' χ at each (x, y)."""
    from scipy.spatial import cKDTree

    k = min(k, len(xs))
    d, i = cKDTree(np.c_[xs, ys]).query(np.c_[x, y], k=k)
    if k == 1:
        d, i = d[:, None], i[:, None]
    w = 1.0 / np.maximum(d, 1e-3)
    return (w * np.asarray(chi_s)[i]).sum(1) / w.sum(1)


def select_stars(shapes, gaia) -> np.ndarray:
    """Indices of catalogue rows that are PSF stars under the frozen rule."""
    from astropy.coordinates import SkyCoord

    g = gaia[np.asarray(gaia["mag"], float) < GAIA_GMAX]
    if not len(g):
        return np.array([], int)
    c = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    gc = SkyCoord(g["ra"], g["dec"], unit="deg")
    idx, d, _ = gc.match_to_catalog_sky(c)
    rows = np.unique(idx[d.arcsec < MATCH_ARCSEC])
    snr = np.asarray(shapes["snr"], float)[rows]
    mag = np.asarray(shapes["mag"], float)[rows]
    return rows[(snr > STAR_MIN_SNR) & (mag >= STAR_MIN_MAG)]


def screen_summary(model, shapes, sargs, psf, extra, label: str) -> dict:
    res = es.shear_screen(model, shapes, sargs, psf, extra)
    sx = res["s_cross"]
    k = int(np.nanargmax(np.abs(sx)))
    gx, gy = res["gx"].ravel(), res["gy"].ravel()
    kf = int(np.nanargmin(np.hypot(gx - FLOOR_CENTRE[0], gy - FLOOR_CENTRE[1])))
    null = np.asarray(res["null_max"], float)
    return {
        "label": label,
        "n_used": int(res["use"].sum()),
        "responsivity": float(res["counts"]["responsivity"]),
        "b_max": float(np.nanmax(np.abs(sx))),
        "b_max_at": [float(gx[k]), float(gy[k])],
        "s_cross_at_floor_centre": float(sx.ravel()[kf]),
        "s_max": float(np.nanmax(res["s"])),
        "null_max_p50_p95_p99": [float(np.percentile(null, q)) for q in (50, 95, 99)],
        "e_used_mean": [
            float(np.mean(res["e"][res["use"]].real)),
            float(np.mean(res["e"][res["use"]].imag)),
        ],
        "r2_adopt": r2_adopt(float(np.nanmax(np.abs(sx)))),
    }


def run(name: str, stars_path: Path, work: Path) -> dict:
    from astropy.table import Table

    spec = ish.field_spec(name)
    (work / name).mkdir(parents=True, exist_ok=True)
    model, shapes, extra, _, _ = ish.ir.load_field(name, work / name, spec)
    sargs = SimpleNamespace(**ish.SHEAR_DEFAULTS)
    sargs.max_radius = spec["max_radius"]
    psf = es.psf_sigma_px(shapes)
    gaia = Table.read(stars_path)
    rows = select_stars(shapes, gaia)
    a = np.asarray(shapes["semimajor_px"], float)
    b = a * (1.0 - np.asarray(shapes["ellipticity"], float))
    chi_all = moments_chi(a, b, shapes["pa_obs"])
    chi_s = chi_all[rows]
    xa, ya = model.to_frame(np.asarray(shapes["ra"]), np.asarray(shapes["dec"]))
    chi_med = complex(np.median(chi_s.real), np.median(chi_s.imag)) if len(rows) else 0j

    base = screen_summary(model, shapes, sargs, psf, extra, "isotropic (committed screen)")
    use0 = es.shear_screen(model, shapes, sargs, psf, extra)["use"]
    from scipy.spatial import cKDTree

    sep = (
        cKDTree(np.c_[xa[rows], ya[rows]]).query(np.c_[xa[use0], ya[use0]])[0]
        if len(rows)
        else np.full(int(use0.sum()), np.inf)
    )
    feasible = r1_map_feasible(sep, len(rows))
    sargs.psf_chi = chi_med
    const = screen_summary(model, shapes, sargs, psf, extra, "constant chi_PSF (median of stars)")
    out = {
        "field": name,
        "provenance": "derived",
        "psf_sigma_px": psf,
        "stars": [
            {
                "row": int(r),
                "ra": float(shapes["ra"][r]),
                "dec": float(shapes["dec"][r]),
                "x": float(xa[r]),
                "y": float(ya[r]),
                "mag": float(shapes["mag"][r]),
                "snr": float(shapes["snr"][r]),
                "chi": [float(c.real), float(c.imag)],
                "abs_chi": float(abs(c)),
                "pa_deg": float(np.mod(0.5 * np.degrees(np.angle(c)), 180.0)),
            }
            for r, c in zip(rows, chi_s, strict=True)
        ],
        "chi_psf_median": [chi_med.real, chi_med.imag],
        "r1": {
            "n_stars": int(len(rows)),
            "median_nearest_star_arcsec": float(np.median(sep)),
            "p90_nearest_star_arcsec": float(np.percentile(sep, 90)),
            "feasible": feasible,
        },
        "screens": [base, const],
        "rules": {
            "gaia_gmax": GAIA_GMAX,
            "match_arcsec": MATCH_ARCSEC,
            "star_min_snr": STAR_MIN_SNR,
            "star_min_mag": STAR_MIN_MAG,
            "r1_max_median_sep_arcsec": R1_MAX_MEDIAN_SEP,
            "r1_min_stars": R1_MIN_STARS,
            "r2_max_bmax": R2_MAX_BMAX,
            "map_k": MAP_K,
        },
    }
    if feasible:
        sargs.psf_chi = idw_map(xa[rows], ya[rows], chi_s, xa, ya)
        out["screens"].append(
            screen_summary(model, shapes, sargs, psf, extra, f"IDW map of {MAP_K} nearest stars")
        )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--field", default="abell2744")
    ap.add_argument("--stars", type=Path, required=True, help="gaia_stars.py ECSV (G < 21)")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    res = run(a.field, a.stars, a.work)
    a.out.write_text(json.dumps(es.json_safe(res), indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("chi_psf_median", "r1")}, indent=1))
    for s in res["screens"]:
        print(s)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
