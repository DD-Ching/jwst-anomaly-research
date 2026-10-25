"""SN Eos (MACS J1931 system 101): predicted images, magnifications and delays from the VENUS
AstroLensPy v2.1 maps.

Frozen plan: ``results/transients/sneos/PLAN.md``. Gate P checks that ``Lens_pot`` is the
potential of ``Alpha_x/y`` (one fitted scalar); gate S requires both catalogued images to be
matched; gate M compares our |μ| with ``Mu_zs5.13``. Every output number is a model_prediction.

    python scripts/sneos_predict.py --maps DIR [--out DIR]

``DIR`` holds ``Alpha_x.fits.zip``, ``Alpha_y.fits.zip``, ``Lens_pot.fits.zip`` and
``Mu_zs5.13.fits.zip`` from Zenodo 10.5281/zenodo.18646525 (md5 checked against the record).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import astropy.constants as const
import astropy.units as u
import lens_consistency as lc
import macs1931_gatek as mg
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.ndimage import map_coordinates

from jwst_anomaly import lensmodel

OUT = Path("results/transients/sneos")
SYSTEM = "101"
Z_SN = 5.133  # catalogue z of 101.1 / 101.2 (PLAN.md rule 2)
P_RADIUS_ARCSEC = 40.0  # PLAN.md rule 1
P_MIN_ALPHA = 1.0  # map units (pixels)
P_MAX_MEDIAN_REL = 0.05  # ASSUMPTION (PLAN.md rule 1)
MATCH_ARCSEC = 0.5  # ASSUMPTION (PLAN.md rule 3)
M_RATIO = (0.67, 1.5)  # ASSUMPTION (PLAN.md rule 5)
TARGET_DT_DAYS = 400.0  # ASSUMPTION (PLAN.md "What changes the plan")
TARGET_ABS_MU = 2.0  # ASSUMPTION
GRID_STEP = 0.25
FILES = ("Alpha_x.fits.zip", "Alpha_y.fits.zip", "Lens_pot.fits.zip", "Mu_zs5.13.fits.zip")
RAD_PER_ARCSEC = np.pi / (180.0 * 3600.0)


def delay_days(phi_arcsec2: np.ndarray, z_l: float, z_s: float, cosmo: FlatLambdaCDM) -> np.ndarray:
    """Fermat potential (arcsec², reduced to the source plane) -> arrival time in days."""
    d_l = cosmo.angular_diameter_distance(z_l)
    d_s = cosmo.angular_diameter_distance(z_s)
    d_ls = cosmo.angular_diameter_distance(z_l, z_s)
    d_dt = (1.0 + z_l) * d_l * d_s / d_ls
    t = (d_dt / const.c).to(u.day).value * RAD_PER_ARCSEC**2
    return t * np.asarray(phi_arcsec2, float)


def fermat(x, y, bx: float, by: float, psi: np.ndarray) -> np.ndarray:
    """|θ − β|²/2 − ψ (all arcsec / arcsec², source-plane reduced)."""
    return 0.5 * ((np.asarray(x) - bx) ** 2 + (np.asarray(y) - by) ** 2) - np.asarray(psi)


def potential_scale(ax, ay, psi, mask) -> dict:
    """PLAN.md rule 1: least-squares s in α ≈ s·∇ψ (per pixel) over ``mask``."""
    dpsi_dj, dpsi_di = np.gradient(np.asarray(psi, np.float64))
    gx, gy = dpsi_di[mask], dpsi_dj[mask]
    a_x, a_y = ax[mask], ay[mask]
    s = float((gx @ a_x + gy @ a_y) / (gx @ gx + gy @ gy))
    rel = np.hypot(s * gx - a_x, s * gy - a_y) / np.hypot(a_x, a_y)
    med = float(np.median(rel))
    return {
        "s": s,
        "n_pixels": int(mask.sum()),
        "median_rel_residual": med,
        "p95_rel_residual": float(np.percentile(rel, 95)),
        "passes": bool(med <= P_MAX_MEDIAN_REL),
    }


def match(pred_x, pred_y, cat_x, cat_y, radius: float = MATCH_ARCSEC) -> dict[int, int]:
    """One-to-one nearest matches (catalogued index -> predicted index) within ``radius``."""
    d = np.hypot(np.subtract.outer(cat_x, pred_x), np.subtract.outer(cat_y, pred_y))
    out: dict[int, int] = {}
    for flat in np.argsort(d, axis=None):
        i, j = np.unravel_index(flat, d.shape)
        if d[i, j] > radius:
            break
        if i not in out and j not in out.values():
            out[int(i)] = int(j)
    return out


def sample(model: lensmodel.MapLensModel, arr: np.ndarray, x, y) -> np.ndarray:
    """Bilinear value of a map on the model grid at frame positions (arcsec)."""
    ra, dec = model.to_sky(np.atleast_1d(x), np.atleast_1d(y))
    i, j = model.wcs.world_to_pixel_values(ra, dec)
    return map_coordinates(arr, [np.asarray(j), np.asarray(i)], order=1, cval=np.nan)


def run(maps: Path, out: Path = OUT) -> dict:
    norm = mg.k_dlds_check()
    if not norm["passes"]:
        raise SystemExit(f"error: k_dlds check fails: {norm}")
    checks = mg.record_checksums()
    shas = {f: mg.md5_ok(maps / f, checks) for f in FILES}
    model_sha = hashlib.sha256((shas[FILES[0]] + shas[FILES[1]]).encode()).hexdigest()
    (ax, hx), (ay, _) = (mg.read_zip_map(maps / f) for f in mg.MAPS)
    psi_map, hp = mg.read_zip_map(maps / FILES[2])
    if psi_map.shape != ax.shape or hp.get("CRVAL1") != hx.get("CRVAL1"):
        raise ValueError("Lens_pot is not on the deflection grid")
    model = mg.build_model(ax, ay, hx, mg.READINGS["P_pixels"], model_sha, "m1931-venus-v2.1")

    images = mg.load_images()
    sn = images[np.asarray(images["system"], str) == SYSTEM]
    cx, cy = model.to_frame(np.asarray(sn["ra"], float), np.asarray(sn["dec"], float))

    # Gate P over pixels within P_RADIUS_ARCSEC of 101.1
    jj, ii = np.indices(ax.shape)
    i0, j0 = model.wcs.world_to_pixel_values(float(sn["ra"][0]), float(sn["dec"][0]))
    r = np.hypot(ii - i0, jj - j0) * model.pixel_arcsec
    mask = (r <= P_RADIUS_ARCSEC) & (np.hypot(ax, ay) >= P_MIN_ALPHA)
    gate_p = potential_scale(ax, ay, psi_map, mask)
    del ax, ay, jj, ii, r, mask

    bt = lensmodel.backtrace_images(model, sn, {})
    bx, by = float(np.mean(bt["beta_x"])), float(np.mean(bt["beta_y"]))
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.compute(model, half, GRID_STEP)
    pred = lensmodel.find_images(model, grid, bx, by, Z_SN)
    px, py = np.asarray(pred["x"], float), np.asarray(pred["y"], float)
    m = match(px, py, np.asarray(cx), np.asarray(cy))
    gate_s = {"matched": {str(sn["image_id"][i]): j for i, j in m.items()}, "passes": len(m) == 2}

    # ψ in arcsec² at D_LS/D_S(5.079), then reduced to z 5.133 (PLAN.md rule 1)
    psi_red = gate_p["s"] * model.pixel_arcsec**2 * mg.dls_ds(Z_SN) / mg.dls_ds(mg.Z_NORM)
    psi_at = sample(model, psi_map, px, py) * psi_red
    cosmo = FlatLambdaCDM(H0=mg.COSMOLOGY[0], Om0=mg.COSMOLOGY[1])
    t = delay_days(fermat(px, py, bx, by, psi_at), mg.Z_LENS, Z_SN, cosmo)
    ref = m.get(0)
    dt = t - t[ref] if ref is not None else np.full(len(t), np.nan)

    mu_map, _ = mg.read_zip_map(maps / FILES[3])
    mu_ours = np.asarray(pred["magnification"], float)
    mu_pub = sample(model, mu_map, px, py)
    ratios = {str(sn["image_id"][i]): float(abs(mu_ours[j]) / abs(mu_pub[j])) for i, j in m.items()}
    gate_m = {
        "ratio_ours_over_map": ratios,
        "passes": len(ratios) == 2 and all(M_RATIO[0] <= v <= M_RATIO[1] for v in ratios.values()),
    }
    cat_of = {j: str(sn["image_id"][i]) for i, j in m.items()}
    rows = []
    for k in np.argsort(dt):
        rows.append(
            {
                "catalogued": cat_of.get(int(k), ""),
                "ra": float(pred["ra"][k]),
                "dec": float(pred["dec"][k]),
                "x_arcsec": float(px[k]),
                "y_arcsec": float(py[k]),
                "mu_ours": float(mu_ours[k]),
                "mu_published_map": float(mu_pub[k]),
                "dt_days_vs_101.1": float(dt[k]),
                "target": bool(
                    int(k) not in cat_of
                    and abs(dt[k]) <= TARGET_DT_DAYS
                    and abs(mu_ours[k]) >= TARGET_ABS_MU
                ),
            }
        )
    summary = {
        "plan": str(out / "PLAN.md"),
        "inputs_sha256": shas,
        "normalisation_check": norm,
        "gate_P": gate_p,
        "source_arcsec": [bx, by],
        "grid_half_width_arcsec": half,
        "gate_S": gate_s,
        "gate_M": gate_m,
        "n_predicted": len(rows),
        "images": rows if gate_s["passes"] and gate_p["passes"] else [],
        "images_unvalidated": [] if gate_s["passes"] and gate_p["passes"] else rows,
        "provenance": "model_prediction",
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    s = run(a.maps, a.out)
    print(json.dumps({k: s[k] for k in ("gate_P", "gate_S", "gate_M", "n_predicted")}, indent=1))
    for row in s["images"] or s["images_unvalidated"]:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
