"""WHL0137 (Sunrise): deflection unit and per-system source redshifts in the Scofield+2025 maps.

Offline: it reads the per-image map values that `scripts/sunrise_scofield.py` stored in
results/sunrise/scofield25_positions.json (raw map units; model_prediction) and needs no download.

1. Unit. The published deflection maps are said to be scaled to z = 6.2 (Scofield+2025 Sect. VI),
   but their finite-difference kappa matches the published z = 6.2 kappa map only after a scale
   s ~ 0.0223" per map unit. Tested reading (hypothesis): the maps are in 0.02" mosaic pixels
   and scaled to a reference redshift z_ref, so s = 0.02" * r(6.2) / r(z_ref) with
   r(z) = D_ds/D_s. Each image gives one z_ref; the reading predicts one value for all.
2. Redshift. With alpha(z) = raw * 0.02" * r(z) / r(z_ref), the source redshift of a system is the
   z that minimises the source-plane rms of its images about their mean. The z = 6.2 system 1 is
   the known case: its subsystems must return z ~ 6.2 before system 2's fit is trusted. The
   linear image-plane rms |A^-1 delta-beta| is reported too but is not used for the fit: it
   diverges for images near critical curves (most of system 1) and near-critical scales.

Cosmology: flat LCDM, h = 0.7, Omega_m = 0.3 (Scofield+2025 Sect. I). Cluster z = 0.566.

Usage: python scripts/sunrise_zfit.py [--positions ...] [--out results/sunrise/scofield25_zfit.json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import brentq

Z_LENS = 0.566
Z_MAPS = 6.2  # the redshift the maps are said to be scaled to (kappa and mu maps agree)
PIXEL = 0.02  # arcsec, the mosaic and map pixel scale
COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)
Z_GRID = np.round(np.arange(0.6, 20.0001, 0.02), 4)


def dratio(zs: float, zl: float = Z_LENS, cosmo=COSMO) -> float:
    """D_ds / D_s for a source at zs behind a lens at zl (0 for zs <= zl)."""
    if zs <= zl:
        return 0.0
    return float(cosmo.angular_diameter_distance_z1z2(zl, zs) / cosmo.angular_diameter_distance(zs))


def zref_for_scale(s: float) -> float:
    """z_ref such that PIXEL * r(Z_MAPS) / r(z_ref) = s (s in arcsec per map unit)."""
    target = PIXEL * dratio(Z_MAPS) / s
    return float(brentq(lambda z: dratio(z) - target, Z_LENS + 1e-3, 1e4))


def image_plane_rms(rows: list[dict], f: float, ref: set[str] | None = None) -> dict:
    """Image-plane offsets of each image about a common source, for alpha = f * raw.

    The source is the mean of the reference images' beta (all images if ref is empty or None);
    each offset is |A^-1 (beta_i - beta_mean)| with A = 1 - f * psi_raw (linear, so approximate
    near critical curves). rms is over the reference images only."""
    use = [r for r in rows if ref and r["image"] in ref] or rows
    bx = np.mean([r["x_arcsec"] - f * r["alpha_x"] for r in use])
    by = np.mean([r["y_arcsec"] - f * r["alpha_y"] for r in use])
    out = {}
    for r in rows:
        d = np.array([r["x_arcsec"] - f * r["alpha_x"] - bx, r["y_arcsec"] - f * r["alpha_y"] - by])
        pxx, pxy, pyy = (f * p for p in r["psi"])
        a = np.array([[1 - pxx, -pxy], [-pxy, 1 - pyy]])
        out[r["image"]] = {
            "dbeta_arcsec": float(np.hypot(*d)),
            "dtheta_linear_arcsec": float(np.hypot(*np.linalg.solve(a, d))),
        }
    names = [r["image"] for r in use]
    return {
        "images": out,
        "rms_theta": float(np.sqrt(np.mean([out[n]["dtheta_linear_arcsec"] ** 2 for n in names]))),
        "rms_beta": float(np.sqrt(np.mean([out[n]["dbeta_arcsec"] ** 2 for n in names]))),
    }


def fit_redshift(rows, zref: float, ref=None, grid=Z_GRID) -> dict:
    """Grid scan of the source-plane rms against source redshift; returns the curves and the best z.

    z_best at the top of the grid means the rms still falls toward z -> infinity."""
    r_ref = dratio(zref)
    fs = [PIXEL * dratio(z) / r_ref for z in grid]
    beta = [image_plane_rms(rows, f, ref)["rms_beta"] for f in fs]
    theta = [image_plane_rms(rows, f, ref)["rms_theta"] for f in fs]
    k = int(np.argmin(beta))
    return {
        "z": [float(z) for z in grid],
        "rms_beta": beta,
        "rms_theta": theta,
        "z_best": float(grid[k]),
        "rms_beta_best": beta[k],
        "at_grid_edge": k in (0, len(grid) - 1),
    }


def interval(fit: dict, level: float) -> list[float] | None:
    """Range of grid redshifts whose source-plane rms is <= level (None if empty)."""
    z = [z for z, v in zip(fit["z"], fit["rms_beta"], strict=True) if v <= level]
    return [min(z), max(z)] if z else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--positions", default="results/sunrise/scofield25_positions.json")
    ap.add_argument("--out", default="results/sunrise/scofield25_zfit.json")
    args = ap.parse_args(argv)
    pos = json.loads(Path(args.positions).read_text())
    rows = pos["rows"]

    # 1. unit: per-image scale from the published kappa map, and the z_ref it implies
    per = {r["image"]: r["kappa_map"] / r["kappa_fd"] for r in rows}
    zref_img = {k: zref_for_scale(v) for k, v in per.items()}
    s_med = float(np.median(list(per.values())))
    zref = zref_for_scale(s_med)
    s16, s84 = (float(v) for v in np.percentile(list(per.values()), [16, 84]))
    s3 = PIXEL * dratio(Z_MAPS) / dratio(3.0)

    # 2. redshifts: system 1 (known case, constraint images) and system 2 (held out, z unknown)
    sys1 = [r for r in rows if r["image"].startswith("1.")]
    cons = {r["image"] for r in sys1 if r["constraint"]}
    sys2 = [r for r in rows if r["image"].startswith("2")]
    ref_rms = image_plane_rms(sys1, s3, cons)["rms_beta"]  # constraint level (system 1 at z = 6.2)
    fits = {}
    c1 = [r for r in sys1 if r["constraint"]]
    for name, members, ref in (
        ("1 (constraints)", [r for r in sys1 if r["constraint"]], cons),
        (
            "1.1 (constraints)",
            [r for r in c1 if r["image"].startswith("1.1")],
            None,
        ),
        ("1.2", [r for r in sys1 if r["image"].startswith("1.2")], None),
        (
            "1.7 (constraints)",
            [r for r in c1 if r["image"].startswith("1.7")],
            None,
        ),
        ("2", sys2, None),
    ):
        f = fit_redshift(members, 3.0, ref)
        f["z_interval_beta_le_constraint_rms"] = interval(f, ref_rms)
        f["z_interval_beta_le_2x_constraint_rms"] = interval(f, 2 * ref_rms)
        fits[name] = f
    at3 = image_plane_rms(sys1, s3, cons)
    result = {
        "provenance": "derived from model_prediction map values (Scofield+2025 M5-G1) "
        "at observed positions",
        "input": args.positions,
        "cosmology": "flat LCDM h=0.7 Om=0.3 (Scofield+2025)",
        "unit": {
            "hypothesis": "maps in 0.02 arcsec pixels, scaled to z_ref (not the stated 6.2)",
            "scale_per_image": per,
            "zref_per_image": zref_img,
            "scale_median": s_med,
            "scale_p16_p84": [s16, s84],
            "zref_median": zref,
            "zref_p16_p84": [zref_for_scale(s84), zref_for_scale(s16)],
            "scale_if_zref_3": s3,
        },
        "system1_at_zref3": {
            "constraint_rms_theta": at3["rms_theta"],
            "constraint_rms_beta": at3["rms_beta"],
            "images": at3["images"],
        },
        "zfit": {
            k: {kk: vv for kk, vv in v.items() if kk not in ("z", "rms_beta", "rms_theta")}
            for k, v in fits.items()
        },
        "zfit_curves": {
            k: {kk: v[kk][::10] for kk in ("z", "rms_beta", "rms_theta")} for k, v in fits.items()
        },
        "scale_at_z_infinity": PIXEL * dratio(1e4) / dratio(3.0),
        "reference_rms_level_arcsec": ref_rms,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n")
    print(
        f"scale median {s_med:.5f} (16-84 {s16:.5f}-{s84:.5f}) -> z_ref {zref:.2f}; "
        f"z_ref=3 gives {s3:.5f}"
    )
    print(
        f"system 1 constraints at z_ref=3: rms theta {at3['rms_theta']:.3f} "
        f"beta {at3['rms_beta']:.3f}"
    )
    for k, v in at3["images"].items():
        print(f"  {k:5s} dbeta {v['dbeta_arcsec']:.3f} dtheta {v['dtheta_linear_arcsec']:.3f}")
    for k, v in fits.items():
        print(
            f"z fit {k}: best {v['z_best']:.2f}{' (grid edge)' if v['at_grid_edge'] else ''} "
            f"beta rms {v['rms_beta_best']:.3f}; "
            f"<= {ref_rms:.3f}: {v['z_interval_beta_le_constraint_rms']}; "
            f"<= 2x: {v['z_interval_beta_le_2x_constraint_rms']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
