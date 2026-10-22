"""MACS0416 D115b vetting follow-up: the χ² cost of a lower member-8515 σ with a fitted source.

Frozen plan: ``results/lens_consistency/macs0416_canucsref/d115_refit/PLAN.md``. For each model
(fiducial; 8515 on z 0.3907 with σ as fitted and with σ 84.4 km/s) the ``validate`` image-plane
χ² is recomputed, and the D115a / D115b terms are replaced by those of a source position fitted
to the image plane (Nelder–Mead, one-to-one image assignment). Writes ``result.json``.

    python scripts/m0416_d115_refit.py
"""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lens_consistency as lc
import numpy as np
from scipy.optimize import linear_sum_assignment, minimize

from jwst_anomaly import lensmodel

MODEL = "macs0416-canucs"
OUT = Path("results/lens_consistency/macs0416_canucsref/d115_refit")
FAMILIES = ("D115a", "D115b")
FLAG = (64.048551, -24.064825)  # D115b extra-image prediction of the #458 run
FLAG_RADIUS = 0.7  # arcsec; as the vetting record's σ scan
MU_MIN = 5.0
DCHI2_MAX = 6.0  # ASSUMPTION: system-51 precedent (330.81 -> 336.44)
GRID_STEP = 0.1
MODELS = {"fiducial": None, "p3907": ["8515=0.3907"], "p3907s84": ["8515=0.3907:84.4"]}


def family_chi2(ox, oy, sig, pred) -> tuple[float, np.ndarray, bool]:
    """χ² of catalogued images (ox, oy, sig) against predicted images (x, y): one-to-one
    assignment when there are enough predicted images, else nearest (``shared`` True)."""
    if not len(pred):
        return np.inf, np.full(len(ox), np.nan), True
    d = np.hypot(ox[:, None] - np.asarray(pred["x"])[None, :], oy[:, None] - np.asarray(pred["y"]))
    if d.shape[1] >= d.shape[0]:
        r, c = linear_sum_assignment((d / sig[:, None]) ** 2)
        dt = d[r, c]
        shared = False
    else:
        dt = d.min(axis=1)
        shared = True
    return float(np.sum((dt / sig) ** 2)), dt, shared


def verdict(dchi2: float, flag_image: bool, all_matched: bool) -> str:
    """Frozen decision rule of the plan."""
    return "ordinary" if dchi2 <= DCHI2_MAX and not flag_image and all_matched else "inconclusive"


def _setup(name):
    files = lc.model_files(MODEL)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model, planes = lc.apply_planes(MODEL, lensmodel.LensModel.from_par(par), MODELS[name])
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    sigma = lc.model_sigpos(MODEL, files, images)
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, GRID_STEP), half, GRID_STEP
    )
    return par, model, images, sigma, grid


def run_model(name: str) -> dict:
    par, model, images, sigma, grid = _setup(name)
    bt = lensmodel.backtrace_images(model, images, par["z_m_limit"])
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    d = np.asarray(ip["dtheta_arcsec"], float)
    ok = np.isfinite(d)
    chi2_validate = float(np.sum((d[ok] / sigma[ok]) ** 2))
    fx, fy = model.to_frame(np.array([FLAG[0]]), np.array([FLAG[1]]))
    systems = np.asarray(bt["system"]).astype(str)
    out = {"model": name, "planes": MODELS[name], "chi2_validate": chi2_validate, "families": {}}
    chi2_fit = chi2_validate
    for fam in FAMILIES:
        sel = np.where(systems == fam)[0]
        obs = bt[sel]
        ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
        sig = sigma[sel]
        z = float(obs["z_used"][0])
        m = ok[sel]  # the same finite terms that chi2_validate summed
        term_validate = float(np.sum((d[sel][m] / sig[m]) ** 2))

        def cost(b, ox=ox, oy=oy, sig=sig, z=z):
            pred = lensmodel.find_images(model, grid, float(b[0]), float(b[1]), z)
            c, _, shared = family_chi2(ox, oy, sig, pred)
            return c + (1e4 if shared else 0.0)

        bxs = np.asarray(obs["beta_x"], float)
        bys = np.asarray(obs["beta_y"], float)
        starts = [(bxs.mean(), bys.mean())] + list(zip(bxs, bys, strict=True))
        best = None
        for s in starts:
            r = minimize(
                cost, np.array(s), method="Nelder-Mead",
                options={"xatol": 1e-4, "fatol": 1e-3, "initial_simplex": None, "maxfev": 300},
            )  # fmt: skip
            if best is None or r.fun < best.fun:
                best = r
        pred = lensmodel.find_images(model, grid, float(best.x[0]), float(best.x[1]), z)
        c, dt, shared = family_chi2(ox, oy, sig, pred)
        df = np.hypot(np.asarray(pred["x"]) - fx[0], np.asarray(pred["y"]) - fy[0])
        mu = np.asarray(pred["magnification"], float)
        near = (df <= FLAG_RADIUS) & (np.abs(mu) >= MU_MIN)
        out["families"][fam] = {
            "z": z,
            "chi2_validate_source": term_validate,
            "shared_validate_source": bool(np.asarray(ip["shared_match"])[sel].any()),
            "beta_mean": [float(bxs.mean()), float(bys.mean())],
            "beta_fit": [float(best.x[0]), float(best.x[1])],
            "chi2_fit_source": c,
            "dtheta_fit_arcsec": {
                str(i): float(v) for i, v in zip(obs["image_id"], dt, strict=True)
            },
            "shared_fit_source": shared,
            "n_predicted": len(pred),
            "flag_image": bool(near.any()),
            "nearest_to_flag_arcsec": float(df.min()) if len(df) else None,
            "nearest_to_flag_abs_mu": float(abs(mu[np.argmin(df)])) if len(df) else None,
        }
        chi2_fit += c - term_validate
    out["chi2_fit"] = chi2_fit
    return out


def main() -> None:
    lc.model_files(MODEL)  # download once before the workers start
    with ProcessPoolExecutor(len(MODELS)) as ex:
        runs = dict(zip(MODELS, ex.map(run_model, MODELS), strict=True))
    low, ref = runs["p3907s84"], runs["p3907"]
    dchi2 = low["chi2_fit"] - ref["chi2_fit"]
    flag_image = low["families"]["D115b"]["flag_image"]
    matched = not any(f["shared_fit_source"] for f in low["families"].values())
    result = {
        "plan": str(OUT / "PLAN.md"),
        "provenance": "model_prediction",
        "grid_step_arcsec": GRID_STEP,
        "dchi2_max_assumption": DCHI2_MAX,
        "runs": runs,
        "dchi2_fit": dchi2,
        "low_sigma_flag_image": flag_image,
        "low_sigma_all_matched": matched,
        "verdict": verdict(dchi2, flag_image, matched),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "runs"}, indent=1))


if __name__ == "__main__":
    main()
