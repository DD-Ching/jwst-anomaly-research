"""MACS0416 D115b vetting: member-8515 σ scan with every family's source fitted.

Frozen plan: ``results/lens_consistency/macs0416_canucsref/sigma_scan/PLAN.md``. 8515 sits on its
z_spec plane 0.3907 with σ = f × σ_fit; for each f every family's source is fitted to the image
plane as in ``m0416_d115_refit.py`` and the summed χ² is compared with f = 1. Writes
``result.json``.

    python scripts/m0416_d115_sigma_scan.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")  # before numpy loads BLAS; spawned workers rerun it

import json  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from functools import lru_cache  # noqa: E402
from pathlib import Path  # noqa: E402

import lens_consistency as lc
import m0416_d115_refit as refit
import numpy as np
from scipy.optimize import minimize

from jwst_anomaly import lensmodel  # noqa: E402

MODEL = refit.MODEL
OUT = Path("results/lens_consistency/macs0416_canucsref/sigma_scan")
SIGMA_FIT = 120.6  # km/s, 8515 in the pinned best.par
FACTORS = (1.00, 0.85, 0.80, 0.75, 0.70)
G1_REF = {"D115a": 2.765, "D115b": 0.538}  # ../d115_refit, σ 84.4
G1_TOL = 0.05
G2_TOL = 1e-6
G3_REF = {"D115a": 0.919, "D115b": 0.386}  # ../d115_refit, σ as fitted
SEARCH_STEP = 0.25  # plan amendment: search grid; final χ² on GRID_STEP
FITS = OUT / "fits.jsonl"
PRIORITY = ("D115a", "D115b", "D34", "D70", "D69", "D15", "D72b", "D74", "D7", "D72a")


def plane_spec(f: float) -> list[str]:
    return [f"8515=0.3907:{f * SIGMA_FIT:.1f}"]


@lru_cache(maxsize=10)
def setup(f: float, step: float = refit.GRID_STEP):
    files = lc.model_files(MODEL)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model, _ = lc.apply_planes(MODEL, lensmodel.LensModel.from_par(par), plane_spec(f))
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    sigma = lc.model_sigpos(MODEL, files, images)
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.cached(model, lc.grid_cache_path(model, half, step), half, step)
    bt = lensmodel.backtrace_images(model, images, par["z_m_limit"])
    return model, grid, sigma, bt


def families() -> list[str]:
    """Family ids in task order; back-traces analytically, so no deflection grid is built."""
    files = lc.model_files(MODEL)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model, _ = lc.apply_planes(MODEL, lensmodel.LensModel.from_par(par), plane_spec(1.0))
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    bt = lensmodel.backtrace_images(model, images, par["z_m_limit"])
    fams = sorted(set(np.asarray(bt["system"]).astype(str)))
    return [f for f in PRIORITY if f in fams] + [f for f in fams if f not in PRIORITY]


def build(f: float) -> float:
    """Build (and cache) both deflection grids of one model; returns f."""
    setup(f)
    setup(f, SEARCH_STEP)
    return f


def fit_family(task: tuple[float, str]) -> dict:
    f, fam = task
    model, grid, sigma, bt = setup(f)
    coarse = setup(f, SEARCH_STEP)[1]
    sel = np.where(np.asarray(bt["system"]).astype(str) == fam)[0]
    obs = bt[sel]
    # validate's residuals are per family (mean back-traced source), so solving only this family
    # gives the same numbers as solving all of them in setup() for every worker, f and grid
    d = np.asarray(lensmodel.imageplane_residuals(model, grid, obs)["dtheta_arcsec"], float)
    ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
    sig = sigma[sel]
    z = float(obs["z_used"][0])
    m = np.isfinite(d)
    term_validate = float(np.sum((d[m] / sig[m]) ** 2))

    def cost(b):
        pred = lensmodel.find_images(model, coarse, float(b[0]), float(b[1]), z)
        c, _, shared = refit.family_chi2(ox, oy, sig, pred)
        return c + (1e4 if shared else 0.0)

    def fine(b):
        pred = lensmodel.find_images(model, grid, float(b[0]), float(b[1]), z)
        c, _, shared = refit.family_chi2(ox, oy, sig, pred)
        return c + (1e4 if shared else 0.0), c, shared, pred

    bxs, bys = np.asarray(obs["beta_x"], float), np.asarray(obs["beta_y"], float)
    starts = [(bxs.mean(), bys.mean())] + list(zip(bxs, bys, strict=True))
    ends = [np.array(s) for s in starts[:1]]  # the mean itself, so G2 holds by construction
    for s in starts:
        r = minimize(
            cost, np.array(s), method="Nelder-Mead",
            options={"xatol": 1e-4, "fatol": 1e-3, "maxfev": 300},
        )  # fmt: skip
        ends.append(r.x)
    _, c, shared, pred, best_x = min((fine(x) + (x,) for x in ends), key=lambda t: t[0])
    out = {
        "family": fam,
        "f": f,
        "n_images": len(sel),
        "chi2_validate_source": term_validate,
        "chi2_fit": c,
        "shared": shared,
        "n_predicted": len(pred),
        "beta_fit": [float(best_x[0]), float(best_x[1])],
    }
    if fam == "D115b":
        fx, fy = model.to_frame(np.array([refit.FLAG[0]]), np.array([refit.FLAG[1]]))
        df = np.hypot(np.asarray(pred["x"]) - fx[0], np.asarray(pred["y"]) - fy[0])
        mu = np.abs(np.asarray(pred["magnification"], float))
        out["flag_image"] = bool(((df <= refit.FLAG_RADIUS) & (mu >= refit.MU_MIN)).any())
        out["nearest_to_flag_arcsec"] = float(df.min()) if len(df) else None
    return out


def summarize(fits: list[dict]) -> dict:
    """Per-f totals against f = 1 (families matched one-to-one in both), plus the losses."""
    by = {(r["f"], r["family"]): r for r in fits}
    fams = sorted({r["family"] for r in fits})
    ref = {fam: by[(1.0, fam)] for fam in fams}
    out = {}
    for f in FACTORS:
        both = [fam for fam in fams if not ref[fam]["shared"] and not by[(f, fam)]["shared"]]
        diff = {fam: by[(f, fam)]["chi2_fit"] - ref[fam]["chi2_fit"] for fam in both}
        top = sorted(diff.items(), key=lambda kv: -abs(kv[1]))[:5]
        out[f"{f:.2f}"] = {
            "sigma_kms": round(f * SIGMA_FIT, 1),
            "chi2_all_matched": float(sum(by[(f, fam)]["chi2_fit"] for fam in both)),
            "dchi2": float(sum(diff.values())),
            "n_families_compared": len(both),
            "lost_match": [
                fam for fam in fams if not ref[fam]["shared"] and by[(f, fam)]["shared"]
            ],
            "flag_image": by[(f, "D115b")]["flag_image"],
            "d115_matched": not (by[(f, "D115a")]["shared"] or by[(f, "D115b")]["shared"]),
            "top5": {k: float(v) for k, v in top},
        }
    return out


def gates(fits: list[dict]) -> dict:
    by = {(r["f"], r["family"]): r for r in fits}
    g1 = all(abs(by[(0.70, fam)]["chi2_fit"] - v) <= G1_TOL for fam, v in G1_REF.items())
    g1 = g1 and not by[(0.70, "D115b")]["flag_image"]
    g3 = all(abs(by[(1.0, fam)]["chi2_fit"] - v) <= G1_TOL for fam, v in G3_REF.items())
    g2_fail = [
        r["family"]
        for r in fits
        if r["f"] == 1.0 and r["chi2_fit"] > r["chi2_validate_source"] + G2_TOL
    ]
    return {"G1": bool(g1), "G2": not g2_fail, "G2_failures": g2_fail, "G3": bool(g3)}


def verdict(scan: dict, gate: dict | None = None) -> tuple[str, str | None]:
    """Frozen rule: f* is the largest f < 1 without the flag image and with D115 matched.

    A failed gate stops the interpretation (plan): the verdict is then ``gates_failed``.
    """
    if gate is not None and not (gate["G1"] and gate["G2"] and gate["G3"]):
        return "gates_failed", None
    cands = [k for k in sorted(scan, reverse=True) if float(k) < 1.0]
    fstar = next((k for k in cands if not scan[k]["flag_image"] and scan[k]["d115_matched"]), None)
    if fstar is None:
        return "inconclusive", None
    s = scan[fstar]
    ok = s["dchi2"] <= refit.DCHI2_MAX and not s["lost_match"]
    return ("ordinary" if ok else "inconclusive"), fstar


def load_fits() -> list[dict]:
    """Finished rows; a last line cut off by a killed run is dropped (and the file trimmed)."""
    if not FITS.exists():
        return []
    lines = FITS.read_text().splitlines()
    try:
        rows = [json.loads(line) for line in lines]
    except json.JSONDecodeError:
        rows = [json.loads(line) for line in lines[:-1]]  # only the last line may be partial
        FITS.write_text("".join(line + "\n" for line in lines[:-1]))
    return rows


def main() -> None:
    fams = families()  # the grids are built in parallel by build(), not here
    OUT.mkdir(parents=True, exist_ok=True)
    fits = load_fits()
    done = {(r["f"], r["family"]) for r in fits}
    tasks = [(f, fam) for fam in fams for f in FACTORS if (f, fam) not in done]
    with ProcessPoolExecutor(4) as ex:
        list(ex.map(build, FACTORS))
        with FITS.open("a") as fh:
            for fut in as_completed([ex.submit(fit_family, t) for t in tasks]):
                r = fut.result()
                fits.append(r)
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                print(
                    f"{len(fits)} / {len(fams) * len(FACTORS)} {r['family']} f={r['f']}", flush=True
                )
    scan = summarize(fits)
    gate = gates(fits)
    v, fstar = verdict(scan, gate)
    result = {
        "plan": str(OUT / "PLAN.md"),
        "provenance": "model_prediction",
        "grid_step_arcsec": refit.GRID_STEP,
        "dchi2_max_assumption": refit.DCHI2_MAX,
        "gates": gate,
        "scan": scan,
        "f_star": fstar,
        "verdict": v,
        "fits": str(FITS),
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "fits"}, indent=1))


if __name__ == "__main__":
    main()
