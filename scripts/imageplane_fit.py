"""Image-plane source fit for map lens models (frozen plan:
``results/lens_consistency/imageplane_fit/PLAN.md``).

``lensmodel.imageplane_residuals`` places each system's source at its mean back-traced
position. Near critical curves that point can form fewer images than are catalogued, so a
model that fits the images (as Lenstool's image-plane optimisation does) can still score badly.
Here the source position of each system is fitted: Nelder-Mead on the sum of squared
image-plane offsets, with catalogued and predicted images paired one-to-one; a catalogued image
left without a partner costs its nearest-prediction offset plus a fixed penalty.

    python scripts/imageplane_fit.py [--steps v1,v2,r]

V1: synthetic sources on the Abell S1063 CATS map; V2: MACS1149 CATS ``arcs.txt``; R: the
D-035 rules of ``s1063_d035`` with fitted sources (run only when V1 and V2 pass). The CATS
files live on ``archive.stsci.edu``: see ``s1063_d035`` for the CA bundle cloud sessions need.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import linear_sum_assignment, minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
import s1063_d035 as d035  # noqa: E402

from jwst_anomaly import lensmodel, schema  # noqa: E402
from jwst_anomaly.paths import repo_root  # noqa: E402

OUT = repo_root() / "results/lens_consistency/imageplane_fit"


@dataclass(frozen=True)
class Params:
    grid_step: float = 0.25  # solver grid step, arcsec (validate's default)
    penalty_arcsec: float = 10.0  # ASSUMPTION: cost of a catalogued image without a partner
    simplex_arcsec: float = 0.5  # initial Nelder-Mead simplex size in the source plane
    xatol: float = 1e-3  # arcsec
    fatol: float = 1e-6  # arcsec^2
    maxfev: int = 200
    start_merge_arcsec: float = 0.05  # traced starts closer than this are one start
    v1_n: int = 10  # synthetic sources
    v1_offset_arcsec: float = 1.0  # start offset from the true source
    v1_max_rms: float = 0.05  # ASSUMPTION (V1): recovered rms per source
    v1_min_pass: int = 9  # of v1_n
    v1_z: float = 2.0
    v1_seed: int = 20261021
    v2_quoted_rms: float = 0.63  # CATS v4.1 MACS1149 quoted rms (D-035 table)
    factor: float = lc.IMAGE_LIST_RMS_FACTOR  # ASSUMPTION (D-035): 1.5 x quoted


def match_cost(ox, oy, px, py, penalty: float) -> tuple[float, np.ndarray, np.ndarray]:
    """Sum of squared offsets of catalogued images ``(ox, oy)`` paired one-to-one with predicted
    ``(px, py)``; unpaired catalogued images take their nearest prediction plus ``penalty``
    (arcsec, added in quadrature). Returns ``(cost, dtheta, paired)`` per catalogued image."""
    ox, oy, px, py = (np.asarray(v, float) for v in (ox, oy, px, py))
    n = len(ox)
    if not len(px):
        return n * penalty**2 * 4.0, np.full(n, np.nan), np.zeros(n, bool)
    d = np.hypot(ox[:, None] - px[None, :], oy[:, None] - py[None, :])
    rows, cols = linear_sum_assignment(d**2)
    dtheta = d.min(axis=1)
    paired = np.zeros(n, bool)
    paired[rows] = True
    dtheta[rows] = d[rows, cols]
    cost = float(np.sum(dtheta**2) + (n - paired.sum()) * penalty**2)
    return cost, dtheta, paired


def starts(bx, by, merge: float) -> list[tuple[float, float]]:
    """The mean traced source, then each image's traced source not within ``merge`` of one
    already listed."""
    out = [(float(np.mean(bx)), float(np.mean(by)))]
    for x, y in zip(bx, by, strict=True):
        if all(np.hypot(x - u, y - v) > merge for u, v in out):
            out.append((float(x), float(y)))
    return out


def fit_source(model, grid, ox, oy, z: float, start_list, p: Params) -> dict:
    """Best source position over Nelder-Mead runs from ``start_list``."""
    cache: dict[tuple[float, float], float] = {}

    def f(b):
        key = (round(float(b[0]), 6), round(float(b[1]), 6))
        if key not in cache:
            pred = lensmodel.find_images(model, grid, key[0], key[1], z)
            cache[key] = match_cost(ox, oy, pred["x"], pred["y"], p.penalty_arcsec)[0]
        return cache[key]

    best = None
    for s in start_list:
        s = np.asarray(s, float)
        simplex = np.array([s, s + [p.simplex_arcsec, 0.0], s + [0.0, p.simplex_arcsec]])
        r = minimize(
            f,
            s,
            method="Nelder-Mead",
            options={
                "initial_simplex": simplex,
                "xatol": p.xatol,
                "fatol": p.fatol,
                "maxfev": p.maxfev,
            },
        )
        if best is None or r.fun < best.fun:
            best = r
    pred = lensmodel.find_images(model, grid, float(best.x[0]), float(best.x[1]), z)
    cost, dtheta, paired = match_cost(ox, oy, pred["x"], pred["y"], p.penalty_arcsec)
    return {
        "beta": [float(best.x[0]), float(best.x[1])],
        "cost": cost,
        "dtheta": dtheta,
        "paired": paired,
        "n_predicted": len(pred),
        "n_eval": len(cache),
        "n_starts": len(start_list),
    }


_WORKER: dict = {}


def _init_worker(name: str, half: float, step: float) -> None:
    """Pool initializer: each worker loads the model and its cached grid by name."""
    model, _, _ = lc.load_model(name)
    grid = lensmodel.DeflectionGrid.cached(model, lc.grid_cache_path(model, half, step), half, step)
    _WORKER.update(model=model, grid=grid)


def _fit_system(job):
    """One system: ``(sel, mean-source dtheta, fit result)``."""
    sel, ox, oy, z, st, p = job
    model, grid = _WORKER["model"], _WORKER["grid"]
    pm = lensmodel.find_images(model, grid, st[0][0], st[0][1], z)
    d_mean = match_cost(ox, oy, pm["x"], pm["y"], p.penalty_arcsec)[1]
    return sel, d_mean, fit_source(model, grid, ox, oy, z, st, p)


def fit_list(name: str, bt: Table, half: float, p: Params, workers: int) -> Table:
    """Per-image fitted residuals for every system with >= 2 traced images (provenance
    ``model_prediction``); ``dtheta_arcsec`` is NaN for an image without a partner. Systems
    are fitted in a process pool of ``workers``."""
    _init_worker(name, half, p.grid_step)
    model = _WORKER["model"]
    systems = np.asarray(bt["system"], str)
    n = len(bt)
    cols = {
        "dtheta_arcsec": np.full(n, np.nan),
        "dtheta_mean_source_arcsec": np.full(n, np.nan),
        "n_predicted": np.zeros(n, int),
        "beta_fit_x": np.full(n, np.nan),
        "beta_fit_y": np.full(n, np.nan),
    }
    ox_all, oy_all = model.to_frame(np.asarray(bt["ra"], float), np.asarray(bt["dec"], float))
    jobs = []
    for s in dict.fromkeys(systems):
        sel = np.where((systems == s) & np.isfinite(bt["beta_x"]))[0]
        if len(sel) < 2:
            continue
        bx, by = np.asarray(bt["beta_x"][sel], float), np.asarray(bt["beta_y"][sel], float)
        st = starts(bx, by, p.start_merge_arcsec)
        jobs.append((sel, ox_all[sel], oy_all[sel], float(bt["z_used"][sel[0]]), st, p))
    t0 = time.time()
    with ProcessPoolExecutor(
        workers, initializer=_init_worker, initargs=(name, half, p.grid_step)
    ) as ex:
        for k, (sel, d_mean, r) in enumerate(ex.map(_fit_system, jobs), 1):
            cols["dtheta_mean_source_arcsec"][sel] = d_mean
            cols["dtheta_arcsec"][sel] = np.where(r["paired"], r["dtheta"], np.nan)
            cols["n_predicted"][sel] = r["n_predicted"]
            cols["beta_fit_x"][sel], cols["beta_fit_y"][sel] = r["beta"]
            print(f"  {name}: {k}/{len(jobs)} systems, {time.time() - t0:.0f}s", flush=True)
    out = Table({"image_id": bt["image_id"], "system": systems, "z_used": bt["z_used"], **cols})
    out.meta.update(model._meta())
    out.meta["provenance"] = schema.Provenance.MODEL_PREDICTION.value
    out.meta["source"] = f"{bt.meta.get('source', 'images')}, image-plane source fit"
    return out


def rms(d) -> float:
    d = np.asarray(d, float)
    ok = np.isfinite(d)
    return float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else float("nan")


def v1_gate(rms_values, max_rms: float, min_pass: int) -> dict:
    """V1: passes when at least ``min_pass`` synthetic sources are recovered to rms <=
    ``max_rms`` with every image paired (NaN counts as a failure)."""
    r = np.asarray(rms_values, float)
    n_ok = int(np.sum(np.isfinite(r) & (r <= max_rms)))
    return {"n": int(len(r)), "n_pass": n_ok, "passes": n_ok >= min_pass}


def grid_for(model, images: Table, step: float):
    """``(grid, half_width)`` covering ``images``; the grid is cached on disk for the workers."""
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.cached(model, lc.grid_cache_path(model, half, step), half, step)
    return grid, half


def step_v1(p: Params) -> dict:
    """Synthetic sources behind the S1063 map: images from a known source, fit started off it."""
    model, _, _ = lc.load_model(d035.MODEL)
    chan = d035.chan_images()
    grid, _ = grid_for(model, chan, p.grid_step)
    rng = np.random.default_rng(p.v1_seed)
    bt = lensmodel.backtrace_images(model, chan, {})
    bx0, by0 = np.asarray(bt["beta_x"], float), np.asarray(bt["beta_y"], float)
    rows = []
    tries = 0
    while len(rows) < p.v1_n and tries < 500:
        tries += 1
        k = rng.integers(len(bt))  # near a real source, so inside the multiple-image region
        b = (bx0[k] + rng.normal(0, 1.0), by0[k] + rng.normal(0, 1.0))
        true = lensmodel.find_images(model, grid, b[0], b[1], p.v1_z)
        if len(true) < 3:
            continue
        ang = rng.uniform(0, 2 * np.pi)
        s0 = (b[0] + p.v1_offset_arcsec * np.cos(ang), b[1] + p.v1_offset_arcsec * np.sin(ang))
        r = fit_source(model, grid, true["x"], true["y"], p.v1_z, [s0], p)
        rr = rms(r["dtheta"]) if r["paired"].all() else float("nan")
        rows.append({"beta_true": [float(b[0]), float(b[1])], "n_images": len(true), "rms": rr})
    return {"sources": rows, "gate": v1_gate([r["rms"] for r in rows], p.v1_max_rms, p.v1_min_pass)}


def step_v2(p: Params, workers: int) -> tuple[dict, Table]:
    """MACS1149 CATS ``arcs.txt`` with its ``params.txt`` redshifts, fitted sources."""
    name = "macs1149-cats"
    model, files, _ = lc.load_model(name)
    images, zml = lc.image_list(name, files, None)
    bt = lensmodel.backtrace_images(model, images, zml)
    bt = bt[np.isfinite(bt["beta_x"])]
    _, half = grid_for(model, bt, p.grid_step)
    fit = fit_list(name, bt, half, p, workers)
    d = np.asarray(fit["dtheta_arcsec"], float)
    limit = p.factor * p.v2_quoted_rms
    gate = {
        "n_images": len(d),
        "n_solved": int(np.isfinite(d).sum()),
        "rms_fit_arcsec": rms(d),
        "rms_mean_source_arcsec": rms(fit["dtheta_mean_source_arcsec"]),
        "limit_arcsec": limit,
        "passes": bool(np.isfinite(d).all() and rms(d) <= limit),
    }
    return gate, fit


def step_r(p: Params, workers: int) -> tuple[dict, Table]:
    """S1063 Chan+2020 list with fitted sources; D-035 rules of ``s1063_d035``."""
    dp = d035.Params()
    model, _, _ = lc.load_model(d035.MODEL)
    chan = d035.chan_images()
    bt = lensmodel.backtrace_images(model, chan, {})
    _, half = grid_for(model, chan, p.grid_step)
    fit = fit_list(d035.MODEL, bt, half, p, workers)
    fit["z_spec"] = chan["z_spec"]
    zs = np.asarray(fit["z_spec"], bool)
    gates = {
        "S": d035.whole_list_gate(fit["dtheta_arcsec"][zs], dp.quoted_rms_arcsec, dp.factor),
        "A": d035.whole_list_gate(fit["dtheta_arcsec"], dp.quoted_rms_arcsec, dp.factor),
    }
    opened = lc.systems_passing_rms(fit[zs], dp.system_max_rms)
    systems = np.asarray(fit["system"], str)
    per_system = {}
    for s in dict.fromkeys(systems):
        sel = systems == s
        per_system[s] = {
            "z_spec": bool(zs[sel][0]),
            "n": int(sel.sum()),
            "n_predicted": int(fit["n_predicted"][sel][0]),
            "rms_fit_arcsec": rms(fit["dtheta_arcsec"][sel]),
            "rms_mean_source_arcsec": rms(fit["dtheta_mean_source_arcsec"][sel]),
            "all_paired": bool(np.isfinite(fit["dtheta_arcsec"][sel]).all()),
        }
    summary = {
        "gates": gates,
        "r2_open_systems_S": opened,
        "outcome": d035.outcome(gates["S"], opened),
        "per_system": per_system,
    }
    return summary, fit


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--steps", default="v1,v2,r")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args(argv)
    steps = args.steps.split(",")
    p = Params()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "summary.json"
    summary = json.loads(path.read_text()) if path.exists() else {}
    summary["params"] = asdict(p)
    t0 = time.time()
    if "v1" in steps:
        summary["v1"] = step_v1(p)
        print("V1", json.dumps(summary["v1"]["gate"]), f"{time.time() - t0:.0f}s", flush=True)
    if "v2" in steps:
        summary["v2"], fit = step_v2(p, args.workers)
        fit.write(OUT / "macs1149_fit.ecsv", format="ascii.ecsv", overwrite=True)
        print("V2", json.dumps(summary["v2"]), f"{time.time() - t0:.0f}s", flush=True)
    if "r" in steps:
        v_ok = summary.get("v1", {}).get("gate", {}).get("passes") and summary.get("v2", {}).get(
            "passes"
        )
        if not v_ok:
            summary["r"] = {"skipped": "V1 or V2 did not pass (frozen plan)"}
        else:
            summary["r"], fit = step_r(p, args.workers)
            fit.write(OUT / "abells1063_fit.ecsv", format="ascii.ecsv", overwrite=True)
        print("R", json.dumps({k: summary["r"].get(k) for k in ("gates", "outcome", "skipped")}))
    path.write_text(json.dumps(summary, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
