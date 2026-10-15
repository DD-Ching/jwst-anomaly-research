"""X4 S1e: joint lens + intrinsic search (SciPy DE), exact continuous refinement, then R1.

Unit file: docs/state/units/X4.md (S1e; diagnostics and plan frozen before the search).
`simulated` / `derived`; no data. Needs lalsuite. Reuses the S1d score, parameter box and objective
(`x4_s1d_de`); one procedure serves validation and the R1 case.

U: DE over the 10 intrinsic parameters, then Nelder-Mead restarts.
P: one DE over 12 parameters (intrinsics + ln M_Lz, ln y), the lens snapped to the nearest cached
grid node; then continuous refinement with exact lens factors (exact-lens 2-D Nelder-Mead
alternated with intrinsic restarts). Final scores: exact lens, pad 64.

  python scripts/x4_s1e_joint.py validate   # -> results/x4/s1e_validate.json
  python scripts/x4_s1e_joint.py case       # -> results/x4/s1e_case.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


e = _load("x4_s1d_de")
m, b = e.m, e.b
TOL = 1e-5
LN_M = (0.0, np.log(1e4))
LN_Y = (np.log(0.01), np.log(5.0))
BOUNDS_P = e.BOUNDS + [LN_M, LN_Y]
DE_U = {"popsize": 12, "maxiter": 600}
DE_P = {"popsize": 12, "maxiter": 400}
MAX_ROUNDS = 10
_G: dict = {}


def node_index(lnm, lny):
    """Nearest frozen grid node in (ln M_Lz, ln y); the grid is uniform in both logs."""
    nm, ny = len(m.GRID_M), len(m.GRID_Y)
    i = int(np.clip(np.rint((lnm - LN_M[0]) / (LN_M[1] - LN_M[0]) * (nm - 1)), 0, nm - 1))
    j = int(np.clip(np.rint((lny - LN_Y[0]) / (LN_Y[1] - LN_Y[0]) * (ny - 1)), 0, ny - 1))
    return i * ny + j


def _obj_p(z, t, grid_f):
    return e.objective(z[:10], t, grid_f[node_index(z[10], z[11])])


def run_de(fun, bounds, cfg, seed, x0, args):
    hist, t1 = [], time.time()

    def cb(intermediate_result):
        hist.append(float(intermediate_result.fun))

    r = differential_evolution(
        fun,
        bounds,
        args=args,
        popsize=cfg["popsize"],
        maxiter=cfg["maxiter"],
        mutation=(0.5, 1.0),
        recombination=0.7,
        init="latinhypercube",
        seed=seed,
        tol=1e-4,
        polish=False,
        x0=x0,
        callback=cb,
    )
    return r, {
        "de_fun": float(r.fun),
        "de_message": str(r.message),
        "de_nit": int(r.nit),
        "de_nfev": int(r.nfev),
        "de_spread": float(np.std(r.population_energies)),
        "de_popsize": int(len(r.population_energies)),
        "de_best_every_10": hist[::10] + ([hist[-1]] if hist else []),
        "de_seconds": round(time.time() - t1, 1),
    }


def nm(t, lensf, x):
    return m.c.nm_restart(
        lambda z: e.objective(z, t, lensf), x, e.STEPS, maxfev=400, max_restarts=8, restart_fev=150
    )


def fit_u(t, seed, x0=None):
    r, info = run_de(e.objective, e.BOUNDS, DE_U, seed, x0, (t, None))
    p = nm(t, None, r.x)
    x, f = (np.asarray(p["x"]), p["f"]) if p["f"] < r.fun else (r.x, float(r.fun))
    info.update(
        {"nm_gain": float(r.fun) - f, "nm_converged": p["converged_restart"], "nm_nfev": p["nfev"]}
    )
    return {"f": f, "x": [float(v) for v in x], "lens": None, "log": [info]}


def refine(t, x, lml, y, f):
    """Continuous refinement with exact lens factors; returns x, lens, f, round log."""
    rounds = []
    for rnd in range(MAX_ROUNDS):
        gain = 0.0

        def lobj(v, x=x):
            if not (LN_M[0] <= v[0] <= LN_M[1] and 0.01 <= v[1] <= 5.0):
                return 1.0
            return e.objective(x, t, b.lens_exact(float(np.exp(v[0])), float(v[1]), +1))

        rl = minimize(
            lobj,
            [lml, y],
            method="Nelder-Mead",
            options={
                "maxfev": 40,
                "xatol": 1e-5,
                "fatol": 1e-8,
                "initial_simplex": [[lml, y], [lml + 0.05, y], [lml, y * 1.05]],
            },
        )
        if rl.fun < f:
            gain += f - float(rl.fun)
            f, lml, y = float(rl.fun), float(rl.x[0]), float(rl.x[1])
        p = nm(t, b.lens_exact(float(np.exp(lml)), y, +1), x)
        if p["f"] < f:
            gain += f - p["f"]
            f, x = p["f"], np.asarray(p["x"])
        rounds.append({"round": rnd, "f": f, "gain": gain, "nm_converged": p["converged_restart"]})
        if gain < TOL:
            break
    return np.asarray(x), (float(np.exp(lml)), float(y)), f, rounds


def fit_p(t, seed, x0=None):
    g = m.grid()
    r, info = run_de(_obj_p, BOUNDS_P, DE_P, seed, x0, (t, g["F"]))
    k = node_index(r.x[10], r.x[11])
    node = g["nodes"][k]
    info.update({"de_node": [float(v) for v in node]})
    x, lens, f, rounds = refine(t, r.x[:10], float(np.log(node[0])), float(node[1]), float(r.fun))
    info.update({"refine_rounds": rounds, "refine_converged": rounds[-1]["gain"] < TOL})
    return {"f": f, "x": [float(v) for v in x], "lens": [float(v) for v in lens], "log": [info]}


# ------------------------------------------------------------------ drivers
def _valid_job(args):
    name, seed = args
    spec = e.injections()[name]
    t, truth_mm = e.injection_target(spec)
    t0 = time.time()
    r = fit_u(t, seed) if spec[1] is None else fit_p(t, seed)
    r = e.finish(t, r, "U" if spec[1] is None else "P")
    r.update({"name": name, "seed": seed, "truth_mismatch": truth_mm})
    r["eps_final"] = 1.0 - r["FF_pad64"]
    r["wall_s"] = round(time.time() - t0, 1)
    return r


def _case_job(args):
    model, seed, x0 = args
    _, t = m.target("D")
    t0 = time.time()
    r = fit_u(t, seed, x0) if model == "U" else fit_p(t, seed, x0)
    r = e.finish(t, r, model)
    r.update({"seed": seed, "wall_s": round(time.time() - t0, 1)})
    return r


def run(jobs, fn, out_name, label, extra=None):
    rows = []
    with Pool(4) as pool:
        for r in pool.imap_unordered(fn, jobs):
            rows.append(r)
            print(
                r.get("name", r["model"]),
                r["seed"],
                r["model"],
                round(r["FF_pad64"], 6),
                r["bounds"],
                r["wall_s"],
                flush=True,
            )
            rec = {"label": label, **(extra or {}), "rows": rows}
            (b.OUT / out_name).write_text(json.dumps(rec, indent=1, default=float) + "\n")


def saved_best():
    """S1c family-D R1 solutions re-scored on the same target (exact lens, pad 64)."""
    _, t = m.target("D")
    out = []
    for z in e.saved_d():
        lf = None if z["lens"] is None else b.lens_exact(*z["lens"], +1)
        hp, hc = e.hpc(*e.to_cart(np.asarray(z["x"])))
        if lf is not None:
            hp, hc = hp * lf, hc * lf
        z["FF_pad64"] = e.stat(t, hp, hc, pad=64)
        out.append(z)
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "validate":
        # P jobs first (longer), so the two waves balance on 4 cores
        jobs = [(n, s) for n in ("V5", "V9", "V6", "V7") for s in (11, 12)]
        run(jobs, _valid_job, "s1e_validate.json", "simulated; S1e tool validation (joint DE)")
    elif argv and argv[0] == "case":
        saved = saved_best()
        x0 = {}
        for model in ("U", "P"):
            best = max((z for z in saved if z["model"] == model), key=lambda z: z["FF_pad64"])
            xx = list(best["x"])
            if model == "P":
                xx += [float(np.log(best["lens"][0])), float(np.log(best["lens"][1]))]
            x0[model] = np.asarray(xx)
        jobs = [("P", s, x0["P"]) for s in (21, 22, 23)] + [("U", s, x0["U"]) for s in (21, 22, 23)]
        run(
            jobs,
            _case_job,
            "s1e_case.json",
            "simulated; S1e R1 case fits (joint DE)",
            {"saved_rescored": saved},
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
