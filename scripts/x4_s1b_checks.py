"""X4 S1b: optimisation and waveform-model checks on four representative S1 cases.

Unit file: docs/state/units/X4.md (S1b; cases, steps and verdict rules frozen there before computing).
Everything here is `simulated` or `derived`; no strain data is used. Needs lalsuite.

  python scripts/x4_s1b_checks.py rescore    # step A: re-score every saved S1 best fit at pad 64
  python scripts/x4_s1b_checks.py fit        # step B: multi-start fits with restart-based convergence
  python scripts/x4_s1b_checks.py validate   # step C: fit-error check on injected ordinary signals
  python scripts/x4_s1b_checks.py xphm       # step D: precessing higher-mode templates (IMRPhenomXPHM)
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


b = _load("x4_s1_degeneracy")
OUT = b.OUT

CASES = {
    "R1": (30.0, 300.0, 0.3),
    "R2": (60.0, 100.0, 1.5),
    "R3": (200.0, 300.0, 1.0),
    "R4": (200.0, 100.0, 1.0),
}
Q_MAX = 15.0
MC_FRAC = 0.60
FIT_PAD, SCORE_PAD = 32, 64
TOL = 1e-5
N_RANDOM_U, N_RANDOM_P = 12, 3


def target(case):
    mtot, ml, y = CASES[case]
    mc0 = b.chirp(mtot, b.Q_INJ)
    h = b.waveform(mc0, b.Q_INJ, b.CHI_INJ, b.CHI_INJ, approx=b.APPROX)
    return mc0, h * b.lens(ml, y, -1), h


def snr_ratio(case):
    _, t, h = target(case)
    w = b.weight()
    return float(np.sqrt(np.sum(np.abs(t) ** 2 * w) / np.sum(np.abs(h) ** 2 * w)))


def set_domain(q_max=Q_MAX, pad=FIT_PAD):
    b.BOUNDS["q_max"] = q_max
    b.OVERLAP_PAD["pad"] = pad


# ----------------------------------------------------------------------------- saved solutions
def saved_solutions(case):
    """All best-fit parameter sets saved by S1 for this case: (model, stage, params, lens)."""
    key = CASES[case]
    out = []

    def k(c):
        return (c["M_total"], c["M_Lz"], c["y"])

    for name in ("s1_cases.json", "s1_robust.json", "s1_refine.json", "s1_refine_wide.json"):
        path = OUT / name
        if not path.exists():
            continue
        for c in json.loads(path.read_text())["cases"]:
            if k(c) != key:
                continue
            out.append(("U", name, c["U"]["params"], None))
            for tag in ("P", "P_6d", "P_grid_pad16"):
                if tag in c and c[tag].get("params"):
                    out.append(
                        ("P", f"{name}:{tag}", c[tag]["params"], (c[tag]["M_Lz"], c[tag]["y"]))
                    )
    path = OUT / "s1_ubound.json"
    if path.exists():
        for c in json.loads(path.read_text())["cases"]:
            if tuple(c["case"]) == key:
                out.append(("U", "s1_ubound.json", c["U_wide"]["params"], None))
    return out


def score(case, model, params, lensp):
    mc0, t, _ = target(case)
    lensf = None if model == "U" else b.lens_exact(lensp[0], lensp[1], +1)
    return 1.0 - b._intrinsic_obj(
        params, t, lensf, mc0, 10.0
    )  # chirp bound not applied when scoring


def _rescore_job(args):
    set_domain(q_max=1e9, pad=SCORE_PAD)
    case, model, stage, params, lensp = args
    return {
        "case": case,
        "model": model,
        "stage": stage,
        "params": params,
        "lens": lensp,
        "FF_pad64": score(case, model, params, lensp),
    }


# ----------------------------------------------------------------------------- optimiser with restarts
def nm_restart(fun, x0, steps, maxfev=400, max_restarts=4, restart_fev=150):
    """Nelder-Mead, then restarts from the end point with a fresh simplex until a restart improves < TOL."""
    x0 = np.asarray(x0, float)
    history, nfev, messages = [], 0, []

    def run(x, fev):
        simplex = np.vstack([x] + [x + np.eye(x.size)[i] * steps[i] for i in range(x.size)])
        r = minimize(
            fun,
            x,
            method="Nelder-Mead",
            options={"maxfev": fev, "xatol": 1e-6, "fatol": 1e-9, "initial_simplex": simplex},
        )
        return r

    r = run(x0, maxfev)
    nfev += r.nfev
    messages.append(r.message)
    best_f, best_x = float(r.fun), r.x
    history.append(best_f)
    converged = False
    for _ in range(max_restarts):
        r = run(best_x, restart_fev)
        nfev += r.nfev
        messages.append(r.message)
        improve = best_f - float(r.fun)
        if r.fun < best_f:
            best_f, best_x = float(r.fun), r.x
        history.append(best_f)
        if improve < TOL:
            converged = True
            break
    return {
        "f": best_f,
        "x": [float(v) for v in best_x],
        "nfev": int(nfev),
        "history": history,
        "messages": [str(m) for m in messages],
        "converged_restart": converged,
    }


U_STEPS = np.array([0.03, 0.5, 0.2, 0.2])


def bounds_hit(x, q_max):
    lmc, q, c1, c2 = x[:4]
    hits = []
    if q >= 0.99 * q_max:
        hits.append("q_max")
    if q <= 1.01:
        hits.append("q_min")
    if max(abs(c1), abs(c2)) >= 0.99 * b.CHI_MAX:
        hits.append("chi")
    if len(x) > 4:
        if x[4] >= 0.99 * np.log(1e4):
            hits.append("M_Lz_max")
        if x[5] >= 0.99 * 5.0 or x[5] <= 0.0101:
            hits.append("y")
    return hits


def fit_u_job(args):
    case, start, q_max, tag = args
    set_domain(q_max=q_max, pad=FIT_PAD)
    mc0, t, _ = target(case)
    res = nm_restart(lambda x: b._intrinsic_obj(x, t, None, mc0, MC_FRAC), start, U_STEPS)
    set_domain(q_max=q_max, pad=SCORE_PAD)
    res["FF_pad64"] = 1.0 - b._intrinsic_obj(res["x"], t, None, mc0, MC_FRAC)
    res.update(
        {
            "case": case,
            "model": "U",
            "start_tag": tag,
            "q_max": q_max,
            "bounds_hit": bounds_hit(res["x"], q_max),
        }
    )
    return res


def fit_p_job(args):
    """Block descent: intrinsic NM (lens fixed) <-> lens NM (intrinsic fixed), then a 6-D polish."""
    case, start, lensp, tag = args
    set_domain(q_max=Q_MAX, pad=FIT_PAD)
    mc0, t, _ = target(case)
    x = np.asarray(start, float)
    lml, y = np.log(lensp[0]), float(lensp[1])
    hist, nfev = [], 0
    f_prev = 1.0
    for _cycle in range(4):
        lf = b.lens_exact(np.exp(lml), y, +1)
        r1 = nm_restart(
            lambda z, lf=lf: b._intrinsic_obj(z, t, lf, mc0, MC_FRAC),
            x,
            U_STEPS,
            maxfev=300,
            max_restarts=2,
        )
        x = np.asarray(r1["x"])
        nfev += r1["nfev"]

        def lens_obj(v, x=x):
            if not (0.0 <= v[0] <= np.log(1e4) and 0.01 <= v[1] <= 5.0):
                return 1.0
            return b._intrinsic_obj(x, t, b.lens_exact(np.exp(v[0]), v[1], +1), mc0, MC_FRAC)

        r2 = minimize(
            lens_obj,
            [lml, y],
            method="Nelder-Mead",
            options={
                "maxfev": 40,
                "xatol": 1e-4,
                "fatol": 1e-7,
                "initial_simplex": [[lml, y], [lml + 0.1, y], [lml, y + 0.05]],
            },
        )
        nfev += r2.nfev
        lml, y = float(r2.x[0]), float(r2.x[1])
        f = float(r2.fun)
        hist.append(f)
        if f_prev - f < TOL:
            break
        f_prev = f
    memo: dict = {}
    pol = nm_restart(
        lambda v: b._p6_obj(v, t, mc0, MC_FRAC, memo),
        [*x, lml, y],
        np.array([0.01, 0.2, 0.05, 0.05, 0.03, 0.02]),
        maxfev=80,
        max_restarts=1,
        restart_fev=40,
    )
    nfev += pol["nfev"]
    xb = pol["x"] if pol["f"] < hist[-1] else [*x, lml, y]
    set_domain(q_max=Q_MAX, pad=SCORE_PAD)
    ff = 1.0 - b._intrinsic_obj(xb[:4], t, b.lens_exact(np.exp(xb[4]), xb[5], +1), mc0, MC_FRAC)
    return {
        "case": case,
        "model": "P",
        "start_tag": tag,
        "x": [float(v) for v in xb],
        "M_Lz": float(np.exp(xb[4])),
        "y": float(xb[5]),
        "FF_pad64": ff,
        "nfev": int(nfev),
        "block_history": hist,
        "polish_converged": pol["converged_restart"],
        "block_converged": len(hist) >= 2 and hist[-2] - hist[-1] < TOL,
        "bounds_hit": bounds_hit(xb, Q_MAX),
    }


# ----------------------------------------------------------------------------- jobs
def u_starts(case, rng):
    mc0 = b.chirp(CASES[case][0], b.Q_INJ)
    starts = [(s[2], f"saved:{s[1]}") for s in saved_solutions(case) if s[0] == "U"]
    for i in range(N_RANDOM_U):
        starts.append(
            (
                [
                    np.log(mc0 * rng.uniform(0.6, 1.4)),
                    rng.uniform(1, 15),
                    *rng.uniform(-0.9, 0.9, 2),
                ],
                f"random{i}",
            )
        )
    return starts


def p_starts(case, rng):
    saved = [s for s in saved_solutions(case) if s[0] == "P"]
    # best three saved (by S1 FF are not comparable across pads: rescore later); keep distinct lenses
    seen, starts = set(), []
    for s in saved:
        key = (round(s[3][0], 1), round(s[3][1], 2))
        if key not in seen:
            seen.add(key)
            starts.append((s[2], s[3], f"saved:{s[1]}"))
    starts = starts[:4]
    mc0 = b.chirp(CASES[case][0], b.Q_INJ)
    for i in range(N_RANDOM_P):
        starts.append(
            (
                [
                    np.log(mc0 * rng.uniform(0.7, 1.3)),
                    rng.uniform(1, 10),
                    *rng.uniform(-0.8, 0.8, 2),
                ],
                (
                    float(np.exp(rng.uniform(np.log(30), np.log(3000)))),
                    float(rng.uniform(0.2, 3.0)),
                ),
                f"random{i}",
            )
        )
    return starts


def main(argv=None) -> int:
    cmd = (argv or sys.argv[1:] or ["rescore"])[0]
    t0 = time.time()
    rng = np.random.default_rng(2026)
    if cmd == "rescore":
        jobs = [(c, s[0], s[1], s[2], s[3]) for c in CASES for s in saved_solutions(c)]
        with Pool(4) as pool:
            rows = pool.map(_rescore_job, jobs)
        best = {}
        for r in rows:
            k = (r["case"], r["model"])
            if k not in best or r["FF_pad64"] > best[k]["FF_pad64"]:
                best[k] = r
        rec = {
            "label": "derived: S1 saved fits re-scored at pad 64",
            "rows": rows,
            "best": {f"{k[0]}:{k[1]}": v for k, v in best.items()},
        }
        (OUT / "s1b_rescore.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
        for k, v in sorted(best.items()):
            print(k, round(v["FF_pad64"], 5), v["stage"])
    elif cmd == "fit":
        ujobs = [(c, s, q, tag) for c in CASES for (s, tag) in u_starts(c, rng) for q in (Q_MAX,)]
        ujobs += [(c, s, 50.0, f"diag50:{tag}") for c in CASES for (s, tag) in u_starts(c, rng)[:6]]
        pjobs = [(c, s, lp, tag) for c in CASES for (s, lp, tag) in p_starts(c, rng)]
        with Pool(4) as pool:
            urows = pool.map(fit_u_job, ujobs)
            print(f"U fits done in {time.time() - t0:.0f} s", flush=True)
            prows = []
            for r in pool.imap_unordered(fit_p_job, pjobs):
                prows.append(r)
                print(r["case"], r["start_tag"], round(r["FF_pad64"], 5), r["nfev"], flush=True)
        rec = {
            "label": "simulated; multi-start fits with restart convergence (S1b step B)",
            "wall_s": round(time.time() - t0, 1),
            "U": urows,
            "P": prows,
        }
        (OUT / "s1b_fits.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
    elif cmd == "validate":
        fits = json.loads((OUT / "s1b_fits.json").read_text())
        jobs = []
        for c in CASES:
            ub = max(
                (r for r in fits["U"] if r["case"] == c and r["q_max"] == Q_MAX),
                key=lambda r: r["FF_pad64"],
            )
            pb = max((r for r in fits["P"] if r["case"] == c), key=lambda r: r["FF_pad64"])
            jobs.append(("U", c, ub["x"], None))
            jobs.append(("P", c, pb["x"][:4], (pb["M_Lz"], pb["y"])))
        with Pool(4) as pool:
            rows = pool.map(_validate_job, jobs)
        (OUT / "s1b_validate.json").write_text(
            json.dumps(
                {"label": "simulated; fit-error check on ordinary injections", "rows": rows},
                indent=1,
                default=float,
            )
            + "\n"
        )
        for r in rows:
            print(r["case"], r["model"], "fit error", r["fit_error"])
    elif cmd == "xphm":
        _load("x4_s1b_xphm").main(argv[1:] if argv else sys.argv[2:])
    return 0


def _validate_job(args):
    """Inject an ordinary signal (a case's best fit), refit from perturbed starts, report 1 - FF."""
    model, case, x, lensp = args
    set_domain(q_max=Q_MAX, pad=FIT_PAD)
    mc0 = b.chirp(CASES[case][0], b.Q_INJ)
    h = b.waveform(float(np.exp(x[0])), *x[1:4])
    rng = np.random.default_rng(99)
    if model == "U":
        tgt = h
        best = 0.0
        for _ in range(3):
            s = [
                x[0] + rng.normal(0, 0.03),
                min(14.0, max(1.0, x[1] * rng.uniform(0.7, 1.3))),
                float(np.clip(x[2] + rng.normal(0, 0.2), -0.9, 0.9)),
                float(np.clip(x[3] + rng.normal(0, 0.2), -0.9, 0.9)),
            ]
            r = nm_restart(lambda z: b._intrinsic_obj(z, tgt, None, mc0, MC_FRAC), s, U_STEPS)
            best = max(best, 1.0 - r["f"])
    else:
        tgt = h * b.lens_exact(lensp[0], lensp[1], +1)
        s = [
            x[0] + 0.02,
            min(14.0, max(1.0, x[1] * 1.2)),
            float(np.clip(x[2] - 0.15, -0.9, 0.9)),
            float(np.clip(x[3] + 0.15, -0.9, 0.9)),
        ]
        r = fit_p_job_target(tgt, mc0, s, (lensp[0] * 1.25, lensp[1] + 0.15))
        best = r
    set_domain(q_max=Q_MAX, pad=SCORE_PAD)
    return {"model": model, "case": case, "truth": x, "lens": lensp, "fit_error": 1.0 - best}


def fit_p_job_target(tgt, mc0, start, lensp):
    """Same block descent as fit_p_job, on an arbitrary target; returns the best FF (pad 32)."""
    x = np.asarray(start, float)
    lml, y = np.log(lensp[0]), float(lensp[1])
    best = 0.0
    f_prev = 1.0
    for _cycle in range(4):
        lf = b.lens_exact(np.exp(lml), y, +1)
        r1 = nm_restart(
            lambda z, lf=lf: b._intrinsic_obj(z, tgt, lf, mc0, MC_FRAC),
            x,
            U_STEPS,
            maxfev=300,
            max_restarts=2,
        )
        x = np.asarray(r1["x"])

        def lens_obj(v, x=x):
            if not (0.0 <= v[0] <= np.log(1e4) and 0.01 <= v[1] <= 5.0):
                return 1.0
            return b._intrinsic_obj(x, tgt, b.lens_exact(np.exp(v[0]), v[1], +1), mc0, MC_FRAC)

        r2 = minimize(
            lens_obj,
            [lml, y],
            method="Nelder-Mead",
            options={
                "maxfev": 40,
                "xatol": 1e-4,
                "fatol": 1e-7,
                "initial_simplex": [[lml, y], [lml + 0.1, y], [lml, y + 0.05]],
            },
        )
        lml, y = float(r2.x[0]), float(r2.x[1])
        best = max(best, 1.0 - float(r2.fun))
        if f_prev - float(r2.fun) < TOL:
            break
        f_prev = float(r2.fun)
    return best


if __name__ == "__main__":
    raise SystemExit(main())
