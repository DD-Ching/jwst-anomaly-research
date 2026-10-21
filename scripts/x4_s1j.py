"""X4 S1j: U polish convergence and spin-bound validation for the S1i method.

Unit file: docs/state/units/X4.md (S1j; settings frozen before any run). `simulated` / `derived`;
no data. Needs lalsuite and dynesty (separate environment). Reuses S1i unchanged (`x4_s1i_dynesty`:
the dynesty search, the S1e polish, the pad-64 score, the S1f verification) and adds one stage: an
extended polish that continues S1e's Nelder-Mead restarts (fresh simplex, 150 evaluations each, the
same steps and the same convergence rule: a restart gaining < 1e-5) from the best parameters, until
convergence, 10,000 extra evaluations or 20 minutes. The S1e restart state cannot be resumed, so
this is "a new polish from the saved best parameters", not a seamless continuation.

  python scripts/x4_s1j.py truth       # -> results/x4/s1j_truth.json (B1, B2 legality and score)
  python scripts/x4_s1j.py validate    # -> results/x4/s1j_validate.json (B1, B2; one run each)
  python scripts/x4_s1j.py refine      # -> results/x4/s1j_refine.json (S1i R1 seeds 31, 32)
  python scripts/x4_s1j.py decide      # -> results/x4/s1j_decision.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from datetime import datetime
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


s = _load("x4_s1i_dynesty")
e, j = s.e, s.j
OUT = s.OUT
TOL = 1e-5  # unchanged S1e restart convergence rule
RESTART_FEV = 150  # unchanged S1e restart size
EXTRA_FEV = 10_000
EXTRA_S = 20 * 60
DEADLINE_UTC = "2026-10-21T14:05:00+00:00"  # hard stop for all computation
EPS_PASS = 5e-4
REL_PASS = 0.1  # epsilon <= 0.1 * (1 - FF_case)
FF_CASE = 0.9151494033695111  # best U on R1 so far (S1i seed 32); for the relative rule only
AGREE = 2e-4
FF_SCREEN = 0.994987  # rho_req > 30 (ordinary counter-example rule)
B_SEEDS = {"B1": 20261022, "B2": 20261023}


def b_spec(name):
    """New spin-bound injections. Everything not fixed below is one uniform box draw (seed above),
    plus psi uniform in [0, pi).
    B1: primary spin magnitude exactly at the cap, 0.99.
    B2: primary spin magnitude 0.975 (near, not at, the cap) with cos tilt -0.97, the anti-aligned
    orientation of R1's U solutions (chi1z / |chi1| = -0.982 and -0.987 in S1i)."""
    r = np.random.default_rng(B_SEEDS[name])
    x = s.LO + r.random(s.NDIM) * (s.HI - s.LO)
    psi = float(r.uniform(0, np.pi))
    if name == "B1":
        x[2] = 0.99
    else:
        x[2], x[3] = 0.975, -0.97
    mc, q, s1, s2, inc, phi = e.to_cart(x)
    return ((mc, q, s1, s2, inc, phi, psi), None), x


def target(name):
    if name in B_SEEDS:
        return e.injection_target(b_spec(name)[0])
    return s.target(name)


def extend(t, x0, f0):
    """S1e restarts continued from (x0, f0) until a restart gains < TOL (converged), the extra
    evaluation budget or the time budget is used, or the deadline passes. Never returns worse."""
    fun = lambda z: e.objective(z, t)  # noqa: E731
    best_x, best_f = np.asarray(x0, float), float(f0)
    used, t0, traj, stop = 0, time.time(), [], "budget_fev"
    deadline = datetime.fromisoformat(DEADLINE_UTC).timestamp()
    while True:
        if used >= EXTRA_FEV:
            stop = "budget_fev"
            break
        if time.time() - t0 >= EXTRA_S:
            stop = "budget_time"
            break
        if time.time() >= deadline:
            stop = "deadline"
            break
        fev = min(RESTART_FEV, EXTRA_FEV - used)
        simplex = np.vstack(
            [best_x] + [best_x + np.eye(s.NDIM)[i] * e.STEPS[i] for i in range(s.NDIM)]
        )
        r = minimize(
            fun,
            best_x,
            method="Nelder-Mead",
            options={"maxfev": fev, "xatol": 1e-6, "fatol": 1e-9, "initial_simplex": simplex},
        )
        used += int(r.nfev)
        gain = best_f - float(r.fun)
        if r.fun < best_f:
            best_f, best_x = float(r.fun), np.asarray(r.x)
        traj.append({"fev": used, "f": best_f, "gain": gain, "s": round(time.time() - t0, 1)})
        if gain < TOL:
            stop = "converged"
            break
    return {
        "x": best_x,
        "f": best_f,
        "extra_fev": used,
        "stop": stop,
        "trajectory": traj,
        "seconds": round(time.time() - t0, 1),
    }


def finish(t, x, f, extra):
    r = e.finish(t, {"f": f, "x": [float(v) for v in x], "lens": None, "log": []}, "U")
    r.update({k: extra[k] for k in ("extra_fev", "stop", "trajectory", "seconds")})
    return r


def truth():
    out = {"label": "derived; S1j new injections (truth legality and score; never a search start)"}
    for n in B_SEEDS:
        (_, _), x = b_spec(n)
        t, mm = target(n)
        mc, q, s1, s2, inc, phi = e.to_cart(x)
        out[n] = {
            "x_box": [float(v) for v in x],
            "inside_box": bool(np.all((x >= s.LO) & (x <= s.HI))),
            "chi1": list(s1),
            "chi1_mag": float(np.linalg.norm(s1)),
            "chi2_mag": float(np.linalg.norm(s2)),
            "mc": mc,
            "q": q,
            "inc": inc,
            "truth_mismatch": mm,
            "objective_at_truth_box": float(e.objective(x, t)),
        }
    return out


def validate_one(name, seed):
    """S1i procedure (dynesty + S1e polish), then the extended polish; one run."""
    s.DEADLINE_UTC = DEADLINE_UTC
    s.target = target  # the S1i run loop builds its target through this name
    r = s.run_one(name, seed)
    t, _ = target(name)
    ext = extend(t, r["x"], r["f"])
    out = finish(t, ext["x"], ext["f"], ext)
    out.update(
        {
            "name": name,
            "seed": seed,
            "s1i_stage": {
                k: r[k]
                for k in (
                    "stop",
                    "ncall",
                    "iterations",
                    "ff_ns_best_pad32",
                    "nm_gain",
                    "nm_converged",
                    "nm_nfev",
                    "FF_pad64",
                    "wall_s",
                )
            },
        }
    )
    out["eps_final"] = 1.0 - out["FF_pad64"]
    return out


def refine_one(seed):
    """Extended polish of the S1i R1 solution of one seed, from that seed's own saved best."""
    rows = {r["seed"]: r for r in json.loads((OUT / "s1i_r1.json").read_text())["rows"]}
    r0 = rows[seed]
    t, _ = s.target("R1")
    x0 = np.asarray(r0["x"])
    f0 = float(e.objective(x0, t))
    ext = extend(t, x0, f0)
    out = finish(t, ext["x"], ext["f"], ext)
    out.update({"seed": seed, "start_FF_pad64_s1i": r0["FF_pad64"], "start_f_pad32": f0})
    return out


def counted_validation(r):
    return r["s1i_stage"]["stop"] in ("dlogz", "maxcall") and r["stop"] != "deadline"


def decide():
    """Frozen S1j decision (X4.md S1j)."""
    v = {r["name"]: r for r in _rows("s1j_validate.json")}
    vres = {}
    for n in B_SEEDS:
        if n not in v:
            vres[n] = {"status": "unfinished"}
            continue
        r = v[n]
        ok = r["eps_final"] <= EPS_PASS and r["eps_final"] <= REL_PASS * (1 - FF_CASE)
        status = "unfinished" if not counted_validation(r) else ("pass" if ok else "fail")
        vres[n] = {
            "status": status,
            "eps_final": r["eps_final"],
            "s1i_stop": r["s1i_stage"]["stop"],
            "extend_stop": r["stop"],
            "extra_fev": r["extra_fev"],
        }
    v_ok = all(x["status"] == "pass" for x in vres.values())
    out = {
        "label": "derived; S1j decision by frozen rules",
        "validation": vres,
        "validation_pass": v_ok,
    }
    rr = {r["seed"]: r for r in _rows("s1j_refine.json")}
    out["refine"] = {
        sd: {
            "start": r["start_FF_pad64_s1i"],
            "FF_pad64": r["FF_pad64"],
            "stop": r["stop"],
            "extra_fev": r["extra_fev"],
            "converged": r["stop"] == "converged",
            "bounds": r["bounds"],
        }
        for sd, r in rr.items()
    }
    out["verified"] = {sd: s.verify(r["x"]) for sd, r in rr.items()}
    if not v_ok:
        out["outcome"] = (
            "spin-bound validation failed or unfinished; R1 not refined"
            " (tool/optimization undetermined)"
        )
        return out
    if set(rr) != {31, 32}:
        out["outcome"] = "unfinished"
        return out
    ffc = {sd: out["verified"][sd]["FF_conservative"] for sd in rr}
    out["gap"] = abs(ffc[31] - ffc[32])
    out["best_conservative"] = max(ffc.values())
    out["counter_example"] = any(
        out["verified"][sd]["FF_conservative"] - out["verified"][sd]["spread"] > FF_SCREEN
        for sd in rr
    )
    both = all(x["converged"] for x in out["refine"].values())
    if both and out["gap"] <= AGREE:
        out["outcome"] = "U cross-start consistency passed this round (not a global-optimum proof)"
    else:
        out["outcome"] = (
            "tool/optimization undetermined ("
            + ("a refinement did not converge" if not both else f"gap {out['gap']:.2e} > {AGREE}")
            + ")"
        )
    return out


def _rows(name):
    p = OUT / name
    return json.loads(p.read_text())["rows"] if p.exists() else []


def main(argv=None) -> int:
    import subprocess

    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else ""
    if cmd == "truth":
        s.save("s1j_truth.json", truth())
    elif cmd == "one_v":
        r = validate_one(argv[1], int(argv[2]))
        s.append_row(
            "s1j_validate.json", r, "simulated; S1j spin-bound validation (S1i + extended polish)"
        )
        print(r["name"], r["eps_final"], r["stop"], r["extra_fev"], flush=True)
    elif cmd == "one_r":
        r = refine_one(int(argv[1]))
        s.append_row(
            "s1j_refine.json", r, "simulated; S1j extended polish of S1i R1 seeds (own starts)"
        )
        print(r["seed"], r["FF_pad64"], r["stop"], r["extra_fev"], flush=True)
    elif cmd in ("validate", "refine"):
        jobs = (
            [["one_v", "B1", "51"], ["one_v", "B2", "52"]]
            if cmd == "validate"
            else [["one_r", "31"], ["one_r", "32"]]
        )
        procs = [subprocess.Popen([sys.executable, __file__, *jb]) for jb in jobs]
        codes = [p.wait() for p in procs]
        s.save(
            f"s1j_{cmd}_exit.json", {"jobs": [[*jb, c] for jb, c in zip(jobs, codes, strict=True)]}
        )
    elif cmd == "decide":
        s.save("s1j_decision.json", decide())
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
