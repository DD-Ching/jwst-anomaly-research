"""X4 S1k: one fixed graded Nelder-Mead polish (step scales 1, 0.1, 0.01) for the unlensed model U.

Unit file: docs/state/units/X4.md (S1k; settings frozen before any run). `simulated` / `derived`;
no data. Needs lalsuite and dynesty (separate environment). Reuses SciPy's Nelder-Mead
(`initial_simplex`, `maxfev`, `xatol`, `fatol`; its status and evaluation count are recorded), the
S1d score (`x4_s1d_de.objective`, pad 32), the pad-64 final score, the S1i search and the S1f
verification, all unchanged.

The graded polish runs three SciPy Nelder-Mead stages in a fixed order, each from the best
legal point so far, with the simplex steps set to S1e's STEPS times 1, 0.1 and 0.01. Every stage
runs, whatever the previous one gained. Periodic angles (spin azimuths, reference phase) are wrapped
into [0, 2 pi) before scoring; other simplex vertices that would leave the box are mirrored inside
it. The best legal point ever scored is kept. Budget: 10,000 evaluations (3,000 / 3,500 / 3,500) or
20 minutes per case.

  python scripts/x4_s1k.py dev        # -> results/x4/s1k_dev.json (B1, B2, V7, V10 saved solutions)
  python scripts/x4_s1k.py truth      # -> results/x4/s1k_truth.json (C1, C2; after the freeze)
  python scripts/x4_s1k.py validate   # -> results/x4/s1k_validate.json (C1, C2; one run each)
  python scripts/x4_s1k.py refine     # -> results/x4/s1k_refine.json (S1i R1 seeds 31, 32)
  python scripts/x4_s1k.py decide     # -> results/x4/s1k_decision.json
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


jj = _load("x4_s1j")
s, e = jj.s, jj.e
OUT = s.OUT
SCALES = (1.0, 0.1, 0.01)
STAGE_FEV = (3000, 3500, 3500)  # sums to 10,000
CASE_S = 20 * 60
XATOL, FATOL = 1e-6, 1e-9  # unchanged S1e Nelder-Mead tolerances
PERIODIC = s.PERIODIC
DEADLINE_UTC = "2026-10-21T22:50:00+00:00"
EPS_PASS = 5e-4
REL_PASS = 0.1
FF_CASE = 0.9151494033695111  # best U on R1 so far (S1i seed 32)
AGREE = 2e-4
FF_SCREEN = 0.994987
C_SPECS = {  # frozen before the dev runs; generated only after the code freeze
    "C1": {
        "seed": 20261031,
        "q": (2.4, 2.8),
        "a1": (0.97, 0.99),
        "cos1": (-1.0, -0.95),
        "search_seed": 71,
    },
    "C2": {
        "seed": 20261032,
        "q": (8.0, 11.0),
        "a1": (0.95, 0.99),
        "cos1": (-1.0, -0.9),
        "search_seed": 72,
    },
}


class _Stop(Exception):
    pass


def wrap(z):
    z = np.array(z, float)
    z[PERIODIC] = np.mod(z[PERIODIC], 2 * np.pi)
    return z


def simplex(x, scale):
    """x plus one vertex per axis at scale * STEPS; periodic axes wrap later, others mirror."""
    lo, hi = s.LO, s.HI
    verts = [x.copy()]
    for i in range(s.NDIM):
        d = scale * e.STEPS[i]
        v = x.copy()
        if i in PERIODIC or x[i] + d <= hi[i]:
            v[i] = x[i] + d
        elif x[i] - d >= lo[i]:
            v[i] = x[i] - d
        else:  # step wider than the room on both sides: half the larger room
            v[i] = x[i] + 0.5 * (
                (hi[i] - x[i]) if hi[i] - x[i] >= x[i] - lo[i] else -(x[i] - lo[i])
            )
        verts.append(v)
    return np.array(verts)


def graded(fun, x0, scales=SCALES, stage_fev=STAGE_FEV, case_s=CASE_S, deadline=None):
    """The fixed graded polish. `fun` scores a legal box point (lower is better).

    `deadline` defaults to the module's DEADLINE_UTC read at call time, so a patched constant
    reaches the polish and the S1i stage alike.
    """
    t0 = time.time()
    t_dead = datetime.fromisoformat(deadline or DEADLINE_UTC).timestamp()
    st = {"n": 0, "best_f": None, "best_x": None}

    def f(z):
        if st["n"] and time.time() - t0 > case_s:  # the start point is always scored
            raise _Stop("time_cap")
        if st["n"] and time.time() > t_dead:
            raise _Stop("deadline")
        z = wrap(z)
        v = float(fun(z))
        st["n"] += 1
        if st["best_f"] is None or v < st["best_f"]:
            st["best_f"], st["best_x"] = v, z
        return v

    x = wrap(x0)
    f0 = f(x)
    stages, stop = [], "completed"
    for sc, fev in zip(scales, stage_fev, strict=True):
        start_f, n_start = st["best_f"], st["n"]
        x = st["best_x"].copy()
        try:
            r = minimize(
                f,
                x,
                method="Nelder-Mead",
                options={
                    "maxfev": fev,
                    "xatol": XATOL,
                    "fatol": FATOL,
                    "initial_simplex": simplex(x, sc),
                },
            )
            status, msg = int(r.status), str(r.message)
        except _Stop as ex:
            status, msg, stop = -1, str(ex), str(ex)
        stages.append(
            {
                "scale": sc,
                "start_f": start_f,
                "end_f": st["best_f"],
                "gain": start_f - st["best_f"],
                "nfev": st["n"] - n_start,
                "scipy_status": status,
                "message": msg,
            }
        )
        if stop != "completed":
            break
    return {
        "x": st["best_x"],
        "f": st["best_f"],
        "f0": f0,
        "nfev": st["n"],
        "stop": stop,
        "stages": stages,
        "seconds": round(time.time() - t0, 1),
        "last_stage_converged": bool(
            stages and stages[-1]["scipy_status"] == 0 and len(stages) == 3
        ),
    }


def _score_row(t, g, extra=None):
    r = e.finish(t, {"f": g["f"], "x": [float(v) for v in g["x"]], "lens": None, "log": []}, "U")
    r.update({k: g[k] for k in ("f0", "nfev", "stop", "stages", "seconds", "last_stage_converged")})
    r.update(extra or {})
    return r


def dev_cases():
    """Saved solutions of the development cases (no new sampling)."""
    sj = {r["name"]: r for r in json.loads((OUT / "s1j_validate.json").read_text())["rows"]}
    si = {r["name"]: r for r in json.loads((OUT / "s1i_validate.json").read_text())["rows"]}
    return {"B1": sj["B1"]["x"], "B2": sj["B2"]["x"], "V7": si["V7"]["x"], "V10": si["V10"]["x"]}


def dev_one(name, x0, control):
    t, _ = jj.target(name)
    sc = (1.0, 1.0, 1.0) if control else SCALES
    g = graded(
        lambda z: e.objective(z, t), np.asarray(x0), scales=sc, deadline="2099-01-01T00:00:00+00:00"
    )
    r = _score_row(t, g, {"name": name, "procedure": "control_scale1x3" if control else "graded"})
    r["eps_final"] = 1.0 - r["FF_pad64"]
    r["eps_start_pad32"] = g["f0"]
    return r


def c_spec(name):
    c = C_SPECS[name]
    r = np.random.default_rng(c["seed"])
    x = s.LO + r.random(s.NDIM) * (s.HI - s.LO)
    x[1] = r.uniform(*c["q"])
    x[2] = r.uniform(*c["a1"])
    x[3] = r.uniform(*c["cos1"])
    psi = float(r.uniform(0, np.pi))
    mc, q, s1, s2, inc, phi = e.to_cart(x)
    return ((mc, q, s1, s2, inc, phi, psi), None), x


def target(name):
    if name in C_SPECS:
        return e.injection_target(c_spec(name)[0])
    return jj.target(name)


def truth():
    out = {
        "label": "derived; S1k new injections (made after the code freeze; never a search start)"
    }
    for n in C_SPECS:
        _, x = c_spec(n)
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


def validate_one(name):
    s.DEADLINE_UTC = DEADLINE_UTC
    s.target = (
        target  # S1i run loop builds its target through this name; non-C names go to jj.target
    )
    r = s.run_one(name, C_SPECS[name]["search_seed"])
    t, _ = target(name)
    g = graded(lambda z: e.objective(z, t), np.asarray(r["x"]))
    out = _score_row(
        t,
        g,
        {
            "name": name,
            "seed": C_SPECS[name]["search_seed"],
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
        },
    )
    out["eps_final"] = 1.0 - out["FF_pad64"]
    return out


def refine_one(seed):
    rows = {r["seed"]: r for r in json.loads((OUT / "s1i_r1.json").read_text())["rows"]}
    t, _ = s.target("R1")
    g = graded(lambda z: e.objective(z, t), np.asarray(rows[seed]["x"]))
    return _score_row(t, g, {"seed": seed, "start_FF_pad64_s1i": rows[seed]["FF_pad64"]})


def _rows(name):
    p = OUT / name
    return json.loads(p.read_text())["rows"] if p.exists() else []


def decide():
    """Frozen S1k decision (X4.md S1k)."""
    v = {}
    # frozen: one run each, so a later appended rerun never counts
    for r in _rows("s1k_validate.json"):
        v.setdefault(r["name"], r)
    vres = {}
    for n in C_SPECS:
        r = v.get(n)
        if r is None or r["s1i_stage"]["stop"] == "deadline" or r["stop"] == "deadline":
            vres[n] = {"status": "unfinished"}
            continue
        ok = r["eps_final"] <= EPS_PASS and r["eps_final"] <= REL_PASS * (1 - FF_CASE)
        vres[n] = {
            "status": "pass" if ok else "fail",
            "eps_final": r["eps_final"],
            "polish_stop": r["stop"],
            "s1i_stop": r["s1i_stage"]["stop"],
        }
    v_ok = all(x["status"] == "pass" for x in vres.values())
    out = {
        "label": "derived; S1k decision by frozen rules",
        "validation": vres,
        "validation_pass": v_ok,
    }
    if not v_ok:
        out["outcome"] = "graded polish failed or did not finish its new validation; R1 not refined"
        return out
    rr = {r["seed"]: r for r in _rows("s1k_refine.json")}
    if set(rr) != {31, 32} or any(r["stop"] == "deadline" for r in rr.values()):
        out["outcome"] = "unfinished"
        return out
    out["verified"] = {sd: s.verify(r["x"]) for sd, r in rr.items()}
    ffc = {sd: out["verified"][sd]["FF_conservative"] for sd in rr}
    out["refine"] = {
        sd: {
            "start": r["start_FF_pad64_s1i"],
            "FF_pad64": r["FF_pad64"],
            "gain": r["FF_pad64"] - r["start_FF_pad64_s1i"],
            "stop": r["stop"],
            "last_stage_converged": r["last_stage_converged"],
        }
        for sd, r in rr.items()
    }
    out["gap"] = abs(ffc[31] - ffc[32])
    out["counter_example"] = any(
        out["verified"][sd]["FF_conservative"] - out["verified"][sd]["spread"] > FF_SCREEN
        for sd in rr
    )
    conv = all(r["last_stage_converged"] for r in rr.values())
    if conv and out["gap"] <= AGREE:
        out["outcome"] = "U cross-start consistency passed this round (not a global-optimum proof)"
    else:
        out["outcome"] = (
            "tool/optimization undetermined ("
            + (
                "a polish was budget-limited in its last stage"
                if not conv
                else f"gap {out['gap']:.2e} > {AGREE}"
            )
            + ")"
        )
    return out


def main(argv=None) -> int:
    import subprocess

    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else ""
    if cmd == "dev_one":
        r = dev_one(argv[1], dev_cases()[argv[1]], argv[2] == "control")
        s.append_row(
            "s1k_dev.json",
            r,
            "simulated; S1k development comparison (saved solutions; no new sampling)",
        )
    elif cmd == "dev":
        jobs = [[n, p] for n in ("B2", "B1", "V7", "V10") for p in ("graded", "control")]
        procs, codes = [], []
        for jb in jobs:  # at most 4 at a time
            procs.append(subprocess.Popen([sys.executable, __file__, "dev_one", *jb]))
            if len(procs) == 4:
                codes += [p.wait() for p in procs]
                procs = []
        codes += [p.wait() for p in procs]
        if any(codes):  # a crashed child leaves its row missing from s1k_dev.json
            print(f"dev: {sum(c != 0 for c in codes)} of {len(jobs)} jobs failed", file=sys.stderr)
            return 1
    elif cmd == "truth":
        s.save("s1k_truth.json", truth())
    elif cmd == "one_v":
        r = validate_one(argv[1])
        s.append_row(
            "s1k_validate.json", r, "simulated; S1k validation (S1i search + graded polish)"
        )
    elif cmd == "one_r":
        r = refine_one(int(argv[1]))
        s.append_row(
            "s1k_refine.json", r, "simulated; S1k graded polish of S1i R1 seeds (own starts)"
        )
    elif cmd in ("validate", "refine"):
        jobs = (
            [["one_v", "C1"], ["one_v", "C2"]]
            if cmd == "validate"
            else [["one_r", "31"], ["one_r", "32"]]
        )
        procs = [subprocess.Popen([sys.executable, __file__, *jb]) for jb in jobs]
        codes = [p.wait() for p in procs]
        s.save(
            f"s1k_{cmd}_exit.json", {"jobs": [[*jb, c] for jb, c in zip(jobs, codes, strict=True)]}
        )
    elif cmd == "decide":
        s.save("s1k_decision.json", decide())
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
