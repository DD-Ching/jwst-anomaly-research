"""X4 S1i: nested sampling (dynesty) as an alternative global search for the unlensed model U.

Unit file: docs/state/units/X4.md (S1i; settings frozen before any run). `simulated` / `derived`;
no data. Needs lalsuite and dynesty (separate environment). Reuses the S1d box, waveform and score
(`x4_s1d_de.objective`, pad 32 during the search) and the S1e Nelder-Mead restart polish and final
pad-64 score (`x4_s1e_joint.nm`, `x4_s1d_de.finish`) unchanged.

Search distribution (not a physical likelihood): ln L(x) = -(RHO^2 / 2) * (1 - FF(x)^2), a strictly
increasing function of FF, so its maximum is the FF maximum. The prior is uniform on the box. The
evidence of this distribution is never used: it is not a Bayes factor and not physical evidence.

  python scripts/x4_s1i_dynesty.py equiv            # -> results/x4/s1i_equiv.json
  python scripts/x4_s1i_dynesty.py cost             # -> results/x4/s1i_cost.json (V6 pilot)
  python scripts/x4_s1i_dynesty.py validate         # -> results/x4/s1i_validate.json (V7, V10)
  python scripts/x4_s1i_dynesty.py r1               # -> results/x4/s1i_r1.json (seeds 31, 32)
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


j = _load("x4_s1e_joint")
e = j.e
OUT = e.OUT
NDIM = 10
LO = np.array([lo for lo, _ in e.BOUNDS])
HI = np.array([hi for _, hi in e.BOUNDS])
PERIODIC = [4, 7, 9]  # spin azimuths and reference phase, [0, 2 pi)
RHO = 20.0  # scale of the search distribution (frozen; see X4.md S1i)
CFG = {"nlive": 400, "bound": "multi", "sample": "rwalk", "dlogz": 0.5, "maxcall": 220_000}
NPROC = 2  # processes per run; two runs at a time on 4 cores
FF_KNOWN_U = 0.9135105524586078  # S1f seed 23 (s1f_case_u.json); comparison only, never a start
V10_SEED = 20261021  # new validation injection, drawn once from the box
DEADLINE_UTC = "2026-10-21T12:40:00+00:00"
_T: dict = {}


def ff_from_obj(f: float) -> float:
    return 1.0 - float(f)


def loglike_from_ff(ff: float) -> float:
    """Strictly increasing in FF on [0, 1]; FF <= 0 (invalid box point) maps to the floor."""
    ff = max(float(ff), 0.0)
    return -0.5 * RHO * RHO * (1.0 - ff * ff)


def prior_transform(u):
    return LO + np.asarray(u, float) * (HI - LO)


def _init(t):
    _T["t"] = t


def loglike(x):
    return loglike_from_ff(ff_from_obj(e.objective(x, _T["t"])))


def v10_spec():
    """New U injection (never used for tuning): one uniform box draw plus a polarisation angle."""
    r = np.random.default_rng(V10_SEED)
    x = LO + r.random(NDIM) * (HI - LO)
    psi = float(r.uniform(0, np.pi))
    mc, q, s1, s2, inc, phi = e.to_cart(x)
    return ((mc, q, s1, s2, inc, phi, psi), None), x


def target(name):
    if name == "R1":
        return j.m.target("D")[1], None
    spec = v10_spec()[0] if name == "V10" else e.injections()[name]
    t, mm = e.injection_target(spec)
    return t, mm


def run_one(name, seed, cfg=CFG):
    """One nested-sampling run, then the S1e polish of its best sample, then the pad-64 score."""
    import dynesty

    t, truth_mm = target(name)
    t0 = time.time()
    trace = []
    with Pool(NPROC, initializer=_init, initargs=(t,)) as pool:
        s = dynesty.NestedSampler(
            loglike,
            prior_transform,
            NDIM,
            nlive=cfg["nlive"],
            bound=cfg["bound"],
            sample=cfg["sample"],
            periodic=PERIODIC,
            rstate=np.random.default_rng(seed),
            pool=pool,
            queue_size=NPROC,
        )
        deadline = datetime.fromisoformat(DEADLINE_UTC).timestamp()
        stop = "dlogz"
        for it, _res in enumerate(s.sample(maxcall=cfg["maxcall"], dlogz=cfg["dlogz"])):
            if it % 500 == 0:
                lmax = float(np.max(s.saved_run["logl"]))
                trace.append(
                    {
                        "it": it,
                        "ncall": int(s.ncall),
                        "logl_max": lmax,
                        "wall_s": round(time.time() - t0, 1),
                    }
                )
            if time.time() > deadline:
                stop = "deadline"
                break
        else:
            # sample() counts its calls after the nlive initial ones (review fix)
            stop = "maxcall" if s.ncall - cfg["nlive"] > cfg["maxcall"] else "dlogz"
        # dead points plus the final live points (the best point may still be live)
        logl = np.concatenate([np.asarray(s.saved_run["logl"]), np.asarray(s.live_logl)])
        xs = np.concatenate([np.asarray(s.saved_run["v"]), np.asarray(s.live_v)])
        n_dead = len(s.saved_run["logl"])
    k = int(np.argmax(logl))
    x_best = xs[k]
    f_ns = float(e.objective(x_best, t))
    p = j.nm(t, None, x_best)
    x, f = (np.asarray(p["x"]), p["f"]) if p["f"] < f_ns else (x_best, f_ns)
    r = {"f": f, "x": [float(v) for v in x], "lens": None, "log": []}
    r = e.finish(t, r, "U")
    r.update(
        {
            "name": name,
            "seed": seed,
            "stop": stop,
            "iterations": int(n_dead),
            "ncall": int(s.ncall),
            "ff_ns_best_pad32": 1.0 - f_ns,
            "nm_gain": f_ns - f,
            "nm_converged": p["converged_restart"],
            "nm_nfev": p["nfev"],
            "trace": trace,
            "truth_mismatch": truth_mm,
            "wall_s": round(time.time() - t0, 1),
        }
    )
    if truth_mm is not None:
        r["eps_final"] = 1.0 - r["FF_pad64"]
    return r


def equiv():
    """Scores of the adapter vs the original objective at fixed points, and the ln L mapping."""
    t, _ = target("R1")
    pts = {
        "s1f_seed23_best": np.asarray(
            json.loads((OUT / "s1f_case_u.json").read_text())["rows"][2]["x"]
        )
    }
    pts["v7_truth"] = e.from_cart(*[e.injections()["V7"][0][k] for k in range(6)])
    corners = {
        "q15_spins099_inc0": [np.log(e.MC0), 15.0, 0.99, 1.0, 0.0, 0.99, -1.0, 0.0, 0.0, 0.0],
        "mc_low_inc_pi": [LO[0], 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np.pi, 6.28],
        "mc_high_spins_tilted": [HI[0], 8.0, 0.99, -1.0, 6.28, 0.99, 1.0, 3.1, 1.5, 3.0],
    }
    pts.update({k: np.asarray(v, float) for k, v in corners.items()})
    rng = np.random.default_rng(7)
    for i in range(5):
        pts[f"random{i}"] = prior_transform(rng.random(NDIM))
    rows = {}
    for k, x in pts.items():
        f = float(e.objective(x, t))
        _init(t)
        rows[k] = {
            "objective": f,
            "FF": 1.0 - f,
            "loglike": loglike(x),
            "loglike_from_FF": loglike_from_ff(1.0 - f),
            "inside_box": bool(np.all((x >= LO) & (x <= HI))),
        }
    u = rng.random((1000, NDIM))
    xs = np.array([prior_transform(v) for v in u])
    ffs = np.linspace(0, 1, 101)
    ll = [loglike_from_ff(v) for v in ffs]
    return {
        "label": "derived; S1i adapter equivalence (R1 target, pad 32 search score)",
        "points": rows,
        "max_abs_diff_loglike": max(
            abs(r["loglike"] - r["loglike_from_FF"]) for r in rows.values()
        ),
        "s1f_seed23_pad64_now": e.stat(t, *e.hpc(*e.to_cart(pts["s1f_seed23_best"])), pad=64),
        "s1f_seed23_pad64_stored": FF_KNOWN_U,
        "prior_inside_box": bool(np.all((xs >= LO) & (xs <= HI))),
        "loglike_strictly_increasing": bool(np.all(np.diff(ll) > 0)),
        "rho": RHO,
    }


FF_TARGET = FF_KNOWN_U - 2e-4  # R1 success per run (pad 64)
EPS_PASS = 5e-4  # validation (S1e rule)


def counted(r):
    """A row counts if sampling ended by dlogz or the call cap; a deadline stop is unfinished."""
    return r["stop"] in ("dlogz", "maxcall")


def verify(x):
    """S1f's checks on one U solution: box legality, pad 32 / 64 / 256, explicit-psi rebuild."""
    t, _ = target("R1")
    x = np.asarray(x)
    mc, q, s1, s2, inc, phi = e.to_cart(x)
    hp, hc = e.hpc(mc, q, s1, s2, inc, phi)
    psi = j.explicit_psi(t, hp, hc)
    hh = np.cos(2 * psi) * hp + np.sin(2 * psi) * hc
    sc = {
        "stat_pad32": e.stat(t, hp, hc, 32),
        "stat_pad64": e.stat(t, hp, hc, 64),
        "stat_pad256": e.stat(t, hp, hc, 256),
        "explicit_psi_ov_pad64": j.m.ov(t, hh, "D", pad=64),
        "explicit_psi_ov_pad256": j.m.ov(t, hh, "D", pad=256),
    }
    return {
        "in_box": bool(np.all((x >= LO) & (x <= HI))),
        "chi_mags": [float(np.linalg.norm(s1)), float(np.linalg.norm(s2))],
        "scores": sc,
        "spread": max(sc.values()) - min(sc.values()),
        "FF_conservative": min(sc.values()),
    }


def _rows(name):
    p = OUT / name
    return json.loads(p.read_text())["rows"] if p.exists() else []


def decide():
    """Frozen S1i decision (X4.md S1i), applied to the stored rows."""
    vrows = {r["name"]: r for r in _rows("s1i_validate.json")}
    v_ok = all(
        n in vrows and counted(vrows[n]) and vrows[n]["eps_final"] <= EPS_PASS
        for n in ("V7", "V10")
    )
    out = {
        "label": "derived; S1i decision by frozen rules",
        "validation": {
            n: {k: r[k] for k in ("stop", "eps_final", "ncall")} for n, r in vrows.items()
        },
        "validation_pass": v_ok,
        "ff_target": FF_TARGET,
        "ff_known_u": FF_KNOWN_U,
    }
    runs = {r["seed"]: r for r in _rows("s1i_r1.json")}
    out["r1"] = {
        s: {
            "stop": r["stop"],
            "FF_pad64": r["FF_pad64"],
            "ncall": r["ncall"],
            "counted": counted(r),
            "success": counted(r) and r["FF_pad64"] >= FF_TARGET,
        }
        for s, r in runs.items()
    }
    out["verified"] = {s: verify(r["x"]) for s, r in runs.items() if r["FF_pad64"] > FF_KNOWN_U}
    n_ok = sum(v["success"] for v in out["r1"].values())
    if not v_ok:
        out["outcome"] = "validation failed or unfinished; R1 not run or not counted"
    elif len(runs) < 2 or not all(v["counted"] for v in out["r1"].values()):
        out["outcome"] = "unfinished"
    elif n_ok == 2:
        # post hoc wording correction (review): the frozen phrase said "recovered the known better
        # U solution"; the runs meet the target without re-finding a specific point
        out["outcome"] = (
            "both runs met the frozen target from prior-only initialisation"
            " (not a global-optimum proof)"
        )
    else:
        out["outcome"] = f"method did not pass this test ({n_ok} of 2 runs reached the target)"
    return out


def save(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=1, default=float) + "\n")


def append_row(name, row, label):
    """Locked read-modify-write with an atomic rename (review fix: two runs may finish together)."""
    import fcntl
    import os

    p = OUT / name
    with open(OUT / f".{name}.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        d = json.loads(p.read_text()) if p.exists() else {"label": label, "rows": []}
        d["rows"].append(row)
        tmp = OUT / f".{name}.tmp"
        tmp.write_text(json.dumps(d, indent=1, default=float) + "\n")
        os.replace(tmp, p)


def run_pair(jobs, out_name, label):
    """Two runs at once as separate processes (each with NPROC pool workers)."""
    import subprocess

    procs = [
        subprocess.Popen([sys.executable, __file__, "one", n, str(s), out_name, label])
        for n, s in jobs
    ]
    codes = [p.wait() for p in procs]
    # review fix: a crashed run (non-zero exit or no row) counts as unfinished
    save(
        out_name.replace(".json", "_exit.json"),
        {"jobs": [[n, s, c] for (n, s), c in zip(jobs, codes, strict=True)]},
    )
    return codes


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else ""
    if cmd == "equiv":
        save("s1i_equiv.json", equiv())
    elif cmd == "cost":  # pilot on V6 (not a validation injection), capped
        cfg = {**CFG, "maxcall": int(argv[1]) if len(argv) > 1 else 20_000}
        r = run_one("V6", 99, cfg)
        save(
            "s1i_cost.json",
            {
                "label": "derived; S1i cost pilot on V6 (not a validation injection)",
                "cfg": cfg,
                "row": r,
            },
        )
        print(json.dumps({k: r[k] for k in ("stop", "iterations", "ncall", "wall_s", "FF_pad64")}))
    elif cmd == "one":  # one run, row appended to argv[3]
        r = run_one(argv[1], int(argv[2]))
        append_row(argv[3], r, argv[4])
        print(
            r["name"],
            r["seed"],
            r["stop"],
            round(r["FF_pad64"], 6),
            r["ncall"],
            r["wall_s"],
            flush=True,
        )
    elif cmd == "decide":
        save("s1i_decision.json", decide())
    elif cmd == "validate":
        run_pair(
            [("V7", 41), ("V10", 41)],
            "s1i_validate.json",
            "simulated; S1i U validation (dynesty + S1e polish)",
        )
    elif cmd == "r1":
        run_pair(
            [("R1", 31), ("R1", 32)],
            "s1i_r1.json",
            "simulated; S1i R1 U searches (dynesty + S1e polish)",
        )
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
