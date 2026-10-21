"""X4 S1h: one bounded extra-seed check of R1's ordinary fits, then the full decision.

Unit file: docs/state/units/X4.md (S1h; plan frozen and committed before any fit). Labels: the R1
target is `simulated`; fitting factors are `derived`. No data is read. Fits need lalsuite; the rule
functions run offline (tests/test_x4_s1h.py).

Only this batch: P (positive point lens) seeds 24-26 with S1g's procedure and budget (S1e joint DE,
1,000 generations, 20 refinement rounds) and U (unlensed) seeds 24-26 with S1f's procedure and
budget (S1e U DE, 600 generations, then Nelder-Mead restarts). Initial populations hold the best
re-scored saved S1c solution of that model, as in S1g (P) and S1f (U). One pool of 4 processes.
Everything is reused from `x4_s1g` / `x4_s1e_joint`; new results go to results/x4/s1h_*.json only.

  python scripts/x4_s1h.py check      # cache / version check before the fits -> s1h_check.json
  python scripts/x4_s1h.py run        # the 6 fits (resumable) -> s1h_r1.json
  python scripts/x4_s1h.py decide     # pooled evaluation and verification -> s1h_decision.json
"""

from __future__ import annotations

import os

# Thread limits before numpy / scipy load their backends (every fit single-threaded).
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import importlib.util  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from datetime import datetime  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g = _load("x4_s1g")
j, e, m, b = g.j, g.e, g.m, g.b
OUT = g.OUT

# ------------------------------------------------------------------ frozen settings (X4.md S1h)
SEEDS = (24, 25, 26)  # the next unused seeds after S1f / S1g's 21-23
MODELS = ("P", "U")  # P jobs first: they are the longer ones
N_PROC = 4
DEADLINE_UTC = "2026-10-21T08:15:00+00:00"  # computation stops; running fits -> "unfinished"
SEED_AGREE = g.SEED_AGREE  # 2e-4
FF_SCREEN = g.FF_SCREEN  # 0.994987
ALLOWED_PASS_WORDING = (
    "under the tested models and this round's reliability standard, this supports planning a "
    "limited real-data trial"
)
UNDETERMINED_WORDING = "tool/optimization undetermined"


def jobs_order():
    return [(mdl, s) for mdl in MODELS for s in SEEDS]


# ------------------------------------------------------------------ rule functions (offline)
def x0_for(model: str, saved: list):
    """Initial-population member, as in the original procedures: the best re-scored saved S1c
    solution of that model (P: with ln M_Lz, ln y appended, as S1g's R1; U: as S1f)."""
    if model == "P":
        best, x0 = g._x0_saved_best_p(saved)
        return best, x0
    best = max((z for z in saved if z["model"] == "U"), key=lambda z: z["FF_pad64"])
    return best, np.asarray(best["x"])


def seed_ffs(rows) -> dict:
    """Finished seeds only: {seed: FF_pad64}. Re-scores and partial states are not seeds."""
    return {int(r["seed"]): float(r["FF_pad64"]) for r in rows if r.get("status") == "finished"}


def top_two_gap(ffs):
    v = sorted((float(x) for x in ffs), reverse=True)
    return (v[0], v[1], v[0] - v[1]) if len(v) >= 2 else (v[0] if v else None, None, None)


def best_of(cands: list):
    """The best candidate by FF_pad64 (seeds, re-scored saved solutions, partial states)."""
    cands = [c for c in cands if c.get("FF_pad64") is not None]
    return max(cands, key=lambda c: c["FF_pad64"]) if cands else None


def validation_statuses(eps_by_fit: dict, ff_case: float) -> dict:
    """Recompute every validation label with the updated FF_case (not carried over)."""
    return {k: g.fit_status(v, ff_case) for k, v in eps_by_fit.items()}


def s1h_outcome(
    tool_p: str,
    tool_u: str,
    p_seeds: list,
    u_seeds: list,
    best_ff: float,
    counterexample_verified: bool,
    feasibility_verified: bool,
    n_unfinished: int,
    ratio: float = g.SNR_RATIO_R1,
) -> dict:
    """Decision order frozen in X4.md S1h."""
    rho = g.rho_req(best_ff)
    intrinsic = rho / ratio
    base = {"best_ordinary_FF": best_ff, "rho_req": rho, "intrinsic_snr_needed": intrinsic}
    if counterexample_verified:
        return {"outcome": 1, "wording": "R1 does not pass the screen", "causes": [], **base}
    if feasibility_verified:
        return {
            "outcome": "feasibility",
            "wording": "R1 fails the screen's intrinsic-SNR feasibility condition",
            "causes": ["verified ordinary solution needs intrinsic SNR > 100 (rho_req <= 30)"],
            **base,
        }
    causes = []
    if tool_p != "pass":
        causes.append(f"P tool validation with the updated FF_case: {tool_p}")
    if tool_u != "pass":
        causes.append(f"U tool validation with the updated FF_case: {tool_u}")
    if n_unfinished:
        causes.append(f"fits unfinished at the stop: {n_unfinished}")
    for name, ffs in (("P", p_seeds), ("U", u_seeds)):
        if not g.seeds_agree(ffs, SEED_AGREE):
            gap = top_two_gap(ffs)[2]
            causes.append(
                f"{name} search reliability: two best seeds differ by {gap:.3g} > 2e-4"
                if gap is not None
                else f"{name} search reliability: fewer than 2 finished seeds"
            )
    if best_ff > FF_SCREEN:
        causes.append("FF above the screen threshold but not verified")
    elif intrinsic > g.INTRINSIC_MAX:
        causes.append("intrinsic SNR needed > 100 but not verified")
    if not causes:
        return {"outcome": 2, "wording": ALLOWED_PASS_WORDING, "causes": [], **base}
    return {"outcome": 3, "wording": UNDETERMINED_WORDING, "causes": causes, **base}


# ------------------------------------------------------------------ fits
def _key(model: str, seed: int) -> str:
    return f"s1h_{model}_R1_s{seed}"


def _job(args):
    """One R1 fit with the original case job (`x4_s1e_joint._case_job`): fit_p (S1g budget via
    `x4_s1g.apply_budget`) or fit_u (S1e / S1f budget, unchanged)."""
    model, seed, x0, run_id = args
    g._RUN.clear()
    g._RUN.update({"key": _key(model, seed), "run_id": run_id})
    t0 = time.time()
    r = j._case_job((model, seed, x0))
    info = r["log"][0]
    info["de_spread_every_10"] = g._RUN.get("spread_every_10", [])
    hook = g._RUN.get("refine_hook", [])
    for rr, hh in zip(info.get("refine_rounds", []), hook, strict=False):
        rr.update(hh)
    info["evaluations"] = {
        "de": info["de_nfev"],
        "nm_intrinsic": int(sum(h["nm_nfev"] for h in hook)),
        "refine_lens": int(sum(h["lens_evals"] for h in hook)) if model == "P" else 0,
    }
    r.update({"name": "R1", "model": model, "status": "finished", "run_id": run_id})
    r["wall_s"] = round(time.time() - t0, 1)
    r["budget"] = (
        {"DE_maxiter": j.DE_P["maxiter"], "refine_round_cap": j.MAX_ROUNDS}
        if model == "P"
        else {"DE_maxiter": j.DE_U["maxiter"], "nm_max_restarts": 8}
    )
    return r


def _partial_row(job) -> dict:
    """A fit stopped at the deadline: its last saved state of this run, re-scored exactly."""
    model, seed, _x0, run_id = job
    row = {"name": "R1", "model": model, "seed": seed, "status": "unfinished", "run_id": run_id}
    st = g.read_state(_key(model, seed), run_id)
    if st is None:
        return row
    x = np.asarray(st["x"][:10])
    lens = None
    if model == "P":
        if st["stage"] == "de":
            node = m.grid()["nodes"][j.node_index(st["x"][10], st["x"][11])]
            lens = [float(node[0]), float(node[1])]
        else:
            lens = st["lens"]
    _, t = m.target("D")
    r = e.finish(t, {"f": float(st["f"]), "x": [float(v) for v in x], "lens": lens}, model)
    row.update({"partial_state": st, "x": r["x"], "lens": lens, "FF_pad64": r["FF_pad64"]})
    row.update({"bounds": r["bounds"], "physical": r["physical"]})
    return row


def _collect(job, ar) -> dict:
    try:
        return ar.get()
    except Exception as exc:  # noqa: BLE001
        model, seed, _x0, run_id = job
        row = {"name": "R1", "model": model, "seed": seed, "status": "error", "run_id": run_id}
        return {**row, "error": repr(exc)}


def _dump(name: str, rec: dict) -> None:
    (OUT / name).write_text(json.dumps(rec, indent=1, default=float) + "\n")


def main_run() -> None:
    g.apply_budget()
    m.grid()  # workers inherit the loaded grid
    run_id = g._run_id()
    path = OUT / "s1h_r1.json"
    saved = j.saved_best()
    if path.exists():  # resume: keep finished rows, run only missing or unfinished seeds
        rec = json.loads(path.read_text())
        rec["rows"] = [r for r in rec["rows"] if r.get("status") == "finished"]
        rec.setdefault("runs", []).append({"run_id": run_id, "start_utc": g._now()})
    else:
        rec = {
            "label": "simulated R1 family-D target; derived FF; S1h extra-seed R1 fits "
            "(P: S1g procedure; U: S1f procedure)",
            "deadline_utc": DEADLINE_UTC,
            "runs": [{"run_id": run_id, "start_utc": g._now()}],
            "saved_rescored": saved,
            "rows": [],
        }
    x0 = {}
    for mdl in MODELS:
        best, x0[mdl] = x0_for(mdl, saved)
        rec.setdefault("x0", {})[mdl] = {"start": best["start"], "FF_pad64": best["FF_pad64"]}
    done = {(r["model"], r["seed"]) for r in rec["rows"]}
    jobs = [(mdl, s, x0[mdl], run_id) for mdl, s in jobs_order() if (mdl, s) not in done]
    _dump("s1h_r1.json", rec)
    for jb in jobs:  # this run's state files start empty
        g._state_path(_key(jb[0], jb[1])).unlink(missing_ok=True)
    deadline = datetime.fromisoformat(DEADLINE_UTC).timestamp()

    def on_row(r):
        rec["rows"].append(r)
        _dump("s1h_r1.json", rec)
        print(g._now(), r["model"], r["seed"], r["status"], r.get("FF_pad64"), flush=True)

    with Pool(N_PROC, initializer=g._init_worker) as pool:
        pending = {(jb[0], jb[1]): (jb, pool.apply_async(_job, (jb,))) for jb in jobs}
        while pending:
            for k, (jb, ar) in list(pending.items()):
                if ar.ready():
                    del pending[k]
                    on_row(_collect(jb, ar))
            if pending and time.time() >= deadline:
                pool.terminate()
                break
            if pending:
                time.sleep(15)
    for jb, _ar in pending.values():
        on_row(_partial_row(jb))
    rec["runs"][-1]["end_utc"] = g._now()
    _dump("s1h_r1.json", rec)


# ------------------------------------------------------------------ verification and decision
def verify(x, lens, label="") -> dict:
    """S1f's three steps for a lensed (P) or unlensed (U, lens None) solution."""
    if lens is not None:
        return g.verify_p(x, lens, label)
    x = np.asarray(x, float)
    _, t = m.target("D")
    mc, q, s1, s2, inc, phi = e.to_cart(x)
    err = None
    try:
        hp, hc = e.hpc(mc, q, s1, s2, inc, phi)
    except Exception as exc:  # noqa: BLE001
        err = repr(exc)
    out = {
        "label": label,
        "x_box": [float(v) for v in x],
        "lens": None,
        "in_box": g.in_box(x, None),
        "chi_mags": [float(np.linalg.norm(s1)), float(np.linalg.norm(s2))],
        "q": q,
        "waveform_error": err,
    }
    if err is not None:
        out.update({"counterexample_verified": False, "feasibility_fail_verified": False})
        return out
    psi = j.explicit_psi(t, hp, hc)
    hh = np.cos(2 * psi) * hp + np.sin(2 * psi) * hc
    scores = {
        "stat_pad32": e.stat(t, hp, hc, 32),
        "stat_pad64": e.stat(t, hp, hc, 64),
        "stat_pad256": e.stat(t, hp, hc, 256),
        "explicit_psi_ov_pad64": m.ov(t, hh, "D", pad=64),
        "explicit_psi_ov_pad256": m.ov(t, hh, "D", pad=256),
    }
    spread = max(scores.values()) - min(scores.values())
    ff_cons = min(scores.values())
    rho = g.rho_req(ff_cons)
    out.update(
        {
            "psi": psi,
            "scores": scores,
            "score_spread": spread,
            "FF_conservative": ff_cons,
            "FF_screen": FF_SCREEN,
            "FF_feasible": g.FF_FEASIBLE,
            **g.screen_flags(ff_cons, spread, g.is_legal(out["in_box"], out["chi_mags"], q, err)),
            "rho_req_observed": rho,
            "intrinsic_snr_needed": rho / g.SNR_RATIO_R1,
        }
    )
    return out


def pooled() -> dict:
    """All R1 ordinary solutions per model: S1h rows, S1g (P) / S1f (U) seeds, saved S1c."""
    h = json.loads((OUT / "s1h_r1.json").read_text())
    s1g = json.loads((OUT / "s1g_r1.json").read_text())
    s1f = json.loads((OUT / "s1f_case_u.json").read_text())
    out = {}
    for mdl, old_rows, old_src in (("P", s1g["rows"], "S1g"), ("U", s1f["rows"], "S1f")):
        new = [r for r in h["rows"] if r["model"] == mdl]
        old = [{**r, "status": r.get("status", "finished")} for r in old_rows]
        seeds = {**seed_ffs(old), **seed_ffs(new)}
        cands = (
            [
                {
                    "src": f"{old_src} seed {r['seed']}",
                    "FF_pad64": r["FF_pad64"],
                    "x": r["x"],
                    "lens": r.get("lens"),
                }
                for r in old
            ]
            + [
                {
                    "src": f"S1h seed {r['seed']} ({r['status']})",
                    "FF_pad64": r.get("FF_pad64"),
                    "x": r.get("x"),
                    "lens": r.get("lens"),
                }
                for r in new
            ]
            + [
                {
                    "src": f"S1c saved {z['start']} (re-scored)",
                    "FF_pad64": z["FF_pad64"],
                    "x": z["x"],
                    "lens": z["lens"],
                }
                for z in s1g["saved_rescored"]
                if z["model"] == mdl
            ]
        )
        out[mdl] = {
            "seeds": seeds,
            "new_status": {str(r["seed"]): r["status"] for r in new},
            "n_unfinished": sum(1 for r in new if r["status"] != "finished"),
            "best": best_of(cands),
            "n_candidates": len(cands),
        }
    return out


def main_decide() -> dict:
    pool_ = pooled()
    ver = {}
    for mdl in ("P", "U"):
        bst = pool_[mdl]["best"]
        ver[mdl] = verify(bst["x"], bst["lens"], label=bst["src"])
    ff_case = {mdl: pool_[mdl]["best"]["FF_pad64"] for mdl in ("P", "U")}
    s1gv = json.loads((OUT / "s1g_validate.json").read_text())
    eps_p = {
        f"{r['name']} s{r['seed']}": r["eps_final"] for rnd in s1gv["rounds"] for r in rnd["rows"]
    }
    s1ev = json.loads((OUT / "s1e_validate.json").read_text())["rows"]
    eps_u = {f"{r['name']} s{r['seed']}": r["eps_final"] for r in s1ev if r["model"] == "U"}
    st_p, st_u = validation_statuses(eps_p, ff_case["P"]), validation_statuses(eps_u, ff_case["U"])
    tool_p = g.tool_verdict(st_p.values()) if len(st_p) == 6 else "unfinished"
    tool_u = g.tool_verdict(st_u.values()) if len(st_u) == 4 else "unfinished"
    best_mdl = max(("P", "U"), key=lambda k: ff_case[k])
    vb = ver[best_mdl]
    cex = any(v.get("counterexample_verified") for v in ver.values())
    feas = any(v.get("feasibility_fail_verified") for v in ver.values())
    outcome = s1h_outcome(
        tool_p,
        tool_u,
        list(pool_["P"]["seeds"].values()),
        list(pool_["U"]["seeds"].values()),
        ff_case[best_mdl],
        cex,
        feas,
        pool_["P"]["n_unfinished"] + pool_["U"]["n_unfinished"],
    )
    out = {
        "label": "derived; S1h pooled decision (rule functions in scripts/x4_s1h.py)",
        "per_model": {
            mdl: {
                "seeds": pool_[mdl]["seeds"],
                "top_two_gap": top_two_gap(pool_[mdl]["seeds"].values()),
                "seeds_agree_2e-4": g.seeds_agree(pool_[mdl]["seeds"].values(), SEED_AGREE),
                "new_status": pool_[mdl]["new_status"],
                "best_src": pool_[mdl]["best"]["src"],
                "FF_case": ff_case[mdl],
                "eps_limit": 0.1 * (1.0 - ff_case[mdl]),
            }
            for mdl in ("P", "U")
        },
        "validation_recomputed": {
            "P_S1g": st_p,
            "U_S1e": st_u,
            "eps_P": eps_p,
            "eps_U": eps_u,
            "tool_P": tool_p,
            "tool_U": tool_u,
        },
        "verify": ver,
        "best_overall_model": best_mdl,
        "best_overall_FF_conservative": vb.get("FF_conservative"),
        **outcome,
    }
    _dump("s1h_decision.json", out)
    return out


# ------------------------------------------------------------------ check before the fits
def main_check() -> dict:
    out: dict = {"label": "derived; S1h cache / version check before the fits", "checks": {}}
    ck = out["checks"]
    s1g = json.loads((OUT / "s1g_r1.json").read_text())
    _, t = m.target("D")
    r21 = next(r for r in s1g["rows"] if r["seed"] == 21)
    hp, hc = e.hpc(*e.to_cart(np.asarray(r21["x"])))
    lf = b.lens_exact(*r21["lens"], +1)
    ff = e.stat(t, hp * lf, hc * lf, pad=64)
    ck["s1g_R1_s21_exact"] = {
        "stored": r21["FF_pad64"],
        "now": ff,
        "abs_diff": abs(ff - r21["FF_pad64"]),
    }
    saved = j.saved_best()
    pairs = zip(saved, s1g["saved_rescored"], strict=True)
    diffs = [abs(a["FF_pad64"] - c["FF_pad64"]) for a, c in pairs]
    ck["saved_rescored_vs_s1g"] = {"n": len(saved), "abs_diff": max(diffs)}
    _, x0 = g._x0_saved_best_p(saved)
    gr = m.grid()
    # a single-point score is not a replay check: at x0 a 4e-16 parameter change (SciPy's internal
    # scaling) moves the objective by ~2e-8. So replay S1g R1 seed 21's DE for 1 generation instead.
    _, info = j.run_de(j._obj_p, j.BOUNDS_P, {"popsize": 12, "maxiter": 1}, 21, x0, (t, gr["F"]))
    stored = r21["log"][0]["de_best_every_10"][0]
    now = info["de_best_every_10"][0]
    ck["grid_de_replay_R1_s21_gen1"] = {"stored": stored, "now": now, "abs_diff": abs(now - stored)}
    ck["grid"] = {
        "nodes": len(gr["nodes"]),
        "cache_files_s+1": len(list(b.lens_cache_dir().glob("s+1_*.npy"))),
    }
    out["pass"] = bool(
        all(c["abs_diff"] <= 1e-12 for c in ck.values() if "abs_diff" in c)
        and ck["grid"]["nodes"] == ck["grid"]["cache_files_s+1"] == 1271
    )
    _dump("s1h_check.json", out)
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    what = argv[0] if argv else ""
    if what == "check":
        print(json.dumps(main_check(), indent=1, default=float))
    elif what == "run":
        main_run()
    elif what == "decide":
        print(json.dumps(main_decide(), indent=1, default=float))
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
