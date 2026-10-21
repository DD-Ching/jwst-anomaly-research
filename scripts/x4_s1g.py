"""X4 S1g: one bounded compute top-up of the S1e positive point-lens (P) search, then R1.

Unit file: docs/state/units/X4.md (S1g; plan frozen and committed before any fit). Labels:
injections and the R1 target are `simulated`; fitting factors are `derived`. No data is read.
Fits need lalsuite; the rule functions and the validation-injection generator run offline
(tests/test_x4_s1g.py).

Reuse, no copy of the procedure: the score, box, injections and `finish` come from `x4_s1d_de`; the
joint DE, the snapped lens grid, the continuous exact-lens refinement, the R1 case job and the saved
re-scoring come from `x4_s1e_joint`. Only two per-fit budgets change (`apply_budget`): DE_P maxiter
400 -> 1000 and the refinement round cap 10 -> 20. Two instrumentation hooks (`_install_hooks`) only
record state: the DE callback (S1e already passes one; the hook returns its value, so the DE is
unchanged) and the refinement's intrinsic restarts (to save each round's x, lens and value).

  python scripts/x4_s1g.py grid        # rebuild the S1c 1,271-node lens-grid cache (3 processes)
  python scripts/x4_s1g.py injections  # validation truths, formal-path check -> s1g_injections.json
  python scripts/x4_s1g.py version     # version check before the fits -> s1g_version_check.json
  python scripts/x4_s1g.py run         # round 1 -> rule -> round 2 -> rule -> R1 -> s1g_*.json
  python scripts/x4_s1g.py summary     # rule functions on the saved rows -> s1g_summary.json
"""

from __future__ import annotations

import os

# Thread limits must be set before numpy / scipy load their BLAS and FFT backends (one fit per
# process). setdefault keeps an explicit shell setting.
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import importlib.util  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from datetime import UTC, datetime  # noqa: E402
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


j = _load("x4_s1e_joint")
e, m, b = j.e, j.m, j.b
OUT = b.OUT

# ------------------------------------------------------------------ frozen settings (X4.md S1g)
S1E_BUDGET = {"DE_P_maxiter": 400, "MAX_ROUNDS": 10}  # as run in S1e (for the record)
DE_P_MAXITER = 1000
MAX_ROUNDS = 20
N_PROC = 3  # one fit per process, OMP/MKL/OPENBLAS threads = 1
DEADLINE_UTC = "2026-10-21T05:15:00+00:00"  # computation stops; running fits -> "unfinished"
EPS_MAX = 5e-4
EPS_CASE_FRAC = 0.1
SEED_AGREE = 2e-4
FF_SCREEN = j.FF_SCREEN  # rho_req > 30  <=>  FF > 0.994987
INTRINSIC_MAX = 100.0
SNR_RATIO_R1 = 0.19338258649590495  # lensed / unlensed SNR of the R1 family-D target (s1c_summary)
# intrinsic SNR needed > 100  <=>  rho_req > 100 * ratio  <=>  FF > 0.98790
FF_FEASIBLE = float(np.sqrt(1.0 - 9.0 / (INTRINSIC_MAX * SNR_RATIO_R1) ** 2))
ROUND_TRIP_TOL = 1e-9  # formal path: |from_cart(to_cart(x)) - x| per box coordinate
TRUTH_OBJ_TOL = 1e-6  # formal path: production objective at the truth (S1d rule)

VAL_SEED = 20261021
N_VAL = 2
VAL_REGION = {  # box coordinates; every range lies inside the S1d / S1e box
    "dln_mc": (-0.15, 0.15),  # ln M_c = ln M_c0 + U(.); M_c0 = R1's 12.74 Msun (detector frame)
    "q": (1.0, 5.0),  # uniform
    "chi_mag": (0.0, 0.99),  # uniform, each spin
    "cos_tilt": (-1.0, 1.0),  # uniform, each spin
    "azimuth": (0.0, 2 * np.pi),  # uniform, each spin
    "cos_inc": (-1.0, 1.0),  # inclination = arccos(U(.)), isotropic
    "phi": (0.0, 2 * np.pi),  # reference phase, uniform
    "psi": (0.0, np.pi),  # polarisation angle, uniform
    "ln_mlz": (np.log(300.0), np.log(3000.0)),  # M_Lz log-uniform, 300-3,000 Msun
    "ln_y": (np.log(0.1), np.log(1.0)),  # y log-uniform, 0.1-1
}
ROUNDS = (
    (("V9", 11), ("VA1", 11), ("VA2", 11)),  # round 1: S1e's failing case + both validation draws
    (("V5", 11), ("V5", 12), ("V9", 12)),  # round 2: the other development fits
)
LABEL = {"V5": "development", "V9": "development", "VA1": "validation", "VA2": "validation"}
R1_SEEDS = (21, 22, 23)


# ------------------------------------------------------------------ rule functions (offline)
def fit_status(eps, ff_case=None) -> str:
    """Pass rule for one injection/seed fit (S1e, unchanged).

    "pass" iff eps <= 5e-4 AND (when the R1 case has been fitted) eps <= 0.1 * (1 - FF_case);
    "unfinished" when the fit has no final score (stopped at the deadline); otherwise "fail"
    (a non-finite eps is a fail).
    """
    if eps is None:
        return "unfinished"
    eps = float(eps)
    if not np.isfinite(eps):
        return "fail"
    ok = eps <= EPS_MAX
    if ff_case is not None:
        ok = ok and eps <= EPS_CASE_FRAC * (1.0 - float(ff_case))
    return "pass" if ok else "fail"


def tool_verdict(statuses) -> str:
    """Whole tool: "pass" iff every fit passes; "fail" if any fails; else "error" if any fit raised;
    otherwise "unfinished"."""
    statuses = list(statuses)
    if not statuses:
        return "unfinished"
    if any(s == "fail" for s in statuses):
        return "fail"
    if all(s == "pass" for s in statuses):
        return "pass"
    if any(s == "error" for s in statuses):
        return "error"
    return "unfinished"


def row_status(row: dict, ff_case=None) -> str:
    """Status of one validation row: the pass rule if it finished; "error" if it raised;
    otherwise "unfinished" (stopped at the deadline)."""
    if row.get("status") == "finished":
        return fit_status(row.get("eps_final"), ff_case)
    if row.get("status") == "error":
        return "error"
    return "unfinished"


def screen_flags(ff_conservative, spread, legal: bool) -> dict:
    """S1f step 3 for one verified solution: FF minus the largest score disagreement, against the
    screen (outcome 1) and against the intrinsic-SNR bound (0.98790 < FF <= 0.994987)."""
    v = float(ff_conservative) - float(spread)
    return {
        "counterexample_verified": bool(legal and v > FF_SCREEN),
        "feasibility_fail_verified": bool(legal and FF_FEASIBLE < v <= FF_SCREEN),
    }


def is_legal(inside: bool, chi_mags, q, waveform_error) -> bool:
    """S1f step 1: inside the box, spins <= 0.99, q <= 15, the waveform generated without error."""
    return bool(inside and max(chi_mags) <= 0.99 + 1e-12 and q <= 15.0 and waveform_error is None)


def formal_path_ok(inside: bool, round_trip_max_diff, objective_at_truth) -> bool:
    """Validation-truth acceptance (frozen plan): inside the box, box round trip within
    ROUND_TRIP_TOL, production objective with the exact lens at the truth <= TRUTH_OBJ_TOL."""
    return bool(
        inside
        and float(round_trip_max_diff) <= ROUND_TRIP_TOL
        and float(objective_at_truth) <= TRUTH_OBJ_TOL
    )


def r1_seeds_to_run(rows) -> list:
    """Resume: R1 seeds without a finished row (finished rows are kept, never re-run)."""
    done = {r["seed"] for r in rows if r.get("status") == "finished"}
    return [s for s in R1_SEEDS if s not in done]


def seeds_agree(ffs, tol=SEED_AGREE) -> bool:
    """Search reliability: the two best seeds agree within tol in FF (needs >= 2 seeds)."""
    v = sorted((float(x) for x in ffs if x is not None), reverse=True)
    return len(v) >= 2 and v[0] - v[1] <= tol


def rho_req(ff: float) -> float:
    return float(3.0 / np.sqrt(max(1.0 - float(ff) ** 2, 1e-24)))


def r1_outcome(
    tool: str,
    p_ffs: list,
    u_ffs: list,
    best_ff: float,
    counterexample_verified: bool,
    p_unfinished: int = 0,
    ratio: float = SNR_RATIO_R1,
    feasibility_verified: bool = False,
) -> dict:
    """R1 outcome under the S1e outcomes and the S1f gate revision (frozen in X4.md S1g).

    1: a verified ordinary solution with FF > 0.994987 (rho_req > 30): R1 does not pass the screen,
       without proof of the tool.
    2: tool passes (with FF_case), both ordinary models reliable (two best seeds within 2e-4),
       rho_req <= 30 and intrinsic SNR <= 100.
    3: otherwise undetermined, with every cause named.
    Pre-registered clarification, checked right after outcome 1 and, like it, constructive: a
    verified solution with 0.98790 < FF <= 0.994987 (rho_req <= 30, intrinsic SNR > 100) is
    reported as failing the screen's intrinsic-SNR feasibility condition.
    """
    rho = rho_req(best_ff)
    intrinsic = rho / ratio
    base = {"best_ordinary_FF": best_ff, "rho_req": rho, "intrinsic_snr_needed": intrinsic}
    if counterexample_verified:
        return {"outcome": 1, "causes": [], **base}
    if feasibility_verified:
        cause = "verified ordinary solution needs intrinsic SNR > 100 (rho_req <= 30)"
        return {"outcome": "feasibility", "causes": [cause], **base}
    causes = []
    if tool != "pass":
        causes.append(f"P tool validation: {tool}")
    if p_unfinished:
        causes.append(f"R1 P fits unfinished at the stop: {p_unfinished}")
    if not seeds_agree(p_ffs):
        causes.append("P search reliability: two best R1 seeds differ by > 2e-4")
    if not seeds_agree(u_ffs):
        causes.append("U search reliability (S1f, not re-run): two best seeds differ by > 2e-4")
    if best_ff > FF_SCREEN:
        causes.append("FF above the screen threshold but not verified")
    elif intrinsic > INTRINSIC_MAX:
        causes.append("intrinsic SNR needed > 100 but not verified")
    if not causes:
        return {"outcome": 2, "causes": [], **base}
    return {"outcome": 3, "causes": causes, **base}


def in_box(x10, lens) -> bool:
    """Inside the S1d box (10 intrinsic coordinates) and the lens box (M_Lz 1-1e4, y 0.01-5)."""
    ok = all(lo <= v <= hi for v, (lo, hi) in zip(x10, e.BOUNDS, strict=True))
    if lens is not None:
        ok = ok and 1.0 <= lens[0] <= 1e4 and 0.01 <= lens[1] <= 5.0
    return bool(ok)


def draw_validation(k: int, seed: int = VAL_SEED) -> list[dict]:
    """The first k validation draws (deterministic): box coordinates, psi and lens.

    Draw order per injection: dln_mc, q, (chi_mag, cos_tilt, azimuth) x 2, cos_inc, phi, psi,
    ln_mlz, ln_y. `injections` keeps the first N_VAL draws that pass the formal-path check.
    """
    rng = np.random.default_rng(seed)
    r = VAL_REGION
    out = []
    for i in range(k):
        x = [np.log(m.MC0) + rng.uniform(*r["dln_mc"]), rng.uniform(*r["q"])]
        for _ in range(2):
            x += [
                rng.uniform(*r["chi_mag"]),
                rng.uniform(*r["cos_tilt"]),
                rng.uniform(*r["azimuth"]),
            ]
        x += [float(np.arccos(rng.uniform(*r["cos_inc"]))), rng.uniform(*r["phi"])]
        psi = rng.uniform(*r["psi"])
        lens = [float(np.exp(rng.uniform(*r["ln_mlz"]))), float(np.exp(rng.uniform(*r["ln_y"])))]
        out.append({"draw": i, "x_box": [float(v) for v in x], "psi": float(psi), "lens": lens})
    return out


def spec_from_draw(z: dict):
    """(physical intrinsics + psi, lens) in the form `x4_s1d_de.injection_target` takes."""
    mc, q, s1, s2, inc, phi = e.to_cart(np.asarray(z["x_box"]))
    return ((mc, q, s1, s2, inc, phi, z["psi"]), tuple(z["lens"]))


def apply_budget() -> None:
    """The only change to the S1e procedure: per-fit DE generations and refinement round cap."""
    j.DE_P = {**j.DE_P, "maxiter": DE_P_MAXITER}
    j.MAX_ROUNDS = MAX_ROUNDS


# ------------------------------------------------------------------ instrumentation (state only)
_RUN: dict = {}
_ORIG: dict = {}


def _state_path(key: str) -> Path:
    env = os.environ.get("X4_S1G_STATE")
    # default: the data root (untracked), never results/
    d = Path(env) if env else b.lens_cache_dir().parent / "s1g_state"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{key}.json"


def _write_state(rec: dict) -> None:
    if "key" not in _RUN:
        return
    p = _state_path(_RUN["key"])
    tmp = p.with_suffix(f".{os.getpid()}.tmp")
    rec = {**rec, "run_id": _RUN.get("run_id"), "utc": _now()}
    tmp.write_text(json.dumps(rec, default=float) + "\n")
    os.replace(tmp, p)


def read_state(key: str, run_id: str):
    """The saved state of one fit, or None when absent or written by another run (stale)."""
    p = _state_path(key)
    if not p.exists():
        return None
    st = json.loads(p.read_text())
    return st if st.get("run_id") == run_id else None


def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _run_id() -> str:
    return f"{_now()}_{os.getpid()}"


def refine_state(p: dict, x, f_in: float, lens, rnd: int) -> dict:
    """Refinement state after one intrinsic restart: refine accepts the restart iff it lowers
    1 - FF below the value it started from (the test in `x4_s1e_joint.refine`)."""
    accepted = float(p["f"]) < float(f_in)
    return {
        "stage": "refine",
        "round": rnd,
        "accepted": accepted,
        "f": float(p["f"]) if accepted else float(f_in),
        "x": [float(v) for v in (p["x"] if accepted else x)],
        "lens": lens,
    }


def _install_hooks() -> None:
    if _ORIG:
        return
    _ORIG.update({"de": j.differential_evolution, "nm": j.nm, "lens_exact": b.lens_exact})

    def de_hooked(*a, callback=None, **kw):
        _RUN["lens_calls_prev"] = _RUN.get("lens_calls", 0)  # refinement counts start after the DE

        def cb(intermediate_result):  # SciPy dispatches on this parameter name
            res = intermediate_result
            nit = int(res.nit)
            if (nit - 1) % 10 == 0:
                _RUN.setdefault("spread_every_10", []).append(
                    float(np.std(res.population_energies))
                )
                _write_state(
                    {"stage": "de", "nit": nit, "f": float(res.fun), "x": [float(v) for v in res.x]}
                )
            return callback(res) if callback is not None else None

        return _ORIG["de"](*a, callback=cb, **kw)

    def lens_hooked(m_lz, y, s):
        _RUN["last_lens"] = [float(m_lz), float(y)]
        _RUN["lens_calls"] = _RUN.get("lens_calls", 0) + 1
        return _ORIG["lens_exact"](m_lz, y, s)

    def nm_hooked(t, lensf, x):
        f_in = float(e.objective(x, t, lensf))  # refine's current value (pure; search unchanged)
        p = _ORIG["nm"](t, lensf, x)
        rr = _RUN.setdefault("refine_hook", [])
        calls = _RUN.get("lens_calls", 0)
        rr.append({"nm_nfev": int(p["nfev"]), "lens_evals": calls - _RUN.get("lens_calls_prev", 0)})
        _RUN["lens_calls_prev"] = calls
        _write_state(refine_state(p, x, f_in, _RUN.get("last_lens"), len(rr) - 1))
        return p

    j.differential_evolution = de_hooked
    j.nm = nm_hooked
    b.lens_exact = lens_hooked


def _init_worker() -> None:
    # thread limits are set at import (module top), before numpy / scipy load
    apply_budget()
    _install_hooks()


# ------------------------------------------------------------------ fits
def _key(kind: str, name: str, seed: int) -> str:
    return f"{kind}_{name}_s{seed}"


def specs() -> dict:
    """Development (S1e V5, V9) and validation (committed VA draws) injection specs. Resolved once
    in the parent (`main_run`) and passed to every job, so one run cannot mix truths."""
    base = e.injections()
    out = {n: base[n] for n in ("V5", "V9")}
    rec = json.loads((OUT / "s1g_injections.json").read_text())
    for z in rec["accepted"]:
        out[z["name"]] = spec_from_draw(z)
    return out


def s1e_trace(name: str, seed: int):
    for r in json.loads((OUT / "s1e_validate.json").read_text())["rows"]:
        if r["name"] == name and r["seed"] == seed:
            return r["log"][0]["de_best_every_10"]
    return None


def _job(args):
    """One fit. args = (kind, name, seed, payload, run_id); payload is the injection spec
    ("valid") or the R1 initial-population member x0 ("r1")."""
    kind, name, seed, payload, run_id = args
    _RUN.clear()
    _RUN.update({"key": _key(kind, name, seed), "run_id": run_id})
    t0 = time.time()
    if kind == "valid":
        t, truth_mm = e.injection_target(payload)
        r = e.finish(t, j.fit_p(t, seed), "P")
        r.update({"name": name, "seed": seed, "set": LABEL[name], "truth_mismatch": truth_mm})
        r["eps_final"] = 1.0 - r["FF_pad64"]
        old = s1e_trace(name, seed)
        if old is not None:  # same seed: the first 400 generations should replay S1e
            new = r["log"][0]["de_best_every_10"]
            n = min(len(new), len(old) - 1, 40)  # S1e's last entry is generation 400 itself
            diff = [abs(a - c) for a, c in zip(new[:n], old[:n], strict=True)]
            r["s1e_replay"] = {
                "compared_entries": n,
                "max_abs_diff": max(diff) if diff else None,
                "first_index_diff_gt_1e-12": next(
                    (i for i, v in enumerate(diff) if v > 1e-12), None
                ),
            }
    else:
        r = j._case_job(("P", seed, payload))
        r.update({"name": "R1", "set": "R1 case"})
    info = r["log"][0]
    info["de_spread_every_10"] = _RUN.get("spread_every_10", [])
    hook = _RUN.get("refine_hook", [])
    for rr, hh in zip(info.get("refine_rounds", []), hook, strict=False):
        rr.update(hh)
    info["evaluations"] = {
        "de": info["de_nfev"],
        "refine_intrinsic": int(sum(h["nm_nfev"] for h in hook)),
        "refine_lens": int(sum(h["lens_evals"] for h in hook)),
    }
    info["refine_n_rounds"] = len(info.get("refine_rounds", []))
    r.update({"kind": kind, "status": "finished", "wall_s": round(time.time() - t0, 1)})
    r.update({"run_id": run_id, "budget": {"DE_P_maxiter": j.DE_P["maxiter"]}})
    r["budget"]["MAX_ROUNDS"] = j.MAX_ROUNDS
    return r


def _base_row(job, status: str) -> dict:
    kind, name, seed, _payload, run_id = job
    row = {"kind": kind, "name": name, "seed": seed, "status": status, "run_id": run_id}
    row["set"] = LABEL[name] if kind == "valid" else "R1 case"
    return row


def collect(job, ar) -> dict:
    """The finished row of one job, or an "error" row if the fit raised (the round goes on)."""
    try:
        return ar.get()
    except Exception as exc:  # noqa: BLE001
        return {**_base_row(job, "error"), "error": repr(exc)}


def _partial_row(job) -> dict:
    """A fit stopped at the deadline: its last saved state of this run, re-scored exactly."""
    kind, name, seed, payload, run_id = job
    row = _base_row(job, "unfinished")
    st = read_state(_key(kind, name, seed), run_id)
    if st is None:
        return row
    x = np.asarray(st["x"][:10])
    if st["stage"] == "de":
        node = m.grid()["nodes"][j.node_index(st["x"][10], st["x"][11])]
        lens = [float(node[0]), float(node[1])]
    else:
        lens = st["lens"]
    t = e.injection_target(payload)[0] if kind == "valid" else m.target("D")[1]
    r = e.finish(t, {"f": float(st["f"]), "x": [float(v) for v in x], "lens": lens}, "P")
    row.update({"partial_state": st, "x": r["x"], "lens": lens, "FF_pad64": r["FF_pad64"]})
    row.update({"bounds": r["bounds"], "physical": r["physical"]})
    return row


def _deadline() -> float:
    return datetime.fromisoformat(DEADLINE_UTC).timestamp()


def run_round(jobs, on_row) -> list:
    rows = []
    for jb in jobs:  # no stale state: this run's job keys start empty
        _state_path(_key(*jb[:3])).unlink(missing_ok=True)
    with Pool(N_PROC, initializer=_init_worker) as pool:
        pending = {i: (jb, pool.apply_async(_job, (jb,))) for i, jb in enumerate(jobs)}
        while pending:
            for i, (jb, ar) in list(pending.items()):
                if ar.ready():
                    rows.append(collect(jb, ar))
                    del pending[i]
                    on_row(rows[-1])
            if pending and time.time() >= _deadline():
                pool.terminate()
                break
            if pending:
                time.sleep(15)
    for jb, _ar in pending.values():
        rows.append(_partial_row(jb))
        on_row(rows[-1])
    return rows


def _dump(name: str, rec: dict) -> None:
    (OUT / name).write_text(json.dumps(rec, indent=1, default=float) + "\n")


def main_run() -> None:
    apply_budget()
    m.grid()  # load the cached grid once; workers inherit it
    run_id = _run_id()
    sp = specs()  # resolved once: every fit of this run uses the same truths
    rec = {
        "label": "simulated injections; derived FF; S1g P-search validation "
        "(S1e procedure, budget top-up)",
        "budget": {"DE_P_maxiter": DE_P_MAXITER, "MAX_ROUNDS": MAX_ROUNDS, "S1e": S1E_BUDGET},
        "deadline_utc": DEADLINE_UTC,
        "run_id": run_id,
        "rounds": [],
    }
    for i, rnd in enumerate(ROUNDS, start=1):
        cur = {"round": i, "start_utc": _now(), "rows": []}
        rec["rounds"].append(cur)

        def on_row(r, cur=cur):
            cur["rows"].append(r)
            _dump("s1g_validate.json", rec)
            print(
                _now(), r.get("name"), r.get("seed"), r.get("status"), r.get("FF_pad64"), flush=True
            )

        run_round([("valid", n, s, sp[n], run_id) for n, s in rnd], on_row)
        cur["end_utc"] = _now()
        cur["statuses"] = [row_status(r) for r in cur["rows"]]
        cur["verdict_eps_only"] = tool_verdict(cur["statuses"])
        _dump("s1g_validate.json", rec)
        print(_now(), f"round {i} verdict", cur["verdict_eps_only"], flush=True)
        if cur["verdict_eps_only"] != "pass":
            print(
                _now(), "stop: the tool did not pass; no further rounds, R1 not fitted", flush=True
            )
            print(json.dumps(summary(), indent=1, default=float))
            return
    main_r1(run_id)


def _x0_saved_best_p(saved):
    best = max((z for z in saved if z["model"] == "P"), key=lambda z: z["FF_pad64"])
    x0 = list(best["x"]) + [float(np.log(best["lens"][0])), float(np.log(best["lens"][1]))]
    return best, np.asarray(x0)


def main_r1(run_id: str) -> None:
    path = OUT / "s1g_r1.json"
    if path.exists():  # resume: keep finished rows; run only the missing or unfinished seeds
        rec = json.loads(path.read_text())
        rec["rows"] = [r for r in rec["rows"] if r.get("status") == "finished"]
        rec.setdefault("resumed", []).append({"utc": _now(), "run_id": run_id})
        saved = rec["saved_rescored"]
        _, x0 = _x0_saved_best_p(saved)
    else:
        saved = j.saved_best()
        best, x0 = _x0_saved_best_p(saved)
        rec = {
            "label": "simulated R1 family-D target; derived FF; S1g R1 P fits (S1e case procedure)",
            "budget": {"DE_P_maxiter": DE_P_MAXITER, "MAX_ROUNDS": MAX_ROUNDS},
            "start_utc": _now(),
            "run_id": run_id,
            "x0_saved_best_P": {"start": best["start"], "FF_pad64": best["FF_pad64"]},
            "saved_rescored": saved,
            "rows": [],
        }

    def on_row(r):
        rec["rows"].append(r)
        _dump("s1g_r1.json", rec)
        print(_now(), "R1", r.get("seed"), r.get("status"), r.get("FF_pad64"), flush=True)

    todo = r1_seeds_to_run(rec["rows"])
    if todo:
        run_round([("r1", "R1", s, x0, run_id) for s in todo], on_row)
    rec["end_utc"] = _now()
    cands = [r for r in rec["rows"] if r.get("FF_pad64") is not None]
    cands += [
        {"FF_pad64": z["FF_pad64"], "x": z["x"], "lens": z["lens"], "seed": z["start"]}
        for z in saved
    ]
    top = max((c for c in cands if c.get("lens") is not None), key=lambda c: c["FF_pad64"])
    rec["verify_best_P"] = verify_p(top["x"], top["lens"], label=str(top.get("seed")))
    _dump("s1g_r1.json", rec)
    print(json.dumps(summary(), indent=1, default=float))


# ------------------------------------------------------------------ verification (S1f 3 steps, P)
def verify_p(x, lens, label="") -> dict:
    """S1f's three-step check, for a lensed solution: legality; an explicit-psi rebuild scored with
    the S1c real-part overlap at pad 64 and 256; FF minus the largest disagreement vs the screen."""
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
        "lens": [float(v) for v in lens],
        "in_box": in_box(x, lens),
        "chi_mags": [float(np.linalg.norm(s1)), float(np.linalg.norm(s2))],
        "q": q,
        "waveform_error": err,
    }
    if err is not None:
        out.update({"counterexample_verified": False, "feasibility_fail_verified": False})
        return out
    lf = b.lens_exact(float(lens[0]), float(lens[1]), +1)
    hp, hc = hp * lf, hc * lf
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
    legal = is_legal(out["in_box"], out["chi_mags"], q, err)
    rho = rho_req(ff_cons)
    out.update(
        {
            "psi": psi,
            "scores": scores,
            "score_spread": spread,
            "FF_conservative": ff_cons,
            "FF_screen": FF_SCREEN,
            "FF_feasible": FF_FEASIBLE,
            **screen_flags(ff_cons, spread, legal),
            "rho_req_observed": rho,
            "intrinsic_snr_needed": rho / SNR_RATIO_R1,
        }
    )
    return out


# ------------------------------------------------------------------ summary
def summary() -> dict:
    v = json.loads((OUT / "s1g_validate.json").read_text())
    r1p = OUT / "s1g_r1.json"
    r1 = json.loads(r1p.read_text()) if r1p.exists() else None
    u = json.loads((OUT / "s1f_case_u.json").read_text())
    u_ffs = [r["FF_pad64"] for r in u["rows"]]
    ff_case = None
    p_ffs, p_unf, best_p = [], 0, None
    if r1 is not None:
        p_ffs = [r["FF_pad64"] for r in r1["rows"] if r["status"] == "finished"]
        p_unf = sum(1 for r in r1["rows"] if r["status"] != "finished")  # unfinished or error
        allp = [r.get("FF_pad64") for r in r1["rows"]] + [
            z["FF_pad64"] for z in r1["saved_rescored"] if z["model"] == "P"
        ]
        best_p = max(x for x in allp if x is not None)
        ff_case = best_p
    fits = []
    for rnd in v["rounds"]:
        for r in rnd["rows"]:
            eps = r.get("eps_final") if r["status"] == "finished" else None
            info = r.get("log", [{}])[0]
            rounds = info.get("refine_rounds", [])
            fits.append(
                {
                    "round": rnd["round"],
                    "name": r["name"],
                    "seed": r["seed"],
                    "set": r["set"],
                    "status_run": r["status"],
                    "eps": eps,
                    "eps_partial": None
                    if eps is not None or "FF_pad64" not in r
                    else 1 - r["FF_pad64"],
                    "rule_eps_only": row_status(r),
                    "rule_with_ff_case": row_status(r, ff_case),
                    "de_message": info.get("de_message"),
                    "de_nit": info.get("de_nit"),
                    "de_spread": info.get("de_spread"),
                    "refine_rounds": len(rounds),
                    "refine_last_gain": rounds[-1]["gain"] if rounds else None,
                    "refine_converged": info.get("refine_converged"),
                    "bounds": r.get("bounds"),
                    "lens": r.get("lens"),
                    "wall_s": r.get("wall_s"),
                }
            )
    ran_all = len(fits) == sum(len(rnd) for rnd in ROUNDS)

    def guarded(v):  # a set that did not run completely can only fail, never pass
        return v if ran_all or v == "fail" else "unfinished"

    out = {
        "label": "derived; S1g summary (rule functions applied to the saved rows)",
        "fits": fits,
        "tool_eps_only": guarded(tool_verdict(f["rule_eps_only"] for f in fits)),
        "FF_case_P": ff_case,
        "tool_final": None,
        "r1": None,
    }
    if r1 is not None:
        out["tool_final"] = guarded(tool_verdict(f["rule_with_ff_case"] for f in fits))
        ver = r1.get("verify_best_P", {})
        if "FF_conservative" in ver:  # step 3 recomputed from the stored verification
            legal = is_legal(ver["in_box"], ver["chi_mags"], ver["q"], ver.get("waveform_error"))
            flags = screen_flags(ver["FF_conservative"], ver["score_spread"], legal)
        else:
            flags = {"counterexample_verified": False, "feasibility_fail_verified": False}
        best_any = max(best_p, max(u_ffs))
        out["r1"] = {
            "P_FF_seeds": {str(r["seed"]): r.get("FF_pad64") for r in r1["rows"]},
            "P_status": {str(r["seed"]): r["status"] for r in r1["rows"]},
            "P_best_incl_saved": best_p,
            "U_FF_seeds_S1f": u_ffs,
            "P_seeds_agree": seeds_agree(p_ffs),
            "U_seeds_agree": seeds_agree(u_ffs),
            "verify_best_P": {
                **{k: ver.get(k) for k in ("label", "FF_conservative", "score_spread")},
                **flags,
            },
            **r1_outcome(
                out["tool_final"],
                p_ffs,
                u_ffs,
                best_any,
                flags["counterexample_verified"],
                p_unf,
                feasibility_verified=flags["feasibility_fail_verified"],
            ),
        }
    _dump("s1g_summary.json", out)
    return out


# ------------------------------------------------------------------ setup commands
def _grid_node(node):
    t0 = time.time()
    b.lens(float(node[0]), float(node[1]), +1)
    return time.time() - t0


def main_grid() -> None:
    nodes = [(mm, y) for mm in m.GRID_M for y in m.GRID_Y]
    t0, cpu = time.time(), 0.0
    with Pool(N_PROC) as pool:
        for dt in pool.imap_unordered(_grid_node, nodes):
            cpu += dt
    print(
        f"{len(nodes)} nodes in {b.lens_cache_dir()}:",
        f"wall {time.time() - t0:.0f} s, cpu {cpu:.0f} s",
    )


def main_injections() -> dict:
    """First N_VAL draws whose truth passes the formal path (box, round trip, objective ~ 0)."""
    acc, rej = [], []
    for z in draw_validation(20):
        spec = spec_from_draw(z)
        rec = dict(z)
        try:
            t, mm = e.injection_target(spec)
            (mc, q, s1, s2, inc, phi, _psi), lens = spec
            xr = e.from_cart(mc, q, s1, s2, inc, phi)
            lf = b.lens_exact(*lens, +1)
            rec.update(
                {
                    "physical": {"mc": mc, "q": q, "chi1": list(s1), "chi2": list(s2), "inc": inc},
                    "in_box": in_box(xr, lens),
                    "round_trip_max_diff": float(np.max(np.abs(xr - np.asarray(z["x_box"])))),
                    "objective_at_truth": float(e.objective(xr, t, lf)),
                    "truth_mismatch_stat": float(mm),
                }
            )
            ok = formal_path_ok(
                rec["in_box"], rec["round_trip_max_diff"], rec["objective_at_truth"]
            )
        except Exception as exc:  # noqa: BLE001
            rec["error"], ok = repr(exc), False
        (acc if ok else rej).append(rec)
        if len(acc) == N_VAL:
            break
    for i, z in enumerate(acc, start=1):
        z["name"] = f"VA{i}"
    out = {
        "label": "simulated; S1g validation injections (generator frozen in X4.md S1g); truths are "
        "never used as starts",
        "rng_seed": VAL_SEED,
        "accept": {"round_trip_tol": ROUND_TRIP_TOL, "objective_at_truth_tol": TRUTH_OBJ_TOL},
        "region": {k: [float(a), float(c)] for k, (a, c) in VAL_REGION.items()},
        "accepted": acc,
        "rejected": rej,
    }
    _dump("s1g_injections.json", out)
    return out


def main_version() -> dict:
    """Exact re-scores of saved solutions, truth checks and a short DE replay against S1e."""
    out: dict = {"label": "derived; S1g version check before the fits", "checks": {}}
    import lal
    import lalsimulation
    import mpmath
    import scipy

    out["versions"] = {
        "lal": lal.__version__,
        "lalsimulation": lalsimulation.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
    }
    ck = out["checks"]
    _, t = m.target("D")
    u = json.loads((OUT / "s1f_case_u.json").read_text())
    r23 = next(r for r in u["rows"] if r["seed"] == 23)
    ff = e.stat(t, *e.hpc(*e.to_cart(np.asarray(r23["x"]))), pad=64)
    ck["s1f_seed23_U"] = {
        "stored": r23["FF_pad64"],
        "now": ff,
        "abs_diff": abs(ff - r23["FF_pad64"]),
    }
    sp = max((z for z in u["saved_rescored"] if z["model"] == "P"), key=lambda z: z["FF_pad64"])
    hp, hc = e.hpc(*e.to_cart(np.asarray(sp["x"])))
    lf = b.lens_exact(*sp["lens"], +1)
    ff = e.stat(t, hp * lf, hc * lf, pad=64)
    ck["s1c_best_P_rescored"] = {
        "stored": sp["FF_pad64"],
        "now": ff,
        "abs_diff": abs(ff - sp["FF_pad64"]),
    }
    inj = e.injections()
    for r in json.loads((OUT / "s1e_validate.json").read_text())["rows"]:
        if r["model"] != "P":
            continue
        tt, _ = e.injection_target(inj[r["name"]])
        hp, hc = e.hpc(*e.to_cart(np.asarray(r["x"])))
        lf = b.lens_exact(*r["lens"], +1)
        eps = 1.0 - e.stat(tt, hp * lf, hc * lf, pad=64)
        key = f"s1e_{r['name']}_s{r['seed']}_final"
        ck[key] = {"stored": r["eps_final"], "now": eps, "abs_diff": abs(eps - r["eps_final"])}
    for name in ("V5", "V9"):
        (mc, q, s1, s2, inc, phi, _psi), lens = inj[name]
        tt, _ = e.injection_target(inj[name])
        xr = e.from_cart(mc, q, s1, s2, inc, phi)
        eps = float(e.objective(xr, tt, b.lens_exact(*lens, +1)))
        ck[f"truth_{name}"] = {"eps": eps, "in_box": in_box(xr, lens)}
    g = m.grid()
    ck["grid"] = {
        "nodes": len(g["nodes"]),
        "cache_files": len(list(b.lens_cache_dir().glob("s+1_*.npy"))),
    }
    tt, _ = e.injection_target(inj["V9"])
    t0 = time.time()
    _, info = j.run_de(j._obj_p, j.BOUNDS_P, {"popsize": 12, "maxiter": 11}, 11, None, (tt, g["F"]))
    old = s1e_trace("V9", 11)[:2]
    new = info["de_best_every_10"][:2]
    ck["de_replay_V9_s11_gen1_gen11"] = {
        "stored": old,
        "now": new,
        "abs_diff": max(abs(a - c) for a, c in zip(old, new, strict=True)),
        "seconds": round(time.time() - t0, 1),
    }
    out["pass"] = bool(
        all(c["abs_diff"] <= 1e-9 for c in ck.values() if "abs_diff" in c)
        and all(
            ck[f"truth_{n}"]["eps"] <= 1e-9 and ck[f"truth_{n}"]["in_box"] for n in ("V5", "V9")
        )
        and ck["grid"]["nodes"] == ck["grid"]["cache_files"] == 1271
    )
    _dump("s1g_version_check.json", out)
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    what = argv[0] if argv else ""
    if what == "grid":
        main_grid()
    elif what == "injections":
        print(json.dumps(main_injections(), indent=1, default=float))
    elif what == "version":
        print(json.dumps(main_version(), indent=1, default=float))
    elif what == "run":
        main_run()
    elif what == "r1":  # resume R1 only, and only after both validation rounds passed
        if summary()["tool_eps_only"] != "pass":
            print("refused: the validation rounds have not passed")
            return 1
        apply_budget()
        m.grid()
        main_r1(_run_id())
    elif what == "summary":
        print(json.dumps(summary(), indent=1, default=float))
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
