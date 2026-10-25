"""W3 / Kepler real-background false-pass test and small search (follow-up to #500).

Plan: ``results/w3_kepler/hard/PLAN.md`` (frozen before any comparison was read and before any
new flux was downloaded). Reuses ``scripts/w3_kepler.py`` unchanged except for the frozen local
shoulder rule (``shoulder_local``, ``screen_star(local=True)``).

Stages (``python scripts/w3_kepler_hard.py <stage>``):

- ``sample``: groups old_bg (seen #500 stars), eb_deep, eb_single, val (25) and search (200);
- ``fetch --group G``: DR25 long-cadence files through the download ledger (3 GB cap);
- ``phaseA``: hard tests on seen backgrounds and real EB light curves; old vs local rule decision;
- ``validate``: the chosen rule on the unread val stars (``gate.json``);
- ``search``: only if ``gate.json`` passes.

Labels: injections simulated; light curves and flags observed; outcomes derived.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w3_kepler as wk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "hard"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "cache" / "w3_kepler_hard"
CAP = 3_000_000_000  # owner limit for this round, bytes
SALT2 = "w3-kepler-hard-2026-10-25"
OLD_BRIGHT = (9137992, 2854411, 2856758, 3834001)  # ../posthoc2.json real_shoulder_like
EB_DEPTH_MIN = 600_000  # ppm: DR25 TCE depth >= 0.6 (0.995 mag at the threshold)
N_VAL, N_SEARCH = 25, 200
N_INJ = 100
N_BRIGHT_PER_STAR = 25
HARD = ("long_dip", "single_point", "quarter_edge", "synth_bright_eclipse")
FAMILIES = tuple(wk.CONTROLS) + HARD
SEEDS = {"A": (5100, 5200, 5290), "val": (5300, 5400)}
GATE_REC, GATE_CTRL = 0.50, 0.10
MIN_BRIGHT_SOURCES, MIN_EB_STARS = 3, 5


def configure() -> None:
    """Point the reused pilot functions at this round's results, cache and ledger cap."""
    wk.OUT, wk.CACHE, wk.LEDGER_CAP = OUT, CACHE, CAP


def usage(stage: str, groups) -> None:
    rec = {"t": round(time.time(), 1), "stage": stage, "groups": sorted(groups)}
    with open(OUT / "source_usage.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


# ------------------------------------------------------------------------------------------ sample


def stage_sample(a) -> None:
    raw = Path(a.pool).read_bytes()
    if hashlib.sha256(raw).hexdigest() != wk.POOL_SHA256:
        raise SystemExit("pool file differs from the pinned sha256")
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    pool = [r for r in rows if r["st_quarters"] == "1" * 17]

    def rank(salt):
        return sorted(
            pool, key=lambda r: hashlib.sha256(f"{salt}:{int(r['kepid'])}".encode()).hexdigest()
        )

    seen = {int(r["kepid"]) for r in rank(wk.SALT)[:101]}  # #500 sample + interface star
    old = json.loads((OUT.parent / "sample.json").read_text())["stars"]
    old_dev = [s["kepid"] for s in old if s["group"] == "dev"]
    old_bg = old_dev + [k for k in OLD_BRIGHT if k not in old_dev]
    tce = list(csv.DictReader(open(OUT / "inputs" / "tce_depth_ge_0.4.csv")))
    eb_deep = sorted({int(r["kepid"]) for r in tce if float(r["tce_depth"]) >= EB_DEPTH_MIN})
    t9 = list(csv.DictReader(open(OUT / "inputs" / "kirk2016_table9.csv")))
    eb_single = sorted({int(r["KIC"]) for r in t9} - set(eb_deep))
    excl = seen | set(eb_deep) | set(eb_single)
    new = [int(r["kepid"]) for r in rank(SALT2) if int(r["kepid"]) not in excl]
    groups = {
        "old_bg": old_bg,
        "eb_deep": eb_deep,
        "eb_single": eb_single,
        "val": new[:N_VAL],
        "search": new[N_VAL : N_VAL + N_SEARCH],
    }
    allk = [k for g in ("val", "search") for k in groups[g]]
    assert len(set(allk)) == len(allk) and not set(allk) & excl
    out = {
        "label": "observed (metadata)",
        "pool_sha256": wk.POOL_SHA256,
        "inputs_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((OUT / "inputs").glob("*.csv"))
        },
        "rule": f"val / search: pool minus the #500 101 and the EB stars, sorted by sha256("
        f"'{SALT2}:<kepid>'); first {N_VAL} val, next {N_SEARCH} search",
        "groups": groups,
    }
    (OUT / "sample_hard.json").write_text(json.dumps(out, indent=1) + "\n")
    print({g: len(v) for g, v in groups.items()})


def group(g: str) -> list[int]:
    return json.loads((OUT / "sample_hard.json").read_text())["groups"][g]


# ------------------------------------------------------------------------------------------- fetch


def _fetch_one(kepid: int) -> dict:
    configure()
    if (CACHE / f"{kepid:09d}.npz").exists():
        return {"kepid": kepid, "cached": True}
    try:
        return wk._fetch_star(kepid)
    except Exception as e:  # recorded, never silently dropped
        return {"kepid": kepid, "error": f"{type(e).__name__}: {e}"}


def stage_fetch(a) -> None:
    usage(f"fetch {a.group}", [a.group])
    with ThreadPoolExecutor(8) as ex:
        info = list(ex.map(_fetch_one, group(a.group)))
    (OUT / f"files_{a.group}.json").write_text(
        json.dumps({"label": "observed", "stars": info}) + "\n"
    )
    used = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP).used()
    print(a.group, len(info), "stars;", sum("error" in x for x in info), "errors; ledger", used)


def have(kepid: int) -> bool:
    return (CACHE / f"{kepid:09d}.npz").exists()


# ------------------------------------------------------------------------------------------ models


def _logu(rng, lo, hi):
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def bright_times(st) -> np.ndarray:
    """Times of uninjected good cadences meeting the #500 shoulder criterion (real brightenings)."""
    g = wk.good_mask(st)
    rel, sig = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], g)
    ok = np.isfinite(rel)
    hit = ok & (rel >= 1 + wk.SHOULDER_MIN) & ((rel - 1) / sig >= wk.SHOULDER_SNR)
    return st["TIME"][hit]


def hard_flux(kind, st, good, rng, bright=None):
    """Flux factor (exposure-averaged where it is a model) and common-mode override of one hard
    control (plan section 2); the #500 controls go to w3_kepler.control_flux unchanged."""
    if kind in wk.CONTROLS:
        return wk.control_flux(kind, st, good, rng)
    t = st["TIME"]
    fin = np.isfinite(t)
    tg = t[good]
    f = np.ones(len(t))
    if kind == "real_bright_eclipse":
        tc = float(rng.choice(bright)) + rng.uniform(-0.5, 0.5)
        dur = _logu(rng, 0.1, 10.0)
        f[fin] = wk.dmag_flux_avg(wk._box(tc - dur / 2, tc + dur / 2), t[fin])
    elif kind == "long_dip":
        dur, t1 = _logu(rng, 0.1, 40.0), rng.uniform(tg.min(), tg.max())
        f[fin] = wk.dmag_flux_avg(wk._box(t1, t1 + dur), t[fin])
    elif kind == "single_point":
        dur, t1 = rng.uniform(0.1, 2.0), rng.uniform(tg.min(), tg.max())
        f[fin] = wk.dmag_flux_avg(wk._box(t1, t1 + dur), t[fin])
        w = wz_frac() * max(dur, 1.0)
        for m in (good & (t >= t1 - w) & (t < t1), good & (t > t1 + dur) & (t <= t1 + dur + w)):
            idx = np.flatnonzero(m)
            if len(idx):
                f[idx[rng.integers(len(idx))]] *= 1.2  # one unflagged +20 % cadence per side
    elif kind == "quarter_edge":
        qs = np.unique(st["QUARTER"][good])
        k = int(rng.integers(len(qs) - 1))  # boundary between consecutive quarters with data
        q_end = t[(st["QUARTER"] == qs[k]) & fin].max()
        q_start = t[(st["QUARTER"] == qs[k + 1]) & fin].min()
        v, dur = int(rng.integers(3)), rng.uniform(0.1, 2.0)
        if v == 0:
            t2 = q_end - rng.uniform(0, 5) * wk.CAD
            t1 = t2 - dur
        elif v == 1:
            t1 = q_start + rng.uniform(0, 5) * wk.CAD
            t2 = t1 + dur
        else:
            t1, t2 = q_end - rng.uniform(0.05, 1.0), q_start + rng.uniform(0.05, 1.0)
        f[fin] = wk.dmag_flux_avg(wk._box(t1, t2), t[fin])
    elif kind == "synth_bright_eclipse":
        amp, ramp, plat = rng.uniform(0.10, 0.20), rng.uniform(0.5, 3.0), rng.uniform(2.0, 10.0)
        t0 = rng.uniform(tg.min(), tg.max())
        dur = _logu(rng, 0.1, 2.0)
        tc = rng.uniform(t0 + ramp + dur / 2, t0 + ramp + plat - dur / 2)

        def model(x):
            x = np.asarray(x)
            up = np.clip((x - t0) / ramp, 0, 1)
            down = np.clip((t0 + 2 * ramp + plat - x) / ramp, 0, 1)
            br = 1 + amp * np.minimum(up, down)
            ecl = np.where(np.abs(x - tc) < dur / 2, 1.5, 0.0)
            return -2.5 * np.log10(br) + ecl

        f[fin] = wk.dmag_flux_avg(model, t[fin])
    else:
        raise ValueError(kind)
    return f, None


def wz_frac() -> float:
    return wk.wz.EDGE_FRAC


# -------------------------------------------------------------------------------------- injections

_CTX: dict = {}


def _init(ctx: dict) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    configure()
    _CTX.clear()
    _CTX.update(ctx)
    _CTX["stars"] = {k: wk.load_star(k) for k in ctx["kids"]}


def _one(job):
    """job = (kind, cell or family, j, seed, carrier or None)."""
    kind, c, j, seed, fixed = job
    rng = np.random.default_rng([seed, j])
    kids = _CTX["kids"]
    kid = fixed if fixed is not None else kids[rng.integers(len(kids))]
    st, frac, base = _CTX["stars"][kid], _CTX["frac"][kid], _CTX["base"][kid]
    good = wk.good_mask(st)
    tg = st["TIME"][good]
    over = None
    if kind == "w3":
        rho, blend, te = wk.CELLS[c]
        t0, u0 = rng.uniform(tg.min(), tg.max()), rng.uniform(0, 2)
        f = wk.w3_flux_avg(st["TIME"], t0, te, u0, rho, blend)
        half = te * np.sqrt(4 - u0**2)
    else:
        f, over = hard_flux(c, st, good, rng, _CTX["bright"].get(kid))
    sap, err = wk.inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, rng)
    res = wk.screen_star(st, sap, err, frac, over, local=True)
    new = [x for x in res if (x["t_start"], x["t_end"]) not in base]
    if kind == "w3":
        new = [x for x in new if x["t_start"] <= t0 + half and x["t_end"] >= t0 - half]
    return {
        "kind": kind,
        "cell": c,
        "j": j,
        "kepid": int(kid),
        "flagged": bool(new),
        "pass_old": any(x["pass"] for x in new),
        "pass_local": any(x["pass_local"] for x in new),
        "edge_insufficient": (not any(x["pass_no_shoulder"] for x in new))
        and any(x["insufficient"] for x in new),
        "local_insufficient": any(
            x["local_insufficient"] and x["pass_no_shoulder"] and not x["pass_local"] for x in new
        ),
    }


def wilson(k, n):
    return wk.wilson(k, n)


def comparison_frac(kids_frac, kids_carrier) -> dict:
    stars = {k: wk.load_star(k) for k in kids_frac if have(k)}
    cf = wk.common_frac(stars)
    return {k: wk.frac_for(cf, k) for k in kids_carrier}


def build_ctx(carriers, frac_from) -> dict:
    carriers = [k for k in carriers if have(k)]
    frac = comparison_frac(frac_from, carriers)
    real, bright = {}, {}
    for k in carriers:
        st = wk.load_star(k)
        real[k] = wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], frac[k], local=True)
        for x in real[k]:
            x["pdc_rel_in_run"] = wk.pdc_in_run(st, x)
        bt = bright_times(st)
        if len(bt):
            bright[k] = bt
    base = {k: {(x["t_start"], x["t_end"]) for x in real[k]} for k in carriers}
    return {"kids": carriers, "frac": frac, "base": base, "real": real, "bright": bright}


def run_jobs(ctx, jobs, cores):
    pc = {k: ctx[k] for k in ("kids", "frac", "base", "bright")}
    with ProcessPoolExecutor(cores, initializer=_init, initargs=(pc,)) as ex:
        return list(ex.map(_one, jobs, chunksize=8))


def family_summary(rows, rule) -> dict:
    out = {}
    for fam in sorted({x["cell"] for x in rows if x["kind"] == "ctl"}):
        r = [x for x in rows if x["kind"] == "ctl" and x["cell"] == fam]
        k = sum(x[f"pass_{rule}"] for x in r)
        per = {}
        for x in r:
            p = per.setdefault(str(x["kepid"]), [0, 0])
            p[0] += x[f"pass_{rule}"]
            p[1] += 1
        out[fam] = {
            "n": len(r),
            "sources": len(per),
            "flagged": sum(x["flagged"] for x in r),
            "passed": k,
            "passed_wilson95": wilson(k, len(r)),
            "sources_with_a_pass": sum(v[0] > 0 for v in per.values()),
            "per_source_passed_n": per if fam == "real_bright_eclipse" else None,
            "edge_insufficient": sum(x["edge_insufficient"] for x in r),
            "local_insufficient": sum(x["local_insufficient"] for x in r),
        }
    return out


def w3_summary(rows, rule) -> list:
    out = []
    for c, (rho, blend, te) in enumerate(wk.CELLS):
        r = [x for x in rows if x["kind"] == "w3" and x["cell"] == c]
        k = sum(x[f"pass_{rule}"] for x in r)
        out.append(
            {
                "rho": rho,
                "blend": blend,
                "tE": te,
                "n": len(r),
                "flagged": sum(x["flagged"] for x in r),
                "recovered": k,
                "recovered_wilson95": wilson(k, len(r)),
                "edge_insufficient": sum(x["edge_insufficient"] for x in r),
                "local_insufficient": sum(x["local_insufficient"] for x in r),
                "sources": len({x["kepid"] for x in r}),
            }
        )
    return out


def eb_summary(real_eb, rule) -> dict:
    """Real EB light curves (uninjected): per flagged window, the single-window outcome without the
    one-off test (window_<rule>) and the full chain (pass / pass_local); per star as well."""
    win_key = "window_old" if rule == "old" else "window_local"
    full_key = "pass" if rule == "old" else "pass_local"
    wins = [(k, x) for k, v in real_eb.items() for x in v]
    stars_flagged = [k for k, v in real_eb.items() if v]
    return {
        "stars_screened": len(real_eb),
        "stars_with_flagged_window": len(stars_flagged),
        "flagged_windows": len(wins),
        "windows_passing_single_window": sum(x[win_key] for _, x in wins),
        "windows_wilson95": wilson(sum(x[win_key] for _, x in wins), len(wins)),
        "stars_with_single_window_pass": sum(
            any(x[win_key] for x in real_eb[k]) for k in stars_flagged
        ),
        "stars_wilson95": wilson(
            sum(any(x[win_key] for x in real_eb[k]) for k in stars_flagged), len(stars_flagged)
        ),
        "stars_with_full_chain_pass": sum(
            any(x[full_key] for x in real_eb[k]) for k in stars_flagged
        ),
    }


def hard_gate(fams, eb, bright_sources) -> dict:
    """Plan section 4: every family passed / n <= 0.10, real EB windows and stars <= 0.10, with the
    minimum sources; returns the per-item decision."""
    g = {f"{k}_le_0.10": v["passed"] / v["n"] <= GATE_CTRL for k, v in fams.items()}
    g["eb_windows_le_0.10"] = (
        eb["flagged_windows"] > 0
        and eb["windows_passing_single_window"] / eb["flagged_windows"] <= GATE_CTRL
    )
    g["eb_stars_le_0.10"] = (
        eb["stars_with_flagged_window"] > 0
        and eb["stars_with_single_window_pass"] / eb["stars_with_flagged_window"] <= GATE_CTRL
    )
    g["eb_min_stars"] = eb["stars_with_flagged_window"] >= MIN_EB_STARS
    g["bright_min_sources"] = bright_sources >= MIN_BRIGHT_SOURCES
    g["pass"] = all(g.values())
    return g


def stage_phase_a(a) -> None:
    usage("phaseA", ["old_bg", "eb_deep", "eb_single"])
    bg = [k for k in group("old_bg") if have(k)]
    ctx = build_ctx(bg, bg)
    s_w3, s_ctl, s_br = SEEDS["A"]
    jobs = [("w3", c, j, s_w3 + c, None) for c in range(len(wk.CELLS)) for j in range(N_INJ)]
    jobs += [("ctl", f, j, s_ctl + i, None) for i, f in enumerate(FAMILIES) for j in range(N_INJ)]
    bright_kids = [k for k in OLD_BRIGHT if k in ctx["bright"]]
    jobs += [
        ("ctl", "real_bright_eclipse", j, s_br + i, k)
        for i, k in enumerate(bright_kids)
        for j in range(N_BRIGHT_PER_STAR)
    ]
    t0 = time.time()
    rows = run_jobs(ctx, jobs, a.cores)
    el = time.time() - t0
    ebk = [k for k in group("eb_deep") + group("eb_single") if have(k)]
    fr = comparison_frac(bg, ebk)
    real_eb = {}
    for k in ebk:
        st = wk.load_star(k)
        real_eb[k] = wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], fr[k], local=True)
    out = {"label": "simulated (injections) and derived from observed light curves (EB)"}
    out["elapsed_s"] = round(el, 1)
    out["bright_sources"] = bright_kids
    for rule in ("old", "local"):
        fams = family_summary(rows, rule)
        eb = eb_summary(real_eb, rule)
        out[rule] = {
            "w3_seen_backgrounds": w3_summary(rows, rule),
            "families": fams,
            "eb": eb,
            "hard_gate": hard_gate(fams, eb, len(bright_kids)),
        }
    out["eb_deep_screened"] = [k for k in group("eb_deep") if have(k)]
    out["decision"] = "old" if out["old"]["hard_gate"]["pass"] else "local"
    runs = {
        "label": "derived from observed light curves (uninjected screens)",
        "eb_runs": {str(k): v for k, v in real_eb.items()},
        "old_bg_real_runs": {str(k): v for k, v in ctx["real"].items()},
    }
    with gzip.open(OUT / "phaseA_runs.json.gz", "wt") as fh:  # every window; > 1 MB as plain JSON
        fh.write(json.dumps(runs) + "\n")
    (OUT / "phaseA.json").write_text(json.dumps(out) + "\n")
    (OUT / "phaseA_rows.json").write_text(json.dumps({"label": "simulated", "rows": rows}) + "\n")
    (OUT / "decision.json").write_text(
        json.dumps({"label": "derived", "rule": out["decision"]}, indent=1) + "\n"
    )
    print(json.dumps({"decision": out["decision"], "old_gate": out["old"]["hard_gate"]}))


def stage_validate(a) -> None:
    rule = json.loads((OUT / "decision.json").read_text())["rule"]
    usage("validate", ["val", "old_bg"])
    val = [k for k in group("val") if have(k)]
    ctx = build_ctx(val, val + group("old_bg"))
    s_w3, s_ctl = SEEDS["val"]
    jobs = [("w3", c, j, s_w3 + c, None) for c in range(len(wk.CELLS)) for j in range(N_INJ)]
    jobs += [("ctl", f, j, s_ctl + i, None) for i, f in enumerate(FAMILIES) for j in range(N_INJ)]
    t0 = time.time()
    rows = run_jobs(ctx, jobs, a.cores)
    el = time.time() - t0
    pa = json.loads((OUT / "phaseA.json").read_text())
    fams = family_summary(rows, rule)
    fams["real_bright_eclipse (phase A, seen backgrounds, conditional)"] = pa[rule]["families"][
        "real_bright_eclipse"
    ]
    eb = pa[rule]["eb"]
    w3 = w3_summary(rows, rule)
    g = hard_gate(fams, eb, len(pa["bright_sources"]))
    g["all_cells_ge_0.50"] = all(x["recovered"] / x["n"] >= GATE_REC for x in w3)
    g["val_stars_available"] = len(val)
    g["pass"] = all(v for k, v in g.items() if k not in ("pass", "val_stars_available"))
    out = {"label": "simulated", "rule": rule, "elapsed_s": round(el, 1), "carriers": len(val)}
    out.update({"w3": w3, "families": fams, "eb": eb, "gate": g})
    out["other_rule"] = {
        "w3": w3_summary(rows, "local" if rule == "old" else "old"),
        "families": family_summary(rows, "local" if rule == "old" else "old"),
    }
    out["val_real_runs"] = {str(k): v for k, v in ctx["real"].items()}
    (OUT / "validation.json").write_text(json.dumps(out) + "\n")
    (OUT / "validation_rows.json").write_text(
        json.dumps({"label": "simulated", "rows": rows}) + "\n"
    )
    (OUT / "gate.json").write_text(
        json.dumps({"label": "derived", "rule": rule, **g}, indent=1) + "\n"
    )
    print(json.dumps({"rule": rule, "gate": g}))


def _screen_search(args):
    kepid, frac = args
    configure()
    st = wk.load_star(kepid)
    res = wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], frac, local=True)
    for x in res:
        x["pdc_rel_in_run"] = wk.pdc_in_run(st, x)
    return kepid, res


def stage_search(a) -> None:
    gate = json.loads((OUT / "gate.json").read_text())
    if not gate["pass"]:
        raise SystemExit("validation gate not passed: the search group is not read")
    rule = gate["rule"]
    usage("search", ["search", "val", "old_bg"])
    order = [k for k in group("search") if have(k)]
    frac = comparison_frac(order + group("val"), order)
    with ProcessPoolExecutor(a.cores) as ex:
        res = dict(ex.map(_screen_search, [(k, frac[k]) for k in order], chunksize=4))
    key = "pass" if rule == "old" else "pass_local"
    surv = [{"kepid": k, "rank": i, **x} for i, k in enumerate(order) for x in res[k] if x[key]]
    out = {"label": "derived from observed light curves", "rule": rule, "stars": len(order)}
    out["missing"] = [k for k in group("search") if not have(k)]
    out["survivors_in_order"] = surv
    out["runs"] = {str(k): v for k, v in res.items() if v}
    (OUT / "search.json").write_text(json.dumps(out, indent=1) + "\n")
    n_runs = sum(len(v) for v in res.values())
    print(len(order), "stars;", n_runs, "flagged runs;", len(surv), "survivors")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("sample")
    p.add_argument("--pool", required=True)
    p.set_defaults(fn=stage_sample)
    p = sub.add_parser("fetch")
    p.add_argument("--group", required=True)
    p.set_defaults(fn=stage_fetch)
    for name, fn in (
        ("phaseA", stage_phase_a),
        ("validate", stage_validate),
        ("search", stage_search),
    ):
        p = sub.add_parser(name)
        p.add_argument("--cores", type=int, default=4)
        p.set_defaults(fn=fn)
    a = ap.parse_args(argv)
    configure()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
