"""W3 / Kepler effective search volume of the #542 sample (injection-recovery, searched stars).

Plan: ``results/w3_kepler/effsearch/PLAN.md`` (frozen before any data were restored or
injected). The #529 / #542 screen is used unchanged (``w3_kepler.screen_star(local=True)``,
survivor ``pass_local``, common mode from ``w3_kepler_expand``). Nothing here changes a threshold
or adds a test.

Stages (``python scripts/w3_kepler_effsearch.py <stage>``):

- ``fetch``: restore the 500 stars' Q1-Q17 files by the names and sha256 in
  ``expand/files_expand.json.gz`` through the download ledger (Q0 not fetched);
- ``check``: re-run the uninjected #542 screen and compare it with ``expand/screen_expand.json.gz``
  and ``expand/runs_expand.json.gz``; per-star sampling windows;
- ``inject``: 2 injections per star and W3 cell, each from the original SAP series;
- ``summarise``: effective search volume per cell, uncertainties, loss stages, conditional limits
  and the expansion estimate.

Labels: injections simulated; light curves observed; recoveries and exposures derived.
"""

from __future__ import annotations

import argparse
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
import w3_kepler_expand as we  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "effsearch"
EXP = ROOT / "results" / "w3_kepler" / "expand"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "cache" / "w3_kepler_effsearch"
CAP = 6_000_000_000  # owner limit for this round, bytes
SEED = 54200  # PLAN section 3: default_rng([SEED, rank, cell, rep])
N_REP = 2  # injections per star and cell
YEAR = 365.25
DIM_REL = we.DIM_REL  # 1 mag
POISSON95 = float(-np.log(0.05))  # classical Poisson 95 % upper limit for zero events (2.996)
POOL_UNUSED = 18240 - 500  # expand/sample_expand.json candidates_after_exclusion minus #542
ORDER = ("common_mode_ok", "sharp_ok", "oneoff_ok", "local_ok")  # as w3_kepler_expand.layers
STAGE_NAMES = {
    "common_mode_ok": "common_mode",
    "sharp_ok": "edge",
    "oneoff_ok": "oneoff",
    "local_ok": "shoulder",
}


def configure() -> None:
    wk.OUT, wk.CACHE, wk.LEDGER_CAP = OUT, CACHE, CAP


def usage(stage: str) -> None:
    with open(OUT / "source_usage.jsonl", "a") as fh:
        fh.write(json.dumps({"t": round(time.time(), 1), "stage": stage}) + "\n")


def stars() -> list[int]:
    return json.loads((EXP / "sample_expand.json").read_text())["stars"]


def file_list() -> dict[int, list[dict]]:
    d = json.loads(gzip.open(EXP / "files_expand.json.gz", "rt").read())
    return {s["kepid"]: s["files"] for s in d["stars"]}


# ------------------------------------------------------------------------------------------- fetch


def columns(b: bytes) -> dict | None:
    """The columns ``w3_kepler._fetch_star`` caches from one ``_llc.fits`` file (None for Q0)."""
    from astropy.io import fits

    with fits.open(io.BytesIO(b)) as h:
        h0, d = h[0].header, h[1].data
        q = int(h0["QUARTER"])
        if not 1 <= q <= 17:
            return None
        arr = {c: np.asarray(d[c], float) for c in wk.COLS}
        arr["QUARTER"] = np.full(len(d), q, float)
        arr["CHANNEL"] = np.full(len(d), int(h0["CHANNEL"]), float)
        return arr


def restore_star(kepid: int, files: list[dict]) -> dict:
    """Fetch the recorded Q1-Q17 files of one star, check each sha256, cache as #542 did."""
    configure()
    kid = f"{kepid:09d}"
    if (CACHE / f"{kid}.npz").exists():
        return {"kepid": kepid, "cached": True}
    lg = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    base = f"{wk.ARCHIVE}/{kid[:4]}/{kid}/"
    parts, got = [], []
    try:
        for f in sorted(files, key=lambda f: f["name"]):
            if not 1 <= f["quarter"] <= 17:
                continue  # Q0 is not fetched (PLAN section 1)
            b = wk.nl.fetch(base + f["name"], lg, est_bytes=f["bytes"], retries=2)
            ok = hashlib.sha256(b).hexdigest() == f["sha256"]
            got.append({"name": f["name"], "bytes": len(b), "sha256_match": ok})
            if not ok:
                return {"kepid": kepid, "error": f"sha256 mismatch {f['name']}", "files": got}
            parts.append(columns(b))
    except Exception as e:  # recorded; the star is not replaced
        return {"kepid": kepid, "error": f"{type(e).__name__}: {e}", "files": got}
    cat = {c: np.concatenate([p[c] for p in parts]) for c in parts[0]}
    o = np.argsort(cat["TIME"], kind="stable")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f"{kid}.npz", **{c: v[o] for c, v in cat.items()})
    return {"kepid": kepid, "files": got}


def stage_fetch(_a) -> None:
    usage("fetch")
    fl = file_list()
    order = stars()
    t0 = time.time()
    with ThreadPoolExecutor(8) as ex:
        info = list(ex.map(lambda k: restore_star(k, fl[k]), order))
    out = {"label": "observed", "elapsed_s": round(time.time() - t0, 1), "stars": info}
    with gzip.open(OUT / "restore_files.json.gz", "wt") as fh:
        fh.write(json.dumps(out) + "\n")
    used = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP).used()
    print(len(info), "stars;", sum("error" in x for x in info), "errors; ledger", used)


# ------------------------------------------------------------------------------------------ check

_CC: dict = {}


def _init(den, num) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    configure()
    _CC["den"], _CC["num"] = den, num


def window(st) -> tuple[float, float]:
    """PLAN section 3: first and last finite TIME of the star's Q1-Q17 cadences (any quality)."""
    t = st["TIME"][np.isfinite(st["TIME"])]
    return float(t.min()), float(t.max())


def uninjected(st) -> list[dict]:
    frac = we.frac_loo(_CC["den"], _CC["num"], st)
    return wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], frac, local=True)


def _check_one(kepid: int) -> dict:
    st = wk.load_star(kepid)
    good = wk.good_mask(st)
    rel, _ = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], good)
    lo, hi = window(st)
    runs = json.loads(json.dumps(uninjected(st)))
    return {
        "kepid": kepid,
        "good_cadences": int(np.isfinite(rel).sum()),
        "t_lo": lo,
        "t_hi": hi,
        "runs": runs,
    }


def counts(kids) -> tuple[np.ndarray, np.ndarray]:
    configure()
    return we.common_counts(kids)


def stage_check(a) -> None:
    usage("check")
    order = stars()
    missing = [k for k in order if not (CACHE / f"{k:09d}.npz").exists()]
    if missing:
        raise SystemExit(f"{len(missing)} stars not restored: {missing[:5]}")
    den, num = counts(order)
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(den, num)) as ex:
        res = list(ex.map(_check_one, order, chunksize=4))
    old_s = json.loads(gzip.open(EXP / "screen_expand.json.gz", "rt").read())["stars"]
    old_s = {r["kepid"]: r for r in old_s}
    old_r = json.loads(gzip.open(EXP / "runs_expand.json.gz", "rt").read())
    bad = []
    for r in res:
        k = r["kepid"]
        if r["good_cadences"] != old_s[k]["good_cadences"]:
            bad.append({"kepid": k, "what": "good_cadences"})
        prev = [
            {x: v for x, v in y.items() if x != "pdc_rel_in_run"} for y in old_r.get(str(k), [])
        ]
        if r["runs"] != prev:
            bad.append({"kepid": k, "what": "runs"})
    out = {
        "label": "derived from observed light curves (uninjected re-screen)",
        "stars": len(res),
        "runs": sum(len(r["runs"]) for r in res),
        "mismatches": bad,
        "match": not bad,
        "per_star": {str(r["kepid"]): [r["good_cadences"], r["t_lo"], r["t_hi"]] for r in res},
        "per_star_columns": ["good_cadences", "t_lo", "t_hi"],
    }
    (OUT / "check.json").write_text(json.dumps(out) + "\n")
    print({k: v for k, v in out.items() if k != "per_star"})


# ----------------------------------------------------------------------------------------- inject


def n_dim_ideal(t_fin, t0, te, u0, rho, blend) -> int:
    """Dim cadences (F <= 1 mag) of the event on a gap-free 29.4-min grid phased to the star's
    nearest finite cadence time (PLAN section 4, loss stage 'short')."""
    k = int(np.searchsorted(t_fin, t0))
    near = [t_fin[i] for i in (k - 1, k) if 0 <= i < len(t_fin)]
    tr = min(near, key=lambda x: abs(x - t0))
    reach = te * (2.0 + rho) + 2 * wk.CAD
    m = int(np.ceil(reach / wk.CAD)) + 1
    grid = t0 + (tr - t0) % wk.CAD + wk.CAD * np.arange(-m, m + 1)
    return int((wk.w3_flux_avg(grid, t0, te, u0, rho, blend) <= DIM_REL).sum())


def depth(x) -> int:
    return next((i for i, t in enumerate(ORDER) if not x[t]), len(ORDER))


def classify(st, sap, err, res, base, f, t0, h, n_ideal) -> dict:
    """PLAN section 4: matched runs, recovery, spurious passes and the loss stage of one event."""
    good = wk.good_mask(st)
    rel, _ = wk.normalise(st, sap, err, good)
    gnorm = np.isfinite(rel)
    t = st["TIME"]
    ts = t[gnorm]
    lo_w, hi_w = t0 - h - wk.CAD / 2, t0 + h + wk.CAD / 2
    new = [x for x in res if (x["t_start"], x["t_end"]) not in base]
    matched = []
    for x in new:
        if x["t_start"] <= t0 + h and x["t_end"] >= t0 - h:
            tin = ts[(ts >= x["t_start"]) & (ts <= x["t_end"])]
            if ((tin >= lo_w) & (tin <= hi_w)).sum() > len(tin) / 2:
                matched.append(x)
    dim = f <= DIM_REL
    exist = np.isfinite(t) & np.isfinite(st["SAP_FLUX"])
    n_exist, n_good = int((dim & exist).sum()), int((dim & gnorm).sum())
    recovered = any(x["pass_local"] for x in matched)
    spurious = any(x["pass_local"] for x in new if x not in matched)
    sub = ""
    if recovered:
        stage = "recovered"
    elif n_ideal < wk.MIN_RUN:
        stage = "short"
    elif n_exist < wk.MIN_RUN:
        stage = "gap"
    elif n_good < wk.MIN_RUN:
        stage = "quality"
    elif not matched:
        stage = "not_flagged"
    else:
        best = max(matched, key=depth)
        stage = STAGE_NAMES[ORDER[depth(best)]]
        if stage == "edge":
            sub = best["edge_class"]
        elif stage == "shoulder":
            sub = "insufficient" if best["local_insufficient"] else "no"
    return {
        "n_ideal": n_ideal,
        "n_exist": n_exist,
        "n_good": n_good,
        "n_new": len(new),
        "n_matched": len(matched),
        "recovered": recovered,
        "spurious": spurious,
        "stage": stage,
        "sub": sub,
    }


EVENT_COLS = ("rank", "kepid", "cell", "rep", "t0", "u0", "h", "n_ideal", "n_exist", "n_good")
EVENT_COLS += ("n_new", "n_matched", "recovered", "spurious", "stage", "sub")


def inject_star(rank: int, kepid: int) -> list[list]:
    """All 40 injections of one star, each from the original SAP series (PLAN section 3)."""
    st = wk.load_star(kepid)
    frac = we.frac_loo(_CC["den"], _CC["num"], st)
    base = {(x["t_start"], x["t_end"]) for x in uninjected(st)}
    t = st["TIME"]
    t_fin = t[np.isfinite(t)]
    lo, hi = window(st)
    rows = []
    for c, (rho, blend, te) in enumerate(wk.CELLS):
        for rep in range(N_REP):
            rng = np.random.default_rng([SEED, rank, c, rep])
            t0, u0 = rng.uniform(lo, hi), rng.uniform(0, 2)
            h = te * np.sqrt(4 - u0**2)
            f = wk.w3_flux_avg(t, t0, te, u0, rho, blend)
            sap, err = wk.inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, rng)
            res = wk.screen_star(st, sap, err, frac, local=True)
            nid = n_dim_ideal(t_fin, t0, te, u0, rho, blend)
            ev = classify(st, sap, err, res, base, f, t0, h, nid)
            ev |= {"rank": rank, "kepid": kepid, "cell": c, "rep": rep}
            ev |= {"t0": round(float(t0), 5), "u0": round(float(u0), 5), "h": round(float(h), 5)}
            rows.append([ev[k] for k in EVENT_COLS])
    return rows


def _inject_job(arg) -> dict:
    rank, kepid = arg
    t0 = time.time()
    rows = inject_star(rank, kepid)
    return {"rank": rank, "kepid": kepid, "rows": rows, "s": round(time.time() - t0, 2)}


def stage_inject(a) -> None:
    check = json.loads((OUT / "check.json").read_text())
    if not check["match"]:
        raise SystemExit("restore check failed: no injection (PLAN section 1)")
    usage("inject")
    order = stars()
    part = CACHE / "events_partial.jsonl"
    done = set()
    if part.exists():
        done = {json.loads(line)["rank"] for line in part.read_text().splitlines()}
    jobs = [(i, k) for i, k in enumerate(order) if i not in done and i < a.limit]
    den, num = counts(order)
    t0 = time.time()
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(den, num)) as ex:
        with open(part, "a") as fh:
            for n, r in enumerate(ex.map(_inject_job, jobs, chunksize=1), 1):
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                if n % 25 == 0:
                    print(n, "stars", round(time.time() - t0), "s", flush=True)
    recs = sorted((json.loads(x) for x in part.read_text().splitlines()), key=lambda r: r["rank"])
    out = {
        "label": "simulated injections on observed light curves; outcomes derived",
        "plan": "results/w3_kepler/effsearch/PLAN.md",
        "seed": SEED,
        "n_rep": N_REP,
        "cells": [list(c) for c in wk.CELLS],
        "stars_done": len(recs),
        "elapsed_s_this_call": round(time.time() - t0, 1),
        "cpu_s_per_star": [r["s"] for r in recs],
        "columns": list(EVENT_COLS),
        "rows": [row for r in recs for row in r["rows"]],
    }
    with gzip.open(OUT / "events.json.gz", "wt") as fh:
        fh.write(json.dumps(out) + "\n")
    print(len(recs), "stars;", len(out["rows"]), "events")


# -------------------------------------------------------------------------------------- summarise


def wilson_p(p: float, n: float, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a proportion p with (effective) trial number n."""
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, float(c - hw)), min(1.0, float(c + hw))


def volume(T: np.ndarray, k: np.ndarray) -> dict:
    """PLAN section 5 for one cell: T per star (yr, sampling window), k recoveries of N_REP."""
    p = k / N_REP
    e = T * p
    E, sum_t, n = float(e.sum()), float(T.sum()), len(T)
    eps = E / sum_t
    n_eff = N_REP * sum_t**2 / float((T**2).sum())
    lo, hi = wilson_p(eps, n_eff)
    m = float(e.mean())
    se = float(e.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    pop = (min(m - 1.96 * se, sum_t * lo / n), max(m + 1.96 * se, sum_t * hi / n))
    return {
        "sampling_star_years": round(sum_t, 2),
        "recovery_rate": round(eps, 4),
        "effective_star_years": round(E, 1),
        "effective_mc95": [round(sum_t * lo, 1), round(sum_t * hi, 1)],
        "n_eff": round(n_eff, 1),
        "per_star_mean_yr": round(m, 4),
        "per_star_se_yr": round(se, 4),
        "per_star_pop95_yr": [round(max(pop[0], 0.0), 4), round(pop[1], 4)],
        "stars_k": [int((k == i).sum()) for i in range(N_REP + 1)],
    }


def limit(e_low: float, T_mean: float) -> dict | None:
    """PLAN section 6: classical Poisson 95 % limit for zero events at the low exposure, and the
    share of event-hosting stars with >= 2 events at that rate."""
    if e_low <= 0:
        return None
    g = POISSON95 / e_low
    lam = g * T_mean
    share2 = 1 - lam * np.exp(-lam) / (1 - np.exp(-lam))
    return {
        "rate95_per_star_year": float(f"{g:.4g}"),
        "lambda_per_star": float(f"{lam:.4g}"),
        "share_hosts_ge2": float(f"{share2:.4g}"),
    }


def stage_summarise(_a) -> None:
    ev = json.loads(gzip.open(OUT / "events.json.gz", "rt").read())
    cols = {c: i for i, c in enumerate(ev["columns"])}
    rows = ev["rows"]
    chk = json.loads((OUT / "check.json").read_text())["per_star"]
    order = stars()
    ranks = sorted({r[cols["rank"]] for r in rows})
    kid = [order[i] for i in ranks]
    T = np.array([(chk[str(k)][2] - chk[str(k)][1]) / YEAR for k in kid])
    t_obs = float(sum(chk[str(k)][0] for k in kid) * wk.CAD / YEAR)
    pos = {r: i for i, r in enumerate(ranks)}
    app = json.loads((EXP / "applicability.json").read_text())
    out_rng = {x["kepid"] for x in app["not_applicable"]}
    eb = set(EB_STARS)
    cells = []
    for c, (rho, blend, te) in enumerate(wk.CELLS):
        rc = [r for r in rows if r[cols["cell"]] == c]
        k = np.zeros(len(ranks))
        for r in rc:
            k[pos[r[cols["rank"]]]] += r[cols["recovered"]]
        v = volume(T, k)
        stages: dict[str, int] = {}
        for r in rc:
            key = r[cols["stage"]] + (f":{r[cols['sub']]}" if r[cols["sub"]] else "")
            stages[key] = stages.get(key, 0) + 1
        cov = [r for r in rc if r[cols["n_good"]] >= wk.MIN_RUN]
        sub = {}
        for name, grp in (("outside_range", out_rng), ("kirk_eb", eb)):
            g = [r for r in rc if r[cols["kepid"]] in grp]
            sub[name] = {"events": len(g), "recovered": sum(r[cols["recovered"]] for r in g)}
        lo_e = v["effective_mc95"][0]
        exp = {}
        for n_new in (500, 5000):
            de = n_new * v["per_star_mean_yr"]
            de_lo, de_hi = (n_new * x for x in v["per_star_pop95_yr"])
            exp[str(n_new)] = {
                "gain_star_years": round(de, 1),
                "gain95": [round(de_lo, 1), round(de_hi, 1)],
                "limit_factor": round(
                    v["effective_star_years"] / (v["effective_star_years"] + de), 3
                )
                if v["effective_star_years"] + de > 0
                else None,
            }
        cells.append(
            {
                "cell": c,
                "rho": rho,
                "blend": blend,
                "tE_d": te,
                "events": len(rc),
                "recovered": int(k.sum()),
                **v,
                "covered_events": len(cov),
                "recovered_given_covered": round(
                    sum(r[cols["recovered"]] for r in cov) / len(cov), 4
                )
                if cov
                else None,
                "stages": dict(sorted(stages.items())),
                "spurious_events": sum(r[cols["spurious"]] for r in rc),
                "subgroups_descriptive": sub,
                "limit_product_level": limit(lo_e, float(T.mean())),
                "expansion": exp,
            }
        )
    out = {
        "label": "derived (simulated injections on observed light curves)",
        "plan": "results/w3_kepler/effsearch/PLAN.md",
        "stars": len(ranks),
        "events": len(rows),
        "observing_star_years": round(t_obs, 1),
        "sampling_star_years": round(float(T.sum()), 1),
        "sampling_window_days_min_max": [
            round(float(T.min() * YEAR), 2),
            round(float(T.max() * YEAR), 2),
        ],
        "pool_unused_after_542": POOL_UNUSED,
        "poisson95_zero_events": round(POISSON95, 4),
        "cells": cells,
    }
    (OUT / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    for x in cells:
        print(
            x["cell"],
            x["rho"],
            x["blend"],
            x["tE_d"],
            x["recovery_rate"],
            x["effective_star_years"],
            x["effective_mc95"],
            x["stages"],
        )


EB_STARS = (
    2717141,
    3954798,
    4563150,
    6060580,
    6311637,
    6699562,
    7877062,
    8218274,
    9772642,
    10334122,
)  # the 500 stars in expand/inputs/kirk2016_catalog.csv

GAP_D = 0.5  # post hoc: a gap is >= 0.5 d between consecutive good cadences
LOSS_GROUPS = {  # post hoc grouping of the frozen loss stages
    "coverage": ("gap", "quality"),
    "tests_reach_data_edge": ("edge:insufficient", "shoulder:insufficient"),
    "short_or_grazing": ("short",),
    "shape_tests": ("edge:slow", "shoulder:no"),
    "oneoff": ("oneoff",),
    "not_flagged": ("not_flagged",),
}


def near_edges(t_good, q_good, lo_n, hi_n) -> bool:
    """Post hoc: does [lo_n, hi_n] reach the start or end of the good data, a quarter boundary or
    a gap of >= GAP_D between consecutive good cadences?"""
    if lo_n <= t_good[0] or hi_n >= t_good[-1]:
        return True
    i, j = np.searchsorted(t_good, [lo_n, hi_n])
    seg_t = t_good[max(i - 1, 0) : j + 1]
    seg_q = q_good[max(i - 1, 0) : j + 1]
    return bool((np.diff(seg_t) >= GAP_D).any() or (np.diff(seg_q) != 0).any())


def _posthoc_star(arg) -> list[list]:
    rank, kepid, evs = arg
    st = wk.load_star(kepid)
    good = wk.good_mask(st)
    rel, _ = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], good)
    g = np.isfinite(rel)
    t_good, q_good = st["TIME"][g], st["QUARTER"][g]
    out = []
    for c, rep, t0, u0, h in evs:
        rho, blend, te = wk.CELLS[c]
        span = max(2 * h, 1.0)
        w = wk.wz.EDGE_FRAC * span + max(1.0, span)  # the farthest reach of the local shoulder
        near = near_edges(t_good, q_good, t0 - h - w, t0 + h + w)
        f = wk.w3_flux_avg(st["TIME"], t0, te, u0, rho, blend)
        dim = (f <= DIM_REL) & g
        qfrac = max(
            (dim & (st["QUARTER"] == q)).sum() / (g & (st["QUARTER"] == q)).sum()
            for q in np.unique(q_good)
        )
        out.append([rank, c, rep, near, round(float(qfrac), 4), bool(u0 >= 2 - rho)])
    return out


def stage_posthoc(a) -> None:
    """Post hoc diagnostics of the loss stages (labelled post hoc; no outcome is changed)."""
    ev = json.loads(gzip.open(OUT / "events.json.gz", "rt").read())
    cols = {c: i for i, c in enumerate(ev["columns"])}
    rows = ev["rows"]
    by_star: dict[int, list] = {}
    for r in rows:
        by_star.setdefault(r[cols["rank"]], []).append(r)
    order = stars()
    jobs = [
        (
            rk,
            order[rk],
            [[r[cols[k]] for k in ("cell", "rep", "t0", "u0", "h")] for r in rs],
        )
        for rk, rs in sorted(by_star.items())
    ]
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(None, None)) as ex:
        diag = [x for part in ex.map(_posthoc_star, jobs, chunksize=4) for x in part]
    key = {(d[0], d[1], d[2]): d for d in diag}
    per_stage: dict[str, dict] = {}
    for r in rows:
        d = key[(r[cols["rank"]], r[cols["cell"]], r[cols["rep"]])]
        stage = r[cols["stage"]] + (f":{r[cols['sub']]}" if r[cols["sub"]] else "")
        te = wk.CELLS[r[cols["cell"]]][2]
        for grp in (stage, f"{stage}|tE={te}"):
            s = per_stage.setdefault(grp, {"events": 0, "near_edges": 0, "quarter_dim_ge_0.5": 0})
            s["events"] += 1
            s["near_edges"] += d[3]
            s["quarter_dim_ge_0.5"] += d[4] >= 0.5
        if stage == "short":
            rho = wk.CELLS[r[cols["cell"]]][0]
            s = per_stage.setdefault(
                f"short|rho={rho}", {"events": 0, "grazing_u0_ge_2_minus_rho": 0}
            )
            s["events"] += 1
            s["grazing_u0_ge_2_minus_rho"] += d[5]
    out = {
        "label": "post hoc diagnostics (derived from the simulated injections; no outcome changed)",
        "definitions": {
            "near_edges": "the window t0 +- (h + 0.25 span + max(1 d, span)), span = max(2h, 1 d), "
            "reaches the start or end of the good data, a quarter boundary or a gap >= 0.5 d "
            "between consecutive good cadences",
            "quarter_dim_ge_0.5": "the noiseless model is >= 1 mag dim on at least half of the "
            "good cadences of some quarter (that quarter's median is then set by the event)",
            "grazing_u0_ge_2_minus_rho": "u0 >= 2 - rho: the source never lies wholly in the umbra",
        },
        "per_stage": dict(sorted(per_stage.items())),
        "loss_groups": {k: list(v) for k, v in LOSS_GROUPS.items()},
    }
    grouped = {}
    for c in range(len(wk.CELLS)):
        rc = [r for r in rows if r[cols["cell"]] == c]
        st_c = [r[cols["stage"]] + (f":{r[cols['sub']]}" if r[cols["sub"]] else "") for r in rc]
        grouped[str(c)] = {g: sum(s in v for s in st_c) for g, v in LOSS_GROUPS.items()}
    all_st = [r[cols["stage"]] + (f":{r[cols['sub']]}" if r[cols["sub"]] else "") for r in rows]
    grouped["all_cells"] = {g: sum(s in v for s in all_st) for g, v in LOSS_GROUPS.items()}
    grouped["all_cells"]["lost"] = sum(s != "recovered" for s in all_st)
    out["grouped_losses"] = grouped
    summ = json.loads((OUT / "summary.json").read_text())
    t_obs = summ["observing_star_years"]
    out["ratios"] = {
        "effective_over_observing_1732": {
            str(x["cell"]): round(x["effective_star_years"] / t_obs, 3) for x in summ["cells"]
        },
        "naive_2.996_over_observing_not_used": float(f"{POISSON95 / t_obs:.4g}"),
    }
    (OUT / "posthoc.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(grouped["all_cells"]), json.dumps(out["ratios"]))
    for k, v in sorted(per_stage.items()):
        if "|" not in k or k.startswith("short|"):
            print(k, v)


EXAMPLE_STAGES = ("recovered", "edge:slow", "edge:insufficient", "shoulder:insufficient")
EXAMPLE_STAGES += ("shoulder:no", "oneoff", "not_flagged", "gap")


def stage_examples(_a) -> None:
    """Look-at figure (not part of the measurement): for each stage in EXAMPLE_STAGES, the first
    event in file order (rank, cell, rep) is injected again with its own seed and plotted."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ev = json.loads(gzip.open(OUT / "events.json.gz", "rt").read())
    cols = {c: i for i, c in enumerate(ev["columns"])}
    pick = {}
    for r in ev["rows"]:
        key = r[cols["stage"]] + (f":{r[cols['sub']]}" if r[cols["sub"]] else "")
        if key in EXAMPLE_STAGES and key not in pick:
            pick[key] = r
    den, num = counts(stars())
    _CC["den"], _CC["num"] = den, num
    fig, axes = plt.subplots(len(pick), 1, figsize=(9, 2.2 * len(pick)))
    meta = []
    for ax, key in zip(np.atleast_1d(axes), [k for k in EXAMPLE_STAGES if k in pick], strict=True):
        r = pick[key]
        rank, kepid, c, rep = (r[cols[x]] for x in ("rank", "kepid", "cell", "rep"))
        rho, blend, te = wk.CELLS[c]
        st = wk.load_star(kepid)
        lo, hi = window(st)
        rng = np.random.default_rng([SEED, rank, c, rep])
        t0, u0 = rng.uniform(lo, hi), rng.uniform(0, 2)
        f = wk.w3_flux_avg(st["TIME"], t0, te, u0, rho, blend)
        sap, err = wk.inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, rng)
        rel, _ = wk.normalise(st, sap, err, wk.good_mask(st))
        h = te * np.sqrt(4 - u0**2)
        w = 3 * max(h, 1.0)
        m = np.isfinite(rel) & (np.abs(st["TIME"] - t0) <= w)
        ax.plot(st["TIME"][m] - t0, rel[m], ".", ms=2, color="k")
        ax.axvspan(-h, h, color="tab:blue", alpha=0.15)
        ax.axhline(DIM_REL, color="tab:red", lw=0.6)
        ax.set_xlim(-w, w)
        ax.set_title(
            f"{key}: KIC {kepid}, rho {rho}, blend {blend}, tE {te} d, u0 {u0:.2f}, rep {rep}",
            fontsize=8,
        )
        ax.tick_params(labelsize=7)
        meta.append({"stage": key, "rank": rank, "kepid": kepid, "cell": c, "rep": rep})
    np.atleast_1d(axes)[-1].set_xlabel("t - t0 (d); shaded: umbra |t - t0| < h", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "examples.png", dpi=90)
    out = {"label": "simulated injections on observed light curves (look-at figure)", "rows": meta}
    (OUT / "examples.json").write_text(json.dumps(out, indent=1) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("fetch")
    p.set_defaults(fn=stage_fetch)
    p = sub.add_parser("check")
    p.add_argument("--cores", type=int, default=4)
    p.set_defaults(fn=stage_check)
    p = sub.add_parser("inject")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--limit", type=int, default=10**9, help="only ranks below this (timing)")
    p.set_defaults(fn=stage_inject)
    p = sub.add_parser("summarise")
    p.set_defaults(fn=stage_summarise)
    p = sub.add_parser("examples")
    p.set_defaults(fn=stage_examples)
    p = sub.add_parser("posthoc")
    p.add_argument("--cores", type=int, default=4)
    p.set_defaults(fn=stage_posthoc)
    a = ap.parse_args(argv)
    configure()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
