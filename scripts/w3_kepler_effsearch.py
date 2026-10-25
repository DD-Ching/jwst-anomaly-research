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
    a = ap.parse_args(argv)
    configure()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
