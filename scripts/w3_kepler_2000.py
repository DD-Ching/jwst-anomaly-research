"""W3 / Kepler fixed-rule search of 2,000 new stars in four batches of 500 (#554).

Plan: ``results/w3_kepler/s2000/PLAN.md`` (frozen before any new flux was read). The #529 / #542
screen (``w3_kepler.screen_star(local=True)``, survivor ``pass_local``) and the #546
injection-recovery measurement (``w3_kepler_effsearch.inject_star``) are used unchanged. Nothing
here changes a threshold or adds a test.

Stages (``python scripts/w3_kepler_2000.py <stage> [--batch b]``):

- ``sample``: metadata-only selection of up to 2,000 never-used pool stars (``sample_2000.json``);
- ``fetch --batch b``: listing + Q1-Q17 ``_llc.fits`` through the download ledger (Q0 not fetched);
- ``screen --batch b``: uninjected screen with the batch's own common-mode set; layers;
- ``inject --batch b``: 2 injections per star and W3 cell, each from the original SAP series;
- ``summarise --batch b``: effective search volume, losses and limits of the batch;
- ``vet --batch b``: frozen diagnostics of the batch's survivors (PLAN section 6);
- ``clean --batch b``: delete the batch's column cache, survivors kept until vetted;
- ``combine``: round totals, and the round plus the #546 batch;
- ``limitcheck``: joint (efficiency-marginalised) limits of #546 from its stored events.

Labels: light curves and flags observed; injections simulated; counts, efficiencies and limits
derived.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w3_kepler as wk  # noqa: E402
import w3_kepler_effsearch as es  # noqa: E402
import w3_kepler_expand as we  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "s2000"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "cache" / "w3_kepler_s2000"
CAP = 18_000_000_000  # owner limit for the whole round, bytes
BATCH_GUARD = 16_000_000_000  # a batch is fetched only if used + BATCH_EST stays below this
BATCH_EST = 3_900_000_000  # bytes per 500-star batch (#542: 3.82 GB with listings and Q0)
SALT = "w3-kepler-s2000-2026-10-26"
N_MAX = 2000
N_BATCH = 500
SEED = 54700  # injection rng key: default_rng([SEED, global rank, cell, rep])
N_VET = 10
Q0_STAMP = "2009131105131"  # the Q0 long-cadence file timestamp (every Q0 file of #542)
CL = 0.95


def configure() -> None:
    wk.OUT, wk.CACHE, wk.LEDGER_CAP = OUT, CACHE, CAP


def usage(stage: str) -> None:
    with open(OUT / "source_usage.jsonl", "a") as fh:
        fh.write(json.dumps({"t": round(time.time(), 1), "stage": stage}) + "\n")


def kirk_ebs() -> set[int]:
    path = ROOT / "results" / "w3_kepler" / "expand" / "inputs" / "kirk2016_catalog.csv"
    return {int(r["KIC"]) for r in csv.DictReader(open(path))}


# ------------------------------------------------------------------------------------------ sample


def used_before() -> dict[str, list[int]]:
    """Every star whose light curve an earlier Kepler W3 round read (#500, #529, #542 / #546)."""
    k = ROOT / "results" / "w3_kepler"
    a = json.loads((k / "sample.json").read_text())
    h = json.loads((k / "hard" / "sample_hard.json").read_text())["groups"]
    e = json.loads((k / "expand" / "sample_expand.json").read_text())
    out = {
        f"pilot_{g}": [s["kepid"] for s in a["stars"] if s["group"] == g] for g in ("dev", "val")
    }
    out["pilot_search"] = [s["kepid"] for s in a["stars"] if s["group"] == "search"]
    out["pilot_interface"] = [a["interface_test_star"]["kepid"]]
    out |= {f"hard_{g}": list(v) for g, v in h.items()}
    out["expand_search_and_effsearch"] = list(e["stars"])
    return {g: sorted(v) for g, v in out.items()}


def select(pool: list[int], used: set[int], n: int = N_MAX) -> list[int]:
    """PLAN section 1: pool minus every used star, sorted by sha256('SALT:<kepid>'), first n."""
    cand = [k for k in pool if k not in used]
    return sorted(cand, key=lambda k: hashlib.sha256(f"{SALT}:{k}".encode()).hexdigest())[:n]


def stage_sample(a) -> None:
    raw = Path(a.pool).read_bytes()
    if hashlib.sha256(raw).hexdigest() != wk.POOL_SHA256:
        raise SystemExit("pool file differs from the pinned sha256")
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    meta = {int(r["kepid"]): r for r in rows if r["st_quarters"] == "1" * 17}
    prev = used_before()
    used = set().union(*map(set, prev.values()))
    chosen = select(list(meta), used)
    eb = kirk_ebs()
    stars = [
        {
            "rank": i,
            "kepid": k,
            "batch": i // N_BATCH + 1,
            "kirk2016_eb": k in eb,
            **{x: float(meta[k][x]) for x in ("kepmag", "teff", "logg")},
        }
        for i, k in enumerate(chosen)
    ]
    out = {
        "label": "observed (metadata)",
        "pool_sha256": wk.POOL_SHA256,
        "pool_all_q1_q17": len(meta),
        "used_before": {g: len(v) for g, v in prev.items()},
        "used_before_distinct": len(used),
        "used_before_in_pool": len(used & set(meta)),
        "candidates": len(meta) - len(used & set(meta)),
        "rule": "pool (all Q1-Q17) minus every star read before, sorted by"
        f" sha256('{SALT}:<kepid>');"
        f" first {N_MAX}; batch b = ranks {N_BATCH}(b-1) .. {N_BATCH}b-1; vetting order = batch,"
        " rank, run start",
        "kirk2016_ebs_selected": sum(s["kirk2016_eb"] for s in stars),
        "stars": stars,
    }
    (OUT / "sample_2000.json").write_text(json.dumps(out, indent=1) + "\n")
    print({k: v for k, v in out.items() if k != "stars"})


def batch_stars(b: int) -> list[tuple[int, int]]:
    s = json.loads((OUT / "sample_2000.json").read_text())["stars"]
    return [(x["rank"], x["kepid"]) for x in s if x["batch"] == b]


def have(k: int) -> bool:
    return (CACHE / f"{k:09d}.npz").exists()


# ------------------------------------------------------------------------------------------- fetch


def is_q0(name: str) -> bool:
    return f"-{Q0_STAMP}_llc" in name


def fetch_star(kepid: int) -> dict:
    """Listing plus every Q1-Q17 long-cadence file of one star; cached as ``w3_kepler`` does."""
    configure()
    kid = f"{kepid:09d}"
    if have(kepid):
        return {"kepid": kepid, "cached": True}
    lg = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    base = f"{wk.ARCHIVE}/{kid[:4]}/{kid}/"
    parts, files = [], []
    try:
        listing = wk.nl.fetch(base, lg, 20000, retries=2).decode()
        names = sorted(set(re.findall(r"kplr\d{9}-\d{13}_llc\.fits", listing)))
        for name in names:
            if is_q0(name):
                continue  # Q0 is not fetched (PLAN section 2)
            b = wk.nl.fetch(base + name, lg, est_bytes=520_000, retries=2)
            arr = es.columns(b)
            files.append(
                {
                    "name": name,
                    "bytes": len(b),
                    "sha256": hashlib.sha256(b).hexdigest(),
                    "used": arr is not None,
                }
            )
            if arr is not None:
                parts.append(arr)
        if not parts:
            raise RuntimeError("no Q1-Q17 file")
    except Exception as e:  # recorded as unread; the star is not replaced
        return {"kepid": kepid, "error": f"{type(e).__name__}: {e}", "files": files}
    cat = {c: np.concatenate([p[c] for p in parts]) for c in parts[0]}
    o = np.argsort(cat["TIME"], kind="stable")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f"{kid}.npz", **{c: v[o] for c, v in cat.items()})
    return {"kepid": kepid, "files": files}


def stage_fetch(a) -> None:
    lg = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    used = lg.used()
    if used + BATCH_EST > BATCH_GUARD:
        raise SystemExit(f"ledger {used} B: batch {a.batch} would pass the guard; not fetched")
    usage(f"fetch {a.batch}")
    t0 = time.time()
    with ThreadPoolExecutor(a.threads) as ex:
        info = list(ex.map(fetch_star, [k for _, k in batch_stars(a.batch)]))
    out = {"label": "observed", "batch": a.batch, "elapsed_s": round(time.time() - t0, 1)}
    out["threads"] = a.threads
    out["ledger_before"], out["ledger_after"] = used, lg.used()
    out["stars"] = info
    with gzip.open(OUT / f"files_b{a.batch}.json.gz", "wt") as fh:
        fh.write(json.dumps(out) + "\n")
    print(a.batch, len(info), "stars;", sum("error" in x for x in info), "errors;", lg.used())


# ------------------------------------------------------------------------------------------ screen

_CC: dict = {}


def _init(den, num) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    configure()
    _CC["den"], _CC["num"] = den, num
    es._CC["den"], es._CC["num"] = den, num


def read_stars(b: int) -> list[tuple[int, int]]:
    """The batch's stars with a cache file: the read stars and the batch's common-mode set."""
    return [(r, k) for r, k in batch_stars(b) if have(k)]


def counts(b: int) -> tuple[np.ndarray, np.ndarray]:
    configure()
    return we.common_counts([k for _, k in read_stars(b)])


def _screen_one(kepid: int) -> dict:
    st = wk.load_star(kepid)
    good = wk.good_mask(st)
    rel, sig = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], good)
    ok = np.isfinite(rel)
    frac = we.frac_loo(_CC["den"], _CC["num"], st)
    runs = wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], frac, local=True)
    for x in runs:
        x["pdc_rel_in_run"] = wk.pdc_in_run(st, x)
    lo, hi = es.window(st)
    return {
        "kepid": kepid,
        "good_cadences": int(ok.sum()),
        "effective_days": round(float(ok.sum() * wk.CAD), 2),
        "t_lo": lo,
        "t_hi": hi,
        "quarters": len({int(q) for q in st["QUARTER"][ok]}),
        "median_sigma_rel": float(np.median(sig[ok])) if ok.any() else None,
        # descriptive only, as #542: 17 quarters and median sigma_rel <= 3.96e-4
        "applicable": len({int(q) for q in st["QUARTER"][ok]}) == 17
        and ok.any()
        and float(np.median(sig[ok])) <= we.SIGMA_MAX,
        "screened": int(ok.sum()) >= 2 * wk.MIN_BRACKET + wk.MIN_RUN,
        "dim_cadences": int((rel[ok] <= we.DIM_REL).sum()),
        "min_rel": round(float(np.nanmin(rel[ok])), 4) if ok.any() else None,
        "runs": runs,
    }


def stage_screen(a) -> None:
    usage(f"screen {a.batch}")
    read = read_stars(a.batch)
    den, num = counts(a.batch)
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(den, num)) as ex:
        res = list(ex.map(_screen_one, [k for _, k in read], chunksize=4))
    lay = we.layers(res, kirk_ebs())
    lay["L1_selected"] = len(batch_stars(a.batch))
    lay["unread"] = [k for _, k in batch_stars(a.batch) if not have(k)]
    rank = dict((k, r) for r, k in read)
    lay["survivors_in_vet_order"] = [
        {"rank": rank[r["kepid"]], "kepid": r["kepid"], "t_start": x["t_start"]}
        for r in sorted(res, key=lambda r: rank[r["kepid"]])
        for x in sorted(r["runs"], key=lambda x: x["t_start"])
        if x["pass_local"]
    ]
    lay["survivor_sources"] = sorted({s["kepid"] for s in lay["survivors_in_vet_order"]})
    (OUT / f"layers_b{a.batch}.json").write_text(
        json.dumps({"label": "derived", "batch": a.batch, **lay}, indent=1) + "\n"
    )
    summ = [{k: v for k, v in r.items() if k != "runs"} | {"n_runs": len(r["runs"])} for r in res]
    with gzip.open(OUT / f"screen_b{a.batch}.json.gz", "wt") as fh:
        fh.write(json.dumps({"label": "derived from observed light curves", "stars": summ}) + "\n")
    with gzip.open(OUT / f"runs_b{a.batch}.json.gz", "wt") as fh:
        fh.write(json.dumps({str(r["kepid"]): r["runs"] for r in res if r["runs"]}) + "\n")
    keys = ("L1_read", "L1_effective_star_years", "L2a_stars_with_dim_cadence", "L2b_flagged_runs")
    keys += ("L2b_flagged_sources", "L4_survivor_runs", "survivor_sources")
    print(json.dumps({k: lay[k] for k in keys}))


# ------------------------------------------------------------------------------------------ inject


def _inject_job(arg) -> dict:
    rank, kepid = arg
    t0 = time.time()
    rows = es.inject_star(rank, kepid, seed=SEED)
    return {"rank": rank, "kepid": kepid, "rows": rows, "s": round(time.time() - t0, 2)}


def stage_inject(a) -> None:
    if not (OUT / f"layers_b{a.batch}.json").exists():
        raise SystemExit("screen first")
    usage(f"inject {a.batch}")
    read = read_stars(a.batch)
    part = CACHE / f"events_b{a.batch}_partial.jsonl"
    done = set()
    if part.exists():
        done = {json.loads(x)["rank"] for x in part.read_text().splitlines()}
    jobs = [(r, k) for r, k in read if r not in done]
    den, num = counts(a.batch)
    t0 = time.time()
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(den, num)) as ex:
        with open(part, "a") as fh:
            for n, r in enumerate(ex.map(_inject_job, jobs, chunksize=1), 1):
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                if n % 100 == 0:
                    print(a.batch, n, "stars", round(time.time() - t0), "s", flush=True)
    recs = sorted((json.loads(x) for x in part.read_text().splitlines()), key=lambda r: r["rank"])
    out = {
        "label": "simulated injections on observed light curves; outcomes derived",
        "plan": "results/w3_kepler/s2000/PLAN.md",
        "batch": a.batch,
        "seed": SEED,
        "n_rep": es.N_REP,
        "cells": [list(c) for c in wk.CELLS],
        "stars_done": len(recs),
        "elapsed_s_this_call": round(time.time() - t0, 1),
        "cores": a.cores,
        "cpu_s_per_star": [r["s"] for r in recs],
        "columns": list(es.EVENT_COLS),
        "rows": [row for r in recs for row in r["rows"]],
    }
    with gzip.open(OUT / f"events_b{a.batch}.json.gz", "wt") as fh:
        fh.write(json.dumps(out) + "\n")
    print(a.batch, len(recs), "stars;", len(out["rows"]), "events")


# --------------------------------------------------------------------------------------- summaries


def _beta_logpdf(x: np.ndarray, a: float, b: float) -> np.ndarray:
    lb = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    return (a - 1) * np.log(x) + (b - 1) * np.log1p(-x) - lb


def joint_limit(sum_t: float, eps: float, n_eff: float, cl: float = CL) -> float | None:
    """PLAN section 7: upper limit on the rate (per star-year) for zero events, with a flat prior on
    the rate and the efficiency integrated over its posterior Beta(K + 1, n - K + 1), K = eps n.
    Solves  int g(e) (1 - exp(-G T e)) / e de = cl * int g(e) / e de  for G. For an exactly known
    efficiency it equals -ln(1 - cl) / (T eps), the classical zero-event limit."""
    if eps <= 0:
        return None
    k = eps * n_eff
    x = np.linspace(1e-7, 1 - 1e-7, 400_001)
    lg = _beta_logpdf(x, k + 1, n_eff - k + 1)
    w = np.exp(lg - lg.max()) / x
    tot = w.sum()

    def frac(g: float) -> float:
        return float((w * -np.expm1(-g * sum_t * x)).sum() / tot)

    lo, hi = 0.0, 50.0 / (sum_t * eps)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if frac(mid) < cl else (lo, mid)
    return 0.5 * (lo + hi)


def per_star_arrays(events: dict, windows: dict) -> tuple[list[int], np.ndarray, np.ndarray]:
    """Per star: window T (yr, sampling) and recoveries per cell (k, shape stars x cells)."""
    cols = {c: i for i, c in enumerate(events["columns"])}
    kids = sorted({r[cols["kepid"]] for r in events["rows"]})
    pos = {k: i for i, k in enumerate(kids)}
    k = np.zeros((len(kids), len(wk.CELLS)))
    for r in events["rows"]:
        k[pos[r[cols["kepid"]]], r[cols["cell"]]] += r[cols["recovered"]]
    T = np.array([(windows[q][1] - windows[q][0]) / es.YEAR for q in kids])
    return kids, T, k


def cell_volumes(T: np.ndarray, k: np.ndarray) -> list[dict]:
    out = []
    for c, (rho, blend, te) in enumerate(wk.CELLS):
        v = es.volume(T, k[:, c])
        eps = float((T * k[:, c] / es.N_REP).sum() / T.sum())
        g_joint = joint_limit(float(T.sum()), eps, v["n_eff"])
        g_old = es.limit(v["effective_mc95"][0], float(T.mean()))
        out.append(
            {
                "cell": c,
                "rho": rho,
                "blend": blend,
                "tE_d": te,
                **v,
                "limit_joint_per_star_year": float(f"{g_joint:.4g}") if g_joint else None,
                "limit_546_style": g_old,
            }
        )
    return out


def load_events(b: int) -> dict:
    return json.loads(gzip.open(OUT / f"events_b{b}.json.gz", "rt").read())


def windows_of(b: int) -> dict:
    s = json.loads(gzip.open(OUT / f"screen_b{b}.json.gz", "rt").read())["stars"]
    return {x["kepid"]: (x["t_lo"], x["t_hi"]) for x in s}


def stage_summarise(a) -> None:
    ev = load_events(a.batch)
    win = windows_of(a.batch)
    kids, T, k = per_star_arrays(ev, win)
    scr = json.loads(gzip.open(OUT / f"screen_b{a.batch}.json.gz", "rt").read())["stars"]
    cols = {c: i for i, c in enumerate(ev["columns"])}
    cells = cell_volumes(T, k)
    for x in cells:
        rc = [r for r in ev["rows"] if r[cols["cell"]] == x["cell"]]
        st: dict[str, int] = {}
        for r in rc:
            key = r[cols["stage"]] + (f":{r[cols['sub']]}" if r[cols["sub"]] else "")
            st[key] = st.get(key, 0) + 1
        x["stages"] = dict(sorted(st.items()))
        x["spurious_events"] = sum(r[cols["spurious"]] for r in rc)
    out = {
        "label": "derived (simulated injections on observed light curves)",
        "batch": a.batch,
        "stars": len(kids),
        "events": len(ev["rows"]),
        "observing_star_years": round(sum(s["good_cadences"] for s in scr) * wk.CAD / es.YEAR, 1),
        "sampling_star_years": round(float(T.sum()), 1),
        "cells": cells,
    }
    (OUT / f"summary_b{a.batch}.json").write_text(json.dumps(out, indent=1) + "\n")
    for x in cells:
        print(
            x["cell"], x["recovery_rate"], x["effective_star_years"], x["limit_joint_per_star_year"]
        )


def stage_combine(_a) -> None:
    """Round totals over the completed batches, and the round plus the #546 batch (same rule, event
    distribution, injection method, recovery definition and common-mode definition)."""
    done = [b for b in (1, 2, 3, 4) if (OUT / f"events_b{b}.json.gz").exists()]
    Ts, ks, t_obs = [], [], 0.0
    for b in done:
        _, T, k = per_star_arrays(load_events(b), windows_of(b))
        Ts.append(T)
        ks.append(k)
        t_obs += json.loads((OUT / f"summary_b{b}.json").read_text())["observing_star_years"]
    T, k = np.concatenate(Ts), np.concatenate(ks)
    eff = ROOT / "results" / "w3_kepler" / "effsearch"
    ev546 = json.loads(gzip.open(eff / "events.json.gz", "rt").read())
    chk = json.loads((eff / "check.json").read_text())["per_star"]
    win546 = {int(q): (v[1], v[2]) for q, v in chk.items()}
    _, T5, k5 = per_star_arrays(ev546, win546)
    t_obs5 = json.loads((eff / "summary.json").read_text())["observing_star_years"]
    out = {
        "label": "derived",
        "batches": done,
        "round": {
            "stars": len(T),
            "observing_star_years": round(t_obs, 1),
            "sampling_star_years": round(float(T.sum()), 1),
            "cells": cell_volumes(T, k),
        },
        "round_plus_546": {
            "stars": len(T) + len(T5),
            "observing_star_years": round(t_obs + t_obs5, 1),
            "sampling_star_years": round(float(T.sum() + T5.sum()), 1),
            "cells": cell_volumes(np.concatenate([T, T5]), np.concatenate([k, k5])),
        },
    }
    (OUT / "combined.json").write_text(json.dumps(out, indent=1) + "\n")
    for name in ("round", "round_plus_546"):
        c = out[name]["cells"]
        e = [x["effective_star_years"] for x in c]
        g = [x["limit_joint_per_star_year"] for x in c]
        print(name, out[name]["stars"], min(e), max(e), min(g), max(g))


def stage_limitcheck(_a) -> None:
    """PLAN section 7: #546's per-cell limits (2.996 / E_low) against the joint limit, from the
    stored #546 summary only (no injection is re-run)."""
    s = json.loads((ROOT / "results" / "w3_kepler" / "effsearch" / "summary.json").read_text())
    rows = []
    for x in s["cells"]:
        # every #546 window is 1,459.49 d (summary.json), so eps = recovered / events exactly
        g = joint_limit(x["sampling_star_years"], x["recovered"] / x["events"], x["n_eff"])
        old = (x["limit_product_level"] or {}).get("rate95_per_star_year")
        if g is None or old is None:  # a cell without recoveries has no limit
            rows.append({"cell": x["cell"], "limit_546": old, "limit_joint": None})
            continue
        rows.append(
            {
                "cell": x["cell"],
                "limit_546": old,
                "limit_joint": float(f"{g:.4g}"),
                "classical_at_E": float(f"{es.POISSON95 / x['effective_star_years']:.4g}"),
                "ratio_546_over_joint": round(old / g, 4),
            }
        )
    out = {
        "label": "derived from the stored #546 summary (no injection re-run)",
        "method": "joint_limit (flat prior on the rate, efficiency Beta posterior integrated)",
        "cells": rows,
    }
    (OUT / "limit_check_546.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        min(r["ratio_546_over_joint"] for r in rows), max(r["ratio_546_over_joint"] for r in rows)
    )


# --------------------------------------------------------------------------------------------- vet


def vet_diagnostics(st: dict, run: dict, others: dict, in_kirk: bool) -> dict:
    """PLAN section 6, V1-V8 for one survivor run (computable conditions; V9 is the pixel check)."""
    g = wk.good_mask(st)
    rel, sig = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], g)
    t, q = st["TIME"], int(run["quarter_start"])
    t0, t1 = run["t_start"], run["t_end"]
    span = max(t1 - t0, 1.0)
    inq = (st["QUARTER"] == q) & np.isfinite(t)
    win = inq & (t >= t0 - 1) & (t <= t1 + 1)
    qual = st["SAP_QUALITY"][win].astype(np.int64)
    v1 = float((qual != 0).mean()) if win.any() else 0.0
    inrun = g & (t >= t0) & (t <= t1)
    outside = g & inq & ((t < t0 - span) | (t > t1 + span))
    cen = {}
    for c in ("MOM_CENTR1", "MOM_CENTR2"):
        m_o, s_o = _robust(st[c][outside])
        m_i, _ = _robust(st[c][inrun])
        cen[c] = None if not np.isfinite(m_i) or not s_o else round(float((m_i - m_o) / s_o), 2)
    cad = st["CADENCENO"][inrun].astype(int)
    comp = []
    ch = int(st["CHANNEL"][inq][0]) if inq.any() else -1
    for so in others.values():
        cq = so["CHANNEL"][so["QUARTER"] == q]
        if not len(cq) or int(cq[0]) != ch:
            continue
        go = wk.good_mask(so)
        ro, _ = wk.normalise(so, so["SAP_FLUX"], so["SAP_FLUX_ERR"], go)
        m = go & np.isin(so["CADENCENO"].astype(int), cad)
        if m.any():
            comp.append(float(np.median(ro[m])))
    v3 = bool(comp) and (np.median(comp) <= 0.9 or sum(x <= we.DIM_REL for x in comp) >= 3)
    tg = t[g]
    qb = t[inq]
    v4 = bool(
        t0 - qb.min() <= 1
        or qb.max() - t1 <= 1
        or (np.diff(tg[(tg >= t0 - 1) & (tg <= t1 + 1)]) >= 0.5).any()
    )
    shallow = (rel <= 10 ** (-0.4 * 0.3)) & g
    dips = sum(
        1 for i, j in wk.wz.runs(shallow) if j - i + 1 >= 3 and not (t0 - span <= t[i] <= t1 + span)
    )
    v5 = in_kirk or dips >= 2
    prof, base = {}, {}
    for side, edge, sgn in (("pre", t0, -1), ("post", t1, 1)):
        d = (edge - t) if sgn < 0 else (t - edge)
        sw = g & inq & (d > 0) & (d <= wk.wz.EDGE_FRAC * span)
        bw = g & inq & (d > wk.wz.EDGE_FRAC * span) & (d <= wk.wz.EDGE_FRAC * span + max(1.0, span))
        base[side] = float(np.median(rel[bw])) if bw.sum() >= 5 else None
        if sw.sum() >= 2 and base[side] is not None:
            order = np.argsort(d[sw])  # moving away from the edge
            rr = rel[sw][order]
            prof[side] = {
                "peak_index_from_edge": int(np.argmax(rr)),
                "excess_area": round(float(np.clip(rr - base[side], 0, None).sum()), 4),
            }
    ratio = None
    if "pre" in prof and "post" in prof and prof["post"]["excess_area"] > 0:
        ratio = round(prof["pre"]["excess_area"] / prof["post"]["excess_area"], 3)
    asym = (ratio is not None and not 1 / 3 <= ratio <= 3) or (
        "pre" in prof and prof["pre"]["peak_index_from_edge"] >= 3
    )
    bright = (rel - 1 >= np.maximum(0.01, 5 * sig)) & g
    n_bright = sum(
        1 for i, j in wk.wz.runs(bright) if j - i + 1 >= 2 and not (t0 - span <= t[i] <= t1 + span)
    )
    v6 = bool(asym or n_bright >= 3)
    bp, bq = base.get("pre"), base.get("post")
    v7 = bp is not None and bq is not None and abs(bp - bq) / (0.5 * (bp + bq)) >= 0.05
    return {
        "V1_flagged_fraction_within_1d": round(v1, 3),
        "V1_quality": v1 >= 0.5,
        "V2_centroid_shift_sigma": cen,
        "V2_centroid": any(x is not None and abs(x) >= 5 for x in cen.values()),
        "V3_comparison_median_rel": comp,
        "V3_common": bool(v3),
        "V4_near_data_edge": v4,
        "V5_other_dips_ge_0.3mag": dips,
        "V5_kirk2016": in_kirk,
        "V5_repeating": bool(v5),
        "V6_shoulder_profiles": prof,
        "V6_excess_ratio_pre_post": ratio,
        "V6_other_bright_events": n_bright,
        "V6_flare_like": v6,
        "V7_baselines": base,
        "V7_baseline_step": bool(v7),
        "V8_pdcsap_rel_in_run": run.get("pdc_rel_in_run"),
    }


def _robust(x: np.ndarray) -> tuple[float, float]:
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return float("nan"), float("nan")
    m = float(np.median(x))
    return m, 1.4826 * float(np.median(np.abs(x - m)))


def pixel_verdict(flux: np.ndarray, in_run: np.ndarray, out_run: np.ndarray) -> dict:
    """PLAN section 6, V9 on a target pixel cube (cadence, y, x): mean in-run image minus mean
    out-of-run image; 'on_target' if the most negative difference pixel is within 1 pixel of the
    out-of-run flux-weighted centroid."""
    img_in = np.nanmean(flux[in_run], axis=0)
    img_out = np.nanmean(flux[out_run], axis=0)
    diff = img_in - img_out
    yy, xx = np.indices(img_out.shape)
    w = np.clip(np.nan_to_num(img_out), 0, None)
    cy, cx = float((w * yy).sum() / w.sum()), float((w * xx).sum() / w.sum())
    iy, ix = np.unravel_index(int(np.nanargmin(diff)), diff.shape)
    dist = float(np.hypot(iy - cy, ix - cx))
    return {
        "diff_min_pixel_yx": [int(iy), int(ix)],
        "centroid_yx": [round(cy, 2), round(cx, 2)],
        "distance_px": round(dist, 2),
        "diff_sum_over_out_sum": round(float(np.nansum(diff) / np.nansum(img_out)), 4),
        "verdict": "on_target" if dist <= 1.0 else "off_target",
    }


def pixel_check(kepid: int, run: dict, llc_name: str) -> dict:
    """V9: fetch the run quarter's official target pixel file through the ledger and apply
    ``pixel_verdict`` (in-run cadences vs the same quarter within one span + 1 d outside)."""
    from astropy.io import fits

    kid = f"{kepid:09d}"
    stamp = llc_name.split("-")[1].split("_")[0]
    url = (
        "https://archive.stsci.edu/missions/kepler/target_pixel_files/"
        f"{kid[:4]}/{kid}/kplr{kid}-{stamp}_lpd-targ.fits.gz"
    )
    lg = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    b = wk.nl.fetch(url, lg, est_bytes=8_000_000, retries=2)
    with fits.open(io.BytesIO(gzip.decompress(b))) as h:
        d = h[1].data
        t, flux, qual = np.asarray(d["TIME"], float), np.asarray(d["FLUX"], float), d["QUALITY"]
    ok = np.isfinite(t) & ((np.asarray(qual).astype(np.int64) & wk.QMASK) == 0)
    t0, t1 = run["t_start"], run["t_end"]
    span = max(t1 - t0, 1.0)
    in_run = ok & (t >= t0) & (t <= t1)
    out_run = ok & (np.abs(t - 0.5 * (t0 + t1)) <= 0.5 * (t1 - t0) + span + 1.0) & ~in_run
    out_run &= (t < t0) | (t > t1)
    res = pixel_verdict(flux, in_run, out_run)
    res |= {"url": url, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
    res |= {"n_in": int(in_run.sum()), "n_out": int(out_run.sum())}
    return res


def vet_class(d: dict, pixel: str | None = None) -> str:
    """PLAN section 6 classification. pixel: None (not checked), 'on_target' or 'off_target'."""
    if d["V1_quality"] or d["V3_common"] or pixel == "off_target":
        return "artifact_supported"
    if d["V5_repeating"] or d["V6_flare_like"]:
        return "ordinary_explanation_supported"
    if pixel == "on_target" and not d["V7_baseline_step"]:
        return "real_dimming_confirmed_unexplained_by_V5_V7"
    return "undetermined"


def stage_vet(a) -> None:
    lay = json.loads((OUT / f"layers_b{a.batch}.json").read_text())
    runs = json.loads(gzip.open(OUT / f"runs_b{a.batch}.json.gz", "rt").read())
    eb = kirk_ebs()
    configure()
    srcs = []
    for s in lay["survivors_in_vet_order"]:
        if s["kepid"] not in srcs:
            srcs.append(s["kepid"])
    others = {k: wk.load_star(k) for _, k in read_stars(a.batch)}
    out = {"label": "derived from observed light curves (vetting diagnostics)", "sources": []}
    for kid in srcs[: a.max]:
        st = others[kid]
        for x in [r for r in runs[str(kid)] if r["pass_local"]]:
            d = vet_diagnostics(st, x, {k: v for k, v in others.items() if k != kid}, kid in eb)
            out["sources"].append(
                {"kepid": kid, "t_start": x["t_start"], **d, "class_without_pixels": vet_class(d)}
            )
    out["not_vetted"] = srcs[a.max :]
    (OUT / f"vet_b{a.batch}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out)[:3000])


def stage_clean(a) -> None:
    lay = json.loads((OUT / f"layers_b{a.batch}.json").read_text())
    keep = set(lay["survivor_sources"]) if not a.all else set()
    n = 0
    for _, k in batch_stars(a.batch):
        p = CACHE / f"{k:09d}.npz"
        if p.exists() and k not in keep:
            p.unlink()
            n += 1
    part = CACHE / f"events_b{a.batch}_partial.jsonl"
    if part.exists():
        part.unlink()
    print(a.batch, "deleted", n, "kept", len(keep))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("sample")
    p.add_argument("--pool", required=True)
    p.set_defaults(fn=stage_sample)
    p = sub.add_parser("fetch")
    p.add_argument("--batch", type=int, required=True, choices=[1, 2, 3, 4])
    p.add_argument("--threads", type=int, default=8, help="concurrent HTTP connections (8-16)")
    p.set_defaults(fn=stage_fetch)
    p = sub.add_parser("summarise")
    p.add_argument("--batch", type=int, required=True, choices=[1, 2, 3, 4])
    p.set_defaults(fn=stage_summarise)
    for name, fn in (("screen", stage_screen), ("inject", stage_inject)):
        p = sub.add_parser(name)
        p.add_argument("--batch", type=int, required=True, choices=[1, 2, 3, 4])
        p.add_argument("--cores", type=int, default=4)
        p.set_defaults(fn=fn)
    p = sub.add_parser("vet")
    p.add_argument("--batch", type=int, required=True, choices=[1, 2, 3, 4])
    p.add_argument("--max", type=int, default=N_VET)
    p.set_defaults(fn=stage_vet)
    p = sub.add_parser("clean")
    p.add_argument("--batch", type=int, required=True, choices=[1, 2, 3, 4])
    p.add_argument("--all", action="store_true", help="also delete survivors (after vetting)")
    p.set_defaults(fn=stage_clean)
    for name, fn in (("combine", stage_combine), ("limitcheck", stage_limitcheck)):
        p = sub.add_parser(name)
        p.set_defaults(fn=fn)
    a = ap.parse_args(argv)
    configure()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
