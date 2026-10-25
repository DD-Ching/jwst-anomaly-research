"""W3 / Kepler fixed-rule limited expansion search (follow-up to #529).

Plan: ``results/w3_kepler/expand/PLAN.md`` (frozen before any new flux was read). The #529 rule is
used unchanged: ``w3_kepler.screen_star(local=True)`` and its ``pass_local``. Nothing here changes a
threshold or adds a method; this module only selects stars, reads them through the ledger, screens
them and stores every layer.

Stages (``python scripts/w3_kepler_expand.py <stage>``): ``sample``, ``fetch``, ``screen``.
Labels: light curves and flags observed; counts derived. No injection is made in this round.
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
OUT = ROOT / "results" / "w3_kepler" / "expand"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "cache" / "w3_kepler_expand"
CAP = 6_000_000_000  # owner limit for this round, bytes
SALT3 = "w3-kepler-expand-2026-10-25"
N_MAX = 500
N_VET = 10
SIGMA_MAX = 2 * 1.98e-4  # applicability: median sigma_rel <= 2 x the #500 maximum (posthoc2.json)
DIM_REL = 10 ** (-0.4 * wk.wz.DIM_MAG)


def configure() -> None:
    wk.OUT, wk.CACHE, wk.LEDGER_CAP = OUT, CACHE, CAP


def usage(stage: str) -> None:
    with open(OUT / "source_usage.jsonl", "a") as fh:
        fh.write(json.dumps({"t": round(time.time(), 1), "stage": stage}) + "\n")


# ------------------------------------------------------------------------------------------ sample


def previous_rounds() -> dict:
    """Every star whose light curve an earlier Kepler W3 round read, by category (deduplicated)."""
    a = json.loads((OUT.parent / "sample.json").read_text())
    b = json.loads((OUT.parent / "hard" / "sample_hard.json").read_text())["groups"]
    s500 = {
        g: {s["kepid"] for s in a["stars"] if s["group"] == g} for g in ("dev", "val", "search")
    }
    cats = {
        "development": sorted(s500["dev"]),
        "validation": sorted(s500["val"] | set(b["val"])),
        "search": sorted(s500["search"] | set(b["search"])),
        "ordinary_event_controls": sorted(set(b["eb_deep"]) | set(b["eb_single"])),
        "interface_check": [a["interface_test_star"]["kepid"]],
    }
    return cats


def stage_sample(a) -> None:
    raw = Path(a.pool).read_bytes()
    if hashlib.sha256(raw).hexdigest() != wk.POOL_SHA256:
        raise SystemExit("pool file differs from the pinned sha256")
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    pool = [int(r["kepid"]) for r in rows if r["st_quarters"] == "1" * 17]
    prev = previous_rounds()
    used = set().union(*map(set, prev.values()))
    cand = [k for k in pool if k not in used]
    rank = sorted(cand, key=lambda k: hashlib.sha256(f"{SALT3}:{k}".encode()).hexdigest())
    stars = rank[:N_MAX]
    assert not set(stars) & used and len(set(stars)) == len(stars)
    out = {
        "label": "observed (metadata)",
        "pool_sha256": wk.POOL_SHA256,
        "pool_all_q1_q17": len(pool),
        "previous_rounds": {k: len(v) for k, v in prev.items()},
        "previous_rounds_distinct": len(used),
        "previous_fgk_pool_stars_read": len(
            set(prev["development"]) | set(prev["validation"]) | set(prev["search"])
        ),
        "candidates_after_exclusion": len(cand),
        "rule": f"pool minus every star read before, sorted by sha256('{SALT3}:<kepid>'); "
        f"first {N_MAX}; vetting order = this order, then run start time",
        "stars": stars,
    }
    (OUT / "sample_expand.json").write_text(json.dumps(out, indent=1) + "\n")
    print({k: v for k, v in out.items() if k != "stars"})


def stars() -> list[int]:
    return json.loads((OUT / "sample_expand.json").read_text())["stars"]


# ------------------------------------------------------------------------------------------- fetch


def _fetch_one(kepid: int) -> dict:
    configure()
    if (CACHE / f"{kepid:09d}.npz").exists():
        return {"kepid": kepid, "cached": True}
    try:
        return wk._fetch_star(kepid)
    except Exception as e:  # recorded; the star is not replaced
        return {"kepid": kepid, "error": f"{type(e).__name__}: {e}"}


def stage_fetch(_a) -> None:
    usage("fetch")
    with ThreadPoolExecutor(8) as ex:
        info = list(ex.map(_fetch_one, stars()))
    with gzip.open(OUT / "files_expand.json.gz", "wt") as fh:
        fh.write(json.dumps({"label": "observed", "stars": info}) + "\n")
    used = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP).used()
    print(len(info), "stars;", sum("error" in x for x in info), "errors; ledger", used)


def have(k: int) -> bool:
    return (CACHE / f"{k:09d}.npz").exists()


# ------------------------------------------------------------------------------------------ screen


def common_counts(kids) -> tuple[np.ndarray, np.ndarray]:
    """Per cadence number, the counts that w3_kepler.common_frac builds (good, dim), as arrays and
    streamed one star at a time (500 stars do not fit in memory as dicts)."""
    den = np.zeros(80000, np.int32)
    num = np.zeros(80000, np.int32)
    for k in kids:
        g, d = own_cadences(wk.load_star(k))
        np.add.at(den, g, 1)
        np.add.at(num, d, 1)
    return den, num


def own_cadences(st) -> tuple[np.ndarray, np.ndarray]:
    good = wk.good_mask(st)
    rel, _ = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], good)
    ok = np.isfinite(rel)
    ks = st["CADENCENO"][ok].astype(int)
    return ks, ks[rel[ok] <= DIM_REL]


def frac_loo(den, num, st) -> dict:
    """Leave-one-out dim fraction at this star's cadences (same values as w3_kepler.frac_for for
    the cadences the star has; only non-zero entries)."""
    g, d = own_cadences(st)
    ks = np.unique(st["CADENCENO"].astype(int))
    n = num[ks] - np.isin(ks, d)
    dd = den[ks] - np.isin(ks, g)
    m = (dd > 0) & (n > 0)
    return {int(k): float(a / b) for k, a, b in zip(ks[m], n[m], dd[m], strict=True)}


_CC: dict = {}


def _init(den, num):
    os.environ["OMP_NUM_THREADS"] = "1"
    configure()
    _CC["den"], _CC["num"] = den, num


def _screen_one(kepid: int) -> dict:
    st = wk.load_star(kepid)
    good = wk.good_mask(st)
    rel, sig = wk.normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], good)
    ok = np.isfinite(rel)
    frac = frac_loo(_CC["den"], _CC["num"], st)
    runs = wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], frac, local=True)
    for x in runs:
        x["pdc_rel_in_run"] = wk.pdc_in_run(st, x)
    quarters = sorted({int(q) for q in st["QUARTER"][ok]})
    med_sig = float(np.median(sig[ok])) if ok.any() else None
    return {
        "kepid": kepid,
        "good_cadences": int(ok.sum()),
        "effective_days": round(float(ok.sum() * wk.CAD), 2),
        "quarters": len(quarters),
        "median_sigma_rel": med_sig,
        "applicable": len(quarters) == 17 and med_sig is not None and med_sig <= SIGMA_MAX,
        "screened": int(ok.sum()) >= 2 * wk.MIN_BRACKET + wk.MIN_RUN,
        "dim_cadences": int((rel[ok] <= DIM_REL).sum()),
        "min_rel": round(float(np.nanmin(rel[ok])), 4) if ok.any() else None,
        "runs": runs,
    }


def layers(res: list[dict], eb: set[int]) -> dict:
    """Plan section 4: counts per layer; sequential tests in the order common mode, bounded fast
    edges, one-off, local shoulder; per flagged source, the first failed test of its best run."""
    runs = [(r["kepid"], x) for r in res for x in r["runs"]]
    order = ("common_mode_ok", "sharp_ok", "oneoff_ok", "local_ok")
    seq, left = {}, runs
    for t in order:
        left = [(k, x) for k, x in left if x[t]]
        seq[f"after_{t}"] = {"runs": len(left), "sources": len({k for k, _ in left})}

    def depth(x):
        return next((i for i, t in enumerate(order) if not x[t]), len(order))

    reasons = {}
    for r in res:
        if r["runs"]:
            best = max(r["runs"], key=depth)
            reasons[str(r["kepid"])] = {
                "first_failed": order[depth(best)]
                if depth(best) < len(order)
                else "none (survivor)",
                "runs": len(r["runs"]),
                "in_kirk2016_eb_catalogue": r["kepid"] in eb,
            }
    surv = [(k, x) for k, x in runs if x["pass_local"]]
    srcs = []
    for k, _ in surv:
        if k not in srcs:
            srcs.append(k)
    return {
        "L1_selected": N_MAX,
        "L1_read": len(res),
        "L1_screened": sum(r["screened"] for r in res),
        "L1_applicable": sum(r["applicable"] for r in res),
        "L1_effective_days": round(sum(r["effective_days"] for r in res), 1),
        "L1_effective_star_years": round(sum(r["effective_days"] for r in res) / 365.25, 1),
        "L2a_stars_with_dim_cadence": sum(r["dim_cadences"] > 0 for r in res),
        "L2a_dim_cadences": sum(r["dim_cadences"] for r in res),
        "L2b_flagged_runs": len(runs),
        "L2b_flagged_sources": len({k for k, _ in runs}),
        "L3_sequential": seq,
        "L3_individual": {t: sum(x[t] for _, x in runs) for t in order},
        "L4_survivor_runs": len(surv),
        "L4_survivor_sources_in_order": srcs,
        "L4_vet_list": srcs[:N_VET],
        "L4_not_vetted": srcs[N_VET:],
        "flagged_source_reasons": reasons,
    }


def stage_screen(a) -> None:
    usage("screen")
    order = stars()
    read = [k for k in order if have(k)]
    den, num = common_counts(read)
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(den, num)) as ex:
        res = list(ex.map(_screen_one, read, chunksize=4))
    eb = {int(r["KIC"]) for r in csv.DictReader(open(OUT / "inputs" / "kirk2016_catalog.csv"))}
    lay = layers(res, eb)
    lay["missing"] = [k for k in order if not have(k)]
    (OUT / "layers.json").write_text(json.dumps({"label": "derived", **lay}, indent=1) + "\n")
    summ = [{k: v for k, v in r.items() if k != "runs"} | {"n_runs": len(r["runs"])} for r in res]
    with gzip.open(OUT / "screen_expand.json.gz", "wt") as fh:
        fh.write(json.dumps({"label": "derived from observed light curves", "stars": summ}) + "\n")
    with gzip.open(OUT / "runs_expand.json.gz", "wt") as fh:
        fh.write(json.dumps({str(r["kepid"]): r["runs"] for r in res if r["runs"]}) + "\n")
    print(json.dumps({k: v for k, v in lay.items() if k != "flagged_source_reasons"}))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("sample")
    p.add_argument("--pool", required=True)
    p.set_defaults(fn=stage_sample)
    p = sub.add_parser("fetch")
    p.set_defaults(fn=stage_fetch)
    p = sub.add_parser("screen")
    p.add_argument("--cores", type=int, default=4)
    p.set_defaults(fn=stage_screen)
    a = ap.parse_args(argv)
    configure()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
