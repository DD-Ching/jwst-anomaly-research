"""S1 M2 search: calibration of the consistency test and the false-candidate rate, then one frozen
full-catalogue re-screen (docs/s1_deform/README.md section 8).

  python scripts/s1_m2_search.py calibrate [--cpu 4]  -> results/s1_m2/calibration_*.json
  python scripts/s1_m2_search.py search [--cpu 4]     -> results/s1_m2/search.json (frozen.json)
  python scripts/s1_m2_search.py efficiency [--cpu 4] -> results/s1_m2/efficiency.json

M2 statistic: max over a log-uniform stretch grid (ratio 1.05 on [0.5, 2]) of the multi-band
cross-correlation with a free per-band scale; consistency: per-band chi2 at a locally refined
stretch and lag.
Injected twins are simulated; recovering them measures search capability only.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import s1_deform as sd  # noqa: E402

from jwst_anomaly import burst_twins as bt  # noqa: E402

OUT = ROOT / "results" / "s1_m2"
P = bt.Params()
GRID_STEP = 1.05
G: dict = {}


def fine_grid(step: float = GRID_STEP) -> tuple[float, ...]:
    n = int(round(math.log(4.0) / math.log(step))) // 2
    return tuple(float(round(step**q, 4)) for q in range(-n, n + 1))


sd.S_GRID = fine_grid()


REFINE = np.arange(-0.03, 0.03001, 0.002)  # ln s offsets around the best grid value


def refined_p(x, ex, y, ey, s_grid_best: float) -> tuple[float, float]:
    """Per-band chi2 p at a locally refined stretch (x, y already at a common level)."""
    best, sb, lb = -1.0, s_grid_best, 0
    for d in REFINE:
        s = s_grid_best * math.exp(d)
        ys = bt.stretch(y, ey, s)[0]
        r, lag = sd.rho_free(x, ys)
        if r > best:
            best, sb, lb = r, s, lag
    ys, eys = bt.stretch(y, ey, sb)
    return sd.chi2_free(x, ex, ys, eys, lb), sb


def m2_only(x, ex, kx, y, ey, ky, with_chi2: bool = False):
    """(m2, s_best[, p_chi2]) for one pair; nan when |dk| > 2."""
    if abs(kx - ky) > sd.MAX_DK_NEW:
        return (np.nan, np.nan, np.nan) if with_chi2 else (np.nan, np.nan)
    x, ex, y, ey = sd.to_common(x, ex, kx, y, ey, ky)
    best, sb = -1.0, 1.0
    for s in sd.S_GRID:
        ys = y if s == 1.0 else bt.stretch(y, ey, s)[0]
        r, _ = sd.rho_free(x, ys)
        if r > best:
            best, sb = r, s
    if not with_chi2:
        return best, sb
    return best, sb, refined_p(x, ex, y, ey, sb)[0]


_scores_orig = sd.scores


def _scores_refined(x, ex, kx, y, ey, ky, with_chi2=False):
    """sd.scores with M2's chi2 p (index 7) replaced by the refined-stretch p."""
    r = _scores_orig(x, ex, kx, y, ey, ky, with_chi2)
    if not with_chi2 or not np.isfinite(r[2]):
        return r
    xc, exc, yc, eyc = sd.to_common(x, ex, kx, y, ey, ky)
    return (*r[:7], refined_p(xc, exc, yc, eyc, r[4])[0])


sd.scores = _scores_refined


def _init(g):
    os.environ["OMP_NUM_THREADS"] = "1"
    G.update(g)
    sd.G.update(g)


def bg_row(i):
    B, allowed, idx = G["B"], G["allowed"], G["idx"]
    out = []
    for j in idx:
        if j <= i or not allowed[i, j]:
            continue
        bi, bj = B[i], B[j]
        m, s = m2_only(bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"])
        out.append((i, j, m, s))
    return out


def setup():
    cat, B, status = sd.st.load()
    elig = np.array([b["npulse"] >= P.min_pulses for b in B])
    ks = np.array([b["k"] for b in B])
    ra, dec, err = (np.array([b[c] for b in B]) for c in ("ra", "dec", "err"))
    sep = bt.angsep_deg(ra[:, None], dec[:, None], ra[None], dec[None])
    incons = bt.position_inconsistent(sep, err[:, None], err[None], P.sys_deg, P.pos_nsigma)
    allowed = (
        incons & elig[:, None] & elig[None] & (np.abs(ks[:, None] - ks[None]) <= sd.MAX_DK_NEW)
    )
    np.fill_diagonal(allowed, False)
    rng = np.random.default_rng(1)  # same split as s1_deform (pre-registered seed)
    perm = rng.permutation(np.flatnonzero(elig))
    cal, val = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2 :])
    return cat, B, elig, ks, sep, incons, allowed, cal, val


def background(B, allowed, idx, cpu):
    g = {"B": B, "allowed": allowed, "idx": idx}
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(g,)) as pool:
        rows = [r for chunk in pool.map(bg_row, list(idx), chunksize=4) for r in chunk]
    return np.array(rows, float)


def pair_chi2(B, i, j):
    bi, bj = B[int(i)], B[int(j)]
    return m2_only(bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"], True)


def run_injections(B, elig_idx, incons, thr_key, seed, per_cell, cpu):
    bf = {(s, e): sd.band_factors(s, e) for s in sd.S_INJ for e in sd.EPEAKS}
    rng = np.random.default_rng(seed)
    jobs = []
    for s0 in sd.S_INJ:
        for e in sd.EPEAKS:
            for ratio in sd.RATIOS:
                for _ in range(per_cell):
                    ia = int(rng.choice(elig_idx))
                    cand = elig_idx[incons[ia, elig_idx]]
                    jobs.append((seed + len(jobs), ia, int(rng.choice(cand)), s0, e, ratio))
    g = {"B": B, "bf": bf}
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(g,)) as pool:
        inj = np.array(list(pool.map(sd.inject_one, jobs, chunksize=8)), float)
    cols = ["a", "b", "s0", "epeak", "ratio", "dk", "elig", "dtd"]
    cols += ["m0", "m1", "m2", "s1", "s2", "p0", "p1", "p2"]
    return {c: inj[:, q] for q, c in enumerate(cols)}


# ------------------------------------------------------------------------------------- calibrate
def calibrate(a):
    t0 = time.time()
    cat, B, elig, ks, sep, incons, allowed, cal, val = setup()
    out = {"grid": list(sd.S_GRID), "far_per_pair": sd.FAR}
    # 1. per-pair threshold on the calibration half's real-pair background (M2 only)
    bg = background(B, allowed, cal, a.cpu)
    v = np.sort(np.nan_to_num(bg[:, 2], nan=-np.inf))[::-1]
    k = max(int(round(sd.FAR * bg.shape[0])) - 1, 0)
    thr = float(v[k])
    out["calibration_half"] = {"pairs": int(bg.shape[0]), "threshold_m2": thr}
    # 2. consistency test on simulated true pairs (independent noise; all eligible templates)
    inj = run_injections(B, np.flatnonzero(elig), incons, None, 777, a.per_cell, a.cpu)
    above = (inj["elig"] > 0) & (np.nan_to_num(inj["m2"], nan=-1) > thr)
    p = inj["p2"][above]
    st = inj["s0"][above]
    out["true_pairs"] = {
        "n_injected": int(inj["a"].size),
        "n_eligible": int((inj["elig"] > 0).sum()),
        "n_above_threshold": int(above.sum()),
        "frac_p_below_1e-3": float((p < 1e-3).mean()),
        "frac_p_below_1e-3_s0_eq_1": float((p[st == 1.0] < 1e-3).mean())
        if (st == 1.0).any()
        else None,
        "frac_p_below_1e-3_s0_ne_1": float((p[st != 1.0] < 1e-3).mean())
        if (st != 1.0).any()
        else None,
        "p_quantiles": {q: float(np.quantile(p, q)) for q in (0.005, 0.01, 0.05, 0.1, 0.5)},
        "frac_p_below_1e-3_by_ratio": {
            str(rt): float((p[inj["ratio"][above] == rt] < 1e-3).mean())
            for rt in sd.RATIOS
            if (inj["ratio"][above] == rt).any()
        },
        "ks_uniform_p": float(__import__("scipy").stats.kstest(p, "uniform").pvalue),
    }
    # declared rule (docs/s1_m2/calibration_plan.md): c = min(1e-3, q05 of refined true-pair p)
    c = float(min(1e-3, np.quantile(p, 0.05)))
    out["consistency_cut"] = c
    # 3. unrelated real pairs (calibration half): chi2 pass rate near and above the threshold
    order = np.argsort(np.nan_to_num(bg[:, 2], nan=-np.inf))[::-1]
    top = order[:2000]
    g = {"B": B}
    _init(g)
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(g,)) as pool:
        jobs = [(int(bg[q, 0]), int(bg[q, 1])) for q in top]
        res = list(pool.map(_chi2_job, jobs, chunksize=10))
    pv = np.array([r[2] for r in res])
    dtd = np.array([abs(B[int(bg[q, 0])]["mjd"] - B[int(bg[q, 1])]["mjd"]) for q in top])
    passed = (pv >= c) & (dtd >= 1.0)
    out["unrelated_pass_at_1e-3"] = int(((pv >= 1e-3) & (dtd >= 1.0)).sum())
    out["unrelated_real_pairs"] = {
        "n_top": int(top.size),
        "m2_range": [float(bg[top[-1], 2]), float(bg[top[0], 2])],
        "pass_consistency_and_delay": int(passed.sum()),
        "pass_top_100": int(passed[:100].sum()),
        "pass_rate_upper95": float(
            3.0 / top.size
            if passed.sum() == 0
            else (passed.sum() + 1.645 * math.sqrt(passed.sum())) / top.size
        ),
        "distinct_bursts_in_top": int(
            len(set(bg[top, 0].astype(int)) | set(bg[top, 1].astype(int)))
        ),
    }
    out["runtime_s"] = round(time.time() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def _chi2_job(ij):
    i, j = ij
    return pair_chi2(G["B"], i, j)


# ---------------------------------------------------------------------------------------- search
def search(a):
    t0 = time.time()
    fz = json.loads((OUT / "frozen.json").read_text())
    assert tuple(fz["grid"]) == tuple(sd.S_GRID)
    cat, B, elig, ks, sep, incons, allowed, cal, val = setup()
    n = len(B)
    iu = np.triu_indices(n, 1)
    e2 = elig[:, None] & elig[None]
    cover = {
        "catalogue_bursts": int(len(cat)),
        "bursts_reduced": n,
        "eligible_bursts_ge2_pulses": int(elig.sum()),
        "eligible_pairs": int(e2[iu].sum()),
        "excluded_position_consistent": int((e2 & ~incons)[iu].sum()),
        "excluded_dk_gt2": int(
            (e2 & incons & (np.abs(ks[:, None] - ks[None]) > sd.MAX_DK_NEW))[iu].sum()
        ),
        "pairs_scored": int(allowed[iu].sum()),
        "pairs_cal_cal": int(allowed[np.ix_(cal, cal)][np.triu_indices(cal.size, 1)].sum()),
        "pairs_val_val": int(allowed[np.ix_(val, val)][np.triu_indices(val.size, 1)].sum()),
    }
    cover["pairs_cross_half"] = (
        cover["pairs_scored"] - cover["pairs_cal_cal"] - cover["pairs_val_val"]
    )
    bg = background(B, allowed, np.flatnonzero(elig), a.cpu)
    sel = np.flatnonzero(np.nan_to_num(bg[:, 2], nan=-np.inf) > fz["threshold_m2"])
    g = {"B": B}
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(g,)) as pool:
        jobs = [(int(bg[q, 0]), int(bg[q, 1])) for q in sel]
        res = list(pool.map(_chi2_job, jobs, chunksize=4))
    flagged, cands = [], []
    for q, r in zip(sel, res, strict=True):
        i, j = int(bg[q, 0]), int(bg[q, 1])
        dtd = abs(B[i]["mjd"] - B[j]["mjd"])
        row = {
            "burst1": B[i]["name"],
            "burst2": B[j]["name"],
            "m2": round(float(bg[q, 2]), 4),
            "s_best": float(bg[q, 3]),
            "chi2_p": float(f"{r[2]:.3g}"),
            "delay_days": round(float(dtd), 3),
            "sep_deg": round(float(sep[i, j]), 1),
            "npulse": [int(B[i]["npulse"]), int(B[j]["npulse"])],
            "half": ["cal" if i in set(cal) else "val", "cal" if j in set(cal) else "val"],
        }
        flagged.append(row)
        if r[2] >= fz["consistency_cut"] and dtd >= fz["min_delay_days"]:
            cands.append(row)
    # de-duplicate: connected components of candidate pairs that share a burst
    comp, parent = [], {}

    def find(x):
        while parent.setdefault(x, x) != x:
            x = parent[x]
        return x

    for c in cands:
        parent[find(c["burst1"])] = find(c["burst2"])
    groups: dict = {}
    for c in cands:
        groups.setdefault(find(c["burst1"]), []).append(c)
    comp = [sorted({b for c in g_ for b in (c["burst1"], c["burst2"])}) for g_ in groups.values()]
    out = {
        "frozen": fz,
        "coverage": cover,
        "pairs_scored": int(bg.shape[0]),
        "flagged_above_threshold": len(flagged),
        "expected_flagged_at_far": round(sd.FAR * bg.shape[0], 1),
        "flagged": sorted(flagged, key=lambda r: -r["m2"]),
        "candidates": cands,
        "candidate_groups": comp,
        "runtime_s": round(time.time() - t0, 1),
    }
    (OUT / "search.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: out[k]
                for k in (
                    "coverage",
                    "flagged_above_threshold",
                    "expected_flagged_at_far",
                    "candidates",
                    "candidate_groups",
                )
            },
            indent=1,
        )
    )


# ------------------------------------------------------------------------------------ efficiency
def efficiency(a):
    fz = json.loads((OUT / "frozen.json").read_text())
    cat, B, elig, ks, sep, incons, allowed, cal, val = setup()
    inj = run_injections(
        B, np.flatnonzero(elig), incons, None, fz["efficiency_seed"], a.per_cell, a.cpu
    )
    el = inj["elig"] > 0
    rec = el & (np.nan_to_num(inj["m2"], nan=-1) > fz["threshold_m2"])
    rec &= (np.nan_to_num(inj["p2"]) >= fz["consistency_cut"]) & (
        inj["dtd"] >= fz["min_delay_days"]
    )
    rows = []
    for s0 in (None,) + sd.S_INJ:
        for ratio in (None,) + sd.RATIOS:
            m = np.ones(el.size, bool)
            if s0 is not None:
                m &= inj["s0"] == s0
            if ratio is not None:
                m &= inj["ratio"] == ratio
            n, ne, nr = int(m.sum()), int((m & el).sum()), int((m & rec).sum())
            rows.append(
                {
                    "s0": s0 if s0 is not None else "all",
                    "ratio": ratio if ratio is not None else "all",
                    "n": n,
                    "selected": ne,
                    "recovered": nr,
                    "selection_rate": round(ne / n, 3) if n else None,
                    "recovery_given_selected": round(nr / ne, 3) if ne else None,
                    "total_recovery": round(nr / n, 3) if n else None,
                }
            )
    out = {"frozen": fz, "n": int(el.size), "table": rows}
    (OUT / "efficiency.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps([r for r in rows if r["ratio"] == "all"], indent=0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=("calibrate", "search", "efficiency"))
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--per-cell", type=int, default=12)
    ap.add_argument("--out", default="calibration.json")
    a = ap.parse_args(argv)
    {"calibrate": calibrate, "search": search, "efficiency": efficiency}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
