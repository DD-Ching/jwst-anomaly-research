"""S1 deformed-twin recovery test (pre-registered in docs/s1_deform/preregistration.md).

Does the S1 screen (scripts/s1_twins.py, D-071) miss same-source bursts that arrive time-stretched
by s with the matching energy shift E -> E/s? Compares the existing statistic (M0) with two minimal
changes (M1: stretch grid in the flagging statistic; M2: M1 plus a free per-band scale) at the same
per-pair false-alarm rate, calibrated on the real-pair background of one half of the eligible
bursts; injections and the hold-out trial use the other half.

  python scripts/s1_deform.py [--cpu 4] [--per-cell 12] [--grid-step 1.05 --out ...]

Injected twins are simulated; recovering them measures search capability only, never evidence of
new physics.
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
from scipy import stats
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import s1_twins as st  # noqa: E402

from jwst_anomaly import burst_twins as bt  # noqa: E402

OUT = ROOT / "results" / "s1_deform"
P = bt.Params()
S_GRID_PREREG = (0.5, 0.63, 0.79, 1.0, 1.26, 1.59, 2.0)  # M1/M2 search grid (pre-registered)
S_GRID = S_GRID_PREREG
S_INJ = (0.5, 0.6, 0.71, 0.85, 1.0, 1.2, 1.41, 1.7, 2.0)  # injected relative factors
EPEAKS = (100.0, 200.0, 400.0)  # keV, ASSUMPTION (typical GBM values)
RATIOS = (1.0, 0.5, 0.3, 0.2, 0.1)
FAR = 1e-4  # per calibration pair, every method
MAX_DK_NEW = 2
G: dict = {}


# ---- physics
def band_photon(e, epeak, alpha=-1.0, beta=-2.3):
    """Band function photon spectrum (arbitrary normalisation)."""
    e0 = epeak / (2 + alpha)
    eb = (alpha - beta) * e0
    if e < eb:
        return (e / 100) ** alpha * math.exp(-e / e0)
    return (
        ((alpha - beta) * e0 / 100) ** (alpha - beta) * math.exp(beta - alpha) * (e / 100) ** beta
    )


def band_factors(s: float, epeak: float) -> np.ndarray:
    """Per-band photon-flux factors of a copy whose energies are E/s and whose times are s t.

    The copy's photon spectrum at observed energy E is the source spectrum at s E (a common 1/s
    rate dilution is irrelevant because brightness is set separately by the ratio). Returned
    relative to s = 1 for the two S1 bands (50-300 keV; 10-50 + 300-1000 keV).
    """

    def flux(s_):
        f = lambda e: band_photon(s_ * e, epeak)  # noqa: E731
        b = quad(f, 50, 300, limit=200)[0]
        r = quad(f, 10, 50, limit=200)[0] + quad(f, 300, 1000, limit=200)[0]
        return np.array([b, r])

    return flux(s) / flux(1.0)


# ---- statistics
def to_common(x, ex, kx, y, ey, ky):
    if kx < ky:
        x, ex = bt.rebin_factor(x, ex, 2 ** (ky - kx))
    elif ky < kx:
        y, ey = bt.rebin_factor(y, ey, 2 ** (kx - ky))
    return x, ex, y, ey


def rho_free(x, y):
    """Max over lag of sqrt(sum_b max(c_b, 0)^2 / ||x_b||^2) / ||y|| (y_b ~ a_b x_b, a_b >= 0)."""
    nx, ny = x.shape[1], y.shape[1]
    L = 1 << int(np.ceil(np.log2(nx + ny)))
    fx = np.fft.rfft(x, L, axis=1)
    fy = np.fft.rfft(y, L, axis=1)
    cc = np.fft.irfft(np.conj(fx) * fy, L, axis=1)  # per band
    nxb = np.maximum((x * x).sum(1), 1e-30)[:, None]
    ny2 = max((y * y).sum(), 1e-30)
    r2 = (np.maximum(cc, 0) ** 2 / nxb).sum(0) / ny2
    k = int(np.argmax(r2))
    lag = k if k < L - nx else k - L
    return float(math.sqrt(r2[k])), lag


def chi2_free(x, ex, y, ey, lag, sys_frac=P.chi2_sys_frac):
    """twin_chi2 with a separate scale per band (M2's vetting)."""
    nx, ny = x.shape[1], y.shape[1]
    lo, hi = max(0, lag), min(ny, nx + lag)
    if hi - lo < 3:
        return 0.0
    ys, eys = y[:, lo:hi], ey[:, lo:hi]
    xs, exs = x[:, lo - lag : hi - lag], ex[:, lo - lag : hi - lag]
    sel = (xs.sum(0) > 2 * np.sqrt((exs**2).sum(0))) | (ys.sum(0) > 2 * np.sqrt((eys**2).sum(0)))
    if sel.sum() < 3:
        return 0.0
    xs, exs, ys, eys = xs[:, sel], exs[:, sel], ys[:, sel], eys[:, sel]
    chi2 = 0.0
    for b in range(xs.shape[0]):
        a = (xs[b] * ys[b]).sum() / max((xs[b] ** 2).sum(), 1e-30)
        for _ in range(3):
            w = 1.0 / (eys[b] ** 2 + a * a * exs[b] ** 2 + (sys_frac * a * xs[b]) ** 2)
            a = (w * xs[b] * ys[b]).sum() / max((w * xs[b] ** 2).sum(), 1e-30)
        w = 1.0 / (eys[b] ** 2 + a * a * exs[b] ** 2 + (sys_frac * a * xs[b]) ** 2)
        chi2 += float((w * (ys[b] - a * xs[b]) ** 2).sum())
    dof = xs.size - xs.shape[0]
    return float(stats.chi2.sf(chi2, dof))


def scores(x, ex, kx, y, ey, ky, with_chi2=False):
    """(m0, m1, m2, s1_best, s2_best[, p0, p1, p2]).

    m0 is nan when |dk| > 1 (as in S1); everything is nan when |dk| > 2.
    """
    dk = abs(kx - ky)
    if dk > MAX_DK_NEW:
        return (np.nan,) * (8 if with_chi2 else 5)
    x, ex, y, ey = to_common(x, ex, kx, y, ey, ky)
    m0, lag0 = bt.xcorr_max(x, y)
    m1 = m2 = -1.0
    s1b = s2b = 1.0
    lag1 = lag2 = 0
    for s in S_GRID:
        ys, eys = (y, ey) if s == 1.0 else bt.stretch(y, ey, s)
        r1, l1 = bt.xcorr_max(x, ys)
        r2, l2 = rho_free(x, ys)
        if r1 > m1:
            m1, s1b, lag1 = r1, s, l1
        if r2 > m2:
            m2, s2b, lag2 = r2, s, l2
    if dk > P.max_dk:
        m0 = np.nan
    if not with_chi2:
        return m0, m1, m2, s1b, s2b
    p0 = bt.twin_chi2(x, ex, y, ey, lag0)[2] if np.isfinite(m0) else np.nan
    ys, eys = (y, ey) if s1b == 1.0 else bt.stretch(y, ey, s1b)
    p1 = bt.twin_chi2(x, ex, ys, eys, lag1)[2]
    ys, eys = (y, ey) if s2b == 1.0 else bt.stretch(y, ey, s2b)
    p2 = chi2_free(x, ex, ys, eys, lag2)
    return m0, m1, m2, s1b, s2b, p0, p1, p2


# ---- background
def _init(g):
    os.environ["OMP_NUM_THREADS"] = "1"
    G.update(g)


def bg_row(i):
    B, allowed, idx = G["B"], G["allowed"], G["idx"]
    out = []
    for j in idx:
        if j <= i or not allowed[i, j]:
            continue
        r = scores(B[i]["flux"], B[i]["errs"], B[i]["k"], B[j]["flux"], B[j]["errs"], B[j]["k"])
        out.append((i, j, *r))
    return out


# ---- injection
def inject_one(args):
    """Deformed copy of burst a's template in burst b's noise slot, scored by all three methods.

    Follows st.injection_one (window/k emulation from a T90 proxy). The copy is stretched by s0 in
    time and its band fluxes are multiplied by band_factors(s0, E_peak); overall brightness ratio
    and gain jitter are as in D-071.
    """
    seed, a, b, s0, epeak, ratio = args
    rng = np.random.default_rng(seed)
    B = G["B"]
    A = B[a]
    fa, ea = A["flux"], A["errs"]
    na = fa.shape[1]
    templ = bt.template(fa)
    # pad so a stretched copy fits, stretch, apply the band factors and the brightness ratio
    pad = int(np.ceil(na * max(s0, 1.0))) - na + 2
    tp = np.pad(templ, ((0, 0), (pad // 2, pad - pad // 2)))
    ts, _ = bt.stretch(tp, np.ones_like(tp), s0)
    ts = ts * G["bf"][(s0, epeak)][:, None]
    peak_t = templ.sum(0).max()
    ts = ts * ratio * peak_t / max(ts.sum(0).max(), 1e-30)
    gain = 1.0 + P.inject_gain_jitter * rng.normal(size=(ts.shape[0], 1))
    ts = ts * gain
    n_c = ts.shape[1]
    # noise: b's background (a's bin width) + Poisson-like source variance scaling with the flux
    eb = B[b]["errs"] * np.sqrt(B[b]["dt"] / A["dt"])
    eb = eb[:, np.arange(n_c) % eb.shape[1]]
    u_a = np.percentile(ea, 20, axis=1, keepdims=True)
    u_b = np.percentile(eb, 20, axis=1, keepdims=True)
    src_var = np.maximum(ea**2 - u_a**2, 0)
    bright = templ > 0
    coef = np.array(
        [
            np.median(src_var[q][bright[q]] / templ[q][bright[q]]) if bright[q].any() else 0.0
            for q in range(2)
        ]
    )[:, None]
    var_c = u_b**2 + coef * np.maximum(ts, 0)
    fi = ts + rng.normal(size=ts.shape) * np.sqrt(var_c)
    ei = np.sqrt(var_c)
    fa2 = templ + rng.normal(size=templ.shape) * ea
    # copy's emulated window / k (as st.injection_one), positions relative to the padded grid
    t90a = max(A["t90"], P.base_dt)
    ia0, ia1 = bt.duration_bins(fa2)
    ic0, ic1 = bt.duration_bins(fi)
    t90c = t90a * max(ic1 - ic0 + 1, 1) / max(ia1 - ia0 + 1, 1)
    lo_c, hi_c, dt_c = bt.window_for(0.0, t90c, P)
    k_raw = int(round(np.log2(dt_c / P.base_dt)))
    k_c = max(k_raw, A["k"])
    span = int(round((hi_c - lo_c) / A["dt"]))
    c0 = max(0, ic0 - int(round((0.0 - lo_c) / A["dt"])))
    c1 = min(n_c, max(c0 + 3, c0 + span))
    yc, eyc = bt.rebin_factor(fi[:, c0:c1], ei[:, c0:c1], 2 ** (k_c - A["k"]))
    dk = abs(k_raw - A["k"])
    npk_c = len(bt.find_pulses(yc, eyc, P)) if yc.shape[1] >= 3 else 0
    npk_a = len(bt.find_pulses(fa2, ea, P))
    elig = npk_c >= P.min_pulses and npk_a >= P.min_pulses
    sc = scores(fa2, ea, A["k"], yc, eyc, k_c, with_chi2=True)
    dtd = abs(A["mjd"] - B[b]["mjd"])
    return (a, b, s0, epeak, ratio, int(dk), bool(elig), float(dtd), *[float(v) for v in sc])


# ---- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--per-cell", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument(
        "--grid-step", type=float, default=None, help="post hoc: log-uniform s grid with this ratio"
    )
    ap.add_argument(
        "--skip-holdout", action="store_true", help="post hoc runs: injections and calibration only"
    )
    ap.add_argument("--out", default="summary.json")
    a = ap.parse_args(argv)
    global S_GRID
    if a.grid_step:
        n = int(round(math.log(4.0) / math.log(a.grid_step))) // 2
        S_GRID = tuple(float(round(a.grid_step**q, 4)) for q in range(-n, n + 1))
    t0 = time.time()
    cat, B, status = st.load()
    elig = np.array([b["npulse"] >= P.min_pulses for b in B])
    ra, dec, err = (np.array([b[c] for b in B]) for c in ("ra", "dec", "err"))
    sep = bt.angsep_deg(ra[:, None], dec[:, None], ra[None], dec[None])
    incons = bt.position_inconsistent(sep, err[:, None], err[None], P.sys_deg, P.pos_nsigma)
    allowed = incons & elig[:, None] & elig[None]
    np.fill_diagonal(allowed, False)
    rng = np.random.default_rng(1)  # pre-registered split seed
    eidx = np.flatnonzero(elig)
    perm = rng.permutation(eidx)
    cal, val = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2 :])
    out = {"n_eligible": int(elig.sum()), "n_cal": int(cal.size), "n_val": int(val.size)}

    def background(idx):
        g = {"B": B, "allowed": allowed, "idx": idx}
        with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(g,)) as pool:
            rows = [r for chunk in pool.map(bg_row, list(idx), chunksize=4) for r in chunk]
        return np.array(rows, float)

    # ---- calibration half: thresholds at equal per-pair FAR
    t1 = time.time()
    bg = background(cal)
    n_pairs = bg.shape[0]
    thr = {}
    for name, col in (("M0", 2), ("M1", 3), ("M2", 4)):
        v = np.nan_to_num(bg[:, col], nan=-np.inf)
        thr[name] = float(np.sort(v)[::-1][max(int(round(FAR * n_pairs)) - 1, 0)])
    out["calibration"] = {
        "pairs": n_pairs,
        "pairs_scored_M0": int(np.isfinite(bg[:, 2]).sum()),
        "far_per_pair": FAR,
        "thresholds": thr,
        "runtime_s": round(time.time() - t1, 1),
    }
    print("calibration", out["calibration"], flush=True)

    # ---- injections on the validation half
    bf = {(s, e): band_factors(s, e) for s in S_INJ for e in EPEAKS}
    out["band_factors"] = {
        f"s={s},Ep={e}": [round(float(x), 4) for x in v] for (s, e), v in bf.items()
    }
    jobs = []
    rng = np.random.default_rng(a.seed)
    for s0 in S_INJ:
        for e in EPEAKS:
            for ratio in RATIOS:
                for _ in range(a.per_cell):
                    ia = int(rng.choice(val))
                    cand = val[incons[ia, val]]
                    ib = int(rng.choice(cand))
                    jobs.append((a.seed + len(jobs), ia, ib, s0, e, ratio))
    g = {"B": B, "bf": bf}
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(g,)) as pool:
        inj = np.array(list(pool.map(inject_one, jobs, chunksize=8)), float)
    cols = [
        "a",
        "b",
        "s0",
        "epeak",
        "ratio",
        "dk",
        "elig",
        "dtd",
        "m0",
        "m1",
        "m2",
        "s1",
        "s2",
        "p0",
        "p1",
        "p2",
    ]
    inj_c = {c: inj[:, q] for q, c in enumerate(cols)}

    def recovered(meth):
        m = {"M0": "m0", "M1": "m1", "M2": "m2"}[meth]
        p = {"M0": "p0", "M1": "p1", "M2": "p2"}[meth]
        sc = np.nan_to_num(inj_c[m], nan=-np.inf)
        return (
            (inj_c["elig"] > 0)
            & (sc > thr[meth])
            & (np.nan_to_num(inj_c[p]) >= P.chi2_p_min)
            & (inj_c["dtd"] >= 1.0)
        )

    rec = {m: recovered(m) for m in ("M0", "M1", "M2")}
    table = []
    for s0 in S_INJ:
        for ratio in (None,) + RATIOS:
            sel = inj_c["s0"] == s0
            if ratio is not None:
                sel &= inj_c["ratio"] == ratio
            n = int(sel.sum())
            row = {
                "s0": s0,
                "ratio": ratio if ratio is not None else "all",
                "n": n,
                "eligible": round(float(inj_c["elig"][sel].mean()), 3),
            }
            for m in ("M0", "M1", "M2"):
                k = int(rec[m][sel].sum())
                row[m] = k
            table.append(row)
    out["injections"] = {"n": int(inj.shape[0]), "per_cell": a.per_cell, "table": table}
    # stage diagnosis for M0 on deformed injections (ratio 1): where are they lost?
    diag = {}
    for s0 in S_INJ:
        sel = (inj_c["s0"] == s0) & (inj_c["ratio"] == 1.0)
        diag[str(s0)] = {
            "n": int(sel.sum()),
            "lost_eligibility": int((sel & (inj_c["elig"] == 0)).sum()),
            "lost_dk_gt1": int((sel & (inj_c["elig"] > 0) & ~np.isfinite(inj_c["m0"])).sum()),
            "below_M0_threshold": int(
                (
                    sel
                    & (inj_c["elig"] > 0)
                    & np.isfinite(inj_c["m0"])
                    & (inj_c["m0"] <= thr["M0"])
                ).sum()
            ),
            "fail_chi2_M0": int(
                (
                    sel
                    & (inj_c["elig"] > 0)
                    & (np.nan_to_num(inj_c["m0"], nan=-1) > thr["M0"])
                    & (np.nan_to_num(inj_c["p0"]) < P.chi2_p_min)
                ).sum()
            ),
            "recovered_M0": int((sel & rec["M0"]).sum()),
            "median_m0": float(np.nanmedian(inj_c["m0"][sel & (inj_c["elig"] > 0)]))
            if (sel & (inj_c["elig"] > 0)).any()
            else None,
        }
    out["m0_stage_diagnosis_ratio1"] = diag
    np.save(OUT / "injections.npy", inj) if os.environ.get("S1D_SAVE_NPY") else None

    # ---- false pairs in the calibration half that pass each method's full chain (chi2 at best s)
    def pair_scores(i, j):
        bi, bj = B[i], B[j]
        return scores(bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"], True)

    def chain_false(idx_bg):
        n_flag, n_pass = {}, {}
        for name, col in (("M0", 2), ("M1", 3), ("M2", 4)):
            sel = np.flatnonzero(np.nan_to_num(idx_bg[:, col], nan=-np.inf) > thr[name])
            n_flag[name] = int(sel.size)
            k = 0
            for q in sel:
                i, j = int(idx_bg[q, 0]), int(idx_bg[q, 1])
                sc = scores(
                    B[i]["flux"],
                    B[i]["errs"],
                    B[i]["k"],
                    B[j]["flux"],
                    B[j]["errs"],
                    B[j]["k"],
                    True,
                )
                pv = {"M0": sc[5], "M1": sc[6], "M2": sc[7]}[name]
                k += int(pv >= P.chi2_p_min and abs(B[i]["mjd"] - B[j]["mjd"]) >= 1.0)
            n_pass[name] = k
        return n_flag, n_pass

    (
        out["calibration"]["flagged_by_threshold"],
        out["calibration"]["false_pairs_passing_chain"],
    ) = chain_false(bg)
    if a.skip_holdout:
        out["runtime_s"] = round(time.time() - t0, 1)
        out["s_grid"] = list(S_GRID)
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / a.out).write_text(json.dumps(out, indent=1, default=float) + "\n")
        return 0

    # ---- hold-out trial: real validation-half pairs, fixed methods
    t1 = time.time()
    bv = background(val)
    flags = {}
    cand = []
    for name, col in (("M0", 2), ("M1", 3), ("M2", 4)):
        sel = np.flatnonzero(np.nan_to_num(bv[:, col], nan=-np.inf) > thr[name])
        flags[name] = int(sel.size)
        for q in sel:
            i, j = int(bv[q, 0]), int(bv[q, 1])
            sc = scores(
                B[i]["flux"], B[i]["errs"], B[i]["k"], B[j]["flux"], B[j]["errs"], B[j]["k"], True
            )
            p = {"M0": sc[5], "M1": sc[6], "M2": sc[7]}[name]
            dtd = abs(B[i]["mjd"] - B[j]["mjd"])
            if p >= P.chi2_p_min and dtd >= 1.0:
                cand.append(
                    {
                        "method": name,
                        "burst1": B[i]["name"],
                        "burst2": B[j]["name"],
                        "score": round(float(bv[q, col]), 4),
                        "s_best": float(sc[3] if name != "M2" else sc[4]),
                        "chi2_p": float(f"{p:.3g}"),
                        "delay_days": round(float(dtd), 2),
                        "sep_deg": round(float(sep[i, j]), 1),
                    }
                )
    out["holdout"] = {
        "pairs": int(bv.shape[0]),
        "flagged_by_threshold": flags,
        "chi2_and_delay_pass": cand,
        "runtime_s": round(time.time() - t1, 1),
    }
    out["runtime_s"] = round(time.time() - t0, 1)
    out["provenance"] = (
        "derived (real GBM bcat light curves from results/s1_twins) and simulated (injections); "
        "scripts/s1_deform.py"
    )
    OUT.mkdir(parents=True, exist_ok=True)
    out["s_grid"] = list(S_GRID)
    (OUT / a.out).write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in ("calibration", "holdout")}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
