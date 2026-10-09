"""S1 burst twins, stage 2: screen, surrogate null, injection-recovery, vetting and limit.

Reads ``results/s1_twins/lc_*.ecsv.gz`` and ``catalogue.ecsv.gz`` (from scripts/s1_ingest.py).
Writes
``bursts.ecsv.gz`` (per-burst eligibility), ``pairs_top.ecsv`` (top real pairs and every flagged
pair),
``null.json`` (surrogate-null summary), ``injections.ecsv.gz`` and ``summary.json``. Plots go to
--plots.

  python scripts/s1_twins.py [--null 100] [--inject 5000] [--inject-sens 1000] [--plots <dir>]

Every number is a screen statistic; a match is an anomaly, never evidence of new physics.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import glob
import io
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from jwst_anomaly import burst_twins as bt  # noqa: E402
from jwst_anomaly.schema import Provenance  # noqa: E402

OUT = ROOT / "results" / "s1_twins"
P = bt.Params()
# FFT length: native (<= 256 bins) + coarsened (<= 128), or two natives, fit without wrap
L = 512
G: dict = {}  # pool globals (fork) / set by initializer


def write_ecsv_gz(t: Table, path: Path) -> None:
    import gzip

    buf = io.StringIO()
    t.write(buf, format="ascii.ecsv")
    with open(path, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
        gz.write(buf.getvalue().encode())


def write_json(obj, path: Path) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float) + "\n")


# -------------------------------------------------------------------------------------------------
# load


def load():
    cat = Table.read(OUT / "catalogue.ecsv.gz", format="ascii.ecsv")
    lc = vstack(
        [Table.read(f, format="ascii.ecsv") for f in sorted(glob.glob(str(OUT / "lc_*.ecsv.gz")))]
    )
    status = {str(r["trigger_name"]): str(r["status"]) for r in lc}
    lc = lc[lc["status"] == "ok"]
    idx = {str(n): i for i, n in enumerate(cat["trigger_name"])}
    bursts = []
    for r in lc:
        name = str(r["trigger_name"])
        c = cat[idx[name]]
        flux, err = bt.decode((r["u1"], r["u2"]), (r["f1"], r["f2"]), (r["e1"], r["e2"]))
        valid = np.isfinite(flux).all(0)
        f, e = bt.clean(flux, err)
        npk = len(bt.find_pulses(f, e, P)) if valid.sum() >= 3 else 0
        bursts.append(
            {
                "name": name,
                "ra": float(c["ra"]),
                "dec": float(c["dec"]),
                "err": float(c["error_radius"]),
                "mjd": float(c["trigger_time"]),
                "t90": float(c["t90"]),
                "fluence": float(c["fluence"]),
                "t_lo": float(r["t_lo"]),
                "dt": float(r["dt"]),
                "k": int(round(np.log2(float(r["dt"]) / P.base_dt))),
                "flux": f,
                "errs": e,
                "valid_frac": float(valid.mean()),
                "npulse": npk,
                "snr": bt.snr_total(f, e),
            }
        )
    return cat, bursts, status


# ------------------------------------------------------------------------------------------------
# score


def level_bank(lcs, ks, K):
    """Rows usable at level K: natives at K, level K-1 coarsened by 2. (ids, coarse, F, norm, n)."""
    ids, coarse, F, nrm, n = [], [], [], [], []
    for i, (fe, k) in enumerate(zip(lcs, ks, strict=True)):
        if k == K:
            x = fe[0]
            co = False
        elif k == K - 1 and P.max_dk >= 1:
            x, _ = bt.rebin_factor(fe[0], fe[1], 2)
            co = True
        else:
            continue
        ids.append(i)
        coarse.append(co)
        F.append(np.fft.rfft(x, L, axis=1))
        nrm.append(np.sqrt((x * x).sum()))
        n.append(x.shape[1])
    if not ids:
        return None
    return np.array(ids), np.array(coarse), np.array(F), np.array(nrm), np.array(n)


def score_all(lcs, ks, allowed):
    """Max multi-band normalised cross-correlation over lag, s = 1, for every allowed pair.

    ``allowed`` is a boolean (N, N) matrix (position cut and eligibility). Returns dict (i, j) ->
    (rho, lag, K).
    """
    out_i, out_j, out_r, out_l, out_k = [], [], [], [], []
    for K in sorted(set(ks)):
        bank = level_bank(lcs, ks, K)
        if bank is None:
            continue
        ids, coarse, F, nrm, n = bank
        for a in range(len(ids)):
            b = np.arange(a + 1, len(ids))
            b = b[~(coarse[a] & coarse[b])]  # coarse-coarse pairs were scored natively at K-1
            b = b[allowed[ids[a], ids[b]]]
            if b.size == 0 or nrm[a] <= 0:
                continue
            cc = np.fft.irfft((np.conj(F[a])[None] * F[b]).sum(1), L, axis=1)
            km = cc.argmax(1)
            rho = cc[np.arange(b.size), km] / np.maximum(nrm[a] * nrm[b], 1e-30)
            lag = np.where(km < L - n[a], km, km - L)
            ia, ib = ids[a], ids[b]
            out_i.append(np.full(b.size, ia))
            out_j.append(ib)
            out_r.append(rho)
            out_l.append(lag)
            out_k.append(np.full(b.size, K))
    if not out_i:
        e = np.array([], int)
        return e, e, np.array([]), e, e
    i, j = np.concatenate(out_i), np.concatenate(out_j)
    sw = i > j
    i, j = np.where(sw, j, i), np.where(sw, i, j)
    return i, j, np.concatenate(out_r), np.concatenate(out_l), np.concatenate(out_k)


def best_per_pair(i, j, r, lag, K):
    """Keep one score per pair (the best, if a pair was scored at two levels)."""
    key = i.astype(np.int64) * 100000 + j
    order = np.lexsort((-r, key))
    keep = np.r_[True, key[order][1:] != key[order][:-1]]
    o = order[keep]
    return i[o], j[o], r[o], lag[o], K[o]


def stretch_scores(lcs, ks, pairs, s_grid):
    """Secondary: max rho over the s grid (j stretched by s, at the coarser level of the pair)."""
    best = np.zeros(len(pairs))
    best_s = np.ones(len(pairs))
    for q, (i, j) in enumerate(pairs):
        xi, ei = lcs[i]
        xj, ej = lcs[j]
        if ks[i] < ks[j]:
            xi, ei = bt.rebin_factor(xi, ei, 2 ** (ks[j] - ks[i]))
        elif ks[j] < ks[i]:
            xj, ej = bt.rebin_factor(xj, ej, 2 ** (ks[i] - ks[j]))
        for s in s_grid:
            ys, _ = bt.stretch(xj, ej, s)
            r, _ = bt.xcorr_max(xi, ys)
            if r > best[q]:
                best[q], best_s[q] = r, s
    return best, best_s


# -------------------------------------------------------------------------------------------------
# null


def _init(g):
    os.environ["OMP_NUM_THREADS"] = "1"
    G.update(g)


def null_realisation(seed):
    rng = np.random.default_rng(seed)
    lcs = [bt.pulse_shuffle(f, e, rng, P) for f, e in G["lcs"]]
    i, j, r, lag, K = score_all(lcs, G["ks"], G["allowed"])
    i, j, r, lag, K = best_per_pair(i, j, r, lag, K)
    hist, _ = np.histogram(r, bins=G["bins"])
    top = np.sort(r)[-20:][::-1]
    return float(r.max()) if r.size else 0.0, hist, top, int(r.size)


# --------------------------------------------------------------------------------------------
# injection


def pair_rho(x, ex, kx, y, ey, ky):
    """s = 1 rho/lag between two light curves, at the coarser of the two levels."""
    if kx < ky:
        x, ex = bt.rebin_factor(x, ex, 2 ** (ky - kx))
    elif ky < kx:
        y, ey = bt.rebin_factor(y, ey, 2 ** (kx - ky))
    r, lag = bt.xcorr_max(x, y)
    return r, lag, x, ex, y, ey


def injection_one(args):
    """One synthetic twin through the chain: eligibility, threshold, chi2, re-trigger veto.

    The copy gets its own stored window and resolution k, emulated from its T90: A's catalogue T90
    scaled by the ratio of the cumulative-fluence durations (copy / A'), then ``bt.window_for``.
    A copy whose k differs
    from A's by more than ``max_dk`` is lost (not scored), as a real pair would be.
    """
    seed, a, b, ratio, jitter = args
    rng = np.random.default_rng(seed)
    B = G["bursts"]
    A = B[a]
    fa, ea = A["flux"], A["errs"]
    eb = B[b]["errs"] * np.sqrt(B[b]["dt"] / A["dt"])  # B's noise at A's bin width
    nb = eb.shape[1]
    na = fa.shape[1]
    eb = eb[:, np.arange(na) % nb] if nb < na else eb[:, :na]
    fa2, ea2, fi, ei = bt.inject_twin(bt.template(fa), ea, eb, ratio, rng, jitter)
    # --- emulate the copy's catalogue T90, stored window and k
    t90a = max(A["t90"], P.base_dt)
    t90s_a = A["t_lo"] + P.pad_frac * t90a + P.pad_s  # window_for: t_lo = t90_start - pad
    ia0, ia1 = bt.duration_bins(fa2)
    ic0, ic1 = bt.duration_bins(fi)
    t90c = t90a * max(ic1 - ic0 + 1, 1) / max(ia1 - ia0 + 1, 1)
    t90s_c = t90s_a + (ic0 - ia0) * A["dt"]
    lo_c, hi_c, dt_c = bt.window_for(t90s_c, t90c, P)
    k_raw = int(round(np.log2(dt_c / P.base_dt)))
    lost_k = abs(k_raw - A["k"]) > P.max_dk
    k_c = max(k_raw, A["k"])  # a finer copy is compared at A's level anyway
    j0 = int(np.clip(np.floor((lo_c - A["t_lo"]) / A["dt"]), 0, na - 1))
    j1 = int(np.clip(np.ceil((hi_c - A["t_lo"]) / A["dt"]), j0 + 1, na))
    yc, eyc = bt.rebin_factor(fi[:, j0:j1], ei[:, j0:j1], 2 ** (k_c - A["k"]))
    if yc.shape[1] < 3:
        yc, eyc, k_c, lost_k = fi, ei, A["k"], True
    npk = len(bt.find_pulses(yc, eyc, P))
    eligible = (
        (not lost_k) and npk >= P.min_pulses and len(bt.find_pulses(fa2, ea2, P)) >= P.min_pulses
    )
    r, lag, x, ex, y, ey = pair_rho(fa2, ea2, A["k"], yc, eyc, k_c)
    chi2, dof, pval = bt.twin_chi2(x, ex, y, ey, lag)
    dtd = abs(A["mjd"] - B[b]["mjd"])
    flagged = eligible and r > G["rho_star"]
    passed = flagged and pval >= P.chi2_p_min and dtd >= P.retrigger_days
    passed_deep = (
        eligible and r > P.rho_prescreen_deep and pval >= P.chi2_p_min and dtd >= P.retrigger_days
    )
    return (
        a,
        b,
        ratio,
        jitter,
        bt.snr_total(fi, ei),
        npk,
        int(k_raw - A["k"]),
        bool(lost_k),
        bool(eligible),
        float(r),
        float(pval),
        bool(flagged),
        bool(passed),
        bool(passed_deep),
    )


# ----------------------------------------------------------------------------------------------
# plotting


def overlay(ax, B, i, j, lag, K, title):
    xi, ei = B[i]["flux"], B[i]["errs"]
    xj, ej = B[j]["flux"], B[j]["errs"]
    r, lag, x, ex, y, ey = pair_rho(xi, ei, B[i]["k"], xj, ej, B[j]["k"])
    dt = P.base_dt * 2 ** max(B[i]["k"], B[j]["k"])
    tx = np.arange(x.shape[1]) * dt
    ty = (np.arange(y.shape[1]) - lag) * dt
    sx, sy = x.sum(0), y.sum(0)
    a = (sx * np.interp(tx, ty, sy, left=0, right=0)).sum() / max((sx * sx).sum(), 1e-30)
    ax.step(tx, a * sx, where="post", lw=0.8, label=B[i]["name"])
    ax.step(ty, sy, where="post", lw=0.8, alpha=0.8, label=B[j]["name"])
    ax.set_title(title, fontsize=7)
    ax.legend(fontsize=5, loc="upper right")
    ax.tick_params(labelsize=5)


def contact_sheet(B, rows, path, ncol=4):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = max(1, len(rows))
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.2 * nrow), squeeze=False)
    for ax, r in zip(axes.flat, rows, strict=False):
        overlay(ax, B, r["i"], r["j"], r["lag"], r["K"], r["title"])
    for ax in axes.flat[len(rows) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# main


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--null", type=int, default=100, help="surrogate catalogues")
    ap.add_argument(
        "--inject", type=int, default=5000, help="injections at the primary gain jitter"
    )
    ap.add_argument(
        "--inject-sens", type=int, default=1000, help="injections per sensitivity jitter"
    )
    ap.add_argument("--jitter-grid", type=float, nargs="*", default=[0.05, 0.10, 0.20])
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--plots", type=Path, default=None)
    ap.add_argument("--seed", type=int, default=20261009)
    a = ap.parse_args(argv)
    t0 = time.time()
    cat, B, status = load()
    N = len(B)
    ks = np.array([b["k"] for b in B])
    elig = np.array([b["npulse"] >= P.min_pulses for b in B])
    ra, dec, err = (np.array([b[c] for b in B]) for c in ("ra", "dec", "err"))
    sep = bt.angsep_deg(ra[:, None], dec[:, None], ra[None], dec[None])
    incons = bt.position_inconsistent(sep, err[:, None], err[None], P.sys_deg, P.pos_nsigma)
    incons_tail = bt.position_inconsistent(
        sep, err[:, None], err[None], P.sys_tail_deg, P.pos_nsigma
    )
    allowed = incons & elig[:, None] & elig[None] & (np.abs(ks[:, None] - ks[None]) <= P.max_dk)
    np.fill_diagonal(allowed, False)
    n_pairs_allowed = int(np.triu(allowed, 1).sum())
    n_elig_pairs = int(elig.sum() * (elig.sum() - 1) // 2)
    print(
        f"bursts: catalogue {len(cat)}, bcat reduced {N}, eligible (>= {P.min_pulses} pulses) "
        f"{elig.sum()}; "
        f"eligible pairs {n_elig_pairs}, after position + resolution cut {n_pairs_allowed}",
        flush=True,
    )

    bt_tab = Table(
        {
            "trigger_name": [b["name"] for b in B],
            "dt": [b["dt"] for b in B],
            "nbin": [b["flux"].shape[1] for b in B],
            "valid_frac": np.round([b["valid_frac"] for b in B], 3),
            "npulse": [b["npulse"] for b in B],
            "snr_peak": np.round([b["snr"] for b in B], 2),
            "eligible": elig,
        }
    )
    bt_tab.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": "results/s1_twins/lc_*.ecsv.gz via jwst_anomaly.burst_twins.find_pulses "
        f"(prominence {P.pulse_prominence_sigma} sigma, ASSUMPTION); scripts/s1_twins.py",
    }
    write_ecsv_gz(bt_tab, OUT / "bursts.ecsv.gz")

    lcs = [(b["flux"], b["errs"]) for b in B]
    # ---- real data, s = 1
    t1 = time.time()
    i, j, r, lag, K = best_per_pair(*score_all(lcs, ks, allowed))
    print(
        f"real: {r.size} pairs scored in {time.time() - t1:.0f} s; max rho {r.max():.4f}",
        flush=True,
    )

    # ---- surrogate null
    bins = np.linspace(-0.2, 1.0, 1201)
    g = {"lcs": lcs, "ks": ks, "allowed": allowed, "bins": bins, "bursts": B}
    G.update(g)
    t1 = time.time()
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(g,)) as pool:
        res = list(pool.map(null_realisation, [a.seed + 1000 + q for q in range(a.null)]))
    maxes = np.array([x[0] for x in res])
    hist = np.sum([x[1] for x in res], axis=0)
    tops = np.concatenate([x[2] for x in res])
    rho_star = float(
        np.quantile(maxes, 0.95)
    )  # family-wise 5 % over all pairs of one catalogue (ASSUMPTION)
    dt_null = time.time() - t1
    print(
        f"null: {a.null} surrogate catalogues in {dt_null:.0f} s ({a.null * r.size / dt_null:.3g} "
        f"pairs/s); "
        f"per-catalogue max rho median {np.median(maxes):.4f}, 95 % {rho_star:.4f}",
        flush=True,
    )
    cdf_real = np.array([(r > x).sum() for x in (0.8, 0.9, 0.95, rho_star)])
    cdf_null = np.array(
        [
            (np.cumsum(hist[::-1])[::-1][np.searchsorted(bins, x)]) / a.null
            for x in (0.8, 0.9, 0.95, rho_star)
        ]
    )
    real_hist, _ = np.histogram(r, bins=bins)

    flagged = np.flatnonzero(r > rho_star)
    print(f"flagged real pairs (rho > rho*): {flagged.size}", flush=True)

    # ---- secondary: stretch grid on the top real pairs (and flags)
    top_idx = np.argsort(r)[::-1][:300]
    sel = np.unique(np.r_[top_idx, flagged])
    s_best, s_arg = stretch_scores(lcs, ks, list(zip(i[sel], j[sel], strict=True)), P.stretch_grid)

    # ---- vetting chain for top/flagged rows
    rows = []
    for q, m in enumerate(sel):
        bi, bj = B[i[m]], B[j[m]]
        _, lg, x, ex, y, ey = pair_rho(
            bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"]
        )
        chi2, dof, pval = bt.twin_chi2(x, ex, y, ey, lg)
        p_s = pval
        if s_arg[q] != 1.0:  # chi2 of the stretched comparison (secondary s grid)
            ys, eys = bt.stretch(y, ey, float(s_arg[q]))
            _, lgs = bt.xcorr_max(x, ys)
            p_s = bt.twin_chi2(x, ex, ys, eys, lgs)[2]
        dtd = abs(bi["mjd"] - bj["mjd"])
        rows.append(
            {
                "burst1": bi["name"],
                "burst2": bj["name"],
                "sep_deg": round(float(sep[i[m], j[m]]), 2),
                "pos_nsigma": round(
                    float(
                        sep[i[m], j[m]]
                        / np.sqrt(err[i[m]] ** 2 + err[j[m]] ** 2 + 2 * P.sys_deg**2)
                    ),
                    2,
                ),
                "incons_tail": bool(incons_tail[i[m], j[m]]),
                "delay_days": round(dtd, 4),
                "dt_s": P.base_dt * 2 ** int(K[m]),
                "lag_bins": int(lag[m]),
                "rho": round(float(r[m]), 4),
                "rho_s": round(float(s_best[q]), 4),
                "s_best": float(s_arg[q]),
                "chi2": round(float(chi2), 1),
                "dof": dof,
                "chi2_p": float(f"{pval:.3g}"),
                "chi2_p_s": float(f"{p_s:.3g}"),
                "flag": bool(r[m] > rho_star),
                "v_chi2": bool(pval >= P.chi2_p_min),
                "v_retrigger": bool(dtd >= P.retrigger_days),
                "t90_1": bi["t90"],
                "t90_2": bj["t90"],
                "npulse_1": bi["npulse"],
                "npulse_2": bj["npulse"],
                "snr_1": round(bi["snr"], 1),
                "snr_2": round(bj["snr"], 1),
                "_i": int(i[m]),
                "_j": int(j[m]),
                "_K": int(K[m]),
            }
        )
    rows.sort(key=lambda d: -d["rho"])
    pt = Table([{k: v for k, v in d.items() if not k.startswith("_")} for d in rows])
    pt.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": "scripts/s1_twins.py over results/s1_twins/lc_*.ecsv.gz; "
        "top 300 real pairs by rho (s = 1) and "
        f"every pair above rho* = {rho_star:.4f} (95 % of the per-catalogue surrogate maximum, "
        f"ASSUMPTION)",
    }
    pt.write(OUT / "pairs_top.ecsv", format="ascii.ecsv", overwrite=True)
    survivors = [d for d in rows if d["flag"] and d["v_chi2"] and d["v_retrigger"]]
    print(
        f"automated survivors (flag, chi2 p >= {P.chi2_p_min}, delay >= {P.retrigger_days} d): "
        f"{len(survivors)}"
    )

    # ---- secondary deep chain: every real pair above rho_prescreen_deep gets the chi2 and re-
    # trigger tests
    deep_idx = np.flatnonzero(r > P.rho_prescreen_deep)
    deep_survivors = []
    deep_chi2_pass = 0
    for m in deep_idx:
        bi, bj = B[i[m]], B[j[m]]
        _, lg, x, ex, y, ey = pair_rho(
            bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"]
        )
        pval = bt.twin_chi2(x, ex, y, ey, lg)[2]
        if pval >= P.chi2_p_min:
            deep_chi2_pass += 1
            if abs(bi["mjd"] - bj["mjd"]) >= P.retrigger_days:
                deep_survivors.append(
                    (bi["name"], bj["name"], round(float(r[m]), 4), float(f"{pval:.3g}"))
                )
    print(
        f"deep chain: {deep_idx.size} pairs with rho > {P.rho_prescreen_deep}; chi2 pass "
        f"{deep_chi2_pass}; "
        f"survivors {len(deep_survivors)}",
        flush=True,
    )

    # ---- positive controls: same-source re-triggers (< 3 d apart, < 2 deg), scored with the
    # position cut ignored
    mjd = np.array([b["mjd"] for b in B])
    o = np.argsort(mjd)
    ctrl = []
    for step in (1, 2, 3):
        for u, v in zip(o[:-step], o[step:], strict=True):
            if mjd[v] - mjd[u] < 3 and sep[u, v] < 2:
                rr, lg, x, ex, y, ey = pair_rho(
                    B[u]["flux"], B[u]["errs"], B[u]["k"], B[v]["flux"], B[v]["errs"], B[v]["k"]
                )
                ctrl.append(
                    {
                        "burst1": B[u]["name"],
                        "burst2": B[v]["name"],
                        "delay_s": round(float((mjd[v] - mjd[u]) * 86400), 1),
                        "sep_deg": round(float(sep[u, v]), 3),
                        "removed_by_position_cut": bool(not incons[u, v]),
                        "rho": round(rr, 4),
                        "chi2_p": float(f"{bt.twin_chi2(x, ex, y, ey, lg)[2]:.3g}"),
                        "npulse": [B[u]["npulse"], B[v]["npulse"]],
                    }
                )
    print(
        f"positive controls (re-triggers): {len(ctrl)}; all removed by position cut: "
        f"{all(c['removed_by_position_cut'] for c in ctrl)}",
        flush=True,
    )

    # ---- injections through the whole chain (primary gain jitter) + jitter sensitivity runs
    G["rho_star"] = rho_star
    rng = np.random.default_rng(a.seed)
    elig_idx = np.flatnonzero(elig)
    ratios = (1.0, 0.5, 0.3, 0.2, 0.1)
    jitters = [P.inject_gain_jitter] + [x for x in a.jitter_grid if x != P.inject_gain_jitter]
    jobs = []
    for jit in jitters:
        n_j = a.inject if jit == P.inject_gain_jitter else a.inject_sens
        for q in range(n_j):
            ia = int(rng.choice(elig_idx))
            ib = int(rng.choice(np.flatnonzero(incons[ia])))
            jobs.append((a.seed + 50000 + len(jobs), ia, ib, ratios[q % len(ratios)], float(jit)))
    gi = dict(g, rho_star=rho_star)
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(gi,)) as pool:
        inj = list(pool.map(injection_one, jobs, chunksize=20))
    it = Table(
        rows=inj,
        names=(
            "a",
            "b",
            "ratio",
            "gain_jitter",
            "snr_copy",
            "npulse_copy",
            "dk_copy",
            "lost_k",
            "eligible",
            "rho",
            "chi2_p",
            "flagged",
            "recovered",
            "recovered_deep",
        ),
    )
    it["a"] = [B[x]["name"] for x in it["a"]]
    it["b"] = [B[x]["name"] for x in it["b"]]
    it["snr_copy"] = np.round(it["snr_copy"], 2)
    it["rho"] = np.round(it["rho"], 4)
    it.meta = {
        "provenance": str(Provenance.SIMULATED),
        "source": "jwst_anomaly.burst_twins.inject_twin: smoothed template of eligible burst a, "
        "re-noised as a and as a copy in slot b (b's background noise, independent realisations, "
        "per-band gain jitter as listed); the copy gets its own window and k from a T90 proxy; "
        f"chain = >= {P.min_pulses} pulses, |dk| <= {P.max_dk}, rho > rho* = {rho_star:.4f}, "
        f"chi2 p >= {P.chi2_p_min}, delay >= {P.retrigger_days} d; recovered_deep is post hoc "
        f"(illustrative); scripts/s1_twins.py seed {a.seed}",
    }
    write_ecsv_gz(it, OUT / "injections.ecsv.gz")

    def efficiencies(tab):
        out = {}
        for rt in ratios:
            m = tab["ratio"] == rt
            out[str(rt)] = {
                "n": int(m.sum()),
                "lost_k": float(tab["lost_k"][m].mean()),
                "eligible": float(tab["eligible"][m].mean()),
                "flagged": float(tab["flagged"][m].mean()),
                "recovered": float(tab["recovered"][m].mean()),
                "recovered_deep_posthoc": float(tab["recovered_deep"][m].mean()),
            }
        return out

    prim = it[it["gain_jitter"] == P.inject_gain_jitter]
    eff = efficiencies(prim)
    eff_mean = float(np.mean([eff[str(x)]["recovered"] for x in ratios]))
    eff_mean_deep = float(np.mean([eff[str(x)]["recovered_deep_posthoc"] for x in ratios]))
    n_surv = len(survivors)
    n_elig = int(elig.sum())
    mu95 = bt.poisson_upper_limit(n_surv)
    mu95_deep = bt.poisson_upper_limit(len(deep_survivors))

    def pair_limit(mu, e):
        return mu / (n_elig * e) if e > 0 else None

    sens = {}
    for jit in jitters:
        e_j = efficiencies(it[it["gain_jitter"] == jit])
        e1 = e_j["1.0"]["recovered"]
        em = float(np.mean([e_j[str(x)]["recovered"] for x in ratios]))
        sens[f"{jit:.2f}"] = {
            "n": int((it["gain_jitter"] == jit).sum()),
            "eff_ratio1": e1,
            "eff_mean_ratio_0.1_1": em,
            "f95_twin_pairs_per_eligible_burst_ratio1": pair_limit(mu95, e1),
            "f95_twin_pairs_per_eligible_burst_mean_ratio_0.1_1": pair_limit(mu95, em),
        }
    f1 = pair_limit(mu95, eff["1.0"]["recovered"])
    fm = pair_limit(mu95, eff_mean)
    limits = {
        "definition": "f95 = mu95 / (N_eligible * efficiency): 95 % upper limit on twin PAIRS "
        "per eligible "
        "burst; the fraction of eligible bursts that have a twin is 2 * f95",
        "chain": "primary (rho > rho*, chi2, re-trigger veto)",
        "gain_jitter": P.inject_gain_jitter,
        "survivors_k": n_surv,
        "poisson_mu95": mu95,
        "n_eligible": n_elig,
        "f95_twin_pairs_per_eligible_burst_ratio1": f1,
        "f95_fraction_of_bursts_with_twin_ratio1": 2 * f1 if f1 else None,
        "f95_twin_pairs_per_eligible_burst_mean_ratio_0.1_1": fm,
        "f95_fraction_of_bursts_with_twin_mean_ratio_0.1_1": 2 * fm if fm else None,
        "gain_jitter_sensitivity": sens,
        "deep_posthoc_illustrative": {
            "note": "rho > 0.90 was chosen after seeing that it gives k = 0; "
            "illustrative, NOT a limit",
            "survivors_k": len(deep_survivors),
            "f95_twin_pairs_per_eligible_burst_ratio1": pair_limit(
                mu95_deep, eff["1.0"]["recovered_deep_posthoc"]
            ),
            "f95_twin_pairs_per_eligible_burst_mean_ratio_0.1_1": pair_limit(
                mu95_deep, eff_mean_deep
            ),
        },
    }
    import astropy
    import scipy

    summary = {
        "argv": sys.argv,
        "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(a).items()},
        "seed": a.seed,
        "n_injections": len(it),
        "n_injections_primary": len(prim),
        "versions": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "astropy": astropy.__version__,
        },
        "catalogue_bursts": len(cat),
        "bcat_status": {
            k: int(v)
            for k, v in zip(*np.unique(list(status.values()), return_counts=True), strict=True)
        },
        "bursts_reduced": N,
        "eligible_bursts": n_elig,
        "eligible_pairs": n_elig_pairs,
        "pairs_scored_s1": int(r.size),
        "params": {k: getattr(P, k) for k in P.__dataclass_fields__},
        "rho_star": rho_star,
        "real_max_rho": float(r.max()),
        "real_pairs_above": dict(
            zip(("0.8", "0.9", "0.95", "rho_star"), cdf_real.tolist(), strict=True)
        ),
        "null_pairs_above_per_catalogue": dict(
            zip(("0.8", "0.9", "0.95", "rho_star"), cdf_null.tolist(), strict=True)
        ),
        "flagged_real": int(flagged.size),
        "automated_survivors": [
            (d["burst1"], d["burst2"], d["rho"], d["chi2_p"]) for d in survivors
        ],
        "deep_chain_posthoc_illustrative": {
            "rho_prescreen": P.rho_prescreen_deep,
            "pairs": int(deep_idx.size),
            "chi2_pass": deep_chi2_pass,
            "survivors": deep_survivors,
        },
        "positive_controls_retrigger": ctrl,
        "injection_efficiency": eff,
        "limits": limits,
        "runtime_s": round(time.time() - t0),
    }
    write_json(summary, OUT / "summary.json")
    write_json(
        {
            "provenance": str(Provenance.SIMULATED),
            "source": "pulse-shuffled surrogate catalogues "
            "(jwst_anomaly.burst_twins.pulse_shuffle), same pair set",
            "n_catalogues": a.null,
            "per_catalogue_max_rho": np.round(np.sort(maxes), 4).tolist(),
            "rho_star_95": rho_star,
            "top_surrogate_rho": np.round(np.sort(tops)[::-1][:50], 4).tolist(),
            "hist_bins": [float(bins[0]), float(bins[-1]), len(bins) - 1],
            "hist_null_mean_counts_rho_gt_0.5": {
                f"{bins[q]:.3f}": round(float(hist[q] / a.null), 3)
                for q in range(len(hist))
                if bins[q] >= 0.5 and hist[q]
            },
            "hist_real_counts_rho_gt_0.5": {
                f"{bins[q]:.3f}": int(real_hist[q])
                for q in range(len(real_hist))
                if bins[q] >= 0.5 and real_hist[q]
            },
        },
        OUT / "null.json",
    )
    print(
        json.dumps(
            {k: summary[k] for k in ("rho_star", "flagged_real", "injection_efficiency", "limits")},
            indent=1,
        )
    )

    if a.plots:
        a.plots.mkdir(parents=True, exist_ok=True)
        sheet = [
            {
                "i": d["_i"],
                "j": d["_j"],
                "lag": d["lag_bins"],
                "K": d["_K"],
                "title": f"rho={d['rho']} p={d['chi2_p']}",
            }
            for d in rows[:24]
        ]
        contact_sheet(B, sheet, a.plots / "top_pairs.png")
        fl = [d for d in rows if d["flag"]][:40]
        if fl:
            contact_sheet(
                B,
                [
                    {
                        "i": d["_i"],
                        "j": d["_j"],
                        "lag": d["lag_bins"],
                        "K": d["_K"],
                        "title": f"rho={d['rho']} p={d['chi2_p']}",
                    }
                    for d in fl
                ],
                a.plots / "flagged.png",
            )
        # injected examples near threshold
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
        ax[0].hist(maxes, bins=30, label="surrogate per-catalogue max")
        ax[0].axvline(rho_star, color="k", ls="--", label="rho*")
        ax[0].axvline(r.max(), color="r", label="real max")
        ax[0].legend(fontsize=7)
        snr = np.asarray(prim["snr_copy"])
        for rt in ratios:
            m = prim["ratio"] == rt
            ax[1].scatter(snr[m], np.asarray(prim["rho"])[m], s=3, label=f"ratio {rt}")
        ax[1].axhline(rho_star, color="k", ls="--")
        ax[1].set_xscale("log")
        ax[1].set_xlabel("copy peak S/N")
        ax[1].set_ylabel("rho")
        ax[1].legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(a.plots / "null_and_injections.png", dpi=110)


if __name__ == "__main__":
    main()
