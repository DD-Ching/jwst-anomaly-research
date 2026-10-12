"""S1 round 3: does catalogue E_peak raise full recovery at equal false-candidate burden?

Pre-declared in docs/s1_energy/plan.md. Compares A (frozen M2 time-shape screen) with B (staged:
E-testable pairs need M2 > thr_B and |z| <= 2, thr_B calibrated for equal flagged burden; other
pairs use A) on the same simulated same-source pairs (validation half) and real unrelated pairs.

  python scripts/s1_energy.py errors              -> results/s1_energy/error_checks.json
  python scripts/s1_energy.py compare [--cpu 4]   -> results/s1_energy/compare.json
  python scripts/s1_energy.py screen  [--cpu 4]   -> results/s1_energy/screen.json (only if B holds)

Catalogue-level simulation only: injected twins are simulated light curves with a modelled
catalogue E_peak; recovering them measures search capability, never evidence of new physics.
Convention: the program stretches b by s_prog to match a; physically b = a's source with
t -> s_phys t, E -> E / s_phys, so s_prog = 1 / s_phys and E_b = E_a * s_prog.
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
import s1_m2_search as ms  # noqa: E402

sd = ms.sd
OUT = ROOT / "results" / "s1_energy"
SPEC = OUT / "gbm_spectral.ecsv.gz"
SIG_MAX, E_MIN, E_MAX = 0.5, 10.0, 1.0e4  # quality rule (plan section 1)
SYS_TEST = 0.20  # test-side systematic in ln E_peak (plan section 2), fixed before validation
Z_CUT = 2.0
SCEN = {
    "optimistic": (0.10, 0.025),
    "central": (0.20, 0.05),
    "pessimistic": (0.29, 0.075),
    "literature": (0.35, 0.075),
}
MIN_GAIN = 0.02
N_DONOR = 20


# ------------------------------------------------------------------------------- catalogue
def _col(t, c):
    v = t[c]
    return np.array(v.filled(np.nan) if hasattr(v, "filled") else v, float)


def lnerr(t, w, m):
    """ln E_peak and mean log error of fit (w, m); nan where the fit is absent."""
    e = _col(t, f"{w}_{m}_epeak")
    p, n = np.abs(_col(t, f"{w}_{m}_epeak_pos_err")), np.abs(_col(t, f"{w}_{m}_epeak_neg_err"))
    ok = np.isfinite(e) & (e > 0) & np.isfinite(p) & np.isfinite(n) & (p > 0)
    e1 = np.where(ok, e, 1.0)
    sp = np.log((e1 + p) / e1)
    sn = np.abs(np.log(np.maximum(e1 - n, 1e-3) / e1))
    return np.where(ok, np.log(e1), np.nan), np.where(ok, 0.5 * (sp + sn), np.nan)


def quality(lne, sig):
    return np.isfinite(lne) & (sig <= SIG_MAX) & (lne >= math.log(E_MIN)) & (lne <= math.log(E_MAX))


def load_spec():
    from astropy.table import Table

    t = Table.read(SPEC, format="ascii.ecsv")
    names = [str(x) for x in t["trigger_name"]]
    lne, sig = lnerr(t, "flnc", "comp")
    lpk, spk = lnerr(t, "pflx", "comp")
    q = quality(lne, sig)
    resid = lpk - lne
    wq = q & quality(lpk, spk)
    return {
        "t": t,
        "names": names,
        "pos": {n: i for i, n in enumerate(names)},
        "lne": lne,
        "sig": sig,
        "q": q,
        "fluence": _col(t, "fluence"),
        "t90": _col(t, "t90"),
        "best": np.array([str(x).strip().lower() for x in t["flnc_best_fitting_model"]]),
        "win": resid[wq],
        "meta": dict(t.meta),
    }


def robust_sd(x):
    x = np.asarray(x, float)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 3), round(c + h, 3)]


def errors(a):
    S = load_spec()
    t = S["t"]
    out = {"catalogue": S["meta"], "rows": len(S["names"])}
    lb, sb = lnerr(t, "flnc", "band")
    best_ok = ((S["best"] == "flnc_band") & np.isfinite(lb)) | (
        (S["best"] == "flnc_comp") & np.isfinite(S["lne"])
    )
    out["availability"] = {
        "round2_best_model_any": round(float(best_ok.mean()), 3),
        "comp_any_fit": round(float(np.isfinite(S["lne"]).mean()), 3),
        "comp_quality_rule": round(float(S["q"].mean()), 3),
    }
    fl, t90 = S["fluence"], S["t90"]
    qs = np.nanpercentile(fl, [0, 25, 50, 75, 100])
    out["availability_by_fluence_quartile"] = [
        {
            "fluence_range": [float(f"{lo:.2g}"), float(f"{hi:.2g}")],
            "comp_quality": round(float(S["q"][(fl >= lo) & (fl <= hi)].mean()), 3),
            "round2": round(float(best_ok[(fl >= lo) & (fl <= hi)].mean()), 3),
        }
        for lo, hi in zip(qs[:-1], qs[1:], strict=True)
    ]
    out["availability_by_type"] = {
        "short_t90_lt_2s": round(float(S["q"][t90 < 2].mean()), 3),
        "long": round(float(S["q"][t90 >= 2].mean()), 3),
        "by_best_model": {
            k: round(float(S["q"][S["best"] == k].mean()), 3)
            for k in sorted(set(S["best"]))
            if (S["best"] == k).sum() > 20
        },
    }
    m = np.isfinite(lb) & S["q"] & (sb <= SIG_MAX)
    d = lb[m] - S["lne"][m]
    out["model_choice_band_minus_comp"] = {
        "n": int(m.sum()),
        "median": round(float(np.median(d)), 3),
        "robust_sd": round(robust_sd(d), 3),
        "median_stat": round(float(np.median(np.hypot(sb[m], S["sig"][m]))), 3),
        "note": "same photons: stat errors are correlated, so the excess is a lower bound",
    }
    w = S["win"]
    lpk, spk = lnerr(t, "pflx", "comp")
    wq = S["q"] & quality(lpk, spk)
    ms_ = float(np.median(np.hypot(S["sig"][wq], spk[wq])))
    out["window_peak_minus_fluence_comp"] = {
        "n": int(w.size),
        "median": round(float(np.median(w)), 3),
        "robust_sd": round(robust_sd(w), 3),
        "median_stat": round(ms_, 3),
        "excess_sd": round(math.sqrt(max(robust_sd(w) ** 2 - ms_**2, 0)), 3),
        "pull_robust_sd": round(robust_sd((w - np.median(w)) / np.hypot(S["sig"][wq], spk[wq])), 2),
    }
    out["sigma_ln_quantiles_quality"] = {
        str(p): round(float(np.quantile(S["sig"][S["q"]], p)), 3) for p in (0.1, 0.5, 0.9)
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "error_checks.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


# ------------------------------------------------------------------------------- compare
def zval(lna, sa, lnb, sb, s_prog, sig_lns):
    var = sa * sa + sb * sb + 2 * SYS_TEST * SYS_TEST + sig_lns * sig_lns
    return (lna - lnb + np.log(s_prog)) / np.sqrt(var)


def burst_epeak(S, B):
    """Per reduced burst: (ln E, sigma, quality) of the flnc Comptonized fit."""
    lne = np.full(len(B), np.nan)
    sig = np.full(len(B), np.nan)
    q = np.zeros(len(B), bool)
    for i, b in enumerate(B):
        k = S["pos"].get(b["name"])
        if k is not None:
            lne[i], sig[i], q[i] = S["lne"][k], S["sig"][k], S["q"][k]
    return lne, sig, q


def ekey(lne):
    """Band-factor cache key: template E_peak rounded to 5 % (200 keV when absent)."""
    e = math.exp(lne) if np.isfinite(lne) else 200.0
    return float(round(math.exp(round(math.log(e) / 0.05) * 0.05), 1))


def injections(B, tmpl, slots, incons, lne_b, seed, per_cell, cpu):
    rng = np.random.default_rng(seed)
    jobs, keys = [], set()
    for s0 in sd.S_INJ:
        for ratio in sd.RATIOS:
            for _ in range(per_cell):
                ia = int(rng.choice(tmpl))
                cand = slots[incons[ia, slots]]
                ek = ekey(lne_b[ia])
                keys.add((s0, ek))
                jobs.append((seed + len(jobs), ia, int(rng.choice(cand)), s0, ek, ratio))
    bf = {k: sd.band_factors(*k) for k in keys}
    g = {"B": B, "bf": bf}
    with cf.ProcessPoolExecutor(cpu, initializer=ms._init, initargs=(g,)) as pool:
        inj = np.array(list(pool.map(sd.inject_one, jobs, chunksize=8)), float)
    cols = ["a", "b", "s0", "epeak", "ratio", "dk", "elig", "dtd"]
    cols += ["m0", "m1", "m2", "s1", "s2", "p0", "p1", "p2"]
    return {c: inj[:, q] for q, c in enumerate(cols)}


def copy_measurement(S, inj, B, lne_b, sig_b, rng, scen):
    """Simulated catalogue E_peak of each injected copy: (ln E_meas, sigma, available)."""
    sd_w, bias = SCEN[scen]
    w0 = S["win"] - np.median(S["win"])
    w_scale = sd_w / robust_sd(S["win"])
    fl_all = np.log10(np.where(S["fluence"] > 0, S["fluence"], np.nan))
    has_fit = np.isfinite(S["lne"]) & np.isfinite(fl_all)
    lf_fit, le_fit = fl_all[has_fit], S["lne"][has_fit] / math.log(10)
    q_fit, s_fit = S["q"][has_fit], S["sig"][has_fit]
    ok_fl = np.isfinite(fl_all)
    n = inj["a"].size
    lnm, sg, av = np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, bool)
    for r in range(n):
        ia = int(inj["a"][r])
        if not np.isfinite(lne_b[ia]):
            continue
        s_phys, ratio = inj["s0"][r], inj["ratio"][r]
        lf_c = math.log10(max(B[ia]["fluence"], 1e-12) * ratio * s_phys)
        ltrue = lne_b[ia] - math.log(s_phys)  # template's catalogue value as the source (+eps_a)
        # 1. any Comp fit at this fluence (nearest 50 in log fluence among all bursts)
        nn = np.argsort(np.abs(fl_all[ok_fl] - lf_c))[:50]
        if rng.random() > has_fit[ok_fl][nn].mean():
            continue
        # 2. donor among the 20 nearest fitted bursts in (log fluence, log E)
        dist = np.hypot(lf_fit - lf_c, le_fit - ltrue / math.log(10))
        k = int(rng.choice(np.argsort(dist)[:N_DONOR]))
        if not q_fit[k]:
            continue
        sc = float(s_fit[k])
        eps_a = rng.normal(0, sig_b[ia])
        w = bias + w_scale * float(rng.choice(w0))
        val = ltrue - eps_a + rng.normal(0, sc) + w
        if not (math.log(E_MIN) <= ltrue <= math.log(E_MAX)):
            continue
        lnm[r], sg[r], av[r] = val, sc, True
    return lnm, sg, av


def background_z(B, bg, lne_b, sig_b, q_b, sig_lns):
    i, j = bg[:, 0].astype(int), bg[:, 1].astype(int)
    test = q_b[i] & q_b[j]
    z = np.full(bg.shape[0], np.nan)
    z[test] = zval(
        lne_b[i[test]], sig_b[i[test]], lne_b[j[test]], sig_b[j[test]], bg[test, 3], sig_lns
    )
    return test, z


def thr_for(m2, zpass, n_target):
    """Largest threshold with count(M2 > thr and z pass) <= n_target."""
    v = np.sort(m2[zpass])[::-1]
    return float(v[n_target]) if n_target < v.size else float(v[-1]) - 1e-9


def chi2_for(B, pairs, cpu):
    g = {"B": B}
    with cf.ProcessPoolExecutor(cpu, initializer=ms._init, initargs=(g,)) as pool:
        return np.array([r[2] for r in pool.map(ms._chi2_job, pairs, chunksize=4)], float)


def compare(a):
    t0 = time.time()
    fz = json.loads((ROOT / "results/s1_m2/frozen.json").read_text())
    thr_a, c_cut = fz["threshold_m2"], fz["consistency_cut"]
    S = load_spec()
    cat, B, elig, ks, sep, incons, allowed, cal, val = ms.setup()
    lne_b, sig_b, q_b = burst_epeak(S, B)
    out = {
        "plan": "docs/s1_energy/plan.md",
        "catalogue_sha": S["meta"].get("votable_sha256_16"),
        "thr_A": thr_a,
        "sys_test": SYS_TEST,
        "eligible_with_quality_epeak": {
            "cal": round(float(q_b[cal].mean()), 3),
            "val": round(float(q_b[val].mean()), 3),
        },
    }
    # 1. sigma_lns from calibration-half injections
    ic = injections(B, cal, cal, incons, lne_b, 31, a.per_cell_cal, a.cpu)
    sel = (ic["elig"] > 0) & (np.nan_to_num(ic["m2"], nan=-1) > 0.9)
    dl = np.log(ic["s2"][sel] * ic["s0"][sel])
    sig_lns = math.ceil(robust_sd(dl) * 100) / 100
    out["sigma_lns"] = {
        "n": int(sel.sum()),
        "robust_sd": round(robust_sd(dl), 4),
        "used": sig_lns,
        "frac_abs_gt_0.1": round(float((np.abs(dl) > 0.1).mean()), 3),
        "frac_abs_gt_0.3": round(float((np.abs(dl) > 0.3).mean()), 3),
    }
    # 2. thr_B on the calibration half (E-testable pairs, flagged level)
    bgc = ms.background(B, allowed, cal, a.cpu)
    m2c = np.nan_to_num(bgc[:, 2], nan=-np.inf)
    testc, zc = background_z(B, bgc, lne_b, sig_b, q_b, sig_lns)
    zpc = testc & (np.abs(np.nan_to_num(zc, nan=99)) <= Z_CUT)
    n_a = int((testc & (m2c > thr_a)).sum())
    thr_b = thr_for(m2c, zpc, n_a)
    lo_n, hi_n = max(n_a - int(round(math.sqrt(n_a))), 0), n_a + int(round(math.sqrt(n_a)))
    thr_b_var = {"n_minus_sqrt": thr_for(m2c, zpc, lo_n), "n_plus_sqrt": thr_for(m2c, zpc, hi_n)}
    out["calibration"] = {
        "pairs": int(bgc.shape[0]),
        "e_testable": int(testc.sum()),
        "flagged_A_testable": n_a,
        "flagged_A_all": int((m2c > thr_a).sum()),
        "thr_B": round(thr_b, 4),
        "thr_B_sensitivity": {k: round(v, 4) for k, v in thr_b_var.items()},
        "z_pass_rate_testable_m2_gt_0.9": round(
            float(zpc[testc & (m2c > 0.9)].mean()) if (testc & (m2c > 0.9)).any() else float("nan"),
            3,
        ),
    }
    # z pass vs M2 (dependence of the two checks among unrelated pairs)
    edges = [0.80, 0.85, 0.90, 0.93, 0.95, 0.96, thr_a, 1.0]
    out["calibration"]["z_pass_by_m2"] = [
        {
            "m2": [round(lo, 4), round(hi, 4)],
            "n_testable": int((testc & (m2c > lo) & (m2c <= hi)).sum()),
            "z_pass": round(float(zpc[testc & (m2c > lo) & (m2c <= hi)].mean()), 3)
            if (testc & (m2c > lo) & (m2c <= hi)).any()
            else None,
        }
        for lo, hi in zip(edges[:-1], edges[1:], strict=True)
    ]
    # 3. validation-half burden (val-val E-testable real pairs)
    bgv = ms.background(B, allowed, val, a.cpu)
    m2v = np.nan_to_num(bgv[:, 2], nan=-np.inf)
    testv, zv = background_z(B, bgv, lne_b, sig_b, q_b, sig_lns)
    zpv = testv & (np.abs(np.nan_to_num(zv, nan=99)) <= Z_CUT)
    fa = testv & (m2v > thr_a)
    fb = zpv & (m2v > thr_b)
    need = np.flatnonzero(fa | fb)
    pv = chi2_for(B, [(int(bgv[q, 0]), int(bgv[q, 1])) for q in need], a.cpu)
    pmap = dict(zip(need.tolist(), pv.tolist(), strict=True))
    dtd = {q: abs(B[int(bgv[q, 0])]["mjd"] - B[int(bgv[q, 1])]["mjd"]) for q in need}

    def cand(mask):
        return int(sum(1 for q in np.flatnonzero(mask) if pmap[q] >= c_cut and dtd[q] >= 1.0))

    out["validation_burden"] = {
        "pairs": int(bgv.shape[0]),
        "e_testable": int(testv.sum()),
        "flagged_A": int(fa.sum()),
        "flagged_B": int(fb.sum()),
        "candidates_A": cand(fa),
        "candidates_B": cand(fb),
        "note": "E-testable val-val pairs; untestable pairs use A in both methods",
    }
    # 4. validation injections and full recovery
    iv = injections(B, val, val, incons, lne_b, 4343, a.per_cell, a.cpu)
    n = iv["a"].size
    ok = (iv["elig"] > 0) & (iv["dtd"] >= 1.0) & (np.nan_to_num(iv["p2"]) >= c_cut)
    m2i = np.nan_to_num(iv["m2"], nan=-1)
    rec_a = ok & (m2i > thr_a)
    tmpl_q = q_b[iv["a"].astype(int)]
    clusters = iv["a"].astype(int)
    uc = np.unique(clusters)
    rng_b = np.random.default_rng(7)
    boots = [rng_b.choice(uc, uc.size) for _ in range(a.nboot)]
    idx_of = {c: np.flatnonzero(clusters == c) for c in uc}
    res = {}
    for scen in SCEN:
        lnm, sgm, av = copy_measurement(S, iv, B, lne_b, sig_b, np.random.default_rng(11), scen)
        testable = tmpl_q & av
        z = np.full(n, np.nan)
        ia = iv["a"].astype(int)
        z[testable] = zval(
            lne_b[ia[testable]],
            sig_b[ia[testable]],
            lnm[testable],
            sgm[testable],
            iv["s2"][testable],
            sig_lns,
        )
        zp = testable & (np.abs(np.nan_to_num(z, nan=99)) <= Z_CUT)

        def rec_b(thr, zp=zp, testable=testable):
            return np.where(testable, ok & zp & (m2i > thr), rec_a)

        rb = rec_b(thr_b)
        d_boot = []
        for bs in boots:
            ix = np.concatenate([idx_of[c] for c in bs])
            d_boot.append(rb[ix].mean() - rec_a[ix].mean())
        lo, hi = np.percentile(d_boot, [2.5, 97.5])
        r = {
            "copy_epeak_available": round(float(av.mean()), 3),
            "e_testable_injections": round(float(testable.mean()), 3),
            "true_pair_z_pass_given_testable": round(float(zp[testable].mean()), 3),
            "R_A": round(float(rec_a.mean()), 4),
            "R_B": round(float(rb.mean()), 4),
            "dR": round(float(rb.mean() - rec_a.mean()), 4),
            "dR_95": [round(float(lo), 4), round(float(hi), 4)],
            "dR_thr_B_sensitivity": {
                k: round(float(rec_b(v).mean() - rec_a.mean()), 4) for k, v in thr_b_var.items()
            },
            "testable": {
                "n": int(testable.sum()),
                "R_A": round(float(rec_a[testable].mean()), 4),
                "R_B": round(float(rb[testable].mean()), 4),
            },
            "untestable": {
                "n": int((~testable).sum()),
                "R_A": round(float(rec_a[~testable].mean()), 4),
                "R_B": round(float(rb[~testable].mean()), 4),
            },
            "by_ratio": {
                str(rt): {
                    "R_A": round(float(rec_a[iv["ratio"] == rt].mean()), 3),
                    "R_B": round(float(rb[iv["ratio"] == rt].mean()), 3),
                }
                for rt in sd.RATIOS
            },
            "gained": int((rb & ~rec_a).sum()),
            "lost": int((rec_a & ~rb).sum()),
        }
        r["holds"] = bool(r["dR"] >= MIN_GAIN and lo > 0)
        res[scen] = r
    out["validation_injections"] = {"n": int(n), "templates": int(uc.size), "scenarios": res}
    vb = out["validation_burden"]
    burden_ok = vb["flagged_B"] <= vb["flagged_A"] + 1.96 * math.sqrt(max(vb["flagged_A"], 1))
    out["decision"] = {
        "burden_ok": bool(burden_ok),
        "holds_central": res["central"]["holds"],
        "holds_pessimistic": res["pessimistic"]["holds"],
        "holds_literature": res["literature"]["holds"],
        "improvement_established": bool(
            burden_ok and all(res[k]["holds"] for k in ("central", "pessimistic", "literature"))
        ),
    }
    out["runtime_s"] = round(time.time() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "compare.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=("errors", "compare"))
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--per-cell", type=int, default=40)
    ap.add_argument("--per-cell-cal", type=int, default=12)
    ap.add_argument("--nboot", type=int, default=1000)
    a = ap.parse_args(argv)
    {"errors": errors, "compare": compare}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
