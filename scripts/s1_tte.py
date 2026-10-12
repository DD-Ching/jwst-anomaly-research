"""S1 TTE round 2 (docs/s1_blind/plan_round2.md): robustness of new multi-pulse short bursts,
controls, joint time-energy (E_peak) constraint, and a limited candidate check.

  python scripts/s1_tte.py fetch        -> TTE cache (scratch), results/s1_blind/tte/manifest.json
  python scripts/s1_tte.py robust       -> results/s1_blind/tte/robustness.json, reduced npz
  python scripts/s1_tte.py epeak        -> results/s1_blind/tte/epeak_test.json
  python scripts/s1_tte.py candidates   -> results/s1_blind/tte/candidates.json

TTE files are cached under $S1_TTE_SCRATCH during the round and are not committed; the manifest
records file names, versions, sizes and sha256 prefixes. Download cap 1 GB per round.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import s1_blind_spots as bs  # noqa: E402

from jwst_anomaly import burst_twins as bt  # noqa: E402

OUT = ROOT / "results" / "s1_blind" / "tte"
CACHE = Path(os.environ.get("S1_TTE_SCRATCH", "/tmp/s1_tte_cache"))
CAP_BYTES = 1.0e9
NEW = ["bn240905186", "bn101129726", "bn210528586", "bn081113230", "bn250330911"]
BANDS = (((50.0, 300.0),), ((8.0, 50.0), (300.0, 900.0)))
RES = (0.002, 0.004, 0.008, 0.016, 0.032, 0.064)
OFFSETS = (0.0, 0.25, 0.5, 0.75)
BKG = {
    "nominal": {"pre": (-20.0, -2.0), "post": (5.0, 30.0), "deg": 1},
    "alt": {"pre": (-40.0, -5.0), "post": (10.0, 60.0), "deg": 0},
}


def load_cat():
    from astropy.table import Table

    cat = Table.read(ROOT / "results/s1_twins/catalogue.ecsv.gz", format="ascii.ecsv")
    bur = Table.read(ROOT / "results/s1_twins/bursts.ecsv.gz", format="ascii.ecsv")
    return cat, bur, {str(n): i for i, n in enumerate(cat["trigger_name"])}


# ------------------------------------------------------------------------------------------ fetch
def fetch(a):
    plan = json.loads(Path(a.plan).read_text())
    masks, sizes = plan["masks"], plan["sizes"]
    rng = np.random.default_rng(7)
    singles = list(rng.choice(plan["single"], size=8, replace=False))
    elig = list(rng.permutation(plan["elig_short"]))
    want = list(NEW) + singles
    budget = CAP_BYTES * 0.97 - sum(sizes[n] for n in want)
    for n in elig:
        if sizes[n] <= budget:
            want.append(n)
            budget -= sizes[n]
    man = {
        "cap_bytes": CAP_BYTES,
        "files": {},
        "groups": {
            "new": NEW,
            "single_pilot": singles,
            "elig_short": [n for n in want if n in elig],
        },
    }
    total = 0
    for name in want:
        yr = "20" + name[2:4]
        for d, m in zip(bs.DETS[:12], masks[name][:12], strict=True):
            if m != "1":
                continue
            for v in ("v00", "v01", "v02"):
                fn = f"glg_tte_{d}_{name}_{v}.fit"
                dst = CACHE / name / fn
                dst.parent.mkdir(parents=True, exist_ok=True)
                try:
                    if not dst.exists():
                        urllib.request.urlretrieve(f"{bs.TRIG}/{yr}/{name}/current/{fn}", dst)
                    b = dst.read_bytes()
                    total += len(b)
                    man["files"].setdefault(name, []).append(
                        {"file": fn, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()[:16]}
                    )
                    break
                except Exception:
                    continue
        if total > CAP_BYTES:
            raise SystemExit(f"download cap exceeded at {name}: {total} bytes")
    man["total_bytes"] = total
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(
        {k: (len(v) if isinstance(v, list) else v) for k, v in man["groups"].items()}, total / 1e6
    )


# ------------------------------------------------------------------------------------------ events
def load_events(name, dets=None):
    """Per detector: (t rel. trigger, energy keV, gti list). Energies are channel mid-points."""
    from astropy.io import fits

    out = {}
    for f in sorted((CACHE / name).glob("glg_tte_*.fit")):
        d = f.name.split("_")[2]
        if dets is not None and d not in dets:
            continue
        with fits.open(f) as h:
            trig = float(h[0].header["TRIGTIME"])
            eb = h["EBOUNDS"].data
            emid = 0.5 * (np.asarray(eb["E_MIN"], float) + np.asarray(eb["E_MAX"], float))
            ev = h["EVENTS"].data
            t = np.asarray(ev["TIME"], float) - trig
            e = emid[np.asarray(ev["PHA"], int)]
            gti = [
                (float(s) - trig, float(x) - trig)
                for s, x in zip(h["GTI"].data["START"], h["GTI"].data["STOP"], strict=True)
            ]
        out[d] = (t, e, gti)
    return out


def in_gti(lo, hi, gti):
    return any(s <= lo and hi <= x for s, x in gti)


def common_gti(evs):
    """Intersection of the detectors' GTIs (list of (start, stop))."""
    gtis = [g for _, _, g in evs.values()]
    if not gtis:
        return []
    cur = gtis[0]
    for g in gtis[1:]:
        cur = [(max(a, c), min(b, d)) for a, b in cur for c, d in g if min(b, d) > max(a, c)]
    return cur


def band_events(evs):
    res = []
    for bands in BANDS:
        ts = []
        for t, e, _ in evs.values():
            m = np.zeros(t.size, bool)
            for lo, hi in bands:
                m |= (e >= lo) & (e < hi)
            ts.append(t[m])
        res.append(np.sort(np.concatenate(ts)) if ts else np.array([]))
    return res


def curve(bt_events, t90s, t90, dt, offset, bkg="nominal", gti=None):
    """Background-subtracted two-band rates on a grid of width dt, start shifted by offset*dt.

    Background: polynomial fit to 1-s bins of the two side windows, using only bins that lie
    inside the GTIs (no-data intervals must not count as zero background).
    """
    p = dataclasses.replace(bt.Params(), base_dt=0.002, max_bins=256, pad_frac=0.5, pad_s=0.1)
    lo, hi, _ = bt.window_for(t90s, t90, p)
    lo -= offset * dt
    n = max(4, int(math.ceil((hi - lo) / dt)))
    edges = lo + dt * np.arange(n + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    b = BKG[bkg]
    end = t90s + t90
    flux, err = np.zeros((2, n)), np.zeros((2, n))
    for k, tk in enumerate(bt_events):
        c = np.histogram(tk, edges)[0].astype(float)
        xs, ys = [], []
        for a_, b_ in (b["pre"], (end + b["post"][0], end + b["post"][1])):
            h1, e1 = np.histogram(tk, np.arange(a_, b_ + 1e-9, 1.0))
            for q in range(h1.size):
                if gti is None or in_gti(e1[q], e1[q + 1], gti):
                    xs.append(0.5 * (e1[q] + e1[q + 1]))
                    ys.append(h1[q])
        coef = np.polyfit(xs, ys, b["deg"])
        bk = np.polyval(coef, mid) * dt
        flux[k] = (c - bk) / dt
        err[k] = np.sqrt(np.maximum(c, 1.0)) / dt
    return flux, err


def npulse_multires(bt_events, t90s, t90, prominence=5.0, offsets=OFFSETS, bkg="nominal", gti=None):
    """Max over resolutions (2-64 ms) and bin offsets of the number of resolved pulses."""
    p = dataclasses.replace(bt.Params(), pulse_prominence_sigma=prominence)
    best = 0
    for dt in RES:
        for off in offsets:
            f, e = curve(bt_events, t90s, t90, dt, off, bkg, gti)
            best = max(best, len(bt.find_pulses(f, e, p)))
    return best


def exposure_ok(evs, t90s, t90, min_side_s=10.0):
    """Burst window inside the common GTI, and >= min_side_s of GTI-covered background per side."""
    g = common_gti(evs)
    end = t90s + t90
    if not in_gti(t90s - 0.2, end + 0.2, g):
        return False

    def covered(lo, hi):
        return sum(max(0.0, min(hi, b) - max(lo, a)) for a, b in g)

    return covered(-40.0, -2.0) >= min_side_s and covered(end + 5.0, end + 60.0) >= min_side_s


# -------------------------------------------------------------------------------- simulation (neg)
def sim_single_pulse(rng, rate_peak, bkg_rates, t90s, t90):
    """Band event lists: one FRED pulse (rise:decay 1:3, FWHM ~ T90/2) on flat backgrounds."""
    out = []
    lo, hi = -45.0, t90s + t90 + 65.0
    tpk = t90s + 0.3 * t90
    tr, td = max(t90 / 8, 0.002), max(3 * t90 / 8, 0.006)
    for k in range(2):
        nb = rng.poisson(bkg_rates[k] * (hi - lo))
        tb = rng.uniform(lo, hi, nb)
        area = rate_peak[k] * (tr + td)
        ns = rng.poisson(area)
        side = rng.random(ns) < tr / (tr + td)
        ts = np.where(side, tpk - rng.exponential(tr, ns), tpk + rng.exponential(td, ns))
        out.append(np.sort(np.r_[tb, ts]))
    return out


def rates_for(name, cat, idx):
    """Peak net rate (8-ms bins) and background rate per band, from the real events."""
    c = cat[idx[name]]
    evs = load_events(name)
    be = band_events(evs)
    f, _ = curve(be, float(c["t90_start"]), float(c["t90"]), 0.008, 0.0)
    peak = np.maximum(f.max(1), 1.0)
    bkg = np.array([np.sum((tk > -20) & (tk < -2)) / 18.0 for tk in be])
    return peak, bkg


# ------------------------------------------------------------------------------------------ robust
def robust(a):
    cat, bur, idx = load_cat()
    man = json.loads((OUT / "manifest.json").read_text())
    out = {"prominence_sigma": 5.0, "events": {}, "controls": {}}
    rng = np.random.default_rng(11)
    # negative control: matched to the 5 new events and the 8 single-pulse pilot bursts
    neg = []
    for name in man["groups"]["new"] + man["groups"]["single_pilot"]:
        c = cat[idx[name]]
        peak, bkg = rates_for(name, cat, idx)
        for _ in range(a.nsim):
            ev = sim_single_pulse(rng, peak, bkg, float(c["t90_start"]), float(c["t90"]))
            neg.append(npulse_multires(ev, float(c["t90_start"]), float(c["t90"])))
    neg = np.array(neg)
    prom = 5.0
    frac = float((neg >= 2).mean())
    out["controls"]["negative_single_pulse"] = {"n": int(neg.size), "false_ge2_at_5sigma": frac}
    if frac > 0.05:  # pre-declared: raise the threshold until <= 5 %
        for prom in (6.0, 7.0, 8.0, 10.0):
            rng2 = np.random.default_rng(11)
            neg2 = []
            for name in man["groups"]["new"] + man["groups"]["single_pilot"]:
                c = cat[idx[name]]
                peak, bkg = rates_for(name, cat, idx)
                for _ in range(a.nsim):
                    ev = sim_single_pulse(rng2, peak, bkg, float(c["t90_start"]), float(c["t90"]))
                    neg2.append(npulse_multires(ev, float(c["t90_start"]), float(c["t90"]), prom))
            f2 = float((np.array(neg2) >= 2).mean())
            out["controls"].setdefault("threshold_scan", []).append(
                {"prominence": prom, "false_ge2": f2}
            )
            if f2 <= 0.05:
                break
    out["prominence_sigma"] = prom
    # positive control: eligible short bursts
    pos = []
    for name in man["groups"]["elig_short"]:
        c = cat[idx[name]]
        evs = load_events(name)
        pos.append(
            npulse_multires(
                band_events(evs), float(c["t90_start"]), float(c["t90"]), prom, gti=common_gti(evs)
            )
        )
    pos = np.array(pos)
    k, n = int((pos >= 2).sum()), int(pos.size)
    out["controls"]["positive_bcat_multi_pulse_short"] = {
        "n": n,
        "recovered_ge2": k,
        "wilson95": wilson(k, n),
    }
    # the five new events: nominal + 12 one-factor variants
    for name in man["groups"]["new"]:
        c = cat[idx[name]]
        t90s, t90 = float(c["t90_start"]), float(c["t90"])
        evs = load_events(name)
        be = band_events(evs)
        g = common_gti(evs)
        nom = npulse_multires(be, t90s, t90, prom, gti=g)
        var = {}
        for off in OFFSETS:
            var[f"offset_{off}"] = npulse_multires(be, t90s, t90, prom, offsets=(off,), gti=g)
        var["bkg_alt"] = npulse_multires(be, t90s, t90, prom, bkg="alt", gti=g)
        var["bkg_nominal_single_offset"] = npulse_multires(
            be, t90s, t90, prom, offsets=(0.0,), gti=g
        )
        net = {d: float(np.sum((t > t90s) & (t < t90s + t90))) for d, (t, _, _) in evs.items()}
        top2 = sorted(net, key=net.get, reverse=True)[:2]
        var["dets_top2"] = npulse_multires(
            band_events({d: evs[d] for d in top2}), t90s, t90, prom, gti=g
        )
        var["dets_all_mask"] = nom
        # remaining variants up to 12: offsets x background alt
        for off in OFFSETS[1:]:
            var[f"offset_{off}_bkg_alt"] = npulse_multires(
                be, t90s, t90, prom, offsets=(off,), bkg="alt", gti=g
            )
        var["dets_top2_bkg_alt"] = npulse_multires(
            band_events({d: evs[d] for d in top2}), t90s, t90, prom, bkg="alt", gti=g
        )
        vals = list(var.values())[:12]
        ok = exposure_ok(evs, t90s, t90)
        out["events"][name] = {
            "t90": t90,
            "nominal": nom,
            "variants": var,
            "n_variants_ge2": int(sum(v >= 2 for v in vals)),
            "exposure_ok": ok,
            "robust": bool(ok and nom >= 2 and sum(v >= 2 for v in vals) >= 10),
        }
        # save reduced light curves (nominal, 8 ms) for reuse
        f, e = curve(be, t90s, t90, 0.008, 0.0, gti=g)
        np.savez_compressed(
            OUT / f"lc_{name}.npz", flux=f, err=e, dt=0.008, t90_start=t90s, t90=t90, dets=list(evs)
        )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "robustness.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 3), round(c + h, 3)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=("fetch", "robust", "epeak", "candidates"))
    ap.add_argument("--plan", default=str(CACHE.parent / "tte_plan.json"))
    ap.add_argument("--nsim", type=int, default=10)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    if a.cmd == "fetch":
        fetch(a)
    elif a.cmd == "robust":
        robust(a)
    elif a.cmd == "epeak":
        epeak(a)
    else:
        candidates(a)
    return 0


SYS_LN = 0.10  # ASSUMPTION: per-burst systematic in ln E_peak (cross-detector / model scatter)


def catalogue_epeak():
    """ln E_peak and its 1-sigma (mean of the asymmetric log errors) from fermigbrst's best-fitting
    fluence model (Band or Comptonized); bursts with other best models get no E_peak."""
    import io
    import urllib.parse

    from astropy.io.votable import parse_single_table

    cols = ["trigger_name", "flnc_best_fitting_model"] + [
        f"flnc_{m}_epeak{x}" for m in ("band", "comp") for x in ("", "_pos_err", "_neg_err")
    ]
    q = f"SELECT {', '.join(cols)} FROM fermigbrst"
    data = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": q, "MAXREC": "20000"}
    ).encode()
    with urllib.request.urlopen(urllib.request.Request(bs.TAP, data=data), timeout=300) as r:
        raw = r.read()
    t = parse_single_table(io.BytesIO(raw)).to_table()
    out = {}
    for row in t:
        m = str(row["flnc_best_fitting_model"]).strip().lower()
        key = "band" if "band" in m else ("comp" if "comp" in m else None)
        if key is None:
            continue
        e = float(row[f"flnc_{key}_epeak"])
        ep, en = float(row[f"flnc_{key}_epeak_pos_err"]), float(row[f"flnc_{key}_epeak_neg_err"])
        if not (np.isfinite(e) and e > 0 and np.isfinite(ep) and np.isfinite(en)):
            continue
        sig = 0.5 * (math.log((e + abs(ep)) / e) + abs(math.log(max(e - abs(en), 1e-3) / e)))
        out[str(row["trigger_name"])] = (math.log(e), sig, key)
    return out, hashlib.sha256(raw).hexdigest()[:16], len(t)


def zstat(lna, sa, lnb, sb, s, sig_lns=0.01, sys_ln=SYS_LN):
    """z for 'b is a's source stretched in time by 1/s' (pair convention of s1_m2_search.m2_only:
    y = burst b stretched by s matches x = burst a, so E_b = E_a * s): ln E_a - ln E_b + ln s."""
    var = sa * sa + sb * sb + 2 * sys_ln * sys_ln + sig_lns * sig_lns
    return (lna - lnb + math.log(s)) / math.sqrt(var)


def epeak(a):
    import s1_m2_search as ms

    sd = ms.sd
    ep, sha, nrow = catalogue_epeak()
    cat, B, elig, ks, sep, incons, allowed, cal, val = ms.setup()
    names = [b["name"] for b in B]
    pos = {n: i for i, n in enumerate(names)}

    def refined_s(i, j):
        bi, bj = B[i], B[j]
        x, ex, y, ey = sd.to_common(
            bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"]
        )
        m, s_grid = ms.m2_only(bi["flux"], bi["errs"], bi["k"], bj["flux"], bj["errs"], bj["k"])
        return ms.refined_p(x, ex, y, ey, s_grid)[1]

    # null 1: the 64 full-catalogue pairs above the M2 threshold
    srch = json.loads((ROOT / "results/s1_m2/search.json").read_text())
    null = [(pos[r["burst1"]], pos[r["burst2"]], "flagged64") for r in srch["flagged"]]
    # null 2: the 2,000 most similar calibration-half pairs
    bg = ms.background(B, allowed, cal, a.cpu)
    order = np.argsort(np.nan_to_num(bg[:, 2], nan=-np.inf))[::-1][:2000]
    null += [(int(bg[q, 0]), int(bg[q, 1]), "cal_top2000") for q in order]
    rows = []
    for i, j, src in null:
        ni, nj = names[i], names[j]
        if ni not in ep or nj not in ep:
            rows.append({"src": src, "a": ni, "b": nj, "has_epeak": False})
            continue
        s = refined_s(i, j)
        z = zstat(ep[ni][0], ep[ni][1], ep[nj][0], ep[nj][1], s)
        rows.append(
            {
                "src": src,
                "a": ni,
                "b": nj,
                "has_epeak": True,
                "s": round(s, 4),
                "z": round(z, 3),
                "pass": abs(z) <= 2.0,
            }
        )
    # true pairs (simulated) from each null pair's first member
    rng = np.random.default_rng(99)
    sims = {}
    for sys_true in (0.10, 0.20):
        res = []
        for r in rows:
            if not r["has_epeak"]:
                continue
            lna, sa, _ = ep[r["a"]]
            s0 = float(rng.choice(sd.S_INJ))
            ratio = float(rng.choice(sd.RATIOS))
            sb = min(1.0, sa / math.sqrt(ratio))
            lna_m = lna + rng.normal(0, sa) + rng.normal(0, sys_true)
            lnb_m = lna - math.log(s0) + rng.normal(0, sb) + rng.normal(0, sys_true)
            s_hat = (1.0 / s0) * math.exp(rng.normal(0, 0.01))
            res.append(abs(zstat(lna_m, sa, lnb_m, sb, s_hat)) <= 2.0)
        sims[str(sys_true)] = {
            "n": len(res),
            "pass": int(sum(res)),
            "rate": round(float(np.mean(res)), 3),
            "wilson95": wilson(int(sum(res)), len(res)),
        }
    out = {
        "catalogue": {"rows": nrow, "with_epeak_best_model": len(ep), "sha256": sha},
        "sys_ln_assumed": SYS_LN,
        "true_pairs": sims,
    }
    for src in ("flagged64", "cal_top2000"):
        rr = [r for r in rows if r["src"] == src]
        ok = [r for r in rr if r["has_epeak"]]
        k = sum(r["pass"] for r in ok)
        out[src] = {
            "pairs": len(rr),
            "with_epeak": len(ok),
            "pass_abs_z_le_2": k,
            "rate": round(k / max(len(ok), 1), 3),
            "wilson95": wilson(k, len(ok)),
            "median_abs_z": round(float(np.median([abs(r["z"]) for r in ok])), 2) if ok else None,
        }
    out["rows_flagged64"] = [r for r in rows if r["src"] == "flagged64"]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "epeak_test.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows_flagged64"}, indent=1))


def tte_pair_curves(ev_i, ev_j, ci, cj):
    """Both bursts from their own events at one common bin width (no rebinning of reduced data)."""
    p = dataclasses.replace(bt.Params(), base_dt=0.002, max_bins=256, pad_frac=0.5, pad_s=0.1)
    dti = bt.window_for(ci[0], ci[1], p)[2]
    dtj = bt.window_for(cj[0], cj[1], p)[2]
    dt = max(dti, dtj)
    fi, ei = curve(band_events(ev_i), ci[0], ci[1], dt, 0.0, gti=common_gti(ev_i))
    fj, ej = curve(band_events(ev_j), cj[0], cj[1], dt, 0.0, gti=common_gti(ev_j))
    return fi, ei, fj, ej, dt


def candidates(a):
    import s1_m2_search as ms

    sd = ms.sd
    cat, bur, idx = load_cat()
    rob = json.loads((OUT / "robustness.json").read_text())
    man = json.loads((OUT / "manifest.json").read_text())
    robust = [n for n, v in rob["events"].items() if v["robust"]]
    thr = json.loads((ROOT / "results/s1_m2/frozen.json").read_text())["threshold_m2"]
    ep, _, _ = catalogue_epeak()
    P = bt.Params()
    names = [str(n) for n in cat["trigger_name"]]
    ra, dec, err = (np.array(cat[c], float) for c in ("ra", "dec", "error_radius"))
    t90 = np.array(cat["t90"], float)
    elig = {str(r["trigger_name"]): bool(r["eligible"]) for r in bur}

    def compat(i, j):
        sep = float(bt.angsep_deg(ra[i], dec[i], ra[j], dec[j]))
        incons = bool(bt.position_inconsistent(sep, err[i], err[j], P.sys_deg, P.pos_nsigma))
        rr = t90[j] / max(t90[i], 1e-3)
        return incons and 0.5 <= rr <= 2.0

    counts = {}
    for n in robust:
        i = idx[n]
        js = [j for j in range(len(names)) if j != i and compat(i, j)]
        counts[n] = {
            "compatible_all_catalogue": len(js),
            "compatible_bcat_eligible": sum(elig.get(names[j], False) for j in js),
            "compatible_with_tte_downloaded": sum(
                names[j] in man["groups"]["elig_short"] for j in js
            ),
        }
    pool = man["groups"]["elig_short"]

    def score(ni, nj):
        ci = (float(cat[idx[ni]]["t90_start"]), float(cat[idx[ni]]["t90"]))
        cj = (float(cat[idx[nj]]["t90_start"]), float(cat[idx[nj]]["t90"]))
        fi, ei, fj, ej, dt = tte_pair_curves(load_events(ni), load_events(nj), ci, cj)
        best, sb = -1.0, 1.0
        for s_ in sd.S_GRID:
            ys = fj if s_ == 1.0 else bt.stretch(fj, ej, s_)[0]
            r, _ = sd.rho_free(fi, ys)
            if r > best:
                best, sb = r, s_
        p_, s_ref = ms.refined_p(fi, ei, fj, ej, sb)
        z = (
            zstat(ep[ni][0], ep[ni][1], ep[nj][0], ep[nj][1], s_ref)
            if ni in ep and nj in ep
            else None
        )
        dtd = abs(float(cat[idx[ni]]["trigger_time"]) - float(cat[idx[nj]]["trigger_time"]))
        return {
            "a": ni,
            "b": nj,
            "dt_s": dt,
            "m2": round(best, 4),
            "s": round(s_ref, 4),
            "chi2_p": float(f"{p_:.3g}"),
            "epeak_z": None if z is None else round(z, 2),
            "delay_days": round(dtd, 2),
        }

    # background: compatible pairs among the downloaded eligible short bursts (unrelated)
    bgr = []
    for x in range(len(pool)):
        for y in range(x + 1, len(pool)):
            if compat(idx[pool[x]], idx[pool[y]]):
                bgr.append(score(pool[x], pool[y]))
    pairs = []
    for n in robust:
        for m in pool + [r for r in robust if r != n]:
            if compat(idx[n], idx[m]) and (m not in robust or n < m):
                pairs.append(score(n, m))
    bm = np.array([r["m2"] for r in bgr]) if bgr else np.array([np.nan])

    def flag(r):
        return (
            r["m2"] > thr
            and r["chi2_p"] >= 1e-3
            and r["delay_days"] >= 1.0
            and (r["epeak_z"] is None or abs(r["epeak_z"]) <= 2.0)
        )

    out = {
        "robust_events": robust,
        "compatible_counts": counts,
        "threshold_m2_frozen_bcat": thr,
        "background_tte_pairs": len(bgr),
        "background_m2_quantiles": {q: float(np.nanquantile(bm, q)) for q in (0.5, 0.9, 0.99)},
        "background_passing_chain": [r for r in bgr if flag(r)],
        "scored_new_pairs": len(pairs),
        "new_pairs": sorted(pairs, key=lambda r: -r["m2"]),
        "candidates": [r for r in pairs if flag(r)],
    }
    (OUT / "candidates.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "new_pairs"}, indent=1))
    print("top new pairs:", out["new_pairs"][:6])


if __name__ == "__main__":
    raise SystemExit(main())
