"""S1 blind spots: can single-pulse bursts be searched (Step A), and do short bursts gain pulses in
TTE (Step B)? Plan: docs/s1_blind/plan.md.

  python scripts/s1_blind_spots.py single-pulse [--cpu 4]  -> results/s1_blind/single_pulse.json
  python scripts/s1_blind_spots.py tte-pilot [--n 30]      -> results/s1_blind/tte_pilot.json

TTE files are streamed one burst at a time into a scratch directory and deleted after reduction.

Injected twins are simulated; recovering them measures search capability only.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import dataclasses
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import s1_m2_search as ms  # noqa: E402  (sets the frozen fine grid and the refined chi2)

from jwst_anomaly import burst_twins as bt  # noqa: E402

sd = ms.sd
OUT = ROOT / "results" / "s1_blind"
MIN_SNR = 10.0


def single_pulse(a):
    t0 = time.time()
    fz = json.loads((ROOT / "results" / "s1_m2" / "frozen.json").read_text())
    sd.P = dataclasses.replace(sd.P, min_pulses=1)  # injections: copies need >= 1 resolved pulse
    cat, B, status = sd.st.load()
    P = sd.P
    npk = np.array([b["npulse"] for b in B])
    snr = np.array([b["snr"] for b in B])
    ks = np.array([b["k"] for b in B])
    pool = (npk == 1) & (snr >= MIN_SNR)
    ra, dec, err = (np.array([b[c] for b in B]) for c in ("ra", "dec", "err"))
    sep = bt.angsep_deg(ra[:, None], dec[:, None], ra[None], dec[None])
    incons = bt.position_inconsistent(sep, err[:, None], err[None], P.sys_deg, P.pos_nsigma)
    allowed = incons & pool[:, None] & pool[None] & (np.abs(ks[:, None] - ks[None]) <= 2)
    np.fill_diagonal(allowed, False)
    rng = np.random.default_rng(2)
    perm = rng.permutation(np.flatnonzero(pool))
    cal, val = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2 :])
    iu = np.triu_indices(len(B), 1)
    out = {
        "pool_npulse1_snr_ge10": int(pool.sum()),
        "npulse1_all": int((npk == 1).sum()),
        "pairs_all_single_pulse_allowed": int(allowed[iu].sum()),
        "n_cal": int(cal.size),
        "n_val": int(val.size),
    }
    bg = ms.background(B, allowed, cal, a.cpu)
    v = np.sort(np.nan_to_num(bg[:, 2], nan=-np.inf))[::-1]
    thr = float(v[max(int(round(sd.FAR * bg.shape[0])) - 1, 0)])
    out["calibration"] = {
        "pairs": int(bg.shape[0]),
        "threshold_m2": thr,
        "threshold_multi_pulse_frozen": fz["threshold_m2"],
    }
    order = np.argsort(np.nan_to_num(bg[:, 2], nan=-np.inf))[::-1][:2000]
    ms._init({"B": B})
    with cf.ProcessPoolExecutor(a.cpu, initializer=ms._init, initargs=({"B": B},)) as ex:
        jobs = [(int(bg[q, 0]), int(bg[q, 1])) for q in order]
        res = list(ex.map(ms._chi2_job, jobs, chunksize=10))
    pv = np.array([r[2] for r in res])
    dtd = np.array([abs(B[int(bg[q, 0])]["mjd"] - B[int(bg[q, 1])]["mjd"]) for q in order])
    passed = (pv >= fz["consistency_cut"]) & (dtd >= 1.0)
    rate = float(passed.mean())
    above = np.nan_to_num(bg[order, 2], nan=-np.inf) > thr  # pairs actually above the threshold
    n_flag_full = sd.FAR * out["pairs_all_single_pulse_allowed"]
    out["unrelated"] = {
        "n_top": int(order.size),
        "m2_range": [float(bg[order[-1], 2]), float(bg[order[0], 2])],
        "pass_rate": rate,
        "pass_top_100": int(passed[:100].sum()),
        "n_above_threshold": int(above.sum()),
        "pass_above_threshold": int(passed[above].sum()),
        "expected_flagged_full_single_pulse": round(n_flag_full, 1),
        "expected_false_candidates_full": round(n_flag_full * rate, 1),
    }
    inj = ms.run_injections(B, val, incons, None, 4343, a.per_cell, a.cpu)
    el = inj["elig"] > 0
    rec = el & (np.nan_to_num(inj["m2"], nan=-1) > thr)
    rec &= (np.nan_to_num(inj["p2"]) >= fz["consistency_cut"]) & (inj["dtd"] >= 1.0)
    d = inj["s0"] != 1.0
    out["injections"] = {
        "n": int(el.size),
        "deformed": {
            "n": int(d.sum()),
            "selected": int((d & el).sum()),
            "recovered": int((d & rec).sum()),
        },
        "undeformed": {
            "n": int((~d).sum()),
            "selected": int((~d & el).sum()),
            "recovered": int((~d & rec).sum()),
        },
    }
    tot = out["injections"]["deformed"]["recovered"] / max(out["injections"]["deformed"]["n"], 1)
    out["decision"] = {
        "a_expected_false_le_3": out["unrelated"]["expected_false_candidates_full"] <= 3,
        "b_deformed_total_recovery_ge_0.10": tot >= 0.10,
        "deformed_total_recovery": round(tot, 3),
    }
    out["runtime_s"] = round(time.time() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "single_pulse.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
TRIG = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers"
DETS = ["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb", "b0", "b1"]


def bcat_masks(names):
    """bcat_detector_mask per trigger from HEASARC TAP (fermigbrst), one query."""
    import io
    import urllib.parse
    import urllib.request

    from astropy.io.votable import parse_single_table

    q = "SELECT trigger_name, bcat_detector_mask FROM fermigbrst WHERE trigger_name IN ({})".format(
        ",".join(f"'{n}'" for n in names)
    )
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": q}).encode()
    with urllib.request.urlopen(urllib.request.Request(TAP, data=data), timeout=120) as r:
        t = parse_single_table(io.BytesIO(r.read())).to_table()
    return {str(a): str(b) for a, b in zip(t["trigger_name"], t["bcat_detector_mask"], strict=True)}


def tte_lightcurve(files, trig_t90_start, t90, base_dt=0.002, max_bins=256):
    """Background-subtracted two-band count-rate curves from NaI TTE (50-300 keV; 8-50 + 300-900).

    Window and resolution as bt.window_for with base_dt = 2 ms. Background: linear fit to
    [-20, -2] s and [t90_start + t90 + 5, + 30] s. Errors: Poisson (counts) / dt.
    """
    from astropy.io import fits

    # short-burst window (pilot setting): pad = 0.5 T90 + 0.1 s, so 256 bins reach 2-16 ms
    p = dataclasses.replace(
        bt.Params(), base_dt=base_dt, max_bins=max_bins, pad_frac=0.5, pad_s=0.1
    )
    lo, hi, dt = bt.window_for(trig_t90_start, t90, p)
    n = int(round((hi - lo) / dt))
    edges = lo + dt * np.arange(n + 1)
    cnt = np.zeros((2, n))
    bkg = np.zeros((2, n))
    for f in files:
        with fits.open(f) as h:
            eb = h["EBOUNDS"].data
            ev = h["EVENTS"].data
            t = np.asarray(ev["TIME"], float) - float(h[0].header["TRIGTIME"])
            emid = 0.5 * (np.asarray(eb["E_MIN"], float) + np.asarray(eb["E_MAX"], float))
            e = emid[np.asarray(ev["PHA"], int)]
        sel_b = [(e >= 50) & (e < 300), ((e >= 8) & (e < 50)) | ((e >= 300) & (e < 900))]
        for k, sb in enumerate(sel_b):
            tk = t[sb]
            cnt[k] += np.histogram(tk, edges)[0]
            # background: counts per second in side windows, linear in time
            side = [(-20.0, -2.0), (trig_t90_start + t90 + 5, trig_t90_start + t90 + 30)]
            xs, ys = [], []
            for a, b in side:
                h1, e1 = np.histogram(tk, np.arange(a, b + 1e-9, 1.0))
                xs += list(0.5 * (e1[1:] + e1[:-1]))
                ys += list(h1)
            if len(xs) >= 4:
                c = np.polyfit(xs, ys, 1)
                bkg[k] += np.polyval(c, 0.5 * (edges[1:] + edges[:-1])) * dt
    flux = (cnt - bkg) / dt
    err = np.sqrt(np.maximum(cnt, 1.0)) / dt
    return flux, err, dt


def tte_pilot(a):
    import shutil
    import urllib.request

    from astropy.table import Table

    t0 = time.time()
    cat = Table.read(ROOT / "results" / "s1_twins" / "catalogue.ecsv.gz", format="ascii.ecsv")
    bur = Table.read(ROOT / "results" / "s1_twins" / "bursts.ecsv.gz", format="ascii.ecsv")
    idx = {str(n): i for i, n in enumerate(cat["trigger_name"])}
    want_eligible = bool(
        a.control
    )  # control: short bursts that bcat already resolves (>= 2 pulses)
    sel = [
        r
        for r in bur
        if bool(r["eligible"]) == want_eligible
        and float(cat[idx[str(r["trigger_name"])]]["t90"]) < 2.0
        and (want_eligible or float(r["snr_peak"]) > 20.0)
    ]
    rng = np.random.default_rng(5)
    pick = [sel[i] for i in rng.choice(len(sel), size=min(a.n, len(sel)), replace=False)]
    names = [str(r["trigger_name"]) for r in pick]
    masks = bcat_masks(names)
    scratch = Path(os.environ.get("S1_TTE_SCRATCH", "/tmp/s1_tte"))
    rows, nbytes = [], 0
    for r in pick:
        name = str(r["trigger_name"])
        c = cat[idx[name]]
        year = "20" + name[2:4]
        mask = masks.get(name, "")
        dets = [d for d, m in zip(DETS[:12], mask[:12], strict=True) if m == "1"]
        d0 = time.time()
        bdir = scratch / name
        bdir.mkdir(parents=True, exist_ok=True)
        files, ok = [], True
        for d in dets:
            got = False
            for v in ("v00", "v01", "v02"):
                url = f"{TRIG}/{year}/{name}/current/glg_tte_{d}_{name}_{v}.fit"
                dst = bdir / Path(url).name
                try:
                    urllib.request.urlretrieve(url, dst)
                    nbytes += dst.stat().st_size
                    files.append(dst)
                    got = True
                    break
                except Exception:
                    continue
            ok &= got
        row = {
            "trigger_name": name,
            "t90": float(c["t90"]),
            "bcat_npulse": int(r["npulse"]),
            "bcat_dt": float(r["dt"]),
            "dets": dets,
            "files_ok": bool(ok and files),
        }
        if files:
            try:
                flux, err, dt = tte_lightcurve(files, float(c["t90_start"]), float(c["t90"]))
                row["tte_dt"] = dt
                row["tte_npulse"] = int(len(bt.find_pulses(flux, err, bt.Params())))
                row["tte_snr"] = round(bt.snr_total(flux, err), 1)
                # multi-resolution count (fix after the control run): rebin by 2, 4, ... up to
                # 64 ms and keep the largest count of 5-sigma pulses
                best, m = row["tte_npulse"], 2
                while dt * m <= 0.064 + 1e-9 and flux.shape[1] // m >= 4:
                    f2, e2 = bt.rebin_factor(flux, err, m)
                    best = max(best, int(len(bt.find_pulses(f2, e2, bt.Params()))))
                    m *= 2
                row["tte_npulse_multires"] = best
            except Exception as exc:  # noqa: BLE001 (record and continue)
                row["error"] = repr(exc)[:200]
        row["seconds"] = round(time.time() - d0, 1)
        shutil.rmtree(bdir, ignore_errors=True)
        rows.append(row)
        print(row, flush=True)
    ok = [x for x in rows if "tte_npulse" in x]
    gain = sum(1 for x in ok if x["tte_npulse_multires"] >= 2)
    gain_single = sum(1 for x in ok if x["tte_npulse"] >= 2)
    out = {
        "selection": (
            "control: eligible (bcat >= 2 pulses) bursts with T90 < 2 s; seed 5"
            if a.control
            else "excluded (npulse < 2) bursts with T90 < 2 s and bcat peak S/N > 20; seed 5"
        ),
        "pool": len(sel),
        "n": len(rows),
        "reduced": len(ok),
        "gain_ge2_pulses_multires": gain,
        "gain_ge2_pulses_single_resolution": gain_single,
        "gain_fraction": round(gain / max(len(ok), 1), 3),
        "mb_streamed": round(nbytes / 1e6, 1),
        "seconds_total": round(time.time() - t0, 1),
        "rows": rows,
        "decision_continue": bool(len(ok) and gain / len(ok) >= 0.25),
    }
    per = out["seconds_total"] / max(len(rows), 1)
    out["extrapolated_hours_all_699_short"] = round(per * 699 / 3600, 2)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / ("tte_control.json" if a.control else "tte_pilot.json")).write_text(
        json.dumps(out, indent=1) + "\n"
    )
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=("single-pulse", "tte-pilot"))
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--control", action="store_true", help="eligible short bursts (pipeline check)")
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--per-cell", type=int, default=8)
    a = ap.parse_args(argv)
    {"single-pulse": single_pulse, "tte-pilot": tte_pilot}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
