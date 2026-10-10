"""E-NF1b: NF-H04 ghost pairs at short lags in IceTracks-DR2 (docs/neutrino_frontier/README.md).

1. Streams the 14 IceTracks-DR2 v3.1 season files (events + good-run lists, ~210 MB) from Dataverse,
   checks each against ``data/manifests/nf_icetracks_dr2.ecsv`` and keeps only a reduced ``.npz``
   under the data root (raw files are never written: cloud-disk decision).
2. Known case: the TXS 0506+056 2014-15 box (T0 = MJD 57020, 185 d) must show an on-source excess.
3. Forecast: the most optimistic reachable R_g per lag bin (every event astrophysical,
   max(3 sqrt(B), 3) / N with B the observed pairs in the bin) for northern energy cuts; it decides
   which bins can improve on the ICECAT-1 limits (E-NF1).
4. Test (pre-registered in the README before the run): northern tracks (Dec > -5 deg) with
   log10(E/GeV) >= 4.0 and >= 4.5, lag bins 0-10 s, 10-100 s, 100 s-1 h; wide pairs (D-074
   class 1), pairs of one (run, event) excluded; null = uptime-aware +-3 d jitter (hour angle
   kept); empirical p, family-wise 3 sigma over 6 cells; injected ghosts (one per event with
   probability R) and 95 % CLs limits on R per event above the cut.

Writes ``results/nf/ghost_pairs_dr2.json``.

  python scripts/nf_ghost_dr2.py [--n 5000] [--cpu 4] [--mjd-min 55694.4]

``--mjd-min 55694.4`` (IC86 seasons only; post hoc robustness check, not part of the
pre-registered family) writes ``results/nf/ghost_pairs_dr2_ic86.json``: the log10 E proxy rate
above the cuts differs by ~18x between IC59 and IC86 (CHANGELOG 2026-10-10).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly import neutrino_frontier as nf  # noqa: E402

OUT = ROOT / "results" / "nf" / "ghost_pairs_dr2.json"
MANIFEST = ROOT / "data" / "manifests" / "nf_icetracks_dr2.ecsv"
URL = "https://dataverse.harvard.edu/api/access/datafile/{}"
DAY = en.DAY
P = en.Params()
#: ASSUMPTIONs (pre-registered): sample, cuts, bins, R grid, trials
DEC_MIN = -5.0
CUTS = (4.0, 4.5)
LAG_EDGES = (0.0, 10.0, 100.0, 3600.0)
ALL_EDGES = (0.0, 10.0, 100.0, 3600.0, DAY, 7 * DAY, 30 * DAY, 180 * DAY)
FORECAST_CUTS = (3.5, 4.0, 4.5, 5.0, 5.5)
R_GRID = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
INJ_TRIALS = 200
DETECT_P = 1.35e-3
TXS = {"ra": 77.3582, "dec": 5.6931, "t0": 57020.0, "dt": 185.0, "log10e": 3.5, "r_deg": 1.0}


def fetch() -> tuple[np.ndarray, np.ndarray]:
    """(events, uptime) arrays; cached as one reduced .npz under the data root."""
    import requests
    from astropy.table import Table

    from jwst_anomaly.paths import data_root

    cache = data_root() / "nf_icetracks_dr2" / "dr2_v3.1_reduced.npz"
    if cache.exists():
        d = np.load(cache)
        return d["ev"], d["up"]
    man = Table.read(MANIFEST, format="ascii.ecsv")

    def get(row):
        r = requests.get(URL.format(int(row["dataverse_file_id"])), timeout=900)
        r.raise_for_status()
        got = hashlib.sha256(r.content).hexdigest()
        if got != row["sha256"]:
            raise SystemExit(f"{row['file']}: sha256 {got} != pinned {row['sha256']}")
        return str(row["file"]), nf.read_icetracks_tab(r.text.splitlines())

    with cf.ThreadPoolExecutor(8) as ex:
        parts = dict(ex.map(get, man))
    ev = nf.dedupe_events(np.vstack([v for k, v in parts.items() if k.startswith("events/")]))
    up = np.vstack([v for k, v in parts.items() if k.startswith("uptime/")])
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, ev=ev, up=up)
    return ev, up


def txs_check(ev: np.ndarray) -> dict:
    """On-source count in the TXS box against the same-Dec-band off-time expectation."""
    t, le, err = ev[:, 3], ev[:, 4], ev[:, 5]
    sep = en.sep_deg(ev[:, 6], ev[:, 7], TXS["ra"], TXS["dec"])
    box = np.abs(t - TXS["t0"]) <= TXS["dt"] / 2
    hi = le >= TXS["log10e"]
    near = hi & (sep < np.maximum(TXS["r_deg"], err))
    band = hi & (np.abs(ev[:, 7] - TXS["dec"]) < 3.0)
    frac = (near & ~box).sum() / max((band & ~box).sum(), 1)
    on, exp = int((near & box).sum()), float(frac * (band & box).sum())
    return {
        **TXS,
        "on": on,
        "expected": round(exp, 3),
        "p_poisson": float(stats.poisson.sf(on - 1, exp)),
    }


def forecast(ev: np.ndarray) -> list[dict]:
    north = ev[:, 7] > DEC_MIN
    out = []
    for c in FORECAST_CUTS:
        tt = np.sort(ev[north & (ev[:, 4] >= c), 3])
        n = len(tt)
        row = {"log10e_min": c, "n": n, "bins": []}
        for lo, hi in zip(ALL_EDGES[:-1], ALL_EDGES[1:], strict=True):
            b = int(
                (
                    np.searchsorted(tt, tt + hi / DAY, "right")
                    - np.searchsorted(tt, tt + lo / DAY, "right")
                ).sum()
            )
            row["bins"].append({"pairs": b, "rg_floor": max(3 * np.sqrt(b), 3.0) / max(n, 1)})
        out.append(row)
    return out


def sample_of(ev: np.ndarray) -> tuple[en.Sample, np.ndarray]:
    t = ev[:, 3]
    s = en.Sample("dr2", t, ev[:, 6], ev[:, 7], ev[:, 5], en.mjd_year(t))
    group = np.round(ev[:, 0]).astype(np.int64) * 10**8 + np.round(ev[:, 1]).astype(np.int64)
    return s, group


_W: dict = {}


def _init(samples, up):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W.update(samples=samples, up=up)


def _null_chunk(seeds):
    up = _W["up"]
    out = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        out.append(
            [
                nf.wide_pair_counts(
                    nf.jitter_uptime(s, up[:, 0], up[:, 1], P.jitter_s / DAY, rng), g, LAG_EDGES, P
                )
                for s, g in _W["samples"]
            ]
        )
    return np.array(out)


def empirical_p(obs, null):
    return (1.0 + (null >= obs[None] - 1e-9).sum(axis=0)) / (null.shape[0] + 1.0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--mjd-min", type=float, default=None, help="post hoc: events from this MJD on")
    a = ap.parse_args(argv)
    out_path = OUT if a.mjd_min is None else OUT.with_name("ghost_pairs_dr2_ic86.json")
    nb = len(LAG_EDGES) - 1
    thresh = DETECT_P / (nb * len(CUTS))
    if 1.0 / (a.n + 1) > thresh:
        raise SystemExit(f"--n {a.n} cannot reach the per-cell threshold {thresh:.2e}")
    t0 = time.time()
    ev, up = fetch()
    txs = txs_check(ev)
    fc = forecast(ev)
    north = ev[:, 7] > DEC_MIN
    if a.mjd_min is not None:
        north &= ev[:, 3] >= a.mjd_min
    samples = [sample_of(ev[north & (ev[:, 4] >= c)]) for c in CUTS]
    obs = np.array([nf.wide_pair_counts(s, g, LAG_EDGES, P) for s, g in samples])
    t1 = time.time()
    seeds = np.arange(a.n) + 9_700_000
    chunks = [seeds[i : i + 250] for i in range(0, a.n, 250)]
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init, initargs=(samples, up)) as ex:
        null = np.concatenate(list(ex.map(_null_chunk, chunks)))  # (n, cuts, bins)
    t_null = time.time() - t1
    pe = empirical_p(obs, null)
    mu, sd = null.mean(axis=0), null.std(axis=0)
    rng = np.random.default_rng(9_800_000)
    cells = []
    labels = en.lag_labels(en.Params(lag_edges=LAG_EDGES))
    for q, (c, (s, g)) in enumerate(zip(CUTS, samples, strict=True)):
        ones = np.ones(len(s.mjd))
        for k in range(nb):
            lo, hi = LAG_EDGES[k], LAG_EDGES[k + 1]
            det, cls = {}, {}
            p_b = (null[:, q, k] <= obs[q, k]).mean()
            for r in R_GRID:
                d = ab = 0
                for _ in range(INJ_TRIALS):
                    si, _w = nf.inject_ghosts(s, ones, r, lo, hi, rng)
                    gi = np.concatenate([g, -1 - np.arange(len(si.mjd) - len(s.mjd))])
                    n0, n1 = len(s.mjd), len(si.mjd)
                    excess = (
                        nf.wide_pair_counts(si, gi, LAG_EDGES, P)[k]
                        - (n1 * (n1 - 1) / (n0 * (n0 - 1)) - 1.0) * mu[q, k]
                        - obs[q, k]
                    )
                    o = null[rng.integers(len(null)), q, k] + excess
                    d += empirical_p(np.array([o]), null[:, q, k][:, None])[0] <= thresh
                    ab += o <= obs[q, k]
                det[r] = d / INJ_TRIALS
                cls[r] = (ab / INJ_TRIALS) / max(p_b, 1e-12)
            ok = [cls[r] <= 0.05 for r in R_GRID]
            cells.append(
                {
                    "log10e_min": c,
                    "n_events": len(s.mjd),
                    "lag": labels[k],
                    "obs": int(obs[q, k]),
                    "null_mean": round(float(mu[q, k]), 3),
                    "null_sd": round(float(sd[q, k]), 3),
                    "p_empirical": float(pe[q, k]),
                    "p_bonferroni": min(1.0, float(pe[q, k]) * nb * len(CUTS)),
                    "r50": next((r for r in R_GRID if det[r] >= 0.5), None),
                    "r_ul95_cls": next((r for i, r in enumerate(R_GRID) if all(ok[i:])), None),
                    "detect_frac": {str(r): v for r, v in det.items()},
                }
            )
    out = {
        "test": (
            "E-NF1b IceTracks-DR2 short-lag ghost pairs (NF-H04), wide separation, D-074 classes"
        ),
        "provenance": (
            "derived (observed IceTracks-DR2 v3.1); injections simulated; limits model_prediction"
        ),
        "limit_meaning": (
            "R = ghosts per event above the cut; per astrophysical neutrino R_g = R / f_astro"
        ),
        "mjd_min_post_hoc": a.mjd_min,
        "n_events_total": len(ev),
        "livetime_days_sum": round(float((up[:, 1] - up[:, 0]).sum()), 1),
        "txs_known_case": txs,
        "forecast_optimistic": fc,
        "n_scrambles": a.n,
        "per_cell_threshold": thresh,
        "min_p_bonferroni": min(1.0, float(pe.min()) * pe.size),
        "r_grid": list(R_GRID),
        "inj_trials": INJ_TRIALS,
        "speed": {"null_s": round(t_null, 1), "cpu": a.cpu},
        "cells": cells,
        "runtime_s": round(time.time() - t0),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(
        json.dumps(
            {k: v for k, v in out.items() if k not in ("cells", "forecast_optimistic")}, indent=1
        )
    )
    for row in fc:
        print(
            row["log10e_min"], row["n"], [f"{b['pairs']}/{b['rg_floor']:.2g}" for b in row["bins"]]
        )
    for cell in cells:
        print({k: v for k, v in cell.items() if k != "detect_frac"})
    return 0


if __name__ == "__main__":
    sys.exit(main())
