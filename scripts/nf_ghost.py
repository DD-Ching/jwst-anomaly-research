"""E-NF1: ICECAT-1 ghost-pair test and limit for NF-H04 (docs/neutrino_frontier/README.md).

Counts ICECAT-1 v4 wide-separation pairs (D-074 classes) in seven lag bins, plain and
signalness-weighted (14 cells). Nulls: D-074 ``jit`` for bins <= 7 d, cyclic +-1 yr jitter for
7-180 d. Per-cell p: empirical (1 + #null >= obs) / (N + 1); trials: the pooled-rank global p
(D-074) and Bonferroni over 14 cells. Sensitivity: injected ghosts
(``neutrino_frontier.inject_ghosts``) per bin and R_g; R_g,50 and a 95 % CLs upper limit on R_g.

Writes ``results/nf/ghost_pairs.json``.

  python scripts/nf_ghost.py [--n 20000] [--cpu 4]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as E  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly import neutrino_frontier as nf  # noqa: E402

OUT = ROOT / "results" / "nf" / "ghost_pairs.json"
GP = nf.GhostParams()
P = en.Params()
NB = nf.n_bins(GP)
STATS = ("count", "signal_weighted")
#: ASSUMPTION: R_g grid and trials per point for the sensitivity and the limit.
RG_GRID = (0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
INJ_TRIALS = 200


ICECAT_FILE = "icecat1_gold_bronze_tracks.csv"


def fetch_icecat() -> Path:
    """Only the ICECAT-1 file, pinned by the D-074 manifest (the other E1 catalogues are not needed
    here; the live HEASARC GBM table no longer matches its pin, see CHANGELOG 2026-10-10)."""
    import hashlib

    import requests
    from astropy.table import Table

    from jwst_anomaly.paths import data_root

    path = data_root() / "e1_events" / ICECAT_FILE
    if not path.exists():
        r = requests.get(E.SOURCES[ICECAT_FILE], timeout=600)
        r.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(r.content)
    pin = {row["file"]: row["sha256"] for row in Table.read(E.MANIFEST, format="ascii.ecsv")}
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != pin[ICECAT_FILE]:
        raise SystemExit(f"{ICECAT_FILE}: sha256 {got} != pinned {pin[ICECAT_FILE]}")
    return path


def load() -> tuple[en.Sample, np.ndarray, list[str]]:
    """ICECAT-1 v4 (CR_VETO dropped, D-074 loader) and its SIGNAL column (signalness)."""
    import pandas as pd

    df = pd.read_csv(fetch_icecat())
    t = en.icecat_events(df)
    df = df.drop_duplicates("NAME")
    names = [str(x) for x in t["name"]]
    w = df.set_index("NAME").loc[names, "SIGNAL"].to_numpy(float)
    return en.Sample.from_table(t), w, names


def stats_of(s: en.Sample, w: np.ndarray) -> np.ndarray:
    return nf.wide_pair_stats(s, w, GP, P)


_W: dict = {}


def _init(s, w):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W.update(s=s, w=w)


def _null_chunk(seeds):
    s, w = _W["s"], _W["w"]
    out = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        short = stats_of(en.scramble_jit(s, rng, P), w)
        long_ = stats_of(nf.scramble_cyclic(s, rng, GP.long_jitter_days), w)
        out.append(np.vstack([short[: GP.long_from_bin], long_[GP.long_from_bin :]]))
    return np.stack(out)


def run_null(s, w, n, cpu, base_seed=9_100_000):
    seeds = np.arange(n) + base_seed
    chunks = [seeds[i : i + 500] for i in range(0, n, 500)]
    if cpu <= 1:
        _init(s, w)
        return np.concatenate([_null_chunk(c) for c in chunks])
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(s, w)) as ex:
        return np.concatenate(list(ex.map(_null_chunk, chunks)))


def pair_ratio(w: np.ndarray, wi: np.ndarray) -> np.ndarray:
    """Ratio of random-pair totals after / before injection: n(n-1) for counts and
    (sum w)^2 - sum w^2 for the weighted statistic (shape (2,))."""
    n0, n1 = len(w), len(wi)
    c = n1 * (n1 - 1) / (n0 * (n0 - 1))
    q = (wi.sum() ** 2 - (wi**2).sum()) / (w.sum() ** 2 - (w**2).sum())
    return np.array([c, q])


def empirical_p(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    return (1.0 + (null >= obs[None] - 1e-9).sum(axis=0)) / (null.shape[0] + 1.0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20_000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    t0 = time.time()
    s, w, names = load()
    obs = stats_of(s, w)
    t1 = time.time()
    null = run_null(s, w, a.n, a.cpu)
    t_null = time.time() - t1
    ncell = obs.size
    pe = empirical_p(obs, null)
    mask = np.ones(obs.shape, bool)
    pmin, pglob = en.global_p(obs, null, mask)
    mu, sd = null.mean(axis=0), null.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 0, (obs - mu) / sd, 0.0)
    thresh = E.DETECT_P / ncell  # per-cell p for family-wise 3 sigma
    labels = [*en.lag_labels(en.Params(lag_edges=GP.lag_edges))]
    rng = np.random.default_rng(9_500_000)
    cells = []
    for k in range(NB):
        lo, hi = GP.lag_edges[k], GP.lag_edges[k + 1]
        det = {st: {} for st in STATS}
        above = {st: {} for st in STATS}
        for rg in RG_GRID:
            d = np.zeros(2)
            ab = np.zeros(2)
            for _ in range(INJ_TRIALS):
                si, wi = nf.inject_ghosts(s, w, rg, lo, hi, rng)
                # remove the random-pair increase from the extra events, so only parent-ghost
                # dependence is tested against the null of the observed catalogue
                excess = stats_of(si, wi)[k] - (pair_ratio(w, wi) - 1.0) * mu[k] - obs[k]
                # under H(R_g): a background realisation (a null draw) plus the ghost excess
                o = null[rng.integers(len(null)), k] + excess
                d += empirical_p(o[None], null[:, k][:, None, :])[0] <= thresh
                ab += o <= obs[k]  # as low as observed under H(R_g)
            p_b = (null[:, k] <= obs[k]).mean(axis=0)  # as low as observed with R_g = 0
            for q, st in enumerate(STATS):
                det[st][rg] = d[q] / INJ_TRIALS
                # CLs = P(<= obs | R_g) / P(<= obs | 0): no exclusion from a low fluctuation alone
                above[st][rg] = (ab[q] / INJ_TRIALS) / max(p_b[q], 1e-12)
        for q, st in enumerate(STATS):
            r50 = next((rg for rg in RG_GRID if det[st][rg] >= 0.5), None)
            ul95 = next((rg for rg in RG_GRID if above[st][rg] <= 0.05), None)
            cells.append(
                {
                    "lag": labels[k],
                    "stat": st,
                    "null": "jit" if k < GP.long_from_bin else "cyclic_1yr",
                    "obs": round(float(obs[k, q]), 3),
                    "null_mean": round(float(mu[k, q]), 3),
                    "null_sd": round(float(sd[k, q]), 3),
                    "z": round(float(z[k, q]), 2),
                    "p_empirical": float(pe[k, q]),
                    "p_bonferroni": min(1.0, float(pe[k, q]) * ncell),
                    "rg50": r50,
                    "rg_ul95": ul95,
                    "detect_frac": {str(r): v for r, v in det[st].items()},
                }
            )
    out = {
        "test": "E-NF1 ICECAT-1 ghost pairs (NF-H04), wide separation, D-074 classes",
        "provenance": (
            "derived (observed ICECAT-1 v4); injections simulated; limits model_prediction"
        ),
        "n_events": len(names),
        "n_scrambles": a.n,
        "n_cells": ncell,
        "per_cell_threshold": thresh,
        "min_cell_p_pooled": pmin,
        "global_p_pooled": pglob,
        "min_p_bonferroni": min(1.0, float(pe.min()) * ncell),
        "reproduces_d074": {"1h-1d_wide": int(obs[3, 0]), "1d-7d_wide": int(obs[4, 0])},
        "rg_grid": list(RG_GRID),
        "inj_trials": INJ_TRIALS,
        "speed": {
            "null_s": round(t_null, 1),
            "scrambles_per_s": round(a.n / t_null, 1),
            "cpu": a.cpu,
        },
        "cells": cells,
        "runtime_s": round(time.time() - t0),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    for c in cells:
        print({k: v for k, v in c.items() if k != "detect_frac"})
    return 0


if __name__ == "__main__":
    sys.exit(main())
