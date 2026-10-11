"""E1: near-antipodal lag channels (D-074 next step).

Counts pairs whose separation is within ``max(ANTI_MIN_DEG, n_sigma * sigma_comb)`` of 180 deg, per
lag bin, in the six localized E1 channels, against the D-074 ``jit`` null (Dec and hour angle kept;
GBM in whole orbits). Trials: the 30 cells pooled (``en.global_p``) and an analytic Bonferroni tail.
Sensitivity: synthetic antipodal pairs (a B event moved to t_A +- lag and to A's antipode, with
the B localization error as scatter) injected per cell; ``n50`` is the smallest injected count
detected at family-wise 3 sigma in at least half of the trials.

Writes ``results/e1_events/antipodal.json``.

  python scripts/e1_antipodal.py [--n 1000] [--cpu 4]
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

OUT = ROOT / "results" / "e1_events" / "antipodal.json"
#: Localized channels (GW events have no positions).
CHANNELS = tuple((a, b) for a, b in E.CHANNELS if "GW" not in (a, b))
#: ASSUMPTION: a pair is antipodal when 180 - sep <= max(ANTI_MIN_DEG, n_sigma * sigma_comb).
ANTI_MIN_DEG = 10.0
#: ASSUMPTION: injected counts tried per cell, and trials per count.
INJ_N = (3, 10, 30, 100)
INJ_TRIALS = 20
P = en.Params()
NLAG = len(P.lag_edges) - 1


def antipodal(sep, sig_a, sig_b, p: en.Params = P) -> np.ndarray:
    sc = np.hypot(sig_a, sig_b)
    return (180.0 - sep) <= np.maximum(ANTI_MIN_DEG, p.n_sigma * sc)


def count_channel(a: en.Sample, b: en.Sample, same: bool, p: en.Params = P) -> np.ndarray:
    """Antipodal pair counts per lag bin (bin 0 includes lag 0)."""
    e = np.asarray(p.lag_edges)
    i, j = en.pairs_within(a.mjd, b.mjd, e[-1], same)
    lag = np.abs(b.mjd[j] - a.mjd[i]) * en.DAY
    sep = en.sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
    m = antipodal(sep, a.sigma[i], b.sigma[j], p) & (lag <= e[-1])
    k = np.clip(np.searchsorted(e, lag[m], side="right") - 1, 0, NLAG - 1)
    return np.bincount(k, minlength=NLAG)


def count_all(s: dict[str, en.Sample]) -> np.ndarray:
    """Counts, shape (n_channels, n_lag)."""
    return np.stack([count_channel(s[a], s[b], a == b) for a, b in CHANNELS])


def inject_antipodal(a, b, n, lag_lo, lag_hi, same, rng, p: en.Params = P) -> en.Sample:
    """Copy of ``b`` with up to n events moved to t_A +- lag (log-uniform) and to the antipode of a
    distinct anchor A within ``p.inject_local_days``, scattered by B's own 1-sigma error
    (provenance: simulated)."""
    mjd, ra, dec = b.mjd.copy(), b.ra.copy(), b.dec.copy()
    used: set[int] = set()
    moved = 0
    for j in rng.permutation(len(b.mjd)):
        if moved >= n:
            break
        near = np.flatnonzero(np.abs(a.mjd - b.mjd[j]) <= p.inject_local_days)
        near = (
            [i for i in near if i not in used and i != j]
            if same
            else [i for i in near if i not in used]
        )
        if not near:
            continue
        i = int(rng.choice(near))
        lag = np.exp(rng.uniform(np.log(max(lag_lo, 1e-3)), np.log(lag_hi)))
        mjd[j] = a.mjd[i] + rng.choice((-1.0, 1.0)) * lag / en.DAY
        s = b.sigma[j] if np.isfinite(b.sigma[j]) else 1.0
        dra, ddec = rng.normal(0, s, 2)
        dec[j] = np.clip(-a.dec[i] + ddec, -90, 90)
        ra[j] = np.mod(a.ra[i] + 180.0 + dra / max(np.cos(np.radians(dec[j])), 0.05), 360)
        used.update((i, int(j)))
        moved += 1
    return en.Sample(b.cat, mjd, ra, dec, b.sigma, b.year)


_W: dict = {}


def _init(s):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W["s"] = s


def _null_chunk(seeds):
    s = _W["s"]
    return np.stack(
        [
            count_all({k: en.scramble_jit(v, np.random.default_rng(sd), P) for k, v in s.items()})
            for sd in seeds
        ]
    )


def run_null(s, n, cpu, base_seed=7_000_000):
    seeds = np.arange(n) + base_seed
    chunks = [seeds[i : i + 50] for i in range(0, n, 50)]
    if cpu <= 1:
        _init(s)
        return np.concatenate([_null_chunk(c) for c in chunks])
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(s,)) as ex:
        return np.concatenate(list(ex.map(_null_chunk, chunks)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    t0 = time.time()
    ev = E.load(E.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items() if k != "GW"}
    obs = count_all(s)
    null = run_null(s, a.n, a.cpu)
    mask = np.ones(obs.shape, bool)
    pmin, pglob = en.global_p(obs, null, mask)
    pa = en.analytic_p(obs, null)
    ncell = int(mask.sum())
    mu, sd = null.mean(axis=0), null.std(axis=0)
    labels = en.lag_labels(P)
    cells = []
    rng = np.random.default_rng(7_500_000)
    for c, (ca, cb) in enumerate(CHANNELS):
        for k in range(NLAG):
            n50 = None
            for ninj in INJ_N:
                det = 0
                for _ in range(INJ_TRIALS):
                    sb = inject_antipodal(
                        s[ca], s[cb], ninj, P.lag_edges[k], P.lag_edges[k + 1], ca == cb, rng
                    )
                    same = ca == cb  # a same-catalogue channel counts sb with itself
                    o = count_channel(sb if same else s[ca], sb, same)[k]
                    p1 = en.analytic_p(np.array([o]), null[:, c, k][:, None])[0]
                    det += p1 * ncell <= E.DETECT_P
                if det >= INJ_TRIALS / 2:
                    n50 = ninj
                    break
            cells.append(
                {
                    "channel": f"{ca}-{cb}",
                    "lag": labels[k],
                    "obs": int(obs[c, k]),
                    "null_mean": round(float(mu[c, k]), 2),
                    "null_sd": round(float(sd[c, k]), 2),
                    "z": round(float(en.z_score(obs[c, k], null[:, c, k])), 2),
                    "analytic_p": float(pa[c, k]),
                    "ul95_extra_pairs": round(
                        float(en.upper_limit_95(obs[c, k : k + 1], null[:, c, k : k + 1])[0]), 1
                    ),
                    "n50_injected": n50,
                }
            )
    out = {
        "test": "E1 antipodal channels (jit null)",
        "anti_min_deg": ANTI_MIN_DEG,
        "n_sigma": P.n_sigma,
        "n_scrambles": a.n,
        "n_cells": ncell,
        "min_cell_p_pooled": pmin,
        "global_p_pooled": pglob,
        "min_analytic_p_bonferroni": min(1.0, float(pa.min()) * ncell),
        "cells": cells,
        "runtime_s": round(time.time() - t0),
    }
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    for r in cells:
        print(r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
