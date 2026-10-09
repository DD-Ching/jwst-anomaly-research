"""E1: signed-lag ("which event comes first") asymmetry in the cross-catalogue channels.

For each cross channel (A, B), lag bin and separation class, D = N(t_B > t_A) - N(t_B < t_A).
Ordinary event pairs and the D-074 ``jit`` null are symmetric in sign, so a causal A -> B link at
a lag no ordinary path explains would show as D > 0 (or B -> A as D < 0). Two-sided: |D - null
mean| ranked with ``en.global_p`` over all cells (trials pooled with the scrambles), and an
analytic Gaussian tail x Bonferroni. Sensitivity: B events moved to t_A + lag (one sign only) of
an A anchor in the tested class; ``n50`` is the smallest injected count detected at family-wise
3 sigma in at least half the trials.

Writes ``results/e1_events/signed_lag.json``.

  python scripts/e1_signed_lag.py [--n 2000] [--cpu 4]
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
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as E  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402

OUT = ROOT / "results" / "e1_events" / "signed_lag.json"
CHANNELS = tuple((a, b) for a, b in E.CHANNELS if a != b)
P = en.Params()
NLAG = len(P.lag_edges) - 1
CLASSES = ("same", "wide", "all")
#: ASSUMPTION: injected counts tried per cell, and trials per count.
INJ_N = (3, 10, 30, 100)
INJ_TRIALS = 20


def signed_counts(a: en.Sample, b: en.Sample, p: en.Params = P) -> np.ndarray:
    """D = N(t_B > t_A) - N(t_B < t_A), shape (n_lag, 3 classes); bin 0 includes lag 0 (sign 0
    pairs count for neither side)."""
    e = np.asarray(p.lag_edges)
    i, j = en.pairs_within(a.mjd, b.mjd, e[-1], False)
    dt = (b.mjd[j] - a.mjd[i]) * en.DAY
    lag = np.abs(dt)
    keep = lag <= e[-1]
    i, j, dt, lag = i[keep], j[keep], dt[keep], lag[keep]
    sep = en.sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
    cls = en.sep_class(sep, a.sigma[i], b.sigma[j], p)
    k = np.clip(np.searchsorted(e, lag, side="right") - 1, 0, NLAG - 1)
    sgn = np.sign(dt).astype(np.int64)
    out = np.zeros((NLAG, 3), np.int64)
    for c, m in enumerate((cls == 0, cls == 1, np.ones(len(cls), bool))):
        out[:, c] = np.bincount(k[m], weights=sgn[m], minlength=NLAG).astype(np.int64)
    return out


def cell_mask(ch=CHANNELS) -> np.ndarray:
    """Pre-registered cells: same/wide for localized channels, all for GW channels."""
    m = np.zeros((len(ch), NLAG, 3), bool)
    for c, (a, b) in enumerate(ch):
        m[c, :, 2 if "GW" in (a, b) else slice(0, 2)] = True
    return m


def inject_after(a, b, n, lo, hi, cls_want, rng, p: en.Params = P) -> en.Sample:
    """Copy of ``b`` with up to n events moved to t_A + lag (lag log-uniform in [lo, hi], B always
    after A) of distinct anchors within ``p.inject_local_days`` whose pair class is ``cls_want``
    (2 = any). Positions kept (provenance: simulated)."""
    mjd = b.mjd.copy()
    used: set[int] = set()
    moved = 0
    for j in rng.permutation(len(b.mjd)):
        if moved >= n:
            break
        near = np.flatnonzero(np.abs(a.mjd - b.mjd[j]) <= p.inject_local_days)
        for i in rng.permutation(near)[:20]:
            if i in used:
                continue
            if cls_want != 2:
                sep = en.sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
                c = en.sep_class(np.atleast_1d(sep), a.sigma[i : i + 1], b.sigma[j : j + 1], p)
                if c[0] != cls_want:
                    continue
            lag = np.exp(rng.uniform(np.log(max(lo, 1e-3)), np.log(hi)))
            mjd[j] = a.mjd[i] + lag / en.DAY
            used.add(int(i))
            moved += 1
            break
    return en.Sample(b.cat, mjd, b.ra, b.dec, b.sigma, b.year)


#: ASSUMPTION: below this null sd (cells of about one pair) the Gaussian tail is not valid; D is
#: then modelled as Skellam(lam, lam), the difference of two equal Poisson counts, with
#: lam = sd^2 / 2 from the null ensemble.
GAUSS_SD_MIN = 1.0


def two_sided_p(d, null) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided p and z per cell (model_prediction): Gaussian tail from the null mean and sd, or
    a Skellam tail where the null sd is below GAUSS_SD_MIN."""
    mu, sd = null.mean(axis=0), null.std(axis=0)
    d = np.asarray(d, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 0, (d - mu) / sd, 0.0)
    x = np.abs(d - np.round(mu))
    lam = np.maximum(sd**2 / 2, 1e-9)
    sk = np.where(x > 0, np.minimum(1.0, 2 * stats.skellam.sf(x - 1, lam, lam)), 1.0)
    return np.where(sd >= GAUSS_SD_MIN, 2 * stats.norm.sf(np.abs(z)), sk), z


_W: dict = {}


def _init(s):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W["s"] = s


def count_all(s) -> np.ndarray:
    return np.stack([signed_counts(s[a], s[b]) for a, b in CHANNELS])


def _null_chunk(seeds):
    s = _W["s"]
    out = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        out.append(count_all({k: en.scramble_jit(v, rng, P) for k, v in s.items()}))
    return np.stack(out)


def run_null(s, n, cpu, base_seed=8_000_000):
    seeds = np.arange(n) + base_seed
    chunks = [seeds[i : i + 50] for i in range(0, n, 50)]
    if cpu <= 1:
        _init(s)
        return np.concatenate([_null_chunk(c) for c in chunks])
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(s,)) as ex:
        return np.concatenate(list(ex.map(_null_chunk, chunks)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    t0 = time.time()
    ev = E.load(E.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items()}
    obs = count_all(s)
    null = run_null(s, a.n, a.cpu)
    mask = cell_mask()
    mu = null.mean(axis=0)
    pmin, pglob = en.global_p(np.abs(obs - mu), np.abs(null - mu), mask)
    pa, z = two_sided_p(obs, null)
    ncell = int(mask.sum())
    labels = en.lag_labels(P)
    rng = np.random.default_rng(8_500_000)
    cells = []
    for c, (ca, cb) in enumerate(CHANNELS):
        for k in range(NLAG):
            for q in np.flatnonzero(mask[c, k]):
                n50 = None
                for ninj in INJ_N:
                    det = 0
                    for _ in range(INJ_TRIALS):
                        sb = inject_after(
                            s[ca], s[cb], ninj, P.lag_edges[k], P.lag_edges[k + 1], q, rng
                        )
                        d = signed_counts(s[ca], sb)[k, q]
                        p1, _ = two_sided_p(np.array([d]), null[:, c, k, q][:, None])
                        det += p1[0] * ncell <= E.DETECT_P
                    if det >= INJ_TRIALS / 2:
                        n50 = ninj
                        break
                cells.append(
                    {
                        "channel": f"{ca}-{cb}",
                        "lag": labels[k],
                        "class": CLASSES[q],
                        "D_obs": int(obs[c, k, q]),
                        "null_mean": round(float(mu[c, k, q]), 2),
                        "null_sd": round(float(null[:, c, k, q].std()), 2),
                        "z": round(float(z[c, k, q]), 2),
                        "p_two_sided": float(pa[c, k, q]),
                        "n50_injected": n50,
                    }
                )
    out = {
        "test": "E1 signed-lag asymmetry D = N(B after A) - N(B before A) (jit null)",
        "n_scrambles": a.n,
        "n_cells": ncell,
        "min_cell_p_pooled": pmin,
        "global_p_pooled": pglob,
        "min_analytic_p_bonferroni": min(1.0, float(pa[mask].min()) * ncell),
        "cells": cells,
        "runtime_s": round(time.time() - t0),
    }
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    for r in sorted(cells, key=lambda r: -abs(r["z"]))[:8]:
        print(r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
