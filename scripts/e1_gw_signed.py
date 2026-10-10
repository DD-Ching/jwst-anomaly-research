"""E1: signed-lag ("which event comes first") GW cells with sky-map classes (D-074 addenda 3-4).

For GW x B (B = GBM, ICECAT, CHIME), lag bin and sky-map class (same / wide / antipodal, as in
``e1_gw_directional``), D = N(t_B > t_GW) - N(t_B < t_GW). Ordinary pairs and the ``jit`` null are
symmetric in sign. GW-GW is left out (D is antisymmetric within one catalogue). The known pair
GW170817 x GRB 170817A is left out of the family, as in addendum 4. Two-sided statistics and the
Skellam tail as in ``e1_signed_lag``; injections put B after the GW event in the tested class.

Writes ``results/e1_events/gw_signed.json``.

  python scripts/e1_gw_signed.py [--n 1000] [--cpu 4]
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
import e1_gw_directional as D  # noqa: E402
import e1_signed_lag as SL  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402

OUT = ROOT / "results" / "e1_events" / "gw_signed.json"
PARTNERS = ("GBM", "ICECAT", "CHIME")
P = D.P
NLAG = D.NLAG
#: ASSUMPTION: injected counts tried per cell, and trials per count.
INJ_N = (3, 10, 30)
INJ_TRIALS = 10


def signed_channel(m, gw, orig, b: en.Sample, skip=(), p: en.Params = P) -> np.ndarray:
    """D per (lag bin, class) for GW x b; only GW events with a map count."""
    e = np.asarray(p.lag_edges)
    i, j = en.pairs_within(gw.mjd, b.mjd, e[-1], False)
    keep = m.has[i]
    for gi, bj in skip:
        keep &= ~((i == gi) & (j == bj))
    i, j = i[keep], j[keep]
    dt = (b.mjd[j] - gw.mjd[i]) * en.DAY
    k = np.clip(np.searchsorted(e, np.abs(dt), side="right") - 1, 0, NLAG - 1)
    c = D.classify_x(m, i, D.rotation(gw, orig), b.ra[j], b.dec[j], b.sigma[j], p)
    sgn = np.sign(dt)
    out = np.zeros((NLAG, 3), np.int64)
    for q in range(3):
        sel = (c >> q) & 1 == 1
        out[:, q] = np.bincount(k[sel], weights=sgn[sel], minlength=NLAG).astype(np.int64)
    return out


def count_all(m, s, orig, skip=None) -> np.ndarray:
    skip = skip or {}
    return np.stack([signed_channel(m, s["GW"], orig, s[b], skip.get(b, ())) for b in PARTNERS])


_W: dict = {}


def _init(m, s, orig, skip):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W.update(m=m, s=s, orig=orig, skip=skip)


def _null_chunk(seeds):
    m, s, orig, skip = _W["m"], _W["s"], _W["orig"], _W["skip"]
    return np.stack(
        [
            count_all(
                m,
                {k: en.scramble_jit(v, np.random.default_rng(sd), P) for k, v in s.items()},
                orig,
                skip,
            )
            for sd in seeds
        ]
    )


def run_null(m, s, orig, n, cpu, skip, base_seed=9_000_000):
    seeds = np.arange(n) + base_seed
    chunks = [seeds[i : i + 25] for i in range(0, n, 25)]
    if cpu <= 1:
        _init(m, s, orig, skip)
        return np.concatenate([_null_chunk(c) for c in chunks])
    with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(m, s, orig, skip)) as ex:
        return np.concatenate(list(ex.map(_null_chunk, chunks)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    t0 = time.time()
    E.MANIFEST = D.EVENTS_MANIFEST
    ev = E.load(E.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items()}
    names = [str(x) for x in ev["GW"]["name"]]
    m = D.GWMaps(names, mjd=s["GW"].mjd, manifest=D.SKYMAP_MANIFEST)
    orig = s["GW"].mjd.copy()
    gi = names.index(E.PC_GW)
    bj = [str(x) for x in ev["GBM"]["name"]].index(E.PC_GBM)
    skip = {"GBM": [(gi, bj)]}
    obs = count_all(m, s, orig, skip)
    pc = signed_channel(m, s["GW"], orig, s["GBM"])[0] - obs[0, 0]
    null = run_null(m, s, orig, a.n, a.cpu, skip)
    mu = null.mean(axis=0)
    mask = np.ones(obs.shape, bool)
    pmin, pglob = en.global_p(np.abs(obs - mu), np.abs(null - mu), mask)
    pa, z = SL.two_sided_p(obs, null)
    ncell = int(mask.sum())
    sd = null.std(axis=0)
    labels = en.lag_labels(P)
    rng = np.random.default_rng(9_500_000)
    cells = []
    for c, b in enumerate(PARTNERS):
        for k in range(NLAG):
            for q, cls in enumerate(D.CLASSES):
                n50 = None
                for ninj in INJ_N:
                    det = 0
                    for _ in range(INJ_TRIALS):
                        sb = D.inject(
                            m,
                            s["GW"],
                            s[b],
                            ninj,
                            P.lag_edges[k],
                            P.lag_edges[k + 1],
                            cls,
                            rng,
                            sign=1,
                        )
                        o = signed_channel(m, s["GW"], orig, sb, skip.get(b, ()))[k, q]
                        p1 = SL.two_sided_p(np.array([o]), null[:, c, k, q][:, None])[0][0]
                        det += p1 * ncell <= E.DETECT_P
                    if det >= INJ_TRIALS / 2:
                        n50 = ninj
                        break
                cells.append(
                    {
                        "channel": f"GW-{b}",
                        "lag": labels[k],
                        "class": cls,
                        "D": int(obs[c, k, q]),
                        "null_mean": round(float(mu[c, k, q]), 2),
                        "null_sd": round(float(sd[c, k, q]), 2),
                        "z": round(float(z[c, k, q]), 2),
                        "analytic_p": float(pa[c, k, q]),
                        "n50_injected": n50,
                    }
                )
    out = {
        "test": "E1 signed-lag GW cells with sky-map classes (jit null)",
        "n_gw_with_map": int(m.has.sum()),
        "positive_control_GW170817_GRB170817A_D_by_lag_class": pc.tolist(),
        "n_scrambles": a.n,
        "n_cells": ncell,
        "min_cell_p_pooled": pmin,
        "global_p_pooled": pglob,
        "min_analytic_p_bonferroni": min(1.0, float(pa.min()) * ncell),
        "cells": cells,
        "runtime_s": round(time.time() - t0),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    for r in sorted(cells, key=lambda r: -abs(r["z"]))[:6]:
        print(r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
