"""Post hoc test of the BAT-GBM long-lag deficit (``results/e1_swift/selfpair/PLAN.md``):
is it the jit null displacing the joint bursts (one GRB seen by both BAT and GBM)?

Variants (the run's jit null on both catalogues, all-separation counts):
``all`` (as the run), ``rm_shared`` (BAT bursts within SHARED_S of a GBM burst dropped) and
``rm_random`` (as many random non-shared BAT bursts dropped). The frozen rule ``outcome`` decides.

  python scripts/e1_swift_selfpair.py [--n 1000] [--cpu 4]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_icetracks_selfpair as sp  # noqa: E402
import e1_swift as sw  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402

OUT = ROOT / "results" / "e1_swift" / "selfpair" / "selfpair.json"
#: ASSUMPTIONs of the frozen plan.
SHARED_S = 100.0  # a BAT burst this close in time to a GBM burst is the same GRB (control window)
EXPLAINED_ABS_Z = 2.0
KEEP_Z = -3.0
CHUNK = 250  # scrambles per pool task
LAGS = ("100s-1h", "1h-1d", "1d-7d")  # lag cells tested by the rule (all-separation class)
ALL = 2  # count_channel class column: every pair
NLAG = len(sw.P.lag_edges) - 1


def outcome(z_rm: dict[str, float], z_ctl_1h1d: float) -> str:
    """Frozen rule: all-separation z of rm_shared per lag in LAGS, and rm_random's 1 h-1 d z."""
    if z_rm["1h-1d"] <= KEEP_Z or z_rm["1d-7d"] <= KEEP_Z:
        return "not explained"
    if all(abs(z_rm[k]) < EXPLAINED_ABS_Z for k in LAGS) and z_ctl_1h1d <= KEEP_Z:
        return "explained"
    return "undetermined"


def _init() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"


def count(b: en.Sample, g: en.Sample) -> np.ndarray:
    return en.count_channel(b, g, sw.P, same=False)[:NLAG, ALL]


def null_chunk(args) -> np.ndarray:
    """All-separation counts per lag bin for ``n`` jit scrambles from ``seed``."""
    b, g, n, seed = args
    rng = np.random.default_rng(seed)
    return np.array(
        [count(en.scramble_jit(b, rng, sw.P), en.scramble_jit(g, rng, sw.P)) for _ in range(n)],
        float,
    )


def summarize(name: str, b: en.Sample, g: en.Sample, null: np.ndarray) -> dict:
    obs = count(b, g)
    return {
        "name": name,
        "n_bat": len(b.mjd),
        "obs_all": obs.tolist(),
        "null_mean": np.round(null.mean(0), 2).tolist(),
        "null_sd": np.round(null.std(0), 2).tolist(),
        "z": np.round(en.z_score(obs, null), 2).tolist(),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    ev = sw.load(sw.fetch(False))
    b = en.Sample.from_table(ev["BAT"])
    g = en.Sample.from_table(ev["GBM"])
    sh = sp.shared_mask(b.mjd, g.mjd, SHARED_S)
    rng = np.random.default_rng(5)
    drop = rng.choice(np.flatnonzero(~sh), size=int(sh.sum()), replace=False)
    ctl = np.ones(len(b.mjd), bool)
    ctl[drop] = False
    variants = {"all": b, "rm_shared": sp.subset(b, ~sh), "rm_random": sp.subset(b, ctl)}
    # same seeds for every variant, chunks of CHUNK scrambles spread over the cores
    seeds = list(range(8_000_000, 8_000_000 + a.n, CHUNK))
    jobs = [
        (v, g, min(CHUNK, 8_000_000 + a.n - s0), s0) for v in variants.values() for s0 in seeds
    ]
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init) as ex:
        chunks = list(ex.map(null_chunk, jobs))
    nchunk = len(seeds)
    res = {
        name: summarize(name, v, g, np.concatenate(chunks[i * nchunk : (i + 1) * nchunk]))
        for i, (name, v) in enumerate(variants.items())
    }
    labels = en.lag_labels(sw.P)
    zl = {k: dict(zip(labels, res[k]["z"], strict=True)) for k in res}
    # pair drops from removing the shared bursts: data, null, and null minus data
    dropped = {}
    for i, lab in enumerate(labels):
        do = res["all"]["obs_all"][i] - res["rm_shared"]["obs_all"][i]
        dn = res["all"]["null_mean"][i] - res["rm_shared"]["null_mean"][i]
        dropped[lab] = {"obs": do, "null": round(dn, 2), "null_minus_obs": round(dn - do, 2)}
    out = {
        "provenance": "derived (post hoc verification; plan results/e1_swift/selfpair/PLAN.md)",
        "n_scrambles": a.n,
        "labels": labels,
        "n_shared_bat": int(sh.sum()),
        "n_gbm": len(g.mjd),
        "variants": res,
        "all_sep_z": zl,
        "rm_shared_drop": dropped,
        "outcome": outcome(zl["rm_shared"], zl["rm_random"]["1h-1d"]),
    }
    print(json.dumps({k: out[k] for k in ("n_shared_bat", "all_sep_z", "outcome")}, indent=1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
