"""Post hoc test of the DR2-ICECAT long-lag deficit (``results/e1_icetracks/selfpair/PLAN.md``):
is it the null displacing the DR2 tracks that are the same neutrinos as ICECAT alerts?

Variants (main null as the run: ICECAT by scramble_jit, DR2 by the uptime jitter):
``all`` (as the run), ``rm_shared`` (DR2 tracks within SHARED_S of an ICECAT alert dropped) and
``rm_random`` (as many random non-shared DR2 tracks dropped). The frozen rule ``outcome`` decides.

  python scripts/e1_icetracks_selfpair.py [--n 1000] [--cpu 4]
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

import e1_icetracks as it  # noqa: E402
import e1_swift as sw  # noqa: E402
import nf_ghost_dr2 as dr2  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly import neutrino_frontier as nf  # noqa: E402

OUT = ROOT / "results" / "e1_icetracks" / "selfpair" / "selfpair.json"
#: ASSUMPTIONs of the frozen plan.
SHARED_S = 10.0  # a DR2 track this close in time to an ICECAT alert is the alert's own event
EXPLAINED_ABS_Z = 2.0
KEEP_Z = -3.0
CHUNK = 250  # scrambles per pool task
LAGS = ("100s-1h", "1h-1d", "1d-7d")  # wide cells tested by the rule


def shared_mask(dr2_mjd: np.ndarray, ic_mjd: np.ndarray, tol_s: float = SHARED_S) -> np.ndarray:
    """True for DR2 times within ``tol_s`` of any ICECAT time."""
    ic = np.sort(ic_mjd)
    if len(ic) < 2:
        return (
            (np.abs(dr2_mjd - ic[0]) * en.DAY <= tol_s) if len(ic) else np.zeros(len(dr2_mjd), bool)
        )
    i = np.clip(np.searchsorted(ic, dr2_mjd), 1, len(ic) - 1)
    d = np.minimum(np.abs(dr2_mjd - ic[i - 1]), np.abs(dr2_mjd - ic[i])) * en.DAY
    return d <= tol_s


def outcome(z_rm: dict[str, float], z_ctl_1d7d: float) -> str:
    """Frozen rule: wide-cell z of rm_shared per lag in LAGS, and rm_random's 1-7 d wide z."""
    if z_rm["1h-1d"] <= KEEP_Z or z_rm["1d-7d"] <= KEEP_Z:
        return "not explained"
    if all(abs(z_rm[k]) < EXPLAINED_ABS_Z for k in LAGS) and z_ctl_1d7d <= KEEP_Z:
        return "explained"
    return "undetermined"


def subset(s: en.Sample, keep: np.ndarray) -> en.Sample:
    return en.Sample(s.cat, s.mjd[keep], s.ra[keep], s.dec[keep], s.sigma[keep], s.year[keep])


NLAG = len(it.P.lag_edges) - 1


def _init() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"


def null_chunk(args) -> np.ndarray:
    """Main-null counts (same, wide) per lag bin for ``n`` scrambles from ``seed``."""
    d, ic, up, n, seed = args
    rng = np.random.default_rng(seed)
    return np.array(
        [
            en.count_channel(
                nf.jitter_uptime(d, up[:, 0], up[:, 1], it.P.jitter_s / en.DAY, rng),
                en.scramble_jit(ic, rng, it.P),
                it.P,
                same=False,
            )[:NLAG, :2]
            for _ in range(n)
        ],
        float,
    )


def summarize(name: str, d: en.Sample, ic: en.Sample, null: np.ndarray) -> dict:
    obs = en.count_channel(d, ic, it.P, same=False)[:NLAG, :2]
    return {
        "name": name,
        "n_dr2": len(d.mjd),
        "obs_same_wide": obs.tolist(),
        "null_mean": np.round(null.mean(0), 2).tolist(),
        "null_sd": np.round(null.std(0), 2).tolist(),
        "z": np.round(en.z_score(obs, null), 2).tolist(),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    raw, up = dr2.fetch()
    ev = sw.load(sw.fetch(False))
    d = en.Sample.from_table(it.dr2_events(raw))
    ic = en.Sample.from_table(ev["ICECAT"])
    sh = shared_mask(d.mjd, ic.mjd)
    rng = np.random.default_rng(5)
    drop = rng.choice(np.flatnonzero(~sh), size=int(sh.sum()), replace=False)
    ctl = np.ones(len(d.mjd), bool)
    ctl[drop] = False
    variants = {"all": d, "rm_shared": subset(d, ~sh), "rm_random": subset(d, ctl)}
    # same seeds for every variant, chunks of CHUNK scrambles spread over the cores
    seeds = list(range(8_000_000, 8_000_000 + a.n, CHUNK))
    jobs = [
        (v, ic, up, min(CHUNK, 8_000_000 + a.n - s0), s0) for v in variants.values() for s0 in seeds
    ]
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init) as ex:
        chunks = list(ex.map(null_chunk, jobs))
    nchunk = len(seeds)
    res = {
        name: summarize(name, v, ic, np.concatenate(chunks[i * nchunk : (i + 1) * nchunk]))
        for i, (name, v) in enumerate(variants.items())
    }
    labels = en.lag_labels(it.P)
    zw = {k: {lab: res[k]["z"][i][1] for i, lab in enumerate(labels)} for k in res}
    # wide-pair drops from removing the shared tracks: data, null, and null minus data
    drop = {}
    for i, lab in enumerate(labels):
        do = res["all"]["obs_same_wide"][i][1] - res["rm_shared"]["obs_same_wide"][i][1]
        dn = res["all"]["null_mean"][i][1] - res["rm_shared"]["null_mean"][i][1]
        drop[lab] = {"obs": do, "null": round(dn, 2), "null_minus_obs": round(dn - do, 2)}
    out = {
        "provenance": "derived (post hoc verification; plan results/e1_icetracks/selfpair/PLAN.md)",
        "n_scrambles": a.n,
        "labels": labels,
        "n_shared_dr2": int(sh.sum()),
        "n_icecat": len(ic.mjd),
        "variants": res,
        "wide_z": zw,
        "rm_shared_wide_drop": drop,
        "outcome": outcome(zw["rm_shared"], zw["rm_random"]["1d-7d"]),
    }
    print(json.dumps({k: out[k] for k in ("n_shared_dr2", "wide_z", "outcome")}, indent=1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
