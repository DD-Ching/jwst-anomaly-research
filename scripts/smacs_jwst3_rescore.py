"""Sensitivity of SMACS galaxy-stratum ranks to the Gaia-position-only stars (no SIMBAD override).

Plan: ``results/smacs_jwst3/PLAN_galaxy_stratum.md``. Re-scores the run's galaxy stratum from its
saved ``features.ecsv`` twice, with the same procedure: as ranked (baseline) and with the
quality-passing ``gaia_position`` stars returned to it. The paired difference bounds what the
missing "a nearer SIMBAD galaxy wins" rule (D-012) can do to the record ranks.

The ECSV round trip sorts ``meta["feature_spec"]`` keys alphabetically, and ``rank.feature_names``
follows that order, so the isolation forest and LOF see permuted features (2822 of 2938 ranks
differed). The spec is therefore restored to column order (the in-memory order), and the script
checks that the baseline reproduces the pipeline's ranks exactly before using it.

Usage: python scripts/smacs_jwst3_rescore.py <run_dir>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly.rank import score_anomalies

RECORDS = ("f200w_921", "f200w_2490", "f277w_743", "f356w_1098", "f200w_389", "f200w_2586")
PREFIX = "jw02736-o001_t001_nircam_"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    d = ap.parse_args().run_dir / "smacs0723_nircam"
    feats = Table.read(d / "features.ecsv")
    spec = feats.meta["feature_spec"]
    feats.meta["feature_spec"] = {c: spec[c] for c in feats.colnames if c in spec}
    ok = {
        str(u): bool(v)
        for u, v in zip(
            *Table.read(d / "quality.ecsv")[["source_uid", "quality_ok"]].columns.values(),
            strict=True,
        )
    }
    pops = Table.read(d / "populations.ecsv")
    basis = {
        str(u): str(b)
        for u, b in zip(pops["source_uid"], np.ma.filled(pops["star_basis"], ""), strict=True)
    }
    ranked = {
        str(u): int(r)
        for u, r in zip(
            *Table.read(d / "scores.ecsv")[["source_uid", "rank"]].columns.values(), strict=True
        )
    }
    added = {u for u, b in basis.items() if b == "gaia_position" and ok[u]}
    uids = np.array([str(u) for u in feats["source_uid"]])
    base = np.isin(uids, list(ranked))
    r0, r1 = (
        {str(u): int(r) for u, r in zip(s["source_uid"], s["rank"], strict=True)}
        for s in (
            score_anomalies(feats[m], random_state=0)
            for m in (base, base | np.isin(uids, list(added)))
        )
    )
    diff = np.array([abs(r0[u] - ranked[u]) for u in ranked])
    if diff.any():
        raise SystemExit(
            f"baseline re-score differs from the pipeline: {int((diff > 0).sum())} ranks"
        )
    print(f"baseline re-score reproduces all {len(diff)} pipeline ranks")
    print(
        "returned gaia_position stars (rank with them):",
        sorted((r1[u], u.removeprefix(PREFIX)) for u in added),
    )
    print("record  pipeline  baseline  with_returned")
    for rec in RECORDS:
        u = PREFIX + rec
        if u not in ranked:
            print(f"{rec}  not in the galaxy stratum")
            continue
        print(f"{rec}  {ranked[u]}  {r0[u]}  {r1[u]}")


if __name__ == "__main__":
    main()
