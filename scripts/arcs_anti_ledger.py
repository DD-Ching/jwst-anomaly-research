"""`arcs` anti-row ledger: the single-component anti class after vetting (post hoc).

Plan: results/arcs/anti-ledger/PLAN.md. Re-reads the stored pooled census table only (no network)
and drops the anti rows that later addenda explained: WHL0137 2194 (`model_sensitive`, D-024
addendum 13) and RXC J0600 581 (`member (model input)` under GLAFIC v3, addendum 19).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_vet_noise as vn  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
CENSUS = REPO / "results/arcs/pooled-census/census.ecsv"
OUT = REPO / "results/arcs/anti-ledger"

# Status of each single-component anti row, from the addenda (derived).
STATUS = {
    "WHL0137 2194": "model_sensitive (D-024 addendum 13)",
    "RXC J0600 581": "member (model input) under GLAFIC v3; arc under Lenstool (addendum 19)",
    "MACS0308 2157": "open: robust to model error (14), MUSE undetermined (16, 17)",
}
REPRODUCE = {"n": 54, "anti": 3, "expected": 0.631, "p_ge_k": 0.0242}
P_MIN = 0.05  # ASSUMPTION threshold, frozen in PLAN.md


def class_rows(t: Table) -> Table:
    """Gate-passed rows that are not two-knot (addendum 9's post hoc `other_gate_passed` class)."""
    sel = np.asarray(t["gate_passed"], bool) & (np.asarray(t["deblend_class"]) != "two-knot")
    return t[sel]


def row_id(t: Table) -> list[str]:
    return [f"{f} {lab}" for f, lab in zip(t["field"], t["label"], strict=True)]


def variant(t: Table, drop_rows=(), drop_probs=True) -> dict:
    """Anti count and Poisson-binomial tail after removing ``drop_rows``.

    ``drop_probs=False`` removes the rows from the count only and keeps every probability.
    """
    ids = row_id(t)
    drop = set(drop_rows)
    anti = np.asarray(t["screen_class"]) == "anti"
    unknown = drop - set(ids)
    if unknown:
        raise ValueError(f"rows not in the class: {sorted(unknown)}")
    gone = np.array([i in drop for i in ids], bool)
    k = int(np.sum(anti & ~gone))
    probs = np.asarray(t["p_anti_row"], float)[~gone if drop_probs else slice(None)]
    return {
        "n": int(probs.size),
        "anti": k,
        "expected": round(float(probs.sum()), 3),
        "p_ge_k": vn.excess_p(probs, k),
    }


def decide(p_primary: float) -> str:
    if not np.isfinite(p_primary):
        return "undetermined"
    if p_primary >= P_MIN:
        return "single-component excess not supported after vetting (post hoc)"
    return "single-component excess remains (post hoc)"


def run(t: Table) -> dict:
    c = class_rows(t)
    rep = variant(c)
    ok = (
        rep["n"] == REPRODUCE["n"]
        and rep["anti"] == REPRODUCE["anti"]
        and rep["expected"] == REPRODUCE["expected"]
        and round(rep["p_ge_k"], 4) == REPRODUCE["p_ge_k"]
    )
    if not ok:
        raise SystemExit(f"stored post hoc split not reproduced: {rep}")
    anti_ids = [i for i, s in zip(row_id(c), c["screen_class"], strict=True) if s == "anti"]
    if sorted(anti_ids) != sorted(STATUS):
        raise SystemExit(f"anti rows differ from the ledger: {anti_ids}")
    vetted = ("WHL0137 2194", "RXC J0600 581")
    out = {
        "label": "post hoc, exploratory; derived from ../pooled-census/census.ecsv",
        "reproduce": rep,
        "primary": variant(c, vetted),
        "drop_2194_only": variant(c, ("WHL0137 2194",)),
        "obs_only": variant(c, vetted, drop_probs=False),
        "rows": {
            i: {"status": STATUS[i], "p_anti_row": round(float(p), 5)}
            for i, p, s in zip(row_id(c), c["p_anti_row"], c["screen_class"], strict=True)
            if s == "anti"
        },
    }
    out["outcome"] = decide(out["primary"]["p_ge_k"])
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--census", type=Path, default=CENSUS)
    p.add_argument("--out", type=Path, default=OUT)
    a = p.parse_args(argv)
    res = run(Table.read(a.census))
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
