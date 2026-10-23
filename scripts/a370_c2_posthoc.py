"""Post hoc diagnosis of the Abell 370 gate C2 failure (CANUCS DR1 reference route).

Plan: ``results/lens_consistency/abell370_buffalo_extra/posthoc_c2/PLAN.md``. Reads only the stored
control tables of the five CANUCS-reference fields; nothing is re-measured and no frozen outcome
changes. Post hoc: the shared-pair counting is a rule proposal, never applied to read a flag.

    python scripts/a370_c2_posthoc.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import m0417_arcref as ar
from astropy.table import Table

LC = Path("results/lens_consistency")
OUT = LC / "abell370_buffalo_extra/posthoc_c2"
# field -> (stored control table of the frozen gate C2, class column that gate read)
FIELDS = {
    "abell370": (LC / "abell370_buffalo_extra/control_abell370.ecsv", "class_local"),
    "macs0416": (LC / "macs0416_canucsref/control_macs0416.ecsv", "class_local"),
    "macs0417": (LC / "macs0417_canucsref/control_macs0417.ecsv", "class"),
    "macs1149": (LC / "macs1149_canucsref/control_macs1149.ecsv", "class_local"),
    "macs1423": (LC / "macs1423_canucsref/control_macs1423.ecsv", "class_local"),
}


# Frozen since 2026-10-23 (D-034 addendum i); kept under this name for the post hoc tests.
merge_shared_pairs = ar.merge_shared_pairs


def _finite(x) -> float | None:
    """JSON has no NaN: a non-finite value is written as null."""
    x = float(x)
    return x if math.isfinite(x) else None


def field_report(path: Path, col: str) -> dict:
    t = Table.read(path, format="ascii.ecsv")
    systems = [str(s) for s in t["system"]]
    classes = [str(c) for c in t[col]]
    arc = [bool(a) for a in t["arc_ref"]]
    pairs = [(str(r), str(c)) for r, c in zip(t["ref_counterpart"], t["counterpart"], strict=True)]
    as_stored = ar.gate_c2(classes, arc)
    shared = ar.gate_c2_shared(classes, arc, systems, pairs)
    merged = {k: shared[k] for k in ("all", "arc_referenced", "passes")}
    ia = [i for i, a in enumerate(arc) if a]
    absent = [
        {
            "system": systems[i],
            "images": str(t["images"][i]),
            "ref_counterpart": pairs[i][0],
            "counterpart": pairs[i][1],
            "ref_a3_arcsec": float(t["ref_a3_arcsec"][i]),
            "own_a3_arcsec": float(t["own_a3_arcsec"][i]),
            "best_ratio": _finite(t["best_ratio"][i]),
        }
        for i in ia
        if classes[i] == "absent"
    ]
    return {
        "table": str(path),
        "class_column": col,
        "gate_C2_as_stored": as_stored,
        "gate_C2_shared_pairs_once": merged,
        "merged_groups": [
            [f"{systems[j]}:{t['images'][j]}" for j in g] for g in shared["merged_groups"]
        ],
        "absent_arc_controls": absent,
    }


def main() -> None:
    out = {"label": "post hoc (not a rule; no flag read)", "fields": {}}
    for f, (p, col) in FIELDS.items():
        out["fields"][f] = field_report(p, col)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "posthoc.json").write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    for f, r in out["fields"].items():
        a, b = r["gate_C2_as_stored"], r["gate_C2_shared_pairs_once"]
        print(
            f"{f}: arc {a['arc_referenced']['n_ok']}/{a['arc_referenced']['n_testable']} "
            f"-> {b['arc_referenced']['n_ok']}/{b['arc_referenced']['n_testable']}, "
            f"pass {a['passes']} -> {b['passes']}, groups {r['merged_groups']}"
        )


if __name__ == "__main__":
    main()
