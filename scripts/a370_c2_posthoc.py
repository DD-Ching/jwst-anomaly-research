"""Post hoc diagnosis of the Abell 370 gate C2 failure (CANUCS DR1 reference route).

Plan: ``results/lens_consistency/abell370_buffalo_extra/posthoc_c2/PLAN.md``. Reads only the stored
control tables of the five CANUCS-reference fields; nothing is re-measured and no frozen outcome
changes. Post hoc: the shared-pair counting is a rule proposal, never applied to read a flag.

    python scripts/a370_c2_posthoc.py
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
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


def merge_shared_pairs(systems, pairs, classes) -> tuple[list[str], list[list[int]]]:
    """Count controls of different systems that share one (reference, counterpart) source pair once.

    Returns the merged class list and the row groups merged (only groups spanning >= 2 systems).
    A merged control is ``absent`` if any member is absent, else the first member's class.
    """
    groups: dict[tuple, list[int]] = defaultdict(list)
    for i, p in enumerate(pairs):
        groups[p].append(i)
    merged, multi, seen = [], [], set()
    for i, p in enumerate(pairs):
        if i in seen:
            continue
        g = groups[p]
        if len({systems[j] for j in g}) < 2:
            merged.append(classes[i])
            seen.add(i)
            continue
        seen.update(g)
        multi.append(g)
        cls = [classes[j] for j in g]
        merged.append("absent" if "absent" in cls else cls[0])
    return merged, multi


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
    # merge only among arc-referenced controls; others keep their rows
    ia = [i for i, a in enumerate(arc) if a]
    m_cls, m_groups = merge_shared_pairs(
        [systems[i] for i in ia], [pairs[i] for i in ia], [classes[i] for i in ia]
    )
    rest = [classes[i] for i, a in enumerate(arc) if not a]
    merged = ar.gate_c2(rest + m_cls, [False] * len(rest) + [True] * len(m_cls))
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
        "merged_groups": [[f"{systems[ia[j]]}:{t['images'][ia[j]]}" for j in g] for g in m_groups],
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
