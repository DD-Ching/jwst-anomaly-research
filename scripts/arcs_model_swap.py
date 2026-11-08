"""Pooled `arcs` anti rows under a second published model of the same field (model-error test).

Plan: results/arcs/anti-model-swap/PLAN.md (frozen before the second-run classes of the anti rows
were read). Stored tables only, no network: each field's deciding-run ``arcs.ecsv`` (as in
``arcs_pooled.FIELDS``) against the field's second run.

    python scripts/arcs_model_swap.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_pooled as ap  # noqa: E402

ROOT = ap.ROOT
OUT = ROOT / "anti-model-swap"
# field -> (deciding run, second run); MACS1423 and RXC J0018 have no second run (plan)
PAIRS = {
    "MACS0417": (ap.FIELDS["MACS0417"], "macs0417-relics/run1-relics-v1/arcs.ecsv"),
    "WHL0137": (ap.FIELDS["WHL0137"], "whl0137-relics/whl0137-relics-glafic/arcs.ecsv"),
    "RXC J0600": (ap.FIELDS["RXC J0600"], "rxc0600-relics/rxc0600-relics-lenstool/arcs.ecsv"),
    "MACS0308": (ap.FIELDS["MACS0308"], "macs0308-relics/macs0308-relics-zitrin/arcs.ecsv"),
}
UNTESTED = ("MACS1423", "RXC J0018")
MATCH_ARCSEC = 0.1  # ASSUMPTION (plan)
CONTROL_MIN = 0.70  # ASSUMPTION (plan; the `arcs` control threshold)
MIN_TESTED = 4  # plan rule 3


def second_class(row, second: Table, match_arcsec: float = MATCH_ARCSEC) -> str:
    """Class of a primary row in the second run: its orientation class, or a reason it has none."""
    hit = second[second["label"] == row["label"]]
    if len(hit) == 0:
        return "unmatched"
    h = hit[0]
    dra = (float(h["ra"]) - float(row["ra"])) * np.cos(np.radians(float(row["dec"])))
    sep = 3600.0 * np.hypot(dra, float(h["dec"]) - float(row["dec"]))
    if sep > match_arcsec:
        return "unmatched"
    if not (bool(h["strong_shear"]) and h["population"] == "background"):
        return "not_strong_background"
    return str(h["orientation_class"])


def control_gate(classes: list[str], minimum: float = CONTROL_MIN) -> dict:
    """Rule 1: retention of matched primary aligned rows as `aligned` in the second run."""
    matched = [c for c in classes if c != "unmatched"]
    r = sum(c == "aligned" for c in matched) / len(matched) if matched else float("nan")
    return {
        "n_matched": len(matched),
        "aligned": sum(c == "aligned" for c in matched),
        "anti": sum(c == "anti" for c in matched),
        "retention": r,
        "passes": bool(matched) and r >= minimum,
    }


def outcome(gate_passes: bool, m: int, k: int, min_tested: int = MIN_TESTED) -> str:
    """Rule 3 (pooled, combined over fields)."""
    if m < min_tested:
        return "undetermined"
    if not gate_passes:
        return "undetermined"
    if k >= m - 1:
        return "persists"
    if k <= m / 2:
        return "model-dependent"
    return "mixed"


def main() -> None:
    anti_rows, control = [], []
    second_stats = {}
    for field, (prim, sec) in PAIRS.items():
        p, pmeta = ap.load_run(prim)
        s_all = Table.read(ROOT / sec)
        s_bg, smeta = ap.load_run(sec)
        if pmeta["assumptions"] != smeta["assumptions"]:
            raise SystemExit(f"error: {field} runs use different assumptions")
        a = pmeta["assumptions"]
        if float(a["sigma_eps"]) != 0.26:
            raise SystemExit(f"error: {field} runs used σ_ε {a['sigma_eps']}, the plan uses 0.26")
        second_stats[field] = ap.field_stats(
            s_bg, float(a["sigma_eps"]), float(a["min_ellipticity"]), float(a["anti_deg"])
        )
        for row in p:
            c = second_class(row, s_all)
            if row["orientation_class"] == "anti":
                anti_rows.append(
                    {
                        "field": field,
                        "label": int(row["label"]),
                        "g_min_primary": float(row["g_min"]),
                        "offset_min_primary": float(row["offset_min"]),
                        "second_class": c,
                    }
                )
            elif row["orientation_class"] == "aligned":
                control.append({"field": field, "label": int(row["label"]), "second_class": c})
    gate = control_gate([c["second_class"] for c in control])
    tested = [r for r in anti_rows if r["second_class"] != "unmatched"]
    m, k = len(tested), sum(r["second_class"] == "anti" for r in tested)
    # secondary: variant B with the second runs replacing the deciding runs in the four fields
    sigma = 0.26
    prim_stats = {}
    for f in ap.VARIANTS["B"]:
        rows, meta = ap.load_run(ap.FIELDS[f])
        a = meta["assumptions"]
        if float(a["sigma_eps"]) != sigma:
            raise SystemExit(f"error: {f} deciding run used σ_ε {a['sigma_eps']}, not {sigma}")
        prim_stats[f] = ap.field_stats(
            rows, sigma, float(a["min_ellipticity"]), float(a["anti_deg"])
        )
    swapped = {f: second_stats.get(f, prim_stats[f]) for f in ap.VARIANTS["B"]}
    summary = {
        "plan": "results/arcs/anti-model-swap/PLAN.md",
        "label": "derived (classes from stored model_prediction tables); exploratory",
        "pairs": PAIRS,
        "untested_fields": list(UNTESTED),
        "assumptions": {
            "match_arcsec": MATCH_ARCSEC,
            "control_min": CONTROL_MIN,
            "min_tested": MIN_TESTED,
        },
        "control": gate,
        "control_by_field": {
            f: control_gate([c["second_class"] for c in control if c["field"] == f]) for f in PAIRS
        },
        "anti_rows": anti_rows,
        "m_tested": m,
        "k_anti_in_second": k,
        "outcome": outcome(gate["passes"], m, k),
        "secondary_variant_B": {
            "sigma_eps": sigma,
            "deciding": ap.pooled(prim_stats, ap.VARIANTS["B"]),
            "second_runs_swapped_in": ap.pooled(swapped, ap.VARIANTS["B"]),
            "per_field_deciding": prim_stats,
            "per_field_second": second_stats,
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    print(
        json.dumps(
            {k_: summary[k_] for k_ in ("control", "m_tested", "k_anti_in_second", "outcome")}
        )
    )


if __name__ == "__main__":
    main()
