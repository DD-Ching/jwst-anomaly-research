"""Corrected headline from run.py (PLAN.md). Usage (repo root): summarise.py
Writes results/dimming/wht_coverage/{summary.json, limits_combined.ecsv}."""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

FIELDS = ("macs0416", "abell2744", "nexus")
HERE = Path("results/dimming/wht_coverage")
COMMITTED = "results/dimming/macs0416_regen/limits_combined.ecsv"

runs = {f: json.loads((HERE / "runs" / f / "summary.json").read_text()) for f in FIELDS}
gates = all(r["gate_R0"] and r["gate_R1"] for r in runs.values())
out = {"gates_passed": gates, "fields": runs}
if gates:
    on = {f: Table.read(HERE / "runs" / f / "limits_veto_on_on.ecsv") for f in FIELDS}
    off = {f: Table.read(HERE / "runs" / f / "limits_veto_off_on.ecsv") for f in FIELDS}
    head = ds.combine_limits(on)
    head.meta["source"] += " (weight-map coverage term on; results/dimming/wht_coverage)"
    head.write(HERE / "limits_combined.ecsv", overwrite=True)
    old = Table.read(COMMITTED)
    swap = ds.combine_limits(off)
    rows = []
    for r, o, s in zip(head, old, swap, strict=True):
        assert (r["t_e_yr"], r["rho"]) == (o["t_e_yr"], o["rho"]) == (s["t_e_yr"], s["rho"])
        rows.append(
            {
                "t_e_yr": float(r["t_e_yr"]),
                "rho": float(r["rho"]),
                "limit_per_source_per_yr": float(r["limit_per_source_per_yr"]),
                "committed_limit_per_source_per_yr": float(o["limit_per_source_per_yr"]),
                "limit_ratio_corrected_over_committed": float(
                    r["limit_per_source_per_yr"] / o["limit_per_source_per_yr"]
                ),
                "area_limit_ratio_corrected_over_committed": float(
                    r["limit_per_deg2_per_yr"] / o["limit_per_deg2_per_yr"]
                ),
                "swap_gain_veto_off_over_on": float(
                    s["exposure_source_yr"] / r["exposure_source_yr"]
                ),
            }
        )
    out["headline"] = rows
    for k in (
        "limit_ratio_corrected_over_committed",
        "area_limit_ratio_corrected_over_committed",
        "swap_gain_veto_off_over_on",
    ):
        v = np.array([x[k] for x in rows])
        out[f"{k}_range"] = [float(v.min()), float(v.max())]
out["provenance"] = (
    "derived (exposures from simulated injections; weight maps and catalogues observed)"
)
(HERE / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: v for k, v in out.items() if k not in ("fields",)}, indent=1))
