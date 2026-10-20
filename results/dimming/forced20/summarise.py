"""Apply the frozen PLAN.md rule to `forced.ecsv` of one field's targets and write
<field>_summary.json and <field>_targets.ecsv (small).
Usage (repo root): summarise.py FIELD SCREEN_OUT
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

DET_EPOCH = {"macs0416": 2, "abell2744": 1}  # the targets' only F200W catalogue epoch


def explanation(row, ne: int, det_epoch: int, band: str = "F200W") -> str:
    """PLAN.md decision rule for one target row of forced.ecsv."""
    if bool(row["survives_artefact_tests"]):
        return "survivor"
    other = [k for k in range(ne) if k != det_epoch]
    finite = [np.isfinite(row[f"{band}_e{k}_fflux"]) for k in other]
    cut = bool(row["cut_edge"] or row["cut_low_weight"] or row["cut_spike"])
    if not bool(row["forced_vanish"]) and any(finite):
        return "catalogue_depth_or_deblend"
    if not any(finite) or cut:
        return "footprint_or_image_quality"
    return "other_forced_test"  # vanish reappears but a later artefact test removes it


def main(field: str, out: Path) -> dict:
    d = DET_EPOCH[field]
    t = Table.read(out / "forced.ecsv")
    ne = max(int(c.split("_e")[1].split("_")[0]) for c in t.colnames if c.startswith("F200W_e")) + 1
    ctl = np.asarray(t["control"], bool)
    depth = {
        f"e{k}": float(np.nanmedian(np.asarray(t[f"F200W_e{k}_ferr"])[ctl])) for k in range(ne)
    }
    rows = []
    for r in t[~ctl]:
        f2, e2 = float(r[f"F200W_e{d}_fflux"]), float(r[f"F200W_e{d}_ferr"])
        snr = {f"e{k}": float(r[f"F200W_e{k}_fflux"] / r[f"F200W_e{k}_ferr"]) for k in range(ne)}
        rows.append(
            {
                "uid": str(r["uid"]),
                "det_epoch_flux_ujy": f2,
                "det_epoch_snr": f2 / e2,
                "snr_by_epoch": snr,
                "expected_snr_other_median": float(
                    np.nanmedian([f2 / r[f"F200W_e{k}_ferr"] for k in range(ne) if k != d])
                ),
                "forced_vanish": bool(r["forced_vanish"]),
                "forced_confirmed": bool(r["forced_confirmed"]),
                "cut": [n for n in ("cut_edge", "cut_low_weight", "cut_spike") if bool(r[n])],
                "chromatic_veto": bool(r["chromatic_veto"]),
                "explanation": explanation(r, ne, d),
            }
        )
    fs = json.loads((out / "forced_summary.json").read_text())
    match = json.loads((out / "targets_match.json").read_text())
    n_ctl = int(ctl.sum())
    summary = {
        "calibrated": fs["calibrated"],
        "n_controls": n_ctl,
        "controls_flagged": fs["controls_flagged"],
        "run_counts": fs["calibrated"] and fs["controls_flagged"] <= 0.01 * n_ctl,
        "noise_scale": fs["noise_scale"],
        "median_control_F200W_ferr_ujy": depth,
        "det_epoch": d,
        "weight_edge_removed": match["weight_edge"],
        "explanations": {
            **{
                k: sum(r["explanation"] == k for r in rows)
                for k in sorted({r["explanation"] for r in rows})
            },
            "weight_edge": len(match["weight_edge"]),
        },
        "targets": rows,
    }
    return summary


if __name__ == "__main__":
    field, out = sys.argv[1], Path(sys.argv[2])
    s = main(field, out)
    here = Path(__file__).parent
    (here / f"{field}_summary.json").write_text(json.dumps(s, indent=1) + "\n")
    t = Table.read(out / "forced.ecsv")
    t[~np.asarray(t["control"], bool)].write(here / f"{field}_targets.ecsv", overwrite=True)
    print(json.dumps({k: v for k, v in s.items() if k != "targets"}, indent=1))
