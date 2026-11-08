"""Gate G for the D-041 model-copy rule on the stored MACS0717 table
(plan: results/lens_consistency/model_copy_rule/PLAN.md)."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import lens_consistency as lc
from astropy.table import Table

# hand classification, docs/fields/macs0717.md (systems of the copy rows, per forced class)
HAND = {
    "confused": ["1"] * 3
    + ["3"] * 2
    + ["7", "12", "17"]
    + ["24"] * 2
    + ["32"] * 2
    + ["34"]
    + ["45"] * 2
    + ["50", "56"]
    + ["62"] * 2
    + ["67"]
    + ["71"] * 2
    + ["72", "75", "76", "79"]
    + ["80"] * 2,
    "absent": ["5", "6", "34"],
}


def gate(table: Table, copies: list[str]) -> dict:
    """Per-class, per-system comparison of rule-flagged copy rows with the hand list."""
    out = {}
    for cls, hand in HAND.items():
        rule = Counter(
            str(r["system"])
            for r, c in zip(table, copies, strict=True)
            if c and r["forced_class"] == cls
        )
        want = Counter(hand)
        out[cls] = {
            "rule_n": sum(rule.values()),
            "hand_n": len(hand),
            "rule_only": sorted((rule - want).elements()),
            "hand_only": sorted((want - rule).elements()),
        }
    out["passes"] = all(not v["rule_only"] and not v["hand_only"] for v in out.values())
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stored", type=Path, default=Path("results/lens_consistency/macs0717-cats"))
    ap.add_argument("--rms", type=float, required=True, help="validate rms_dtheta_arcsec")
    ap.add_argument("--out", type=Path, default=Path("results/lens_consistency/model_copy_rule"))
    a = ap.parse_args()
    table = Table.read(a.stored / "images_predicted.ecsv")
    unpred = json.loads((a.stored / "images.json").read_text())["unpredicted_catalogued_images"]
    files = lc.model_files("macs0717-cats", ("arcs.dat",))
    images = lc.lensmodel.load_lenstool_images(files["arcs.dat"])
    copies = lc.model_copies(table, images, unpred, a.rms)
    rows = [
        {
            "system": str(r["system"]),
            "forced_class": str(r["forced_class"]),
            "copy_of": c,
            "ra": round(float(r["ra"]), 6),
            "dec": round(float(r["dec"]), 6),
        }
        for r, c in zip(table, copies, strict=True)
        if r["forced_class"] in lc.COPY_CLASSES
    ]
    summary = {
        "provenance": "derived",
        "model": "macs0717-cats",
        "rms_arcsec": a.rms,
        "factor": lc.COPY_RMS_FACTOR,
        "detection_condition_applied": False,
        "n_catalogued_images": len(images),
        "n_absent_confused": len(rows),
        "n_copy": sum(bool(r["copy_of"]) for r in rows),
        "gate_G": gate(table, copies),
        "rows": rows,
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
