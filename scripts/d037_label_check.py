"""Check that ``images``' d037_label column (``lens_consistency.d037_columns``) reproduces the
D-037 gate-R repeated-pair labels on the stored MACS1149 CATS forced table
(results/lens_consistency/d037_rules/summary.json)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
from astropy.table import Table


def compare(table: Table, gate_rows: list[dict]) -> dict:
    """Per gate flag row (system, forced_class): the column's label against the gate's
    ``shared_pair``-derived label; ``passes`` when every one agrees."""
    got = {}
    for r in table:
        if str(r["forced_class"]) in lc.COPY_CLASSES and str(r["image_class"]) != "observed":
            got.setdefault((str(r["system"]), str(r["forced_class"])), []).append(
                [str(r["d037_label"]), str(r["d037_pair"])]
            )
    res = {}
    for g in gate_rows:
        key = (g["system"], g["forced_class"])
        want = ["model_topology" if g["shared_pair"] else "", g["shared_pair"]]
        res[f"{key[0]}:{key[1]}"] = {
            "want": want,
            "got": got.get(key),
            "ok": got.get(key) == [want],
        }
    extra = sorted(f"{k[0]}:{k[1]}" for k in got if f"{k[0]}:{k[1]}" not in res)
    return {
        "rows": res,
        "rows_not_in_gate": extra,
        "passes": all(v["ok"] for v in res.values()) and not extra,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", default="macs1149-cats")
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--imageplane", type=Path, required=True)
    ap.add_argument("--gate-summary", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    table = Table.read(a.table, format="ascii.ecsv")
    for c in ("forced_class", "image_class"):
        if hasattr(table[c], "filled"):
            table[c] = [str(x) for x in table[c].filled("")]
    model, files, par = lc.load_model(a.model)
    images, _ = lc.image_list(a.model, files, par)
    lc.apply_frame_offset(a.model, model, images)
    lc.d037_columns(table, images, lc.shared_match_ids(a.imageplane))
    gate_rows = json.loads(a.gate_summary.read_text())["flag_rows"]
    summary = {
        "provenance": "derived (forced classes: stored 2026-10-14 run; labels: d037_columns)",
        "model": a.model,
        "inputs": {k: str(v) for k, v in vars(a).items() if k != "out"},
        "assumptions": {"pair_radius_arcsec": lc.PAIR_RADIUS_ARCSEC},
        "check": compare(table, gate_rows),
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
