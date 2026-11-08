"""Gate U for the D-042 second copy rule, and its diagnostic on stored forced tables
(plan: results/lens_consistency/copy_rule_unique/PLAN.md)."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import lens_consistency as lc
import numpy as np
from astropy.table import Table

# D-042 hand reading on MACS0416 (docs/fields/macs0416.md): copy rows 1, 15, 122 x2 (confused);
# 45 (mu -9) is not a copy (opposite parity)
HAND_U = {"confused": ["1", "15", "122", "122"], "absent": []}
NOT_COPY_U = {"45"}


def image_mu(model, images: Table, table: Table, ids: list[str]) -> dict[str, float]:
    """Signed model magnification at each catalogued image, at its system redshift in ``table``."""
    zsys = {str(r["system"]): float(r["z_sys"]) for r in table}
    out = {}
    for img in images:
        iid, s = str(img["image_id"]), str(img["system"])
        if iid in ids and s in zsys and np.isfinite(zsys[s]):
            m = model.magnification([float(img["ra"])], [float(img["dec"])], zsys[s])
            out[iid] = float(m["magnification"][0])
    return out


def label(table: Table, copies: list[str]) -> dict[str, Counter]:
    return {
        cls: Counter(
            str(r["system"])
            for r, c in zip(table, copies, strict=True)
            if c and r["forced_class"] == cls
        )
        for cls in lc.COPY_CLASSES
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True)
    ap.add_argument(
        "--table", type=Path, required=True, help="images_predicted.ecsv with forced_class"
    )
    ap.add_argument(
        "--unpredicted", type=Path, help="images summary JSON (unpredicted_catalogued_images)"
    )
    ap.add_argument("--unpredicted-ids", help="comma-separated ids when no summary JSON is stored")
    ap.add_argument("--rms", type=float, required=True)
    ap.add_argument(
        "--gate-u", action="store_true", help="compare with the D-042 MACS0416 hand list"
    )
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    table = Table.read(a.table)
    if hasattr(table["forced_class"], "filled"):
        table["forced_class"] = [str(x) for x in table["forced_class"].filled("")]
    if a.unpredicted:
        unpred = json.loads(a.unpredicted.read_text())["unpredicted_catalogued_images"]
    else:
        unpred = [s for s in (a.unpredicted_ids or "").split(",") if s]
    model, files, par = lc.load_model(a.model)
    images, _ = lc.image_list(a.model, files, par)
    lc.apply_frame_offset(a.model, model, images)  # same frame as the `images` run
    mu = image_mu(model, images, table, unpred)
    copies = lc.model_copies_unique(table, images, unpred, a.rms, mu)
    base = lc.model_copies(table, images, unpred, a.rms)
    pairs = lc.copy_partners(table, images, unpred, lc.COPY_UNIQUE_RMS_FACTOR * a.rms)
    rows = [
        {
            "system": str(r["system"]),
            "forced_class": str(r["forced_class"]),
            "magnification": round(float(r["magnification"]), 2),
            "copy_of_unique": c,
            "copy_of_d041": b,
            "partners_arcsec": {k: round(v, 2) for k, v in pairs[i].items()},
            "ra": round(float(r["ra"]), 6),
            "dec": round(float(r["dec"]), 6),
        }
        for i, (r, c, b) in enumerate(zip(table, copies, base, strict=True))
        if r["forced_class"] in lc.COPY_CLASSES
    ]
    got = label(table, copies)
    summary = {
        "provenance": "derived",
        "model": a.model,
        "table": str(a.table),
        "rms_arcsec": a.rms,
        "factor": lc.COPY_UNIQUE_RMS_FACTOR,
        "unpredicted_catalogued_images": unpred,
        "image_mu_model_prediction": {k: round(v, 3) for k, v in mu.items()},
        "forced_classes": dict(Counter(str(x) for x in table["forced_class"] if str(x))),
        "n_absent_confused": len(rows),
        "n_copy_unique": sum(bool(r["copy_of_unique"]) for r in rows),
        "labelled": {k: dict(v) for k, v in got.items()},
        "rows": rows,
    }
    if a.gate_u:
        res = {}
        for cls, hand in HAND_U.items():
            want = Counter(hand)
            res[cls] = {
                "rule_only": sorted((got[cls] - want).elements()),
                "hand_only": sorted((want - got[cls]).elements()),
            }
        labelled_not = sorted(NOT_COPY_U & {s for v in got.values() for s in v})
        res["not_copy_labelled"] = labelled_not
        res["passes"] = not labelled_not and all(
            not v["rule_only"] and not v["hand_only"] for k, v in res.items() if k in HAND_U
        )
        summary["gate_U"] = res
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
