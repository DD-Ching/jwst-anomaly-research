"""PA convention of the deblend two-knot axis ``pa_sep`` against the catalogue axis ``pa_obs``.

Plan (frozen before the run): ``results/arcs/pa-convention/PLAN.md``. Reads the stored census
outputs of S1063 (#273) and MACS1423 (#282) with their ``arcs.ecsv`` tables; no download.
``pa_obs`` is ``observed``; ``pa_sep`` and the offsets are ``derived``.

    python scripts/arcs_pa_convention.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel

ROOT = Path(__file__).resolve().parents[1] / "results" / "arcs"
OUT = ROOT / "pa-convention"
FIELDS = {
    "abells1063-cats": (
        ROOT / "abells1063-cats" / "vetting" / "census.json",
        ROOT / "abells1063-cats" / "arcs.ecsv",
    ),
    "macs1423-canucs": (
        ROOT / "macs1423-canucs" / "vetting" / "census.json",
        ROOT / "macs1423-canucs" / "macs1423-canucs-maps" / "arcs.ecsv",
    ),
}
MIN_ROWS = 8  # PLAN.md: fewer primary rows → undetermined
MAX_MEDIAN = 20.0  # PLAN.md: identity median offset, degrees (ASSUMPTION)
MIN_MARGIN = 15.0  # PLAN.md: margin below both alternatives, degrees (ASSUMPTION)
MIN_ELL = 0.2  # PLAN.md: secondary set ellipticity cut (ASSUMPTION)
CONVENTIONS = {
    "identity": lambda pa: np.mod(pa, 180.0),
    "mirrored": lambda pa: np.mod(180.0 - pa, 180.0),
    "rotated90": lambda pa: np.mod(pa + 90.0, 180.0),
}


def medians(pa_sep, pa_obs) -> dict[str, float]:
    """Median axis offset (deg, 0–90) of ``pa_sep`` from each candidate convention of ``pa_obs``."""
    pa_sep = np.asarray(pa_sep, float)
    pa_obs = np.asarray(pa_obs, float)
    return {
        k: float(np.median(lensmodel.axis_offset_deg(pa_sep, f(pa_obs))))
        for k, f in CONVENTIONS.items()
    }


def outcome(n_rows: int, med: dict[str, float]) -> str:
    """Frozen PLAN.md rule."""
    if n_rows < MIN_ROWS:
        return "undetermined"
    for name, value in med.items():
        others = [v for k, v in med.items() if k != name]
        if value <= MAX_MEDIAN and all(v - value >= MIN_MARGIN for v in others):
            return "convention agrees" if name == "identity" else f"convention differs: {name}"
    return "undetermined"


def load_rows() -> Table:
    rows = []
    for field, (census_path, arcs_path) in FIELDS.items():
        census = json.loads(census_path.read_text())["rows"]
        arcs = Table.read(arcs_path)
        by_label = {str(int(lab)): i for i, lab in enumerate(arcs["label"])}
        for label, r in census.items():
            i = by_label[label]
            rows.append(
                {
                    "field": field,
                    "label": int(label),
                    "deblend_class": r["class"],
                    "gate_passed": bool(r["gate"]["passed"]),
                    "sep_arcsec": float(r["sep_arcsec"]),
                    "pa_sep": float(r["pa_sep"]),
                    "pa_obs": float(arcs["pa_obs"][i]),
                    "ellipticity": float(arcs["ellipticity"][i]),
                    "pa_pred_z2": float(arcs["pa_pred_z2"][i]),
                    "screen_class": str(arcs["orientation_class"][i]),
                }
            )
    t = Table(rows=rows)
    t["offset_identity"] = lensmodel.axis_offset_deg(t["pa_sep"], t["pa_obs"])
    # post hoc (not in PLAN.md): knot axis and catalogue axis against the z_s = 2 prediction
    t["offset_sep_pred_z2"] = lensmodel.axis_offset_deg(t["pa_sep"], t["pa_pred_z2"])
    t["offset_obs_pred_z2"] = lensmodel.axis_offset_deg(t["pa_obs"], t["pa_pred_z2"])
    return t


def summarise(t: Table, mask) -> dict:
    s = t[np.asarray(mask, bool)]
    med = medians(s["pa_sep"], s["pa_obs"]) if len(s) else {}
    per_field = {
        f: {"n": int(np.sum(s["field"] == f)), "medians": medians(*_cols(s, f))}
        for f in FIELDS
        if np.any(s["field"] == f)
    }
    return {"n": len(s), "medians": med, "per_field": per_field}


def _cols(s: Table, field: str):
    m = s["field"] == field
    return s["pa_sep"][m], s["pa_obs"][m]


def main() -> dict:
    t = load_rows()
    primary = t["gate_passed"] & (t["deblend_class"] == "two-knot")
    secondary = t["gate_passed"] & (t["ellipticity"] >= MIN_ELL)
    res = {
        "plan": "results/arcs/pa-convention/PLAN.md",
        "label": "derived",
        "assumptions": {
            "max_median_deg": MAX_MEDIAN,
            "min_margin_deg": MIN_MARGIN,
            "min_rows": MIN_ROWS,
            "secondary_min_ellipticity": MIN_ELL,
        },
        "primary": summarise(t, primary),
        "secondary": summarise(t, secondary),
    }
    res["outcome"] = outcome(res["primary"]["n"], res["primary"]["medians"])
    OUT.mkdir(parents=True, exist_ok=True)
    t.write(OUT / "rows.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))
    return res


if __name__ == "__main__":
    main()
