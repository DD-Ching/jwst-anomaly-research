"""Gate R for D-037's repeated-pair and model-dependent-mu rules as functions on the stored MACS1149
CATS forced table (plan: results/lens_consistency/d037_rules/PLAN.md)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import numpy as np
from astropy.table import Table

# D-037 hand labels of the three MACS1149 forced flags (system, forced_class) -> label
HAND_R = {
    ("16", "absent"): "model_topology",
    ("2", "confused"): "model_topology",
    ("1192", "confused"): "mu_model_dependent",
}


def labels(topo: str, mu_dep: bool | None) -> list[str]:
    out = []
    if topo:
        out.append("model_topology")
    if mu_dep:
        out.append("mu_model_dependent")
    return out


def gate(rows: list[dict]) -> dict:
    """Each hand flag must carry exactly its hand label (no other D-037 label)."""
    res = {}
    for (sys_id, cls), want in HAND_R.items():
        got = [r["labels"] for r in rows if r["system"] == sys_id and r["forced_class"] == cls]
        res[f"{sys_id}:{cls}"] = {"want": [want], "got": got, "ok": got == [[want]]}
    extra = [r for r in rows if r["labels"] and (r["system"], r["forced_class"]) not in HAND_R]
    res["rule_only_rows"] = [f"{r['system']}:{r['forced_class']}" for r in extra]
    res["passes"] = all(v["ok"] for k, v in res.items() if ":" in k) and not extra
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", default="macs1149-cats")
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--imageplane", type=Path, required=True, help="validate's imageplane.ecsv")
    ap.add_argument("--crossmodel", type=Path, required=True, help="D-041 cross-model summary")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    table = Table.read(a.table, format="ascii.ecsv")
    if hasattr(table["forced_class"], "filled"):
        table["forced_class"] = [str(x) for x in table["forced_class"].filled("")]
    ip = Table.read(a.imageplane, format="ascii.ecsv")
    shared = [str(i) for i in ip["image_id"][ip["shared_match"]]]
    model, files, par = lc.load_model(a.model)
    images, _ = lc.image_list(a.model, files, par)
    lc.apply_frame_offset(a.model, model, images)
    topo = lc.repeated_pair_topology(table, images, shared)
    xm_rows = json.loads(a.crossmodel.read_text())["rows"]

    def second_model(sys_id: str, ra: float, dec: float) -> dict:
        """The cross-model row of the same system within 0.05" (its radec is rounded to 1e-6 deg;
        rounding both sides to 1e-5 deg mis-keys half-way cases)."""
        for x in xm_rows:
            if x["system"] == sys_id and (
                abs(x["radec"][0] - ra) * 3600 * np.cos(np.radians(dec)) < 0.05
                and abs(x["radec"][1] - dec) * 3600 < 0.05
            ):
                return x.get("sharon") or {}
        return {}

    rows, diag = [], []
    all_topo = lc.repeated_pair_topology(
        table, images, shared, classes=tuple(sorted({str(x) for x in table["forced_class"]}))
    )
    for r, t, t_all in zip(table, topo, all_topo, strict=True):
        if r["image_class"] == "observed":
            continue
        sh = second_model(str(r["system"]), float(r["ra"]), float(r["dec"]))
        mu2 = sh.get("nearest_extra_mu")
        rec = {
            "system": str(r["system"]),
            "forced_class": str(r["forced_class"]),
            "mu_first": round(float(r["magnification"]), 2),
            "mu_second_nearest": mu2,
            "second_sep_arcsec": sh.get("nearest_extra_sep_arcsec"),
            "mu_ratio": (
                round(
                    max(abs(r["magnification"]), abs(mu2)) / min(abs(r["magnification"]), abs(mu2)),
                    2,
                )
                if mu2
                else None
            ),
            "shared_pair": t_all,
        }
        if r["forced_class"] in lc.COPY_CLASSES:
            md = lc.mu_model_dependent(float(r["magnification"]), mu2)
            rec.update(mu_model_dependent=md, labels=labels(t, md))
            rows.append(rec)
        else:
            diag.append(rec)
    summary = {
        "provenance": "derived (forced classes: stored 2026-10-14 run; mu: model_prediction)",
        "model": a.model,
        "inputs": {k: str(v) for k, v in vars(a).items() if k != "out"},
        "assumptions": {
            "pair_radius_arcsec": lc.PAIR_RADIUS_ARCSEC,
            "mu_model_factor": lc.MU_MODEL_FACTOR,
        },
        "validate_rms_arcsec": round(float((ip["dtheta_arcsec"] ** 2).mean() ** 0.5), 4),
        "n_imageplane": len(ip),
        "shared_match_images": shared,
        "flag_rows": rows,
        "diagnostic_other_rows": diag,
        "gate_R": gate(rows),
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
