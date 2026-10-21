"""MACS0717: D-041 cross-model readings with the pinned ``macs0717-keeton`` and ``macs0717-sharon``.

Every non-``observed`` row of the CATS v4.1 ``images_predicted.ecsv`` is solved in each second
model's maps from its system's magnification-weighted least-squares source
(``m1149_refs.second_model_images``) and read with ``m1149_refs.d041_verdict``. A per-model
known-case gate (the share of catalogued images the model reproduces) must pass before that model's
readings are used. Flagged rows that a second model supports are looked up in the frozen Step 6
record of ``docs/fields/macs0717.md``. No photometry is re-measured. Frozen plan:
``results/lens_consistency/macs0717_crossmodel/PLAN.md``.

Positions and magnifications are ``model_prediction``; catalogued positions are ``observed``.

    python scripts/m0717_crossmodel.py [--out results/lens_consistency/macs0717_crossmodel]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_uncover as au  # noqa: E402
import lens_consistency as lc  # noqa: E402
import m1149_crossmodel as mc  # noqa: E402
import m1149_refs as mr  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

OUT = Path("results/lens_consistency/macs0717_crossmodel")
PRED = Path("results/lens_consistency/macs0717-cats/images_predicted.ecsv")
SECOND = ("macs0717-keeton", "macs0717-sharon")
GRID_STEP = 0.25  # arcsec, as the CATS run
FLAGGED = ("absent", "confused", "inconsistent_reference")
RECORD_MATCH_ARCSEC = 1.0  # Step 6 record lookup radius (ASSUMPTION)
_r = au._r

# Step 6 record of docs/fields/macs0717.md (frozen in PLAN.md):
# (system, ra, dec, verdict, independent)
STEP6 = (
    ("5", 109.379428, 37.746088, "copy of 5.1 (3.1 arcsec)", True),
    ("6", 109.374860, 37.770130, "copy of 6.3 (3.2 arcsec)", True),
    ("7", 109.366048, 37.765120, "position error, near-critical mu (D-036)", True),
    ("7", 109.365638, 37.764574, "position error, near-critical mu (D-036)", True),
    ("34", 109.373742, 37.744904, "copy of 34.3 (2.8 arcsec)", True),
    ("68", 109.395425, 37.739381, "on a critical curve, mu spread > 2x (D-036)", True),
    ("19", 109.409131, 37.754055, "copy of 19.1 (2.3 arcsec)", True),
    ("19", 109.397479, 37.747960, "central image below sensitivity", True),
    ("20", 109.371506, 37.763382, "below sensitivity (predicted S/N 5.8)", True),
    ("67", 109.377758, 37.768356, "counterpart within the position uncertainty", True),
    ("70", 109.403161, 37.746014, "below sensitivity", True),
    ("75", 109.382840, 37.744173, "mu spread > 2x (D-036)", True),
    ("75", 109.386305, 37.742128, "mu spread > 2x (D-036)", True),
    ("3", 109.393233, 37.747986, "CATS-only extra image", False),
    ("19", 109.392674, 37.752049, "CATS-only extra image", False),
    ("45", 109.384136, 37.743892, "CATS-only extra image", False),
    ("68", 109.391572, 37.744043, "CATS-only extra image", False),
    ("68", 109.384180, 37.744954, "CATS-only extra image", False),
    ("71", 109.398568, 37.743334, "CATS-only extra image", False),
)
# Step 6: systems of the 28 confused copy rows, plus the 65 pair (D-041 copy rule)
COPY_SYSTEMS = frozenset("1 3 7 12 17 24 32 34 45 50 56 62 67 71 72 75 76 79 80 65".split())


def combine(readings: list[str]) -> str:
    """Combined D-041 reading over the gated models that test the row (PLAN.md)."""
    for v in ("both models", "catalogued image misplaced by CATS", "model-dependent (CATS only)"):
        if v in readings:
            return v
    return ""


def outcome(forced_class: str, combined: str) -> str:
    """Frozen per-row outcome (PLAN.md)."""
    if forced_class not in FLAGGED:
        return "no_flag"
    if not combined:
        return "untestable"
    if combined == "both models":
        return "lensing_residual_input" if forced_class == "absent" else "second_model_supported"
    return "explained_model_dependent"


def step6_record(system: str, ra: float, dec: float):
    """The Step 6 record entry within ``RECORD_MATCH_ARCSEC`` of the row, or None."""
    best, best_d = None, RECORD_MATCH_ARCSEC
    for rec in STEP6:
        if rec[0] != system:
            continue
        d = float(au._sep_arcsec(np.array([rec[1]]), np.array([rec[2]]), ra, dec)[0])
        if d <= best_d:
            best, best_d = rec, d
    return best


def final(out: str, forced_class: str, system: str, ra: float, dec: float) -> str:
    """Frozen final outcome: a second-model-supported flag needs an independent Step 6 reason."""
    if out not in ("lensing_residual_input", "second_model_supported"):
        return out
    rec = step6_record(system, ra, dec)
    if rec is not None:
        return "explained_step6" if rec[4] else "to_vet_candidate"
    if forced_class == "confused" and system in COPY_SYSTEMS:
        return "explained_step6"
    if system == "65":  # the 65 pair's copies (inconsistent_reference) and flux-ratio flags
        return "explained_step6"
    return "to_vet_candidate"


def post_hoc_arcs_z(rows: list[dict], gates: dict, zcat: dict[str, bool]) -> dict:
    """Post hoc (deviation from PLAN.md): the frozen rules on the systems with an ``arcs.txt``
    redshift only, with the models whose secondary gate passes there."""
    used = [n for n, g in gates.items() if g["secondary_arcs_z"]["passes"]]
    out_rows, counts = [], {}
    for r in rows:
        if not zcat.get(r["system"]):
            continue
        comb = combine([r["models"][n]["d041"] for n in used if n in r["models"]])
        o = outcome(r["forced_class"], comb)
        f = final(o, r["forced_class"], r["system"], *r["radec"])
        counts[f] = counts.get(f, 0) + 1
        out_rows.append(
            {"system": r["system"], "radec": r["radec"], "forced_class": r["forced_class"]}
            | {"combined": comb, "outcome": o, "final": f}
        )
    return {"models_used": used, "outcome": counts, "n_rows": len(out_rows), "rows": out_rows}


def read_row(sh: Table, ids: np.ndarray, ra: float, dec: float, cats_ids: set[str]) -> dict:
    """D-041 reading of one CATS row against a second model's image set ``sh``."""
    d = au._sep_arcsec(sh["ra"], sh["dec"], ra, dec) if len(sh) else np.array([])
    ext, obs = np.flatnonzero(ids == ""), np.flatnonzero(ids != "")
    ke = int(ext[np.argmin(d[ext])]) if len(ext) else -1
    ko = int(obs[np.argmin(d[obs])]) if len(obs) else -1
    sep_e = float(d[ke]) if ke >= 0 else np.nan
    sep_o = float(d[ko]) if ko >= 0 else np.nan
    obs_id = str(ids[ko]) if ko >= 0 else ""
    return {
        "nearest_extra_sep_arcsec": _r(sep_e),
        "nearest_extra_mu": _r(float(sh["magnification"][ke])) if ke >= 0 else None,
        "nearest_catalogued_sep_arcsec": _r(sep_o),
        "nearest_catalogued_id": obs_id,
        "d041": mr.d041_verdict(sep_e, sep_o, obs_id in cats_ids),
    }


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred", type=Path, default=PRED)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    pred = Table.read(args.pred)
    psys = np.asarray(pred["system"]).astype(str)
    pclass = np.asarray(pred["image_class"]).astype(str)
    fc = np.array([au._col(pred, "forced_class", k, "") for k in range(len(pred))]).astype(str)
    tgt = np.flatnonzero(pclass != "observed")
    _, files, par = lc.load_model("macs0717-cats")
    arcs, _ = lc.image_list("macs0717-cats", files, par)
    asys = np.asarray(arcs["system"]).astype(str)
    zcat = {
        str(s): bool(np.isfinite(z) and z > 0)
        for s, z in zip(asys, np.ma.filled(np.ma.asarray(arcs["z"], float), np.nan), strict=True)
    }
    systems = list(dict.fromkeys(psys[tgt]))
    rows = [
        {
            "system": str(psys[k]),
            "radec": [round(float(pred["ra"][k]), 6), round(float(pred["dec"][k]), 6)],
            "image_class": str(pclass[k]),
            "mu_cats": _r(float(pred["magnification"][k])),
            "forced_class": str(fc[k]),
            "models": {},
        }
        for k in tgt
    ]
    gates = {}
    for name in SECOND:
        model, _, _ = lc.load_model(name)
        grid = lensmodel.DeflectionGrid.compute(
            model, lc.grid_half_width(model, arcs[np.isin(asys, systems)]), GRID_STEP
        )
        n_rep = n_test = n_rep_z = n_test_z = 0
        for s in systems:
            cat = arcs[asys == s]
            z = float(pred["z_sys"][np.flatnonzero(psys == s)[0]])
            cat_in = mc.inside_maps(model, cat["ra"], cat["dec"], z)
            if not cat_in.all():
                continue
            sh = mr.second_model_images(model, grid, cat, z)
            ids = mr.split_second(sh, cat)
            rep = len({str(i) for i in ids if i})
            n_rep, n_test = n_rep + rep, n_test + len(cat)
            if zcat.get(s):
                n_rep_z, n_test_z = n_rep_z + rep, n_test_z + len(cat)
            cats_ids = {
                str(m)
                for m in pred["matched_image"][(psys == s) & (pclass == "observed")]
                if str(m)
            }
            for k, row in zip(tgt, rows, strict=True):
                if psys[k] != s:
                    continue
                ra, dec = row["radec"]
                if not mc.inside_maps(model, ra, dec, z)[0]:
                    continue
                row["models"][name] = {"z": z, "n_images": len(sh)} | read_row(
                    sh, ids, float(pred["ra"][k]), float(pred["dec"][k]), cats_ids
                )
        gates[name] = {
            "reproduced": n_rep,
            "tested": n_test,
            "passes": mc.gate(n_rep, n_test),
            "secondary_arcs_z": {
                "reproduced": n_rep_z,
                "tested": n_test_z,
                "passes": mc.gate(n_rep_z, n_test_z),
            },
        }
    used = [n for n in SECOND if gates[n]["passes"]]
    counts: dict[str, int] = {}
    for row in rows:
        row["combined"] = combine([row["models"][n]["d041"] for n in used if n in row["models"]])
        row["outcome"] = outcome(row["forced_class"], row["combined"])
        rec = step6_record(row["system"], *row["radec"])
        row["step6"] = rec[3] if rec else ""
        row["final"] = final(row["outcome"], row["forced_class"], row["system"], *row["radec"])
        counts[row["final"]] = counts.get(row["final"], 0) + 1
    post_hoc = post_hoc_arcs_z(rows, gates, zcat)
    comparison = [
        {
            "system": rec[0],
            "radec": [rec[1], rec[2]],
            "pinned_combined": next(
                (r["combined"] for r in rows if step6_record(r["system"], *r["radec"]) is rec),
                None,
            ),
        }
        for rec in STEP6
        if not rec[4]
    ]
    summary = {
        "provenance": "positions and mu: model_prediction (CATS v4.1, Keeton v4, Sharon v4cor); "
        "catalogued positions: observed; forced classes: derived (CATS images run 2026-10-21)",
        "inputs": {"pred": str(args.pred), "second_models": list(SECOND)},
        "assumptions": {
            "grid_step_arcsec": GRID_STEP,
            "match_arcsec": au.MATCH_ARCSEC,
            "xmodel_loose_arcsec": mr.XMODEL_LOOSE_ARCSEC,
            "gate_min_reproduced": mc.GATE_MIN_REPRODUCED,
            "record_match_arcsec": RECORD_MATCH_ARCSEC,
        },
        "forced_class_counts": {v: int((fc == v).sum()) for v in dict.fromkeys(fc)},
        "image_class_counts": {v: int((pclass == v).sum()) for v in dict.fromkeys(pclass)},
        "gates": gates,
        "models_used": used,
        "outcome": counts if used else "tool not validated",
        "n_rows": len(rows),
        "n_systems": len(systems),
        "combined_counts": {
            v: sum(r["combined"] == v for r in rows)
            for v in dict.fromkeys(r["combined"] for r in rows)
        },
        "comparison_cats_only": comparison,
        "post_hoc_arcs_z": post_hoc,
        "rows": rows,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


if __name__ == "__main__":
    s = main()
    print("gates", json.dumps(s["gates"]), "outcome", s["outcome"])
    print("combined", s["combined_counts"], "comparison", s["comparison_cats_only"])
    for r in s["rows"]:
        if r["forced_class"] in FLAGGED:
            m = {n[9:]: v.get("nearest_extra_sep_arcsec") for n, v in r["models"].items()}
            print(r["system"], r["forced_class"], r["combined"], r["final"], m, r["step6"])
