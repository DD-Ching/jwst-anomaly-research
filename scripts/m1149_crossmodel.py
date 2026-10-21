"""MACS1149: D-041 cross-model rule with the pinned ``macs1149-sharon`` on every CATS extra row.

Every non-``observed`` row of the tracked CATS v4.1 ``images_predicted.ecsv`` (16 rows) is solved in
the Sharon v4cor maps from its system's magnification-weighted least-squares source
(``m1149_refs.second_model_images``) and read with ``m1149_refs.d041_verdict``. A known-case gate
(the share of catalogued images Sharon reproduces) must pass before any verdict is used. No
photometry is re-measured; the forced classes come from the CATS table. Frozen plan:
``results/lens_consistency/macs1149_crossmodel/PLAN.md``.

Positions and magnifications are ``model_prediction``; catalogued positions are ``observed``.

    python scripts/m1149_crossmodel.py [--out results/lens_consistency/macs1149_crossmodel]
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
import m1149_refs as mr  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

OUT = Path("results/lens_consistency/macs1149_crossmodel")
PRED = Path("results/lens_consistency/macs1149-cats/images_predicted.ecsv")
SCRATCH = Path("results/lens_consistency/macs1149-cats/no_reference_refs.json")
SECOND = "macs1149-sharon"
GRID_STEP = 0.25  # arcsec, as the CATS run and m1149_refs
GATE_MIN_REPRODUCED = (
    0.8  # share of catalogued images reproduced within au.MATCH_ARCSEC (ASSUMPTION)
)
SECONDARY_Z = {"16": 2.7}  # D-037 addendum 2026-10-14: z ≈ 2.5–2.9 (pre-declared in PLAN.md)
_r = au._r


def outcome(inside: bool, forced_class: str, d041: str) -> str:
    """Frozen per-row outcome (PLAN.md)."""
    if not inside:
        return "untestable"
    if forced_class == "absent":
        return "lensing_residual_input" if d041 == "both models" else "explained_model_dependent"
    return "no_flag"


def inside_maps(model, ra, dec, z: float) -> np.ndarray:
    """Whether each position lies inside the map model's deflection maps (finite deflection)."""
    x, y = model.to_frame(np.asarray(ra, float), np.asarray(dec, float))
    lm = model.lens_map(np.atleast_1d(x), np.atleast_1d(y), np.full(np.size(x), z))
    return np.isfinite(np.asarray(lm["beta_x"], float))


def gate(n_reproduced: int, n_tested: int) -> bool:
    return n_tested > 0 and n_reproduced / n_tested >= GATE_MIN_REPRODUCED


def solve_row(model, grid, cat: Table, row, z: float, cats_ids: set[str]) -> dict:
    """D-041 reading of one CATS extra row against ``model`` solved at ``z``."""
    ra, dec = float(row["ra"]), float(row["dec"])
    sh = mr.second_model_images(model, grid, cat, z)
    ids = mr.split_second(sh, cat)
    d = au._sep_arcsec(sh["ra"], sh["dec"], ra, dec) if len(sh) else np.array([])
    ext, obs = np.flatnonzero(ids == ""), np.flatnonzero(ids != "")
    ke = int(ext[np.argmin(d[ext])]) if len(ext) else -1
    ko = int(obs[np.argmin(d[obs])]) if len(obs) else -1
    sep_e = float(d[ke]) if ke >= 0 else np.nan
    sep_o = float(d[ko]) if ko >= 0 else np.nan
    obs_id = str(ids[ko]) if ko >= 0 else ""
    return {
        "z": z,
        "n_images": len(sh),
        "n_extra": len(ext),
        "reproduced": sorted(set(str(i) for i in ids if i)),
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
    _, files, par = lc.load_model("macs1149-cats")
    arcs, _ = lc.image_list("macs1149-cats", files, par)
    asys = np.asarray(arcs["system"]).astype(str)
    model, _, _ = lc.load_model(SECOND)
    systems = list(dict.fromkeys(psys[tgt]))
    grid = lensmodel.DeflectionGrid.compute(
        model, lc.grid_half_width(model, arcs[np.isin(asys, systems)]), GRID_STEP
    )
    scratch = {}
    if SCRATCH.exists():
        for r in json.loads(SCRATCH.read_text())["rows"]:
            scratch[(r["system"], round(r["radec"][0], 5), round(r["radec"][1], 5))] = r["sharon"]
    rows, n_rep, n_test = [], 0, 0
    for s in systems:
        cat = arcs[asys == s]
        z = float(pred["z_sys"][np.flatnonzero(psys == s)[0]])
        cat_in = inside_maps(model, cat["ra"], cat["dec"], z)
        # CATS reproduces a catalogued image when one of its observed rows matched it
        cats_ids = {
            str(m) for m in pred["matched_image"][(psys == s) & (pclass == "observed")] if str(m)
        }
        for k in tgt[psys[tgt] == s]:
            r = pred[k]
            inside = bool(cat_in.all() and inside_maps(model, r["ra"], r["dec"], z)[0])
            row = {
                "system": s,
                "radec": [round(float(r["ra"]), 6), round(float(r["dec"]), 6)],
                "image_class": str(r["image_class"]),
                "mu_cats": _r(float(r["magnification"])),
                "forced_class": str(fc[k]),
                "inside_footprint": inside,
                "catalogued": [str(i) for i in cat["image_id"]],
            }
            if inside:
                row["sharon"] = solve_row(model, grid, cat, r, z, cats_ids)
                if s in SECONDARY_Z:
                    row["sharon_secondary"] = solve_row(
                        model, grid, cat, r, SECONDARY_Z[s], cats_ids
                    )
            d041 = row.get("sharon", {}).get("d041", "")
            row["outcome"] = outcome(inside, row["forced_class"], d041)
            key = (s, round(row["radec"][0], 5), round(row["radec"][1], 5))
            if key in scratch:
                old = scratch[key]["d041"]
                row["scratch_d041"] = old
                row["scratch_agrees"] = old == d041
            rows.append(row)
        if cat_in.all():  # known cases: the system's catalogued images, at the CATS z
            n_rep += len(solve_row(model, grid, cat, cat[0], z, cats_ids)["reproduced"])
            n_test += len(cat)
    gate_ok = gate(n_rep, n_test)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    summary = {
        "provenance": "positions and mu: model_prediction (CATS v4.1, Sharon v4cor); catalogued "
        "positions: observed; forced classes: derived (CATS images run 2026-10-14)",
        "inputs": {"pred": str(args.pred), "second_model": SECOND},
        "assumptions": {
            "grid_step_arcsec": GRID_STEP,
            "match_arcsec": au.MATCH_ARCSEC,
            "xmodel_loose_arcsec": mr.XMODEL_LOOSE_ARCSEC,
            "gate_min_reproduced": GATE_MIN_REPRODUCED,
            "secondary_z": SECONDARY_Z,
        },
        "gate": {"reproduced": n_rep, "tested": n_test, "passes": gate_ok},
        "outcome": "tool not validated" if not gate_ok else counts,
        "n_systems": len(systems),
        "d041_counts": {
            v: sum(r.get("sharon", {}).get("d041") == v for r in rows)
            for v in dict.fromkeys(r.get("sharon", {}).get("d041", "") for r in rows)
        },
        "scratch_disagreements": [
            r["system"] for r in rows if "scratch_agrees" in r and not r["scratch_agrees"]
        ],
        "rows": rows,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


if __name__ == "__main__":
    s = main()
    print("gate", s["gate"], "outcome", s["outcome"], "scratch diff", s["scratch_disagreements"])
    for r in s["rows"]:
        sh = r.get("sharon", {})
        print(
            r["system"],
            r["image_class"],
            r["forced_class"],
            r["inside_footprint"],
            sh.get("nearest_extra_sep_arcsec"),
            sh.get("nearest_catalogued_sep_arcsec"),
            sh.get("d041"),
            r["outcome"],
            r.get("sharon_secondary", {}).get("d041", ""),
        )
