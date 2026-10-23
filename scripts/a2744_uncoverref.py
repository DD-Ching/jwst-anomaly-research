"""Abell 2744: gate K of the UNCOVER-reference extra-image test (stages 1 and 1b).

Frozen plans: ``results/lens_consistency/abell2744_uncoverref/PLAN.md`` (stage 1) and ``PLAN_1b.md``
(``--stage 1b``: gate K with the WLS source of ``a2744_uncover``; a deviation after stage 1 failed).
Runs gate K of ``m0417_images`` on the spectroscopic images of the UNCOVER v2.0 multiple-image
table with the v2.0 deflection maps (``abell2744-uncover2-maps``). The extra-image stages follow
in a later cycle only if the gate passes. Writes ``summary.json`` to the plan's directory.

    python scripts/a2744_uncoverref.py [--stage 1b]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_uncover as au
import lens_consistency as lc
import m0417_images as mi
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel
from jwst_anomaly.photometry import fetch_catalog

FIELD = "abell2744"
MODEL = "abell2744-uncover2-maps"
OUT = Path("results/lens_consistency/abell2744_uncoverref")
IMAGES = "uncover2_images.ecsv"
GRID_STEP_K = 0.25  # m0417_images default (PLAN.md rule 1)

lc.MODELS[MODEL] = {
    "files": {k: au.UNCOVER_V2[k] for k in ("alpha_x", "alpha_y")},
    "kind": "maps",
    "z_lens": au.Z_LENS,
    "cosmology": au.COSMOLOGY,
}


def lenstool_id(image_id: str) -> str:
    """UNCOVER id -> an id whose ``lensmodel.image_family`` is ``au.uncover_system(image_id)``:
    "S.KI" -> "S.K" + letter I, "S.I" -> "S" + letter I (PLAN.md Inputs)."""
    idx = int(image_id.split(".")[1][-1])
    if idx == 0:  # "S.K0" / "S.10" cannot be told apart once stored as floats
        raise ValueError(f"UNCOVER id {image_id!r}: image number 0 is ambiguous")
    return f"{au.uncover_system(image_id)}{chr(ord('a') + idx - 1)}"


def image_list(raw: Table) -> Table:
    """The release's image table with ``z_spec`` from ``Redshift_type``; ``None`` rows dropped."""
    kind = np.asarray(raw["Redshift_type"], str)
    keep = kind != "None"
    ids = [repr(float(i)) for i in raw["ID"][keep]]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate UNCOVER image ids")
    t = Table(
        {
            "image_id": [lenstool_id(i) for i in ids],
            "uncover_id": ids,
            "ra": np.asarray(raw["RA"], float)[keep],
            "dec": np.asarray(raw["Dec"], float)[keep],
            "z": np.asarray(raw["z"], float)[keep],
            "z_spec": kind[keep] == "spec",
        }
    )
    t["ra"].unit = t["dec"].unit = "deg"
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{au.UNCOVER_V2['images'][0]} (UNCOVER v2.0 multiple images; z_spec = Redshift_type spec)"
    )
    return t


def read_summary() -> dict:
    p = OUT / "summary.json"
    return json.loads(p.read_text()) if p.exists() else {}


def write_summary(summary: dict) -> None:
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")


def gate_k_wls(model, grid, sp: Table) -> tuple[Table, dict]:
    """PLAN_1b.md rule 1: ``a2744_uncover.imageplane`` (WLS family source) on the spectroscopic
    images ``sp``, decided by the unchanged ``m0417_images.gate_k``; rule 2's reported counts."""
    bt = lensmodel.backtrace_images(model, sp, {})
    ip = au.imageplane(model, grid, bt)
    d = np.asarray(ip["dtheta_arcsec"], float)
    bad = ~(np.isfinite(d) & (d <= mi.MATCH_ARCSEC))
    counts = {}
    for fam in dict.fromkeys(np.asarray(ip["system"], str)):
        m = np.asarray(ip["system"], str) == fam
        counts[fam] = (int(m.sum()), int(np.max(ip["n_predicted"][m])))
    return ip, {
        "gate_K": mi.gate_k(ip),
        "estimator": "wls_source (a2744_uncover.imageplane)",
        "not_reproduced": {
            str(i): float(x) for i, x in zip(ip["image_id"][bad], d[bad], strict=True)
        },
        "shared_match": int(np.sum(ip["shared_match"])),
        "extra_predicted": {f: c for f, c in counts.items() if c[1] > c[0]},
        "fewer_predicted": {f: c for f, c in counts.items() if c[1] < c[0]},
    }


def stage_1b() -> dict:
    """PLAN_1b.md: gate K-WLS on the committed image list, same Solver grid as stage 1."""
    images = mi.load_images(OUT / IMAGES)
    sp = images[images["z_spec"]]
    solver = mi.Solver(MODEL, images, GRID_STEP_K)
    ip, k = gate_k_wls(solver.model, solver.grid, sp)
    ip.meta.update(solver.model._meta())
    ip.meta["provenance"] = "model_prediction"
    ip.meta["source"] = f"{MODEL}, WLS family source (PLAN_1b.md)"
    lc._write(ip, OUT / f"imageplane_wls_{FIELD}_{MODEL}.ecsv")
    k["model_sha256"] = solver.model.sha256
    if k["gate_K"]["passes"]:
        k["next"] = "gate K-WLS passed (deviation): freeze stage 2 with the WLS source"
    else:
        k["stopped"] = "gate K-WLS failed: model not used (PLAN_1b.md rule 1)"
    return k


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", choices=("1", "1b"), default="1")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "1b":
        summary = read_summary()
        summary["stage1b_gate_K_wls"] = _clean(stage_1b())
        summary["status"] = (
            "stage 1 stopped (mean source); stage 1b passed (WLS source, deviation): stage 2 next"
            if summary["stage1b_gate_K_wls"]["gate_K"]["passes"]
            else "stage 1 and stage 1b failed: Abell 2744 extra-image route stopped with this model"
        )
        write_summary(summary)
        print(json.dumps(summary["stage1b_gate_K_wls"]["gate_K"]))
        return 0
    raw = Table.read(fetch_catalog(*au.UNCOVER_V2["images"]))
    image_list(raw).write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    k = mi.run_field(FIELD, OUT, GRID_STEP_K)
    # a stage-1 rerun keeps the later stages' keys (stage 1b adds to the same file)
    summary: dict = {
        x: v for x, v in read_summary().items() if x.startswith("stage1b") or x == "status"
    }
    summary["stage1_gate_K"] = k
    if MODEL not in k["passing_models"]:
        summary["stopped"] = "gate K failed: model not used (PLAN.md rule 1)"
    else:
        summary["next"] = "gate K passed: pin an UNCOVER reference, run the extra-image stages"
    write_summary(_clean(summary))
    print(json.dumps({"gate_K": k[MODEL]["gate_K"], **{x: summary.get(x) for x in ("stopped",)}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
