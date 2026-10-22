"""Abell 370: gate K of the CANUCS DR1 reference extra-image test (stage 1 only).

Frozen plan: ``results/lens_consistency/abell370_canucsref/PLAN.md``. Runs gate K of
``m0417_images`` on the 115 images of the pinned CANUCS Lenstool model (``abell370-canucs``,
D-044). Stages 2–3 of ``m0416_canucsref`` follow in a later cycle only if the gate passes (the
CANUCS DR1 Abell 370 catalogue and mosaic are pinned then). Writes ``summary.json`` to the plan's
directory.

    python scripts/a370_canucsref.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

FIELD = "abell370"
MODEL = "abell370-canucs"
OUT = Path("results/lens_consistency/abell370_canucsref")
IMAGES = "canucs_a370_images.ecsv"
GRID_STEP_K = 0.25  # m0417_images default (PLAN.md rule 1)


def write_images(out: Path = OUT) -> Path:
    """The model's multiple-image file, ``z_spec`` = true (PLAN.md ASSUMPTION: its z are fixed)."""
    _, files, _ = lc.load_model(MODEL)
    src = lensmodel.load_lenstool_images(files["arcs.dat"])
    t = Table({k: src[k] for k in ("image_id", "ra", "dec", "z")})
    t["ra"].unit = t["dec"].unit = "deg"
    t["z_spec"] = np.ones(len(t), bool)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{lensmodel.ABELL370_CANUCS['arcs.dat'][0]} (CANUCS DR1 Lenstool multfile; "
        "z_spec = true: fixed z, PLAN.md ASSUMPTION)"
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args(argv)
    write_images()
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    k = mi.run_field(FIELD, OUT, GRID_STEP_K)
    summary: dict = {"stage1_gate_K": k}
    if MODEL not in k["passing_models"]:
        summary["stopped"] = "gate K failed: model not used (PLAN.md rule 1)"
    else:
        summary["next"] = "gate K passed: pin the CANUCS DR1 catalogue and mosaic, run stages 2-3"
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps({"gate_K": k[MODEL]["gate_K"], **{x: summary.get(x) for x in ("stopped",)}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
