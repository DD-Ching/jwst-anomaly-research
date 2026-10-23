"""Abell 2744: gate K of the UNCOVER-reference extra-image test (stage 1 only).

Frozen plan: ``results/lens_consistency/abell2744_uncoverref/PLAN.md``. Runs gate K of
``m0417_images`` on the spectroscopic images of the UNCOVER v2.0 multiple-image table with the
v2.0 deflection maps (``abell2744-uncover2-maps``). The extra-image stages follow in a later cycle
only if the gate passes. Writes ``summary.json`` to the plan's directory.

    python scripts/a2744_uncoverref.py
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


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    raw = Table.read(fetch_catalog(*au.UNCOVER_V2["images"]))
    image_list(raw).write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    k = mi.run_field(FIELD, OUT, GRID_STEP_K)
    summary: dict = {"stage1_gate_K": k}
    if MODEL not in k["passing_models"]:
        summary["stopped"] = "gate K failed: model not used (PLAN.md rule 1)"
    else:
        summary["next"] = "gate K passed: pin an UNCOVER reference, run the extra-image stages"
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps({"gate_K": k[MODEL]["gate_K"], **{x: summary.get(x) for x in ("stopped",)}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
