"""Abell 370: gate K of the BUFFALO Niemiec+2023 model A deflection maps (stage 1 only).

Frozen plan: ``results/lens_consistency/abell370_buffalo/PLAN.md``. Gate L checks the map
loading (κ from the deflection Hessian vs the published κ map); gate K runs ``m0417_images`` on
the Lagattuta et al. 2019 gold spectroscopic images. Writes ``summary.json`` to the plan's
directory. ``archive.stsci.edu`` needs a CA bundle with Mozilla's roots in cloud sessions
(docs/operations.md).

    python scripts/a370_buffalo.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import m1149_sharon as ms
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly.photometry import fetch_catalog

FIELD = "abell370"
MODEL = "abell370-buffalo-a-maps"
OUT = Path("results/lens_consistency/abell370_buffalo")
IMAGES = "lagattuta19_gold_zspec.ecsv"
GRID_STEP_K = 0.25  # m0417_images default (PLAN.md rule 2)
L_MAX_MEDIAN_REL = ms.G1_MAX_MEDIAN_REL  # gate L (PLAN.md rule 1; ASSUMPTION)
L_MAP_STEP = 37  # κ sub-grid step in map pixels, as m0417_images.map_checks
LAGATTUTA19 = (
    "https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/485/3738/table1.dat",
    "6c7e41483a179b07d43f062790b4e868bc29162019c59abe8ebf95290960cdc5",
)


def parse_table1(path: Path) -> Table:
    """CDS ``table1.dat`` of Lagattuta et al. 2019 by its byte-by-byte ReadMe (1-based columns)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        line = line.ljust(55)
        zspec = line[26:32].strip()
        rows.append(
            (
                line[0:4].strip(),
                float(line[5:14]),
                float(line[15:25]),
                float(zspec) if zspec else np.nan,
                line[45:51].strip(),
                line[52:55].strip(),
            )
        )
    return Table(rows=rows, names=("image_id", "ra", "dec", "zspec", "class", "note"))


def select(t: Table) -> np.ndarray:
    """PLAN.md selection: gold, zspec present, Note not containing 1."""
    cls = np.asarray(t["class"], str)
    notes = np.asarray(t["note"], str)
    return (
        (cls == "gold")
        & np.isfinite(np.asarray(t["zspec"], float))
        & np.array(["1" not in n.split(",") for n in notes])
    )


def write_images(out: Path = OUT) -> Path:
    """The selected images in the ``m0417_images`` list format (z = zspec, fixed: ASSUMPTION)."""
    src = parse_table1(fetch_catalog(*LAGATTUTA19))
    sel = src[select(src)]
    t = Table({"image_id": sel["image_id"], "ra": sel["ra"], "dec": sel["dec"], "z": sel["zspec"]})
    t["ra"].unit = t["dec"].unit = "deg"
    t["z_spec"] = np.ones(len(t), bool)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{LAGATTUTA19[0]} (Lagattuta et al. 2019; gold, zspec present, Note without 1; "
        "z = zspec fixed, PLAN.md ASSUMPTION)"
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def gate_l(check: dict) -> dict:
    rel = check["median_rel_diff"]
    return {**check, "passes": bool(np.isfinite(rel) and rel <= L_MAX_MEDIAN_REL)}


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args(argv)
    write_images()
    model, _, _ = lc.load_model(MODEL)
    check = lc.map_check(model, lc.model_files(MODEL)["kappa_map"], "kappa", None, L_MAP_STEP)
    check["map"] = Path(check["map"]).name  # cache file name, not a local path
    summary: dict = {"gate_L": gate_l(check)}
    if not summary["gate_L"]["passes"]:
        summary["stopped"] = "gate L failed: map loading not verified (PLAN.md rule 1)"
    else:
        mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
        k = mi.run_field(FIELD, OUT, GRID_STEP_K)
        summary["stage1_gate_K"] = k
        if MODEL not in k["passing_models"]:
            summary["stopped"] = "gate K failed: model not used (PLAN.md rule 2)"
        else:
            summary["next"] = (
                "gate K passed: pin the CANUCS DR1 catalogue and mosaic, run stages 2-3"
            )
    summary["rules"] = {
        "L_MAX_MEDIAN_REL": L_MAX_MEDIAN_REL,
        "K_MAX_RMS_ARCSEC": mi.K_MAX_RMS_ARCSEC,
        "K_MIN_MATCHED": mi.K_MIN_MATCHED,
        "MATCH_ARCSEC": mi.MATCH_ARCSEC,
        "GRID_STEP_K": GRID_STEP_K,
    }
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    k = summary.get("stage1_gate_K", {}).get(MODEL, {}).get("gate_K")
    print(json.dumps({"gate_L": summary["gate_L"]["median_rel_diff"], "gate_K": k}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
