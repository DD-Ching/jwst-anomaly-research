"""HST-era MACS1149 models (CATS v4.1, Sharon v4cor) vs the Schuldt et al. 2024 new images (held-out).

Frozen plan: ``results/lens_consistency/m1149_hst_heldout/PLAN.md``. Rules 1-5 of the
``mars_hff_heldout`` stage, unchanged, per model, with that model's own catalogued images as
the "seen" list. No flag, flux or extra-image reading.

    python scripts/m1149_hst_heldout.py
"""

from __future__ import annotations

import json
from pathlib import Path

import lens_consistency as lc
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import mars_m1149_heldout as m1
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/m1149_hst_heldout")
SHARON_GOLD = Path("results/lens_consistency/macs1149_sharon/sharon_gold_images.ecsv")
MODELS = ("macs1149-cats", "macs1149-sharon")


def seen_list(name: str, files: dict, par) -> Table:
    """The model's own catalogued images (PLAN.md Inputs): CATS ``arcs.dat``, Sharon's gold list."""
    if name == "macs1149-sharon":
        return Table.read(SHARON_GOLD)
    return lc.image_list(name, files, par)[0]


def model_block(name: str, images: Table, out: Path) -> dict:
    model, files, par = lc.load_model(name)
    seen = seen_list(name, files, par)
    off = mo.frame_offset(images, seen)
    t = mo.split_seen(mo.shift(images, off), seen)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    r = mo.heldout_residuals(model, grid, t)
    r["dtheta_arcsec"].info.format = ".3f"
    r.write(out / f"heldout_{name}.ecsv", format="ascii.ecsv", overwrite=True)
    g = mo.gate_h(r)
    held = ~np.asarray(t["seen"], bool)
    d = np.asarray(r["dtheta_arcsec"], float) if len(r) else np.array([])
    return {
        "files": {k: Path(v).name for k, v in files.items()},
        "seen_list_rows": len(seen),
        "frame_offset": off,
        "seen": int((~held).sum()),
        "held_out": int(held.sum()),
        "types": {k: int(np.sum(np.asarray(r["type"]) == k)) for k in "ENM"} if len(r) else {},
        "outside_map": int(np.sum(~np.asarray(r["inside"], bool))) if len(r) else 0,
        "gate_H": g,
        "outcome": mo.outcome(g),
        "over_3arcsec": [str(i) for i in r["image_id"][d > 3.0]] if len(r) else [],
    }


def run(out: Path = OUT) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    images = Table.read(m1.IMAGES)
    summary = {
        "plan": str(out / "PLAN.md"),
        "images": str(m1.IMAGES),
        "provenance": "positions and z: observed (Schuldt+2024 Table 2; model image lists); "
        "residuals: model_prediction; thresholds: ASSUMPTION (mars_hff_heldout/PLAN.md)",
        "models": {n: model_block(n, images, out) for n in MODELS},
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    s = run()
    print(json.dumps(_clean(s["models"]), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
