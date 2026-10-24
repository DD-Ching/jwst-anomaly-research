"""HST-era MACS1149 models (CATS v4.1, Sharon v4cor) vs the Schuldt+2024 new images (held-out).

Frozen plan: ``results/lens_consistency/m1149_hst_heldout/PLAN.md``. Rules 1-5 of the
``mars_hff_heldout`` stage, unchanged, per model, with that model's own catalogued images as
the "seen" list. No flag, flux or extra-image reading.

    python scripts/m1149_hst_heldout.py [--out DIR]

``posthoc.json`` (not in the frozen plan) compares the same images across the four MACS1149
models from the stored residual tables.
"""

from __future__ import annotations

import argparse
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
MARS_DIR = Path("results/lens_consistency/mars_m1149_heldout")
TABLES = {  # post hoc comparison: model -> stored per-image residual table
    "cats": OUT / "heldout_macs1149-cats.ecsv",
    "sharon": OUT / "heldout_macs1149-sharon.ecsv",
    "mars": MARS_DIR / "heldout_MACSJ1149.ecsv",
    "canucs_saw_them": MARS_DIR / "posthoc_canucs.ecsv",
}
GROUPS = {
    "202-204": ("202", "203", "204"),
    "201-204": ("201", "202", "203", "204"),
    "401": ("401",),
}


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
    tested = set(np.asarray(r["image_id"], str)) if len(r) else set()
    d = np.asarray(r["dtheta_arcsec"], float) if len(r) else np.array([])
    return {
        "files": {k: Path(v).name for k, v in files.items()},
        "seen_list_rows": len(seen),
        "frame_offset": off,
        "seen": int((~held).sum()),
        "seen_rows": {
            str(i): round(float(d), 3)
            for i, d in zip(t["image_id"][~held], t["sep_mars_arcsec"][~held], strict=True)
        },
        "held_out": int(held.sum()),
        "held_out_untested": [str(i) for i in t["image_id"][held] if str(i) not in tested],
        "types": {k: int(np.sum(np.asarray(r["type"]) == k)) for k in "ENM"} if len(r) else {},
        "outside_map": int(np.sum(~np.asarray(r["inside"], bool))) if len(r) else 0,
        "gate_H": g,
        "outcome": mo.outcome(g),
        "over_3arcsec": [str(i) for i in r["image_id"][d > 3.0]] if len(r) else [],
    }


def posthoc(tables: dict[str, Path] = TABLES) -> dict:
    """Post hoc: rms and max residual per family group and model, on each model's tested images."""
    out: dict = {"label": "post hoc; not in the frozen plan", "groups": {}}
    for gname, fams in GROUPS.items():
        row = {}
        for m, path in tables.items():
            r = Table.read(path)
            d = np.asarray(r["dtheta_arcsec"], float)[np.isin(np.asarray(r["system"], str), fams)]
            d = d[np.isfinite(d)]
            row[m] = (
                {
                    "n": len(d),
                    "rms_arcsec": round(float(np.sqrt(np.mean(d**2))), 3),
                    "max_arcsec": round(float(d.max()), 2),
                }
                if len(d)
                else {"n": 0}
            )
        out["groups"][gname] = row
    return out


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
    if out == OUT:
        (out / "posthoc.json").write_text(json.dumps(posthoc(), indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT)
    s = run(ap.parse_args(argv).out)
    print(json.dumps(_clean(s["models"]), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
