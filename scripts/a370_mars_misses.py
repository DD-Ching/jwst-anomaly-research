"""Abell 370: the MARS held-out images under BUFFALO model A (same anchors, frozen plan).

Frozen plan: ``results/lens_consistency/a370_mars_misses/PLAN.md``. Re-uses the seen / held-out
split and anchors of ``mars_hff_heldout`` (MARS catalogue, D-034 addendum ap) and solves every
Abell 370 held-out image with the BUFFALO model A maps (``abell370-buffalo-a-maps``). Per image:
``classify``; per field: ``mars_hff_heldout.gate_h`` (reported). No flag reading.

    python scripts/a370_mars_misses.py --zip DIR
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import lens_consistency as lc
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import requests
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/a370_mars_misses")
FIELD = "Abell 370"
MODEL = "abell370-buffalo-a-maps"
MARS_TABLE = Path("results/lens_consistency/mars_hff_heldout/heldout_A370.ecsv")
REPRODUCED_ARCSEC = 3.0  # ASSUMPTION (PLAN.md rule 2; gate H's match radius)


def classify(d_mars: float, d_buffalo: float, r: float = REPRODUCED_ARCSEC) -> str:
    """PLAN.md rule 2: which model reproduces the image within ``r`` (NaN / inf = not)."""
    m = bool(np.isfinite(d_mars) and d_mars <= r)
    b = bool(np.isfinite(d_buffalo) and d_buffalo <= r)
    return {(True, True): "both", (True, False): "MARS only", (False, True): "BUFFALO only"}.get(
        (m, b), "neither"
    )


def run(zdir: Path, out: Path = OUT) -> dict:
    meta = mh.fetch_zip(requests.Session(), zdir)
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        cat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
    mars = cat[np.asarray(cat["field"]) == FIELD]
    canucs = Table.read(mo.FIELDS[FIELD][0])
    off = mo.frame_offset(canucs, mars)
    seen = mo.split_seen(mo.shift(canucs, off), mars)["seen"]
    t = mo.split_seen(canucs, mars[:0])  # unshifted positions (PLAN.md rule 1)
    t["seen"] = seen
    model, _, _ = lc.load_model(MODEL)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    r = mo.heldout_residuals(model, grid, t)
    prev = Table.read(MARS_TABLE)
    dm = dict(zip(prev["image_id"], np.asarray(prev["dtheta_arcsec"], float), strict=True))
    r["dtheta_mars_arcsec"] = [dm[i] for i in r["image_id"]]
    r["class"] = [
        classify(a, b) for a, b in zip(r["dtheta_mars_arcsec"], r["dtheta_arcsec"], strict=True)
    ]
    out.mkdir(parents=True, exist_ok=True)
    r.write(out / "residuals.ecsv", format="ascii.ecsv", overwrite=True)
    cls = np.asarray(r["class"], str)
    summary = {
        "plan": str(out / "PLAN.md"),
        "zip": meta,
        "model": MODEL,
        "n_tested": len(r),
        "gate_H_buffalo_reported": mo.gate_h(r),
        "classes": {
            k: int(np.sum(cls == k)) for k in ("both", "MARS only", "BUFFALO only", "neither")
        },
        "mars_misses": {
            str(i): {"class": str(c), "dtheta_buffalo_arcsec": float(d), "n_predicted": int(n)}
            for i, c, d, n in zip(
                r["image_id"], cls, r["dtheta_arcsec"], r["n_predicted"], strict=True
            )
            if c in ("BUFFALO only", "neither")
        },
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the MARS zip")
    a = ap.parse_args(argv)
    try:
        s = run(a.zip)
    finally:
        (a.zip / mh.ZIP).unlink(missing_ok=True)
    print(
        json.dumps(
            _clean({k: s[k] for k in ("gate_H_buffalo_reported", "classes", "mars_misses")}),
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
