"""Vet 78.3c, test 1: knot c-image positions of Abell 2744 system 78 in the CATS v4.1 maps.

Plan (frozen first): ``results/lens_consistency/abell2744_uncoverref/PLAN_vet78.md``. For each knot
(78.1, 78.2, 78.3) the a and b images are back-traced, the ``wls_source`` is solved for its images,
and the image farthest from the a/b midpoint is the c-type image (ASSUMPTION). ``decide`` is the
plan's frozen rule. Everything computed is ``model_prediction``; image positions are ``observed``.

Usage: python scripts/a2744_vet78_cats.py [--out FILE]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_s700 as s700  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MODEL = "abell2744-cats"
HERE = Path("results/lens_consistency/abell2744_uncoverref")
IMAGES = HERE / "uncover2_images_spec.ecsv"
KNOTS = ("78.1", "78.2", "78.3")
Z = 2.3151  # UNCOVER catalogue spec-z of system 78
UNCOVER_78_3C = (3.575198, -30.355756)  # stage-2c predicted 78.3c (UNCOVER v2.0 maps)
TOOL_ARCSEC = 1.5  # D-041 match radius (ASSUMPTION)
SEARCH_ARCSEC = 1.0  # stage-2c search radius (ASSUMPTION)


def sep_arcsec(ra1, dec1, ra2, dec2) -> float:
    """Small-angle separation (arcsec)."""
    dra = (ra1 - ra2) * np.cos(np.deg2rad(0.5 * (dec1 + dec2)))
    return float(np.hypot(dra, dec1 - dec2) * 3600.0)


def c_type(pred: Table, ab: Table) -> int | None:
    """Index of the predicted image farthest from the a/b midpoint (None if fewer than 3)."""
    if len(pred) < 3:
        return None
    ra0, dec0 = float(np.mean(ab["ra"])), float(np.mean(ab["dec"]))
    d = [sep_arcsec(float(r["ra"]), float(r["dec"]), ra0, dec0) for r in pred]
    return int(np.argmax(d))


def decide(knots: dict[str, dict]) -> str:
    """The plan's frozen rule on the per-knot rows (``c_to_catalogued_arcsec``,
    ``c_to_uncover_78_3c_arcsec``; ``c_ra`` None when no c-type image)."""
    for k in ("78.1", "78.2"):
        d = knots[k]
        if d["c_ra"] is None or not d["c_to_own_catalogued_arcsec"] <= TOOL_ARCSEC:
            return "tool not validated"
    t = knots["78.3"]
    if t["c_ra"] is None:
        return "tool not validated"
    if t["c_to_catalogued_arcsec"] <= SEARCH_ARCSEC:
        return "explained: model position error of the shared c image"
    if t["c_to_uncover_78_3c_arcsec"] > SEARCH_ARCSEC:
        return "untestable: an independent model disagrees"
    return "flag survives test 1"


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=HERE / "vet78_cats.json")
    args = ap.parse_args(argv)
    imgs = Table.read(IMAGES)
    ids = np.asarray(imgs["image_id"]).astype(str)
    cat_c = {k: imgs[ids == f"{k}c"] for k in ("78.1", "78.2")}
    model, files, _ = lc.load_model(MODEL)
    try:
        sys78 = imgs[np.char.startswith(ids, "78.")]
        grid = lensmodel.DeflectionGrid.compute(
            model, lc.grid_half_width(model, sys78), s700.GRID_STEP
        )
        knots = {}
        for k in KNOTS:
            ab = imgs[np.isin(ids, [f"{k}a", f"{k}b"])]
            pred, beta = s700.solve(model, ab, grid, Z)
            i = c_type(pred, ab)
            row: dict = {
                "n_predicted": len(pred),
                "source_spread_arcsec": s700._r(float(np.ptp(beta, axis=0).max())),
                "mu": [s700._r(float(m), 2) for m in pred["magnification"]],
                "c_ra": None,
            }
            if i is not None:
                cra, cdec = float(pred["ra"][i]), float(pred["dec"][i])
                d_cat = {
                    c: sep_arcsec(cra, cdec, float(t["ra"][0]), float(t["dec"][0]))
                    for c, t in cat_c.items()
                }
                row.update(
                    c_ra=cra,
                    c_dec=cdec,
                    c_abs_mu=s700._r(abs(float(pred["magnification"][i])), 2),
                    c_to_catalogued_arcsec=s700._r(min(d_cat.values())),
                    c_to_catalogued_each_arcsec={c: s700._r(v) for c, v in d_cat.items()},
                    c_to_uncover_78_3c_arcsec=s700._r(sep_arcsec(cra, cdec, *UNCOVER_78_3C)),
                )
                if k in cat_c:
                    row["c_to_own_catalogued_arcsec"] = row["c_to_catalogued_each_arcsec"][k]
            knots[k] = row
            print(k, row, flush=True)
        out = {
            "model": MODEL,
            "z": Z,
            "grid_step_arcsec": s700.GRID_STEP,
            "grid_half_width_arcsec": float(grid.x[-1]),
            "thresholds_arcsec": {"tool": TOOL_ARCSEC, "search": SEARCH_ARCSEC},
            "knots": knots,
            "outcome": decide(knots),
            "provenance": "model_prediction (positions observed)",
        }
        args.out.write_text(json.dumps(out, indent=1) + "\n")
    finally:
        if os.environ.get("CLAUDE_CODE_REMOTE") == "true":  # cloud-disk rule
            for p in files.values():
                Path(p).unlink(missing_ok=True)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
