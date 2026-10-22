"""MACS1149: extra-image test with CANUCS DR1 reference fluxes, the local-noise rule and a
containing-source counterpart rule (third-field gate for the CANUCS DR1 reference).

Frozen plan: ``results/lens_consistency/macs1149_canucsref/PLAN.md``. Stage 1 runs gate K of
``m0417_images`` on the Sharon gold list with the CANUCS v1 maps; stage 2 builds the predicted and
control tables with ``m0417_extra.run_field``; stage 3 runs ``m0417_arcref.run`` with the CANUCS
DR1 catalogue (``m0417_canucsref.to_cat``), ``contain=True`` and the ``m1423_canucsref``
local-noise rule. Everything is written to the plan's directory.

    python scripts/m1149_canucsref.py --cat <photometry-cat.fits.gz> --image <f200w bgsub-sci> \
        [--cache data/cache/m1149x]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import m0417_arcref as ar
import m0417_canucsref as cr
import m0417_extra as mx
import m0417_images as mi
import m1423_canucsref as m14
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from m0417_canucsref_posthoc import _clean

from jwst_anomaly.provenance import file_sha256

FIELD = "macs1149"
MODEL = "macs1149-canucs-maps"
OUT = Path("results/lens_consistency/macs1149_canucsref")
GOLD = Path("results/lens_consistency/macs1149_sharon/sharon_gold_images.ecsv")
IMAGES = "sharon_gold_images.ecsv"
OBS, TILE = "o008", "t004"
BASE = "https://archive.stsci.edu/hlsps/canucs/dr1/macs1149/clu/"
CAT_NAME = "hlsp_canucs_jwst-hst_multi_macs1149-clu_multi_v1_photometry-cat.fits.gz"
IMAGE_NAME = "hlsp_canucs_jwst_nircam_macs1149-clu-40mas_f200w_v1_bgsub-sci.fits.gz"
FRAME = OUT / "frame"  # stage-2 predicted / control tables and their summary.json
GRID_STEP_K = 0.25  # m0417_images default (not stated in PLAN.md; deviation, see RESULT.md)
CAT_SHA = "08ab67347f2c3dfe4f743cc1b2a9d4b77a1e40eb66ddeb611648bb296adc0739"
IMAGE_SHA = "34d8f2609fc49f0c88fc93d4ddcccc613fee27348204bdbeb86297d4e8339c26"


def write_images(out: Path = OUT) -> Path:
    """The gold list with ``z_spec`` = true (plan: all gold images are spectroscopic)."""
    t = Table.read(GOLD, format="ascii.ecsv")
    t["z_spec"] = np.ones(len(t), bool)
    t.meta["source"] = f"{GOLD} (z_spec column added: gold = spectroscopic, PLAN.md)"
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def register() -> None:
    """Point the field tables of the reused scripts at this field and this directory."""
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    mx.FIELDS[FIELD] = (IMAGES, MODEL, OBS, TILE)
    mx.IMAGES = OUT
    mx.OUT = FRAME


def counterpart_offsets(cat: dict) -> dict:
    """Reported, not a gate: offsets (arcsec, East / North) of gold images from the centroid of
    their rule-1 CANUCS counterpart."""
    t = Table.read(OUT / IMAGES)
    d = []
    for r in t:
        k = ar.counterpart(float(r["ra"]), float(r["dec"]), cat)
        if k < 0:
            continue
        dx = (cat["ra"][k] - r["ra"]) * np.cos(np.radians(r["dec"])) * 3600.0
        d.append((float(dx), float((cat["dec"][k] - r["dec"]) * 3600.0)))
    a = np.asarray(d, float).reshape(-1, 2)
    return {
        "n": len(a),
        "median_dra_arcsec": float(np.median(a[:, 0])) if len(a) else None,
        "median_ddec_arcsec": float(np.median(a[:, 1])) if len(a) else None,
    }


def stage3(cat, depth, inputs, d, w, src_x, src_y, contain: bool, out: Path) -> dict:
    s = ar.run(
        Path("."),
        cat=cat,
        depth=depth,
        cat_inputs=inputs,
        out=out,
        field=FIELD,
        frame=FRAME,
        offset=(0.0, 0.0),
        contain=contain,
    )
    rng = np.random.default_rng(m14.SEED)
    ctrl = Table.read(out / f"control_{FIELD}.ecsv")
    pred = Table.read(out / f"predicted_{FIELD}.ecsv")
    m14.apply_local(ctrl, d, w, src_x, src_y, rng)
    if len(pred):
        m14.apply_local(pred, d, w, src_x, src_y, rng)
    c2 = ar.gate_c2(list(ctrl["class_local"]), list(ctrl["arc_ref"]))
    flags = [
        f"{r['system']}@{r['ra']:.6f},{r['dec']:.6f}"
        for r in pred
        if c2["passes"] and str(r["class_local"]) == "absent"
    ]
    ctrl.write(out / f"control_{FIELD}.ecsv", overwrite=True)
    pred.write(out / f"predicted_{FIELD}.ecsv", overwrite=True)
    return {
        "gate_C2": c2,
        "control_classes": m14._count(ctrl["class_local"]),
        "target_classes": m14._count(pred["class_local"]) if len(pred) else {},
        "flags": flags,
        "without_local_noise_rule": {
            k: s[k] for k in ("gate_C2", "control_classes", "target_classes", "flags")
        },
        "posthoc_absent_controls": s["posthoc_absent_controls"],
        "inputs": s["inputs"],
        "rules": s["rules"],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m1149x"))
    args = ap.parse_args(argv)
    # cached names carry a sha256 prefix, so the files are checked by content
    if file_sha256(args.cat) != CAT_SHA or file_sha256(args.image) != IMAGE_SHA:
        raise ValueError("--cat / --image are not the CANUCS DR1 MACS1149 files named in the plan")
    write_images()
    register()
    summary: dict = {}
    k = mi.run_field(FIELD, OUT, GRID_STEP_K)  # stage 1: gate K
    summary["stage1_gate_K"] = k
    if MODEL not in k["passing_models"]:
        summary["stopped"] = "gate K failed: model not used (PLAN.md rule 1)"
    else:
        FRAME.mkdir(parents=True, exist_ok=True)
        stage2 = mx.run_field(FIELD, args.cache)  # stage 2
        summary["stage2_pipeline_forced"] = stage2
        (FRAME / "summary.json").write_text(
            json.dumps(_clean({FIELD: stage2}), indent=1, allow_nan=False) + "\n"
        )
        with fits.open(args.image) as hl:
            if hl[0].header.get("BUNIT") != "nJy":
                raise ValueError(f"unexpected BUNIT {hl[0].header.get('BUNIT')!r}")
            w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
        if abs(proj_plane_pixel_scales(w.celestial).mean() * 3600.0 - m14.PIX) > 1e-4:
            raise ValueError("image pixel scale is not 0.04 arcsec")
        t = Table.read(args.cat)
        cat = cr.to_cat(t)
        inputs = {
            "canucs_cat": CAT_NAME,
            "canucs_url": BASE + CAT_NAME,
            "canucs_cat_sha256": CAT_SHA,
            "n_rows_finite_flux": int(len(cat["flux"])),
            "flux_column": f"FLUX_KRON_TOTAL_{cr.BAND}",
            "image": IMAGE_NAME,
            "image_url": BASE + IMAGE_NAME,
            "image_sha256": IMAGE_SHA,
        }
        depth = cr.depth_mag(cat["flux"], cat["snr"])
        src_x, src_y = (
            np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["RA"], t["DEC"], unit="deg"))
        )
        cmp_dir = OUT / "without_contain"
        summary["comparison_without_contain"] = stage3(
            cat, depth, inputs, d, w, src_x, src_y, False, cmp_dir
        )
        summary["frozen"] = stage3(cat, depth, inputs, d, w, src_x, src_y, True, OUT)
        summary["frame_check_counterpart_offsets"] = counterpart_offsets(cat)
        summary["rules_local_noise"] = {
            "aperture_radius_arcsec": m14.AP_ARCSEC,
            "n_random": m14.N_RAND,
            "annulus_arcsec": list(m14.R_RAND),
            "source_exclusion_arcsec": m14.SRC_EXCL,
            "seed": m14.SEED,
            "min_pred_over_sigma": m14.MIN_SIGMA,
            "absent_pred_over_sigma": m14.ABSENT_SIGMA,
        }
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k.startswith(("stopped", "frozen"))}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
