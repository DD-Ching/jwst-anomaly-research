"""Abell 370: extra-image test with the BUFFALO model A maps and CANUCS DR1 reference fluxes.

Frozen plan: ``results/lens_consistency/abell370_buffalo_extra/PLAN.md``. Same stages as
``m0416_canucsref``: gate K of ``m0417_images`` on the 104 images of the BUFFALO gold list
(``a370_buffalo.write_buffalo_images``) with ``abell370-buffalo-a-maps``; predicted and control
tables from ``m0417_extra.run_field`` (JWST pipeline level-3 F200W of 1208 ``o002_t001``);
``m0416_canucsref.stage3`` with the CANUCS DR1 Abell 370 catalogue and F200W mosaic. The CANUCS
file checksums are pinned at first download (``CAT_SHA`` / ``IMAGE_SHA``), before any run; the
script refuses to run while they are unset.

    python scripts/a370_buffalo_extra.py --cat <photometry-cat.fits.gz> --image <f200w bgsub-sci> \
        [--cache data/cache/a370x]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a370_buffalo as ab
import m0416_canucsref as m16
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

FIELD = "abell370"
MODEL = ab.MODEL  # abell370-buffalo-a-maps (gate K2 passed, ../abell370_buffalo)
OUT = Path("results/lens_consistency/abell370_buffalo_extra")
IMAGES = ab.IMAGES_K2  # the BUFFALO gold list, written to OUT
OBS, TILE = "o002", "t001"  # TARGPROP ABELL370 (i2d header read 2026-10-22)
BASE = "https://archive.stsci.edu/hlsps/canucs/dr1/a370/clu/"
CAT_NAME = "hlsp_canucs_jwst-hst_multi_a370-clu_multi_v1_photometry-cat.fits.gz"
IMAGE_NAME = "hlsp_canucs_jwst_nircam_a370-clu-40mas_f200w_v1_bgsub-sci.fits.gz"
FRAME = OUT / "frame"
GRID_STEP_K = 0.25  # m0417_images default (PLAN.md rule 1)
CAT_SHA: str | None = (
    "f5622f2867aa3094df6861981ee67f7f1879b8381b348224748c7381436ae17b"  # pinned 2026-10-22
)
IMAGE_SHA: str | None = "81501effc709f3484fce989c221b2e88d2fefbe6fbdd23bcd308e6f95ffd3d7d"


def register() -> None:
    """Point the field tables of the reused scripts at this field and this directory."""
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    mx.FIELDS[FIELD] = (IMAGES, MODEL, OBS, TILE)
    mx.IMAGES = OUT
    mx.OUT = FRAME


def check_inputs(cat: Path, image: Path) -> None:
    if CAT_SHA is None or IMAGE_SHA is None:
        raise SystemExit("error: pin CAT_SHA / IMAGE_SHA at first download (PLAN.md) before a run")
    # cached names carry a sha256 prefix, so the files are checked by content
    if file_sha256(cat) != CAT_SHA or file_sha256(image) != IMAGE_SHA:
        raise SystemExit("error: --cat / --image are not the CANUCS DR1 Abell 370 files pinned")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/a370x"))
    args = ap.parse_args(argv)
    check_inputs(args.cat, args.image)
    ab.write_buffalo_images(OUT)
    register()
    summary: dict = {}
    k = mi.run_field(FIELD, OUT, GRID_STEP_K)  # stage 1: gate K (repeats gate K2)
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
        args3 = (cat, depth, inputs, d, w, src_x, src_y)
        summary["comparison_without_contain"] = m16.stage3(
            *args3, False, OUT / "without_contain", field=FIELD, frame=FRAME
        )
        summary["frozen"] = m16.stage3(*args3, True, OUT, field=FIELD, frame=FRAME)
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k.startswith(("stopped", "frozen"))}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
