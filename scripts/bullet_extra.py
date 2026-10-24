"""Bullet Cluster: extra-image test with the CANUCS DR1 best-fit maps (WLS source) and CANUCS DR1
reference fluxes.

Frozen plan: ``results/lens_consistency/bullet_extra/PLAN.md``. Rule 1 (frame check F) compares the
gold list with the CANUCS catalogue; rule 2 repeats gate K-WLS (``bullet_gatek.run_wls``); rule 3
builds the frame tables with ``m0417_extra.run_field`` and the WLS source
(``a2744_uncover.wls_backtrace``) on GO 4598 ``o001_t001``; rules 4-8 are
``a2744_uncoverref._reference_run`` with the CANUCS catalogue and F200W bCG-subtracted mosaic. The
CANUCS file checksums are pinned (``CAT_SHA`` / ``IMAGE_SHA``) before any run; the script refuses
to run while they are unset.

    python scripts/bullet_extra.py --cat <photometry-cat.fits> --image <f200w bgsub-sci.fits.gz> \
        [--cache data/cache/bulletx]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_uncover as au
import a2744_uncoverref as ur
import astropy.units as u
import bullet_gatek as bg
import m0417_canucsref as cr
import m0417_extra as mx
import m1423_canucsref as m14
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from m0417_canucsref_posthoc import _clean

from jwst_anomaly.provenance import file_sha256

FIELD = "bullet"
MODEL = bg.MODEL  # bullet-canucs-maps (gate K-WLS passed, ../bullet_gatek)
OUT = Path("results/lens_consistency/bullet_extra")
FRAME = OUT / "frame"
GATEK = Path(
    "results/lens_consistency/bullet_gatek"
)  # holds imageplane_bullet_<MODEL>.ecsv (z_used)
GOLD_REL = "../bullet_inputs/bullet_gold_images.ecsv"  # relative to GATEK
OBS, TILE = "o001", "t001"  # TARGPROP BulletCluster; all 135 gold positions in its grid
S3_4598 = "s3://stpubdata/jwst/public/jw04598/L3/t/{0}/jw04598-{0}_{1}_nircam_clear-f200w_{2}"
BASE = "https://www.canfar.net/storage/vault/file/canucs/BULLET_DR1/"
CAT_NAME = "catalogs/hlsp_canucs_jwst-hst_multi_bullet-clu_multi_v1_photometry-cat.fits"
IMAGE_NAME = "images/hlsp_canucs_jwst_nircam_bullet-clu-40mas_f200w_v1_bgsub-sci.fits.gz"
CAT_SHA: str | None = (
    "536b4d56cad8f00a89a79b7bfcf1acbe1cf68ac37fea7af8355623cf2d830f8e"  # pinned 2026-10-24
)
IMAGE_SHA: str | None = "5023e937be3071faed1dc4219139dade69aaf7d590b46dcdecf7d511743f2bac"
# Frozen rule 1 (PLAN.md; ASSUMPTIONs)
F_MATCH_ARCSEC = 0.5
F_MAX_MEDIAN_ARCSEC = 0.1
F_MIN_N = 30


def frame_check_f(gra, gdec, cra, cdec) -> dict:
    """Rule 1: nearest catalogue row within ``F_MATCH_ARCSEC`` of each gold image; passes iff
    ``n >= F_MIN_N`` and |median dx| and |median dy| <= ``F_MAX_MEDIAN_ARCSEC`` (gold - catalogue,
    dx along RA cos dec, arcsec)."""
    cra, cdec = np.asarray(cra, float), np.asarray(cdec, float)
    ok = np.isfinite(cra) & np.isfinite(cdec)
    g = SkyCoord(np.asarray(gra, float) * u.deg, np.asarray(gdec, float) * u.deg)
    c = SkyCoord(cra[ok] * u.deg, cdec[ok] * u.deg)
    idx, sep, _ = g.match_to_catalog_sky(c)
    m = sep.arcsec <= F_MATCH_ARCSEC
    dx, dy = c[idx[m]].spherical_offsets_to(g[m])
    dx, dy = dx.arcsec, dy.arcsec
    n = int(m.sum())
    mdx = float(np.median(dx)) if n else float("nan")
    mdy = float(np.median(dy)) if n else float("nan")
    return {
        "n_gold": len(g),
        "n_matched": n,
        "median_dx_arcsec": mdx,
        "median_dy_arcsec": mdy,
        "median_sep_arcsec": float(np.median(sep.arcsec[m])) if n else float("nan"),
        "rules": {
            "match_arcsec": F_MATCH_ARCSEC,
            "max_abs_median_arcsec": F_MAX_MEDIAN_ARCSEC,
            "min_n": F_MIN_N,
        },
        "passes": bool(
            n >= F_MIN_N and abs(mdx) <= F_MAX_MEDIAN_ARCSEC and abs(mdy) <= F_MAX_MEDIAN_ARCSEC
        ),
    }


def register() -> None:
    """Point ``m0417_extra``'s field table at this field, tile and directory."""
    mx.FIELDS[FIELD] = (GOLD_REL, MODEL, OBS, TILE)
    mx.IMAGES = GATEK
    mx.OUT = FRAME
    mx.S3 = S3_4598


def check_inputs(cat: Path, image: Path) -> None:
    if CAT_SHA is None or IMAGE_SHA is None:
        raise SystemExit("error: pin CAT_SHA / IMAGE_SHA (PLAN.md) before a run")
    if file_sha256(cat) != CAT_SHA or file_sha256(image) != IMAGE_SHA:
        raise SystemExit("error: --cat / --image are not the CANUCS DR1 Bullet files pinned")


def read_image(path: Path) -> tuple:
    """CANUCS F200W mosaic -> (WCS, nJy per pixel, pixel scale); nJy and 0.04 arcsec checked."""
    with fits.open(path) as hl:
        hdu = next(h for h in hl if h.data is not None)
        if hdu.header.get("BUNIT") != "nJy":
            raise ValueError(f"unexpected BUNIT {hdu.header.get('BUNIT')!r}")
        w, d = WCS(hdu.header), np.asarray(hdu.data, float)
    pix = float(proj_plane_pixel_scales(w.celestial).mean() * 3600.0)
    if abs(pix - m14.PIX) > 1e-4:
        raise ValueError(f"image pixel scale {pix} is not 0.04 arcsec")
    return w, d, pix


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/bulletx"))
    args = ap.parse_args(argv)
    check_inputs(args.cat, args.image)
    OUT.mkdir(parents=True, exist_ok=True)
    t = Table.read(args.cat)
    gold = Table.read(GATEK / GOLD_REL)
    summary: dict = {"plan": str(OUT / "PLAN.md")}
    summary["rule1_frame_F"] = f = frame_check_f(gold["ra"], gold["dec"], t["RA"], t["DEC"])
    if not f["passes"]:
        summary["stopped"] = "frame check F failed (PLAN.md rule 1)"
    else:
        k = bg.run_wls(OUT)["gate_K_wls"]  # rule 2
        summary["rule2_gate_K_wls"] = k
        if not k["passes"]:
            summary["stopped"] = "gate K-WLS failed on repeat (PLAN.md rule 2)"
    if "stopped" not in summary:
        register()
        FRAME.mkdir(parents=True, exist_ok=True)
        frame = mx.run_field(FIELD, args.cache, bt_hook=au.wls_backtrace)  # rule 3
        summary["rule3_frame_tables"] = frame
        (FRAME / "summary.json").write_text(
            json.dumps(_clean({FIELD: frame}), indent=1, allow_nan=False) + "\n"
        )
        w, d, pix = read_image(args.image)
        rc = cr.to_cat(t)  # rule 4
        inputs = {
            "canucs_cat": CAT_NAME,
            "canucs_url": BASE + CAT_NAME,
            "canucs_cat_sha256": CAT_SHA,
            "n_rows_finite_flux": int(len(rc["flux"])),
            "flux_column": f"FLUX_KRON_TOTAL_{cr.BAND}",
            "image": IMAGE_NAME,
            "image_url": BASE + IMAGE_NAME,
            "image_sha256": IMAGE_SHA,
            "noise_pixel_scale_arcsec": pix,
        }
        depth = cr.depth_mag(rc["flux"], rc["snr"])
        src_x, src_y = (
            np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["RA"], t["DEC"], unit="deg"))
        )
        for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 8 / rules 5-7
            out.mkdir(parents=True, exist_ok=True)
            summary["frozen" if contain else "without_contain"] = ur._reference_run(
                rc, depth, inputs, out, contain, d, w, src_x, src_y, pix, field=FIELD, frame=FRAME
            )
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(
        json.dumps(
            {k: summary[k] for k in summary if k.startswith(("rule1", "stopped"))}
            | {k: summary[k].get("gate_C2_shared") for k in ("frozen",) if k in summary}
            | {"flags": summary.get("frozen", {}).get("flags")},
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
