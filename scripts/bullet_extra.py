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
# Frozen PLAN_1b.md rules W1 / W2 (ASSUMPTIONs)
W1_RMS_ARCSEC, W1_P99_ARCSEC = 0.01, 0.04
W2_N, W2_MIN_USABLE, W2_EDGE, W2_HALF, W2_ISOL_ARCSEC = 100, 50, 10, 4, 1.0
W2_MAX_MEDIAN_PIX, W2_SAT_FRAC, W2_SAT_PCT = 0.25, 0.5, 99.999


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


def header_wcs(header) -> WCS:
    """PLAN_1b.md W0: header CRPIX, CRVAL, CTYPE and PC; CD removed, CDELT 0.04 arcsec."""
    h = header.copy()
    for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2"):
        h.remove(k, ignore_missing=True)
    h["CDELT1"] = h["CDELT2"] = m14.PIX / 3600.0
    return WCS(h).celestial


def wcs_catalogue_check(w: WCS, x, y, ra, dec) -> dict:
    """PLAN_1b.md W1: ``w`` applied to catalogue X, Y read 0-based and 1-based, against RA, DEC."""
    x, y, ra, dec = (np.asarray(v, float) for v in (x, y, ra, dec))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(ra) & np.isfinite(dec)
    sky = SkyCoord(ra[ok] * u.deg, dec[ok] * u.deg)
    res: dict = {"n_rows": int(ok.sum())}
    for origin in (0, 1):
        sep = w.pixel_to_world(x[ok] - origin, y[ok] - origin).separation(sky).arcsec
        rms, p99 = float(np.sqrt(np.mean(sep**2))), float(np.percentile(sep, 99))
        res[f"origin_{origin}"] = {
            "rms_arcsec": rms,
            "p99_arcsec": p99,
            "passes": bool(rms <= W1_RMS_ARCSEC and p99 <= W1_P99_ARCSEC),
        }
    best = min((0, 1), key=lambda o: res[f"origin_{o}"]["rms_arcsec"])
    res["catalogue_origin"] = best
    res["passes"] = bool(res["origin_0"]["passes"] or res["origin_1"]["passes"])
    return res


def centroid_offsets(d, xs, ys, half: int = W2_HALF) -> np.ndarray:
    """Intensity-weighted centroid minus (x, y) in a (2 half + 1)^2 box at the nearest pixel, after
    subtracting the median of the box border; NaN rows where the box sum is not finite and > 0."""
    out = np.full((len(xs), 2), np.nan)
    k = np.arange(-half, half + 1)
    for i, (x, y) in enumerate(zip(xs, ys, strict=True)):
        ix, iy = int(np.rint(x)), int(np.rint(y))
        box = np.asarray(d[iy - half : iy + half + 1, ix - half : ix + half + 1], float)
        if box.shape != (2 * half + 1, 2 * half + 1):
            continue
        border = np.concatenate([box[0], box[-1], box[1:-1, 0], box[1:-1, -1]])
        b = box - np.median(border)
        s = b.sum()
        if not (np.isfinite(s) and s > 0):
            continue
        out[i] = (ix + (b.sum(0) * k).sum() / s - x, iy + (b.sum(1) * k).sum() / s - y)
    return out


def wcs_pixel_check(w: WCS, d, t: Table) -> dict:
    """PLAN_1b.md W2: centroids of bright isolated unsaturated point sources at ``w`` positions."""
    ra, dec = np.asarray(t["RA"], float), np.asarray(t["DEC"], float)
    ok = np.isfinite(ra) & np.isfinite(dec)
    sky = SkyCoord(ra[ok] * u.deg, dec[ok] * u.deg)
    _, sep2, _ = sky.match_to_catalog_sky(sky, nthneighbor=2)
    isolated = np.zeros(len(t), bool)
    isolated[np.flatnonzero(ok)] = sep2.arcsec > W2_ISOL_ARCSEC
    with np.errstate(invalid="ignore"):
        x, y = (np.asarray(v, float) for v in w.world_to_pixel_values(ra, dec))
    ny, nx = d.shape
    inside = (x >= W2_EDGE) & (x <= nx - 1 - W2_EDGE) & (y >= W2_EDGE) & (y <= ny - 1 - W2_EDGE)
    sat = W2_SAT_FRAC * float(np.nanpercentile(d, W2_SAT_PCT))
    cand = np.flatnonzero(np.asarray(t["FLAG_POINTSRC"], bool) & isolated & inside)
    cand = cand[np.argsort(-np.asarray(t["SNR_DETIMG"], float)[cand], kind="stable")]
    chosen = []
    for i in cand:
        ix, iy = int(np.rint(x[i])), int(np.rint(y[i]))
        pk = np.nanmax(d[iy - W2_HALF : iy + W2_HALF + 1, ix - W2_HALF : ix + W2_HALF + 1])
        if pk < sat:
            chosen.append(i)
        if len(chosen) == W2_N:
            break
    chosen = np.asarray(chosen, int)
    off = centroid_offsets(d, x[chosen], y[chosen])
    good = np.isfinite(off).all(1)
    n = int(good.sum())
    mdx = float(np.median(off[good, 0])) if n else float("nan")
    mdy = float(np.median(off[good, 1])) if n else float("nan")
    return {
        "n_selected": len(chosen),
        "n_usable": n,
        "median_dx_pix": mdx,
        "median_dy_pix": mdy,
        "saturation_guard_njy": sat,
        "passes": bool(
            n >= W2_MIN_USABLE and abs(mdx) <= W2_MAX_MEDIAN_PIX and abs(mdy) <= W2_MAX_MEDIAN_PIX
        ),
    }


def read_image_1b(path: Path, t: Table) -> tuple:
    """CANUCS F200W mosaic with the PLAN_1b.md WCS -> (WCS, nJy per pixel, scale, W1, W2)."""
    with fits.open(path) as hl:
        hdu = next(h for h in hl if h.data is not None)
        if hdu.header.get("BUNIT") != "nJy":
            raise ValueError(f"unexpected BUNIT {hdu.header.get('BUNIT')!r}")
        w, d = header_wcs(hdu.header), np.asarray(hdu.data, float)
    pix = float(proj_plane_pixel_scales(w).mean() * 3600.0)
    if abs(pix - m14.PIX) > 1e-4:
        raise ValueError(f"image pixel scale {pix} is not 0.04 arcsec")
    w1 = wcs_catalogue_check(w, t["X"], t["Y"], t["RA"], t["DEC"])
    w2 = wcs_pixel_check(w, d, t) if w1["passes"] else None
    return w, d, pix, w1, w2


def run_1b(cat: Path, image: Path) -> dict:
    """PLAN_1b.md: rule 1 recomputed, W1 / W2, then PLAN.md rules 4-8 on the stored frame tables."""
    t = Table.read(cat)
    gold = Table.read(GATEK / GOLD_REL)
    summary: dict = {"plan": str(OUT / "PLAN.md"), "deviation": str(OUT / "PLAN_1b.md")}
    summary["rule1_frame_F"] = f = frame_check_f(gold["ra"], gold["dec"], t["RA"], t["DEC"])
    summary["rule2_gate_K_wls"] = "stored: summary_wls.json (not re-run, PLAN_1b.md)"
    summary["rule3_frame_tables"] = "stored: frame/ (not re-run, PLAN_1b.md)"
    if not f["passes"]:
        summary["stopped"] = "frame check F failed on recompute (PLAN_1b.md)"
        return summary
    w, d, pix, w1, w2 = read_image_1b(image, t)
    summary["W1_catalogue_wcs"], summary["W2_pixel_wcs"] = w1, w2
    if not w1["passes"]:
        summary["stopped"] = "W1 failed (PLAN_1b.md)"
        return summary
    if not w2["passes"]:
        summary["stopped"] = "W2 failed (PLAN_1b.md)"
        return summary
    register()
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
        "image_wcs": "PLAN_1b.md W0 (header CRPIX/CRVAL/PC, CDELT 0.04 arcsec)",
        "noise_pixel_scale_arcsec": pix,
    }
    depth = cr.depth_mag(rc["flux"], rc["snr"])
    src_x, src_y = (np.asarray(v, float) for v in w.world_to_pixel_values(t["RA"], t["DEC"]))
    for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 8 / rules 5-7
        out.mkdir(parents=True, exist_ok=True)
        summary["frozen" if contain else "without_contain"] = ur._reference_run(
            rc, depth, inputs, out, contain, d, w, src_x, src_y, pix, field=FIELD, frame=FRAME
        )
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/bulletx"))
    ap.add_argument(
        "--plan-1b", action="store_true", help="PLAN_1b.md: mosaic WCS gates, then rules 4-8"
    )
    args = ap.parse_args(argv)
    check_inputs(args.cat, args.image)
    OUT.mkdir(parents=True, exist_ok=True)
    if args.plan_1b:
        return _write(run_1b(args.cat, args.image))
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
    return _write(summary)


def _write(summary: dict) -> int:
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(
        json.dumps(
            {k: summary[k] for k in summary if k.startswith(("rule1", "W1", "W2", "stopped"))}
            | {k: summary[k].get("gate_C2_shared") for k in ("frozen",) if k in summary}
            | {"flags": summary.get("frozen", {}).get("flags")},
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
