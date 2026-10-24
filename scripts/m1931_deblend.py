"""Own deblended F200W reference catalogue (photutils, cluster light filtered) and its stage V
validation on MACS0417 through the unchanged frozen ``m0417_arcref`` route.

Frozen plan: ``results/lens_consistency/macs1931_counter/PLAN.md`` ("Stage V: frozen rules").

    python scripts/m1931_deblend.py stage-v    # V1 (rules 1-6)
    python scripts/m1931_deblend.py stage-v2   # V2 amendment
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import m0417_arcref as ar
import m0417_canucsref as cr
import m0417_canucsref_posthoc as cp
import m0417_extra as mx
import m1423_canucsref as m14
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from jwst_anomaly.paths import data_root

LC = Path("results/lens_consistency")
FRAME = LC / "macs0417_frame"
OUT = LC / "macs0417_ownref"
OUT_V2 = LC / "macs0417_ownref_v2"
# Frozen rules (PLAN.md "Stage V: frozen rules"; ASSUMPTIONs)
MARGIN_ARCSEC = 15.0
BKG_BOX, BKG_FILTER, BKG_SIGMA = 64, 3, 3.0
KERNEL_FWHM, KERNEL_SIZE = 2.0, 5
DETECT_NSIGMA, NPIXELS = 3.0, 25
NLEVELS, CONTRAST = 32, 0.001
EXT_MIN_SIGMA_PIX = 1.5


def section_box(ra, dec, wcs: WCS, shape: tuple[int, int], margin_pix: float) -> tuple:
    """Rule 1: (y0, y1, x0, x1) bounding every position plus ``margin_pix``, clipped to
    ``shape``."""
    x, y = wcs.celestial.world_to_pixel_values(np.asarray(ra, float), np.asarray(dec, float))
    x0 = max(int(np.floor(np.min(x) - margin_pix)), 0)
    x1 = min(int(np.ceil(np.max(x) + margin_pix)) + 1, shape[1])
    y0 = max(int(np.floor(np.min(y) - margin_pix)), 0)
    y1 = min(int(np.ceil(np.max(y) + margin_pix)) + 1, shape[0])
    if x1 <= x0 or y1 <= y0:
        raise ValueError("positions fall outside the mosaic")
    return y0, y1, x0, x1


def to_njy(a, pixar_sr: float) -> np.ndarray:
    """MJy/sr -> nJy per pixel."""
    return np.asarray(a, float) * pixar_sr * 1e15


def sky_pa(x, y, theta_deg, wcs: WCS) -> np.ndarray:
    """Rule 4: major-axis direction ``theta_deg`` (CCW from +x) at pixel (x, y) as a sky PA East
    of North (deg, 0-180), carried through the WCS by a 1-pixel step."""
    x, y, t = (np.asarray(v, float) for v in (x, y, np.radians(theta_deg)))
    c0 = wcs.pixel_to_world(x, y)
    c1 = wcs.pixel_to_world(x + np.cos(t), y + np.sin(t))
    return np.mod(c0.position_angle(c1).deg, 180.0)


def build(sci: np.ndarray, err: np.ndarray, wcs: WCS, bkg_box: int = BKG_BOX) -> tuple:
    """Rules 1-4 on one section (nJy per pixel): (catalogue in ``m0417_arcref`` form, metadata,
    deblended segmentation image, background-subtracted image with masked pixels NaN, (x, y)
    pixel positions of the kept rows). ``bkg_box`` (pixels) is stage 3c's one change."""
    from astropy.convolution import convolve
    from astropy.stats import SigmaClip
    from photutils.background import Background2D, SExtractorBackground
    from photutils.segmentation import (
        SourceCatalog,
        deblend_sources,
        detect_sources,
        make_2dgaussian_kernel,
    )

    mask = ~np.isfinite(sci) | ~np.isfinite(err) | (err <= 0)
    sci = np.where(mask, 0.0, sci)
    err = np.where(mask, 0.0, err)
    bkg = Background2D(
        sci,
        bkg_box,
        filter_size=BKG_FILTER,
        coverage_mask=mask,
        sigma_clip=SigmaClip(sigma=BKG_SIGMA),
        bkg_estimator=SExtractorBackground(),
        fill_value=0.0,
    )
    data = np.where(mask, 0.0, sci - bkg.background)
    kernel = make_2dgaussian_kernel(KERNEL_FWHM, size=KERNEL_SIZE)
    conv = convolve(data, kernel, mask=mask, normalize_kernel=True)
    thr = DETECT_NSIGMA * bkg.background_rms
    segm = detect_sources(conv, thr, npixels=NPIXELS, mask=mask)
    if segm is None:
        raise ValueError("no sources detected")
    deb = deblend_sources(
        conv,
        segm,
        npixels=NPIXELS,
        nlevels=NLEVELS,
        contrast=CONTRAST,
        mode="exponential",
        progress_bar=False,
    )
    sc = SourceCatalog(data, deb, convolved_data=conv, error=err, mask=mask, wcs=wcs)
    t = sc.to_table(
        columns=[
            "label", "xcentroid", "ycentroid", "sky_centroid", "semimajor_sigma",
            "semiminor_sigma", "orientation", "ellipticity", "kron_flux", "kron_fluxerr",
        ]
    )  # fmt: skip
    scale = float(proj_plane_pixel_scales(wcs.celestial).mean() * 3600.0)
    cat = to_cat(t, wcs, scale)
    if not len(cat["flux"]):
        raise ValueError("no source with a positive finite Kron flux and error")
    xy = (cat.pop("x"), cat.pop("y"))
    meta = {
        "bkg_box_pix": int(bkg_box),
        "n_segments_detected": int(segm.n_labels),
        "n_segments_deblended": int(deb.n_labels),
        "n_rows_kept": int(len(cat["flux"])),
        "pixel_scale_arcsec": scale,
        "masked_fraction": float(mask.mean()),
        "kron_snr_min": float(np.min(cat["snr"])),
        "n_rows_snr_4_6": int(np.sum((cat["snr"] >= 4.0) & (cat["snr"] <= 6.0))),
    }
    return cat, meta, deb, np.where(mask, np.nan, data), xy


def kron_depth(cat: dict) -> float:
    """V2-a: AB magnitude of 5 x the median Kron flux error (Jy) of the catalogue rows."""
    return float(-2.5 * np.log10(5.0 * np.median(cat["flux"] / cat["snr"]) / 3631.0))


def to_cat(t: Table, wcs: WCS, scale: float) -> dict:
    """Rule 4 column mapping (flux in Jy; rows without a positive finite Kron flux and error
    are dropped)."""
    f = np.asarray(t["kron_flux"], float)
    e = np.asarray(t["kron_fluxerr"], float)
    keep = keep_rows(f, e)
    t = t[keep]
    sky = t["sky_centroid"]
    a = np.asarray(t["semimajor_sigma"], float)
    return {
        "label": np.asarray(t["label"]),
        "ra": sky.ra.deg,
        "dec": sky.dec.deg,
        "a": a * scale,
        "b": np.asarray(t["semiminor_sigma"], float) * scale,
        "pa": sky_pa(t["xcentroid"], t["ycentroid"], np.asarray(t["orientation"], float), wcs),
        "flux": f[keep] * 1e-9,
        "snr": f[keep] / e[keep],
        "ext": a > EXT_MIN_SIGMA_PIX,
        "ell": np.asarray(t["ellipticity"], float),
        "scale": scale,
        "x": np.asarray(t["xcentroid"], float),  # pixel positions; ``build`` pops them
        "y": np.asarray(t["ycentroid"], float),
    }


def keep_rows(f, e) -> np.ndarray:
    """Rule 4 row filter: positive finite Kron flux and error."""
    f, e = np.asarray(f, float), np.asarray(e, float)
    return np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)


def read_section(fs, path: str, ra, dec) -> tuple:
    """Rule 1: SCI and ERR sections (nJy per pixel) and the section WCS from one open file."""
    with fs.open(path, block_size=32 * 2**20) as fh, fits.open(fh) as hl:
        h = hl["SCI"].header
        if h.get("BUNIT") != "MJy/sr":
            raise ValueError(f"unexpected BUNIT {h.get('BUNIT')!r}")
        w = WCS(h)
        scale = float(proj_plane_pixel_scales(w.celestial).mean() * 3600.0)
        box = section_box(ra, dec, w, (h["NAXIS2"], h["NAXIS1"]), MARGIN_ARCSEC / scale)
        y0, y1, x0, x1 = box
        pix = float(h["PIXAR_SR"])
        sci = to_njy(hl["SCI"].section[y0:y1, x0:x1], pix)
        err = to_njy(hl["ERR"].section[y0:y1, x0:x1], pix)
    return sci, err, w[y0:y1, x0:x1], box


def _section():
    import s3fs

    _, _, obs, tile = mx.FIELDS["macs0417"]
    uri = mx.S3.format(obs, tile, "i2d.fits")
    pos = vstack(
        [Table.read(FRAME / f"{k}_macs0417.ecsv")["ra", "dec"] for k in ("control", "predicted")]
    )
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    return uri, read_section(fs, uri[5:], pos["ra"], pos["dec"])


def _write(out: Path, s: dict, keys) -> int:
    s = cp._clean(s)  # JSON has no NaN: non-finite -> null
    (out / "summary.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: s[k] for k in keys}, indent=1))
    return 0


def _pairs(ctrl: Table) -> list:
    return [
        (str(r), str(c)) for r, c in zip(ctrl["ref_counterpart"], ctrl["counterpart"], strict=True)
    ]


def stage_v() -> int:
    """Rules 1-6 (V1)."""
    uri, (sci, err, w, box) = _section()
    cat, meta, deb, _, _ = build(sci, err, w)
    OUT.mkdir(parents=True, exist_ok=True)
    cache = data_root() / "cache"  # segmentation and science section for the contact sheet
    cache.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache / "m0417_ownref_section.npz", segm=deb.data, sci=sci, box=box)
    depth = cr.depth_mag(cat["flux"], cat["snr"])
    inputs = {"image": uri, "section_y0_y1_x0_x1": [int(v) for v in box], **meta}
    s = ar.run(Path("."), cat=cat, depth=depth, cat_inputs=inputs, out=OUT)
    ctrl = Table.read(OUT / "control_macs0417.ecsv")
    shared = ar.gate_c2_shared(
        [str(c) for c in ctrl["class"]],
        [bool(a) for a in ctrl["arc_ref"]],
        [str(x) for x in ctrl["system"]],
        _pairs(ctrl),
    )
    s["gate_C2_shared_pairs"] = shared
    s["stage_V_passes"] = bool(shared["passes"])
    s["predicted_classes_note"] = "written, not read (PLAN.md stage V rule 6)"
    return _write(OUT, s, ("gate_C2", "gate_C2_shared_pairs", "stage_V_passes"))


def stage_v2(out_dir: Path = OUT_V2, bkg_box: int = BKG_BOX) -> int:
    """V2-a to V2-c: catalogue 5-sigma Kron depth, containing-source rule, local noise
    (stage 3c reruns it with ``bkg_box`` into its own ``out_dir``)."""
    uri, (sci, err, w, box) = _section()
    cat, meta, _, data, (sx, sy) = build(sci, err, w, bkg_box)
    out_dir.mkdir(parents=True, exist_ok=True)
    depth = kron_depth(cat)
    inputs = {"image": uri, "section_y0_y1_x0_x1": [int(v) for v in box], **meta}
    s = ar.run(Path("."), cat=cat, depth=depth, cat_inputs=inputs, out=out_dir, contain=True)
    rng = np.random.default_rng(m14.SEED)
    ctrl = Table.read(out_dir / "control_macs0417.ecsv")
    pred = Table.read(out_dir / "predicted_macs0417.ecsv")
    m14.apply_local(ctrl, data, w, sx, sy, rng, meta["pixel_scale_arcsec"])
    if len(pred):
        m14.apply_local(pred, data, w, sx, sy, rng, meta["pixel_scale_arcsec"])
    ctrl.write(out_dir / "control_macs0417.ecsv", overwrite=True)
    pred.write(out_dir / "predicted_macs0417.ecsv", overwrite=True)
    cl = [str(c) for c in ctrl["class_local"]]
    shared = ar.gate_c2_shared(
        cl, [bool(a) for a in ctrl["arc_ref"]], [str(x) for x in ctrl["system"]], _pairs(ctrl)
    )
    out = {
        "gate_C2_shared": shared,
        "stage_V2_passes": bool(shared["passes"]),
        "depth_mag_v2a": depth,
        "control_classes": m14._count(ctrl["class_local"]),
        "target_classes": m14._count(pred["class_local"]) if len(pred) else {},
        "predicted_classes_note": "written, not read (PLAN.md stage V2-c)",
        "without_local_noise_rule": {
            k: s[k] for k in ("gate_C2", "control_classes", "target_classes", "flags")
        },
        "posthoc_absent_controls": s["posthoc_absent_controls"],
        "inputs": s["inputs"],
        "rules": s["rules"],
    }
    return _write(out_dir, out, ("gate_C2_shared", "stage_V2_passes", "control_classes"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=["stage-v", "stage-v2"])
    a = ap.parse_args(argv)
    return stage_v() if a.stage == "stage-v" else stage_v2()


if __name__ == "__main__":
    raise SystemExit(main())
