"""Forced two-epoch aperture photometry at fixed sky positions (verifies transient candidates).

Catalog-level comparisons across epochs inherit each catalog's segmentation and deblending, and
those changed between jwst pipeline versions. This script instead takes byte-range cutouts of
both epochs' level-3 ``_i2d.fits`` images (S3; no full downloads) at each candidate position. It
measures the flux in the same circular aperture, minus a local background annulus (median), with
errors from the ``ERR`` extension. Every result is ``derived``.

``dmag`` is m(epoch2) − m(epoch1), with the zero point set by a 3-sigma-clipped median over all
positions. That assumes most candidates are not variable, which holds for catalog-level
candidates (mostly deblending differences). The aperture sits at each candidate's RA/Dec through
the cutout WCS, not at the cutout's centre pixel. A candidate is ``confirmed`` when
|dmag| ≥ ``--min-dmag`` at ≥ ``--min-sigma`` in every band given.

Positions from another frame or catalog (e.g. a literature position) can sit 0.1″ off the source,
and then PSF-wing differences between epochs fake a change (Earendel, 2026-10-08: −0.75 mag at
the literature position, +0.01 mag at the centroid). ``--recentre-arcsec`` moves the aperture to
the epoch-1 centroid in each band and measures epoch 2 at that same sky position.

    python scripts/transient_forced.py --candidates outputs/transients/coincident.ecsv \\
        --band F444W jw02736-o001_t001_nircam_clear-f444w jw06882-o057_t057_nircam_clear-f444w \\
        --band F150W jw02736-o001_t001_nircam_clear-f150w jw06882-o057_t057_nircam_clear-f150w \\
        --config configs/reference_sample.yaml --out outputs/transients/forced
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import cutouts, pipeline


def aperture_flux(
    sci: np.ndarray,
    err: np.ndarray | None,
    x: float,
    y: float,
    r_px: float,
    r_in: float,
    r_out: float,
) -> tuple[float, float]:
    """Background-subtracted sum within ``r_px`` of ``(x, y)`` and its error (same units as
    ``sci`` times pixels). NaN when the aperture has no-data pixels."""
    ny, nx = sci.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x, yy - y)
    ap = r <= r_px
    ann = (r >= r_in) & (r <= r_out) & np.isfinite(sci)
    if not ap.any() or not np.isfinite(sci[ap]).all() or ann.sum() < 10:
        return float("nan"), float("nan")
    bkg = float(np.median(sci[ann]))
    flux = float(np.sum(sci[ap] - bkg))
    if err is None or not np.isfinite(err[ap]).all():
        e = float(np.std(sci[ann]) * np.sqrt(ap.sum()))
    else:
        e = float(np.sqrt(np.sum(err[ap] ** 2)))
    return flux, e


def recentre(sci: np.ndarray, x: float, y: float, half_px: int) -> tuple[float, float]:
    """Intensity-weighted centroid of the positive, background-subtracted pixels within a
    ``(2 half_px + 1)`` box around ``(x, y)``; the input position when the box is empty."""
    if not (np.isfinite(x) and np.isfinite(y)):
        return x, y
    xi, yi = int(round(x)), int(round(y))
    y0, x0 = max(yi - half_px, 0), max(xi - half_px, 0)
    sub = sci[y0 : yi + half_px + 1, x0 : xi + half_px + 1]
    if sub.size == 0 or not np.isfinite(sub).any():
        return x, y
    w = np.clip(np.nan_to_num(sub - np.nanmedian(sci)), 0, None)
    if w.sum() <= 0:
        return x, y
    yy, xx = np.mgrid[0 : sub.shape[0], 0 : sub.shape[1]]
    return float(x0 + (w * xx).sum() / w.sum()), float(y0 + (w * yy).sum() / w.sum())


def measure(
    paths: list[str],
    r_arcsec: float,
    scale_arcsec: float,
    ra: np.ndarray | None = None,
    dec: np.ndarray | None = None,
    recentre_arcsec: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Aperture fluxes and errors at each target's sky position (cutout WCS), or at the cutout
    centre when no position is given (NaN for missing files). With ``recentre_arcsec`` > 0 the
    aperture moves to the centroid within that half-width first. Also returns the RA/Dec used."""
    flux, err = np.full(len(paths), np.nan), np.full(len(paths), np.nan)
    ra_used, dec_used = np.full(len(paths), np.nan), np.full(len(paths), np.nan)
    r_px = r_arcsec / scale_arcsec
    for k, p in enumerate(paths):
        if not p:
            continue
        with fits.open(p) as h:
            sci = np.asarray(h["SCI"].data, float)
            e = np.asarray(h["ERR"].data, float) if "ERR" in h else None
            header = h["SCI"].header
        cy, cx = (sci.shape[0] - 1) / 2, (sci.shape[1] - 1) / 2
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(header)
        if ra is not None and dec is not None:
            cx, cy = (float(v) for v in wcs.world_to_pixel_values(ra[k], dec[k]))
        if recentre_arcsec > 0:
            cx, cy = recentre(sci, cx, cy, max(1, int(round(recentre_arcsec / scale_arcsec))))
        ra_used[k], dec_used[k] = (float(v) for v in wcs.pixel_to_world_values(cx, cy))
        flux[k], err[k] = aperture_flux(sci, e, cx, cy, r_px, 2.5 * r_px, 4.0 * r_px)
    return flux, err, ra_used, dec_used


def compare(
    f1, e1, f2, e2, sys_floor: float = 0.05, min_zp_refs: int = 10
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(dmag, dmag significance, flux-difference significance)`` between two epochs.

    The zero point is a 3-sigma-clipped median of ``dmag`` over all positions with positive flux in
    both epochs, used only with at least ``min_zp_refs`` such positions (otherwise 0: both epochs
    are calibrated MJy/sr). ``sys_floor`` (mag) is added in quadrature: PSF and registration
    differences dominate the ERR-based error for bright sources. The flux-space significance,
    ``(f2 - f1 * 10**(-0.4 zp)) / err``, also works for non-detections (flux <= 0), where a
    magnitude is undefined: that is the test for ``appeared``/``disappeared`` sources.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        both = (f1 > 0) & (f2 > 0)
        dm = np.where(both, -2.5 * np.log10(f2 / f1), np.nan)
        zp = 0.0
        if np.isfinite(dm).sum() >= min_zp_refs:
            zp = float(np.nanmedian(dm))
            for _ in range(3):  # sigma-clipped: real variables do not drag the zero point
                mad = 1.4826 * np.nanmedian(np.abs(dm - zp))
                keep = np.abs(dm - zp) <= 3 * max(mad, 1e-3)
                if keep.any():
                    zp = float(np.nanmedian(dm[keep]))
        sig_dm = np.hypot((2.5 / np.log(10)) * np.hypot(e1 / f1, e2 / f2), sys_floor)
        dm = dm - zp
        f1_scaled = f1 * 10 ** (-0.4 * zp)
        floor = sys_floor * np.log(10) / 2.5 * np.fmax(np.abs(f1_scaled), np.abs(f2))
        sig_flux = (f2 - f1_scaled) / np.sqrt(e1**2 + e2**2 + floor**2)
        return dm, dm / sig_dm, sig_flux


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--candidates", type=Path, required=True, help="table with ra, dec columns")
    ap.add_argument(
        "--band",
        nargs=3,
        action="append",
        required=True,
        metavar=("BAND", "OBS_ID_EPOCH1", "OBS_ID_EPOCH2"),
    )
    ap.add_argument("--config", type=Path, required=True, help="for the cloud key pattern")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--radius-arcsec", type=float, default=0.15)
    ap.add_argument("--min-dmag", type=float, default=0.3)
    ap.add_argument("--min-sigma", type=float, default=5.0)
    ap.add_argument("--sys-floor", type=float, default=0.05, help="mag, added in quadrature")
    ap.add_argument(
        "--recentre-arcsec",
        type=float,
        default=0.0,
        help="move the aperture to the epoch-1 centroid within this half-width, per band, and use "
        "that sky position in epoch 2 too (for positions from another frame or catalog)",
    )
    args = ap.parse_args(argv)

    cfg = pipeline.load_config(args.config)
    cand = Table.read(args.candidates)
    targets = Table(
        {
            "source_uid": [f"c{k:04d}" for k in range(len(cand))],
            "ra": np.asarray(cand["ra"], float),
            "dec": np.asarray(cand["dec"], float),
        }
    )
    out = Table({"source_uid": targets["source_uid"], "ra": targets["ra"], "dec": targets["dec"]})
    confirmed = np.ones(len(targets), bool)
    for band, obs1, obs2 in args.band:
        fl, er = [], []
        ra_b, dec_b = np.asarray(targets["ra"]), np.asarray(targets["dec"])
        for epoch, obs in (("1", obs1), ("2", obs2)):
            uri = pipeline.l3_image_uri(cfg.get("cloud"), obs)
            t = cutouts.make_cutouts(
                uri,
                targets,
                size_arcsec=2.0,
                extensions=["ERR"],
                out_dir=args.out / "cutouts" / f"{band}_e{epoch}",
            )
            scale = float(np.mean(t.meta["pixel_scale_arcsec"]))
            paths = {str(u): str(p) for u, p in zip(t["source_uid"], t["path"], strict=True)}
            f, e, ra_used, dec_used = measure(
                [paths.get(str(u), "") for u in targets["source_uid"]],
                args.radius_arcsec,
                scale,
                ra_b,
                dec_b,
                args.recentre_arcsec if epoch == "1" else 0.0,
            )
            if epoch == "1" and args.recentre_arcsec > 0:
                ok = np.isfinite(ra_used)
                ra_b, dec_b = np.where(ok, ra_used, ra_b), np.where(ok, dec_used, dec_b)
                out[f"{band}_ra"], out[f"{band}_dec"] = ra_b, dec_b
            fl.append(f)
            er.append(e)
        dm, sig, sig_flux = compare(fl[0], er[0], fl[1], er[1], args.sys_floor)
        out[f"{band}_flux1"], out[f"{band}_flux2"] = fl[0], fl[1]
        out[f"{band}_dmag"], out[f"{band}_sigma"] = dm, sig
        out[f"{band}_flux_sigma"] = sig_flux
        with np.errstate(invalid="ignore"):
            changed = (np.abs(dm) >= args.min_dmag) & (np.abs(sig) >= args.min_sigma)
            # a non-detection in one epoch (flux <= 0) is judged in flux space
            gone = (np.fmin(fl[0], fl[1]) <= 0) & (np.abs(sig_flux) >= args.min_sigma)
            confirmed &= changed | gone
    out["confirmed"] = confirmed
    for c in cand.colnames:
        if c not in ("ra", "dec") and c not in out.colnames:
            out[c] = cand[c]
    out.meta.update(
        provenance="derived",
        radius_arcsec=args.radius_arcsec,
        recentre_arcsec=args.recentre_arcsec,
        thresholds={
            "min_dmag": args.min_dmag,
            "min_sigma": args.min_sigma,
            "provenance": "assumption",
        },
        bands=[list(b) for b in args.band],
    )
    args.out.mkdir(parents=True, exist_ok=True)
    out.write(args.out / "forced.ecsv", overwrite=True)
    print(json.dumps({"n": len(out), "confirmed": int(confirmed.sum())}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
