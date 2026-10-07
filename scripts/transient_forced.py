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


def measure(
    paths: list[str],
    r_arcsec: float,
    scale_arcsec: float,
    ra: np.ndarray | None = None,
    dec: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Aperture fluxes and errors at each target's sky position (cutout WCS), or at the cutout
    centre when no position is given (NaN for missing files)."""
    flux, err = np.full(len(paths), np.nan), np.full(len(paths), np.nan)
    r_px = r_arcsec / scale_arcsec
    for k, p in enumerate(paths):
        if not p:
            continue
        with fits.open(p) as h:
            sci = np.asarray(h["SCI"].data, float)
            e = np.asarray(h["ERR"].data, float) if "ERR" in h else None
            header = h["SCI"].header
        cy, cx = (sci.shape[0] - 1) / 2, (sci.shape[1] - 1) / 2
        if ra is not None and dec is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FITSFixedWarning)
                cx, cy = (float(v) for v in WCS(header).world_to_pixel_values(ra[k], dec[k]))
        flux[k], err[k] = aperture_flux(sci, e, cx, cy, r_px, 2.5 * r_px, 4.0 * r_px)
    return flux, err


def compare(f1, e1, f2, e2, sys_floor: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """m2 − m1 with the median zero point removed, and its significance. ``sys_floor`` (mag) is
    added in quadrature: PSF and registration differences between epochs dominate the ERR-based
    error for bright sources."""
    with np.errstate(divide="ignore", invalid="ignore"):
        dm = -2.5 * np.log10(f2 / f1)
        zp = np.nanmedian(dm)
        for _ in range(3):  # sigma-clipped zero point: real variables do not drag it
            mad = 1.4826 * np.nanmedian(np.abs(dm - zp))
            keep = np.abs(dm - zp) <= 3 * max(mad, 1e-3)
            if keep.any():
                zp = np.nanmedian(dm[keep])
        sig_dm = np.hypot((2.5 / np.log(10)) * np.hypot(e1 / f1, e2 / f2), sys_floor)
        dm = dm - zp
        return dm, dm / sig_dm


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
            f, e = measure(
                [paths.get(str(u), "") for u in targets["source_uid"]],
                args.radius_arcsec,
                scale,
                np.asarray(targets["ra"]),
                np.asarray(targets["dec"]),
            )
            fl.append(f)
            er.append(e)
        dm, sig = compare(fl[0], er[0], fl[1], er[1], args.sys_floor)
        out[f"{band}_flux1"], out[f"{band}_flux2"] = fl[0], fl[1]
        out[f"{band}_dmag"], out[f"{band}_sigma"] = dm, sig
        with np.errstate(invalid="ignore"):
            confirmed &= (np.abs(dm) >= args.min_dmag) & (np.abs(sig) >= args.min_sigma)
    out["confirmed"] = confirmed
    for c in cand.colnames:
        if c not in ("ra", "dec") and c not in out.colnames:
            out[c] = cand[c]
    out.meta.update(
        provenance="derived",
        radius_arcsec=args.radius_arcsec,
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
