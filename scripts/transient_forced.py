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
and then PSF-wing differences between epochs fake a change (Earendel, docs/fields/sunrise.md).
``--recentre-arcsec`` moves the aperture to the source centroid in each band, found in the epoch
where the source exists (epoch 2 for ``appeared`` candidates, else epoch 1), and measures both
epochs at that same sky position.

ERR-based errors underestimate the noise of these measurements (1.2–1.5x in Sunrise, D-027).
``--controls`` measures ordinary sources (a random subset of a catalog, e.g. the epoch-1
pipeline ``_cat.ecsv``, within ``--control-mag``) the same way. Per band, the robust std
(1.4826 MAD) of their significances is the noise scale. Candidate significances are divided by it
(never by less than 1) before thresholding, and the raw values are kept as ``*_sigma_raw``.

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


def recentre(
    sci: np.ndarray, x: float, y: float, half_px: int, max_iter: int = 10
) -> tuple[float, float]:
    """Iterated intensity-weighted centroid near ``(x, y)``.

    Each pass takes the positive pixels of a ``(2 half_px + 1)`` box around the current position,
    minus the median of a local ring (``2–3 half_px``), and moves the box to their centroid until
    it moves by < 0.05 px. The input position is returned when it lies outside the image, when a
    box is empty, or when the centroid drifts more than ``2 half_px`` from the input (the box then
    followed a neighbour, not the target)."""
    ny, nx = sci.shape
    if not (np.isfinite(x) and np.isfinite(y)) or not (0 <= x < nx and 0 <= y < ny):
        return x, y
    yy, xx = np.mgrid[0:ny, 0:nx]
    cx, cy = x, y
    for _ in range(max_iter):
        r = np.maximum(np.abs(xx - cx), np.abs(yy - cy))
        box = (r <= half_px) & np.isfinite(sci)
        ring = (r >= 2 * half_px) & (r <= 3 * half_px) & np.isfinite(sci)
        if not box.any():
            return x, y
        bkg = float(np.median(sci[ring])) if ring.sum() >= 8 else float(np.nanmedian(sci))
        w = np.clip(sci[box] - bkg, 0, None)
        if w.sum() <= 0:
            return x, y
        nxc, nyc = float((w * xx[box]).sum() / w.sum()), float((w * yy[box]).sum() / w.sum())
        done = np.hypot(nxc - cx, nyc - cy) < 0.05
        cx, cy = nxc, nyc
        if done:
            break
    if np.hypot(cx - x, cy - y) > 2 * half_px:
        return x, y
    return cx, cy


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
    aperture moves to the centroid found from that half-width first (``recentre``). Also returns
    the RA/Dec used."""
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
            cx, cy = recentre(sci, cx, cy, max(1, int(np.ceil(recentre_arcsec / scale_arcsec))))
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


def robust_std(x: np.ndarray) -> float:
    """1.4826 x median absolute deviation of the finite values (NaN when fewer than 10)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 10:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def select_controls(
    cat: Table,
    n: int,
    mag_range: tuple[float, float],
    avoid_ra: np.ndarray,
    avoid_dec: np.ndarray,
    avoid_arcsec: float = 1.0,
    seed: int = 0,
) -> Table:
    """A reproducible random subset of ``n`` catalog sources (``ra``/``dec`` or the pipeline's
    ``sky_centroid``) as noise controls.

    Uses ``aper_total_abmag`` within ``mag_range`` when the column exists, and drops sources within
    ``avoid_arcsec`` of any candidate, so a control never measures a candidate."""
    if "ra" in cat.colnames:
        ra, dec = np.asarray(cat["ra"], float), np.asarray(cat["dec"], float)
    else:  # a raw pipeline catalog
        ra, dec = cat["sky_centroid"].ra.deg, cat["sky_centroid"].dec.deg
    keep = np.isfinite(ra) & np.isfinite(dec)
    if "aper_total_abmag" in cat.colnames:
        m = np.asarray(np.ma.filled(cat["aper_total_abmag"], np.nan), float)
        with np.errstate(invalid="ignore"):
            keep &= (m >= mag_range[0]) & (m <= mag_range[1])
    if len(avoid_ra):
        from astropy import units as u
        from astropy.coordinates import SkyCoord

        c = SkyCoord(ra * u.deg, dec * u.deg)
        _, d, _ = c.match_to_catalog_sky(SkyCoord(avoid_ra * u.deg, avoid_dec * u.deg))
        keep &= d.arcsec > avoid_arcsec
    idx = np.flatnonzero(keep)
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(idx, size=min(n, len(idx)), replace=False))
    return Table({"ra": ra[idx], "dec": dec[idx]})


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
        help="move the aperture to the source centroid (box half-width) per band, found in epoch 1 "
        "(epoch 2 for 'appeared' candidates), and measure both epochs there",
    )
    ap.add_argument(
        "--controls",
        type=Path,
        default=None,
        help="ordinary sources (ra, dec; e.g. the epoch-1 _cat.ecsv) for noise calibration",
    )
    ap.add_argument("--n-controls", type=int, default=200)
    ap.add_argument(
        "--control-mag", nargs=2, type=float, default=(25.5, 28.0), metavar=("BRIGHT", "FAINT")
    )
    args = ap.parse_args(argv)

    cfg = pipeline.load_config(args.config)
    cand = Table.read(args.candidates)
    n_cand = len(cand)
    ra_all, dec_all = np.asarray(cand["ra"], float), np.asarray(cand["dec"], float)
    uids = [f"c{k:04d}" for k in range(n_cand)]
    kinds = np.asarray(cand["kind"]).astype(str) if "kind" in cand.colnames else np.full(n_cand, "")
    if args.controls is not None:
        ctl = select_controls(
            Table.read(args.controls), args.n_controls, tuple(args.control_mag), ra_all, dec_all
        )
        ra_all = np.concatenate([ra_all, np.asarray(ctl["ra"])])
        dec_all = np.concatenate([dec_all, np.asarray(ctl["dec"])])
        uids += [f"n{k:04d}" for k in range(len(ctl))]
        kinds = np.concatenate([kinds, np.full(len(ctl), "control")])
    is_control = np.arange(len(uids)) >= n_cand
    targets = Table({"source_uid": uids, "ra": ra_all, "dec": dec_all})
    out = Table({"source_uid": targets["source_uid"], "ra": targets["ra"], "dec": targets["dec"]})
    out["control"] = is_control
    confirmed = ~is_control
    noise = {}
    # appeared sources exist only in epoch 2, so they are centroided there; others in epoch 1
    from_e2 = np.char.startswith(kinds, "appeared")
    for band, obs1, obs2 in args.band:
        files, scales = [], []
        for epoch, obs in (("1", obs1), ("2", obs2)):
            uri = pipeline.l3_image_uri(cfg.get("cloud"), obs)
            t = cutouts.make_cutouts(
                uri,
                targets,
                size_arcsec=2.0,
                extensions=["ERR"],
                out_dir=args.out / "cutouts" / f"{band}_e{epoch}",
            )
            scales.append(float(np.mean(t.meta["pixel_scale_arcsec"])))
            paths = {str(u): str(p) for u, p in zip(t["source_uid"], t["path"], strict=True)}
            files.append([paths.get(str(u), "") for u in targets["source_uid"]])
        ra_b, dec_b = np.asarray(targets["ra"]), np.asarray(targets["dec"])
        if args.recentre_arcsec > 0:
            cen = [
                measure(files[k], args.radius_arcsec, scales[k], ra_b, dec_b, args.recentre_arcsec)[
                    2:
                ]
                for k in (0, 1)
            ]
            ra_c = np.where(from_e2, cen[1][0], cen[0][0])
            dec_c = np.where(from_e2, cen[1][1], cen[0][1])
            ok = np.isfinite(ra_c) & np.isfinite(dec_c)
            ra_b, dec_b = np.where(ok, ra_c, ra_b), np.where(ok, dec_c, dec_b)
            out[f"{band}_ra"], out[f"{band}_dec"] = ra_b, dec_b
        fl, er = [], []
        for k in (0, 1):  # both epochs at the same sky position
            f, e, _, _ = measure(files[k], args.radius_arcsec, scales[k], ra_b, dec_b)
            fl.append(f)
            er.append(e)
        dm, sig, sig_flux = compare(fl[0], er[0], fl[1], er[1], args.sys_floor)
        out[f"{band}_flux1"], out[f"{band}_flux2"] = fl[0], fl[1]
        out[f"{band}_dmag"] = dm
        if is_control.any():
            s_dm, s_fl = robust_std(sig[is_control]), robust_std(sig_flux[is_control])
            noise[band] = {
                "sigma": s_dm,
                "flux_sigma": s_fl,
                "n_controls": int(np.isfinite(sig_flux[is_control]).sum()),
            }
            out[f"{band}_sigma_raw"], out[f"{band}_flux_sigma_raw"] = sig, sig_flux
            # never sharpen: a scale below 1 (or undefined) leaves the ERR-based value
            sig = sig / (s_dm if s_dm > 1 else 1.0)
            sig_flux = sig_flux / (s_fl if s_fl > 1 else 1.0)
        out[f"{band}_sigma"], out[f"{band}_flux_sigma"] = sig, sig_flux
        with np.errstate(invalid="ignore"):
            changed = (np.abs(dm) >= args.min_dmag) & (np.abs(sig) >= args.min_sigma)
            # a non-detection in one epoch (flux <= 0) is judged in flux space
            gone = (np.fmin(fl[0], fl[1]) <= 0) & (np.abs(sig_flux) >= args.min_sigma)
            confirmed &= changed | gone
    out["confirmed"] = confirmed
    for c in cand.colnames:
        if c not in ("ra", "dec") and c not in out.colnames:
            if is_control.any():  # controls get masked values in candidate-only columns
                pad = np.ma.masked_all(int(is_control.sum()), dtype=cand[c].dtype)
                out[c] = np.ma.concatenate([np.ma.asarray(cand[c]), pad])
            else:
                out[c] = cand[c]
    if is_control.any():
        out["kind"] = kinds
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
    if noise:
        out.meta["noise_scale"] = {
            **noise,
            "controls": str(args.controls),
            "control_mag": list(args.control_mag),
            "provenance": "derived",
        }
    args.out.mkdir(parents=True, exist_ok=True)
    out.write(args.out / "forced.ecsv", overwrite=True)
    print(
        json.dumps(
            {"n": int((~is_control).sum()), "confirmed": int(confirmed.sum()), "noise": noise}
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
