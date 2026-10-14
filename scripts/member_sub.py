"""Member-galaxy subtraction for counter-image rows buried in cluster-member light (D-029 revisit).

``lens_consistency.py images --forced-image`` leaves predicted images ``confused`` when a member
galaxy within the search radius is many times brighter than the prediction; a median filter does
not remove a compact member core (D-029 addendum 2026-10-14). This script re-measures such rows
after subtracting parametric models of the bright neighbours:

1. neighbours: catalogue sources within the row's search radius plus ``NB_MARGIN_ARCSEC`` whose
   0.2" aperture flux in the forced image exceeds the predicted image's flux (fainter ones cannot
   confuse it), brightest ``MAX_COMPONENTS`` first;
2. fit: one elliptical Sersic per neighbour plus a constant (``astropy.modeling``, weighted by
   ERR, not PSF-convolved) with the whole search disc masked except each neighbour's 0.3" core,
   so the fit cannot absorb the image it is testing;
3. measure: the best 0.2" aperture flux on the residual within the disc (apertures clear of the
   cores), against an empirical noise: the scatter of residual aperture fluxes on a ring around
   the nearest neighbour at the target's distance (core residuals of a model dominate there);
4. calibrate: the same fit with a compact image of the predicted flux injected at the position
   and four points at half the radius; the lowest recovered fraction scales the prediction
   (simulated; ~1 by construction while the disc is masked). Controls: the same masked fit in
   discs rotated about the nearest neighbour, where no image is predicted (``CONTROL_ANGLES``).

Classes (thresholds are ASSUMPTIONs, as ``lens_consistency.forced_class``): ``no_member`` (nothing
to subtract), ``blend`` (a neighbour core within 0.5" of the position), ``untestable`` (recovery
< ``MIN_RECOVERY``, calibrated prediction < 10 sigma, a control disc above a third of it, or an
absence without a control disc), then ``absent`` / ``recovered`` / ``confused`` / ``ambiguous``
from the residual flux over the calibrated prediction. Positions and fluxes of the prediction are
model_prediction; residual fluxes derived; the injection is simulated.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.modeling import fitting, models
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

NB_MARGIN_ARCSEC = 0.5  # neighbours within search + this are modelled; the fit uses pixels
# within search + 2x this
MAX_COMPONENTS = 4  # Sersic components per stamp (brightest neighbours first)
CORE_CLEARANCE_ARCSEC = 0.3  # neighbour cores (kept in the fit; test apertures stay clear)
STAMP_HALF_ARCSEC = 7.0  # room for control discs rotated about a neighbour up to 3.5" away
CONTROL_ANGLES = (90.0, 180.0, 270.0)  # control discs, rotated about the nearest neighbour
RING_N = 24  # empirical-noise apertures around the nearest neighbour
RING_EXCLUDE_ARCSEC = 0.6  # ring apertures this close to the target are skipped
MIN_RING = 8  # fewer valid ring apertures: no empirical noise, untestable
INJECT_FWHM_ARCSEC = 0.1  # Gaussian stand-in for a compact image at 2-3 um (~F277W PSF)
MIN_RECOVERY = 0.5  # the fit must keep at least half of an injected image
MIN_PRED_SNR = 10.0  # calibrated prediction over the empirical noise
SERSIC_N = (0.5, 8.0)
SERSIC_REFF_ARCSEC = (0.02, 6.0)


def neighbour_offsets(cat: Table, ra0: float, dec0: float) -> tuple[np.ndarray, ...]:
    """Catalogue offsets (West, North arcsec, as the stamp's ``xx``/``yy``) and shapes near a
    position: ``(dx, dy, ellip, theta_rad)``."""
    dx = -(np.asarray(cat["ra"], float) - ra0) * np.cos(np.radians(dec0)) * 3600.0
    dy = (np.asarray(cat["dec"], float) - dec0) * 3600.0
    ell = np.clip(np.nan_to_num(np.asarray(cat["ellipticity"], float)), 0.0, 0.8)
    # sky_orientation is a PA east of north; Sersic2D's theta runs from +x (West) to +y (North)
    theta = np.radians(90.0 + np.nan_to_num(np.asarray(cat["pa_obs"], float)))
    return dx, dy, ell, theta


def _sersic_fit(sci, err, xx, yy, comps, fit_mask):
    """Least-squares fit of ``comps`` (list of (x0, y0, amp, ellip, theta)) Sersics plus a
    constant; returns the model evaluated on the whole stamp (NaN where ``sci`` is NaN)."""
    model = models.Const2D(amplitude=float(np.nanmedian(sci[fit_mask])))
    for x0, y0, amp, ell, th in comps:
        s = models.Sersic2D(
            amplitude=max(amp * np.exp(-(2 * 2.5 - 1 / 3)), 1e-9),  # peak -> I(r_eff), n = 2.5
            r_eff=0.3,
            n=2.5,
            x_0=x0,
            y_0=y0,
            ellip=min(ell, 0.7),
            theta=th,
            bounds={
                "n": SERSIC_N,
                "r_eff": SERSIC_REFF_ARCSEC,
                "ellip": (0.0, 0.9),
                "x_0": (x0 - 0.3, x0 + 0.3),
                "y_0": (y0 - 0.3, y0 + 0.3),
                "amplitude": (0.0, None),
            },
        )
        model = model + s
    w = 1.0 / np.where(np.isfinite(err) & (err > 0), err, np.inf)
    good = fit_mask & np.isfinite(sci) & (w > 0)
    fitter = fitting.TRFLSQFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = fitter(model, xx[good], yy[good], sci[good], weights=w[good], maxiter=400)
    out = fit(xx, yy)
    out[~np.isfinite(sci)] = np.nan
    return out


def _gaussian(xx, yy, flux, fwhm=INJECT_FWHM_ARCSEC):
    s = fwhm / 2.3548
    g = np.exp(-(xx**2 + yy**2) / (2 * s**2))
    return g * flux / g.sum()


def _testable(cx, cy, radius, centres) -> bool:
    """A test centre lies within ``radius`` and its aperture clears every neighbour's core."""
    clear = CORE_CLEARANCE_ARCSEC + lc.FORCED_R_AP
    return np.hypot(cx, cy) <= radius and all(np.hypot(cx - x, cy - y) > clear for x, y in centres)


def _best_in_disc(resid, err, xx, yy, radius, centres, step=0.05):
    """Highest residual aperture flux on a grid of testable centres: (flux, dx, dy)."""
    best = (-np.inf, np.nan, np.nan)
    for cx in np.arange(-radius, radius + 1e-9, step):
        for cy in np.arange(-radius, radius + 1e-9, step):
            if not _testable(cx, cy, radius, centres):
                continue
            f, _ = lc.aperture_snr(resid, err, xx, yy, float(cx), float(cy))
            if np.isfinite(f) and f > best[0]:
                best = (f, float(cx), float(cy))
    return best


def _inject_points(radius, centres):
    """Injection positions: the predicted position and four at half the test radius."""
    pts = [(0.0, 0.0)] + [
        (0.5 * radius * np.cos(a), 0.5 * radius * np.sin(a)) for a in np.arange(4) * np.pi / 2
    ]
    return [p for p in pts if _testable(p[0], p[1], radius, centres)]


def _recovery(sci, err, xx, yy, resid, comps, fit_mask, flux, px, py) -> float:
    """Fraction of a compact image of ``flux`` injected at (px, py) left in the residual."""
    if not flux > 0:
        return float("nan")
    inj = sci + _gaussian(xx - px, yy - py, flux)
    resid_inj = inj - _sersic_fit(inj, err, xx, yy, comps, fit_mask)
    f0, _ = lc.aperture_snr(resid, err, xx, yy, px, py)
    f1, _ = lc.aperture_snr(resid_inj, err, xx, yy, px, py)
    return (f1 - f0) / flux


def _control_fluxes(sci, err, xx, yy, comps, mask_r, near) -> list[float]:
    """Best residual flux in discs like the target's at ``CONTROL_ANGLES`` about the nearest
    neighbour's centre (no image predicted there): what the member model leaves under a mask.
    A control disc must lie inside the stamp."""
    half = float(min(np.nanmax(np.abs(xx)), np.nanmax(np.abs(yy))))
    d, a0 = np.hypot(near[0], near[1]), np.arctan2(-near[1], -near[0])
    out = []
    for da in CONTROL_ANGLES:
        a = a0 + np.radians(da)
        cx, cy = near[0] + d * np.cos(a), near[1] + d * np.sin(a)
        if max(abs(cx), abs(cy)) + mask_r + lc.FORCED_ANNULUS[1] > half:
            continue
        r = np.hypot(xx - cx, yy - cy)
        core = np.zeros_like(r, dtype=bool)
        for c in comps:
            core |= np.hypot(xx - c[0], yy - c[1]) <= CORE_CLEARANCE_ARCSEC
        fit_mask = ((r > mask_r) | core) & (r <= mask_r + 2 * NB_MARGIN_ARCSEC)
        resid = sci - _sersic_fit(sci, err, xx, yy, comps, fit_mask)
        centres = [(c[0] - cx, c[1] - cy) for c in comps]
        f, _, _ = _best_in_disc(resid, err, xx - cx, yy - cy, mask_r, centres)
        if np.isfinite(f):
            out.append(float(f))
    return out


def _ring_sigma(resid, err, xx, yy, cx, cy, d):
    """Robust scatter (1.4826 MAD) of residual aperture fluxes on a ring of radius ``d`` around
    (cx, cy), skipping apertures within ``RING_EXCLUDE_ARCSEC`` of the target."""
    vals = []
    for k in range(RING_N):
        a = 2 * np.pi * k / RING_N
        px, py = cx + d * np.cos(a), cy + d * np.sin(a)
        if np.hypot(px, py) < RING_EXCLUDE_ARCSEC:
            continue
        f, _ = lc.aperture_snr(resid, err, xx, yy, px, py)
        if np.isfinite(f):
            vals.append(f)
    if len(vals) < MIN_RING:
        return float("nan"), len(vals)
    v = np.asarray(vals)
    return float(1.4826 * np.median(np.abs(v - np.median(v)))), len(vals)


def member_class(ratio: float, snr: float) -> str:
    """Verdict from residual flux over calibrated prediction and empirical S/N (ASSUMPTIONs)."""
    if not np.isfinite(ratio):
        return "ambiguous"
    if ratio < lc.MIN_FLUX_RATIO:
        return "absent"
    if snr >= 5 and ratio > lc.MAX_FLUX_RATIO:
        return "confused"
    if snr >= 5:
        return "recovered"
    return "ambiguous"


def measure(sci, err, xx, yy, pred_flux: float, search: float, nb: list[tuple]) -> dict:
    """Member-subtracted verdict for one predicted image at the stamp centre.

    ``nb`` lists the candidate neighbours as (dx, dy, ellip, theta) in stamp offsets."""
    out = {"n_components": 0, "member_class": "no_member"}
    comps = []
    for dx, dy, ell, th in nb:
        if np.hypot(dx, dy) > search + NB_MARGIN_ARCSEC:
            continue
        f, _ = lc.aperture_snr(sci, err, xx, yy, dx, dy)
        if np.isfinite(f) and f > pred_flux:
            near = np.hypot(xx - dx, yy - dy) <= 0.1
            amp = float(np.nanmax(sci[near])) if near.any() else f
            comps.append((f, (dx, dy, amp, ell, th)))
    if not comps:
        return out
    comps = [c for _, c in sorted(comps, key=lambda t: -t[0])[:MAX_COMPONENTS]]
    d_core = min(np.hypot(c[0], c[1]) for c in comps)
    out.update(n_components=len(comps), d_core_arcsec=round(float(d_core), 3))
    if d_core < CORE_CLEARANCE_ARCSEC + lc.FORCED_R_AP:
        out["member_class"] = "blend"
        return out
    # the whole search disc is masked, except each neighbour's core, which pins its model
    mask_r = search
    r = np.hypot(xx, yy)
    core = np.zeros_like(r, dtype=bool)
    for c in comps:
        core |= np.hypot(xx - c[0], yy - c[1]) <= CORE_CLEARANCE_ARCSEC
    fit_mask = ((r > mask_r) | core) & (r <= mask_r + 2 * NB_MARGIN_ARCSEC)
    resid = sci - _sersic_fit(sci, err, xx, yy, comps, fit_mask)
    near = min(comps, key=lambda c: np.hypot(c[0], c[1]))
    sigma, n_ring = _ring_sigma(resid, err, xx, yy, near[0], near[1], d_core)
    _, e0 = lc.aperture_snr(sci, err, xx, yy, 0.0, 0.0)
    sigma = max(sigma, e0) if np.isfinite(sigma) else float("nan")
    centres = [c[:2] for c in comps]
    recs = []
    for px, py in _inject_points(mask_r, centres):
        recs.append(_recovery(sci, err, xx, yy, resid, comps, fit_mask, pred_flux, px, py))
    recs = [v for v in recs if np.isfinite(v)]
    rec = float(np.min(recs)) if recs else float("nan")
    best, bx, by = _best_in_disc(resid, err, xx, yy, mask_r, centres)
    f0, _ = lc.aperture_snr(resid, err, xx, yy, 0.0, 0.0)
    calib = pred_flux * rec
    ctrl = _control_fluxes(sci, err, xx, yy, comps, mask_r, near)
    out.update(
        mask_arcsec=round(float(mask_r), 3),
        sigma_emp=float(sigma),
        n_ring=n_ring,
        recovery=float(rec),
        n_inject=len(recs),
        calib_pred_snr=float(calib / sigma) if sigma > 0 else float("nan"),
        resid_flux_at_position=float(f0),
        best_resid_flux=float(best),
        best_dx=bx,
        best_dy=by,
        control_best_fluxes=[round(v, 4) for v in ctrl],
    )
    if not (np.isfinite(rec) and rec >= MIN_RECOVERY and np.isfinite(sigma) and sigma > 0):
        out["member_class"] = "untestable"
        return out
    if calib / sigma < MIN_PRED_SNR:
        out["member_class"] = "untestable"
        return out
    # a control disc with residual flux above a third of the prediction: the member model (or
    # an unrelated source) leaves sources as bright as an "absent" threshold (conservative)
    if ctrl and max(ctrl) > lc.MIN_FLUX_RATIO * calib:
        out["member_class"] = "untestable"
        return out
    ratio = best / calib
    out.update(ratio=float(ratio), best_snr_emp=float(best / sigma))
    out["member_class"] = member_class(ratio, best / sigma)
    if out["member_class"] == "absent" and not ctrl:
        out["member_class"] = "untestable"  # an absence needs an empty control disc to compare
    return out


def run(table_path: Path, catalog: Path, image: str, classes: tuple[str, ...]) -> list[dict]:
    t = Table.read(table_path)
    fc = np.asarray(
        t["forced_class"].filled("") if hasattr(t["forced_class"], "filled") else t["forced_class"]
    ).astype(str)
    rows = t[np.isin(fc, classes)]
    cat = lc.load_shapes(catalog)
    stamp, close = lc.image_stamper(image, STAMP_HALF_ARCSEC)
    results = []
    try:
        for r in rows:
            ra, dec = float(r["ra"]), float(r["dec"])
            sci, err, xx, yy = stamp(ra, dec)
            dx, dy, ell, th = neighbour_offsets(cat, ra, dec)
            close_ = np.hypot(dx, dy) <= lc.MAX_SEARCH_ARCSEC + NB_MARGIN_ARCSEC
            nb = list(zip(dx[close_], dy[close_], ell[close_], th[close_], strict=True))
            res = measure(sci, err, xx, yy, float(r["pred_flux"]), float(r["search_arcsec"]), nb)
            res.update(
                system=str(r["system"]),
                ra=round(ra, 7),
                dec=round(dec, 7),
                magnification=round(float(r["magnification"]), 2),
                pred_flux=float(r["pred_flux"]),
                pred_snr=round(float(r["pred_snr"]), 1),
                forced_class=str(r["forced_class"]),
            )
            results.append(
                {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()}
            )
    finally:
        close()
    return results


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--table", type=Path, required=True, help="images_predicted.ecsv of an images run"
    )
    ap.add_argument(
        "--catalog", type=Path, required=True, help="pipeline _cat.ecsv used by that run"
    )
    ap.add_argument("--forced-image", required=True, help="the run's forced _i2d (path or S3)")
    ap.add_argument(
        "--classes", default="confused", help="forced classes to re-measure (comma list)"
    )
    ap.add_argument("--out", type=Path, required=True, help="output JSON")
    args = ap.parse_args(argv)
    res = run(args.table, args.catalog, args.forced_image, tuple(args.classes.split(",")))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "table": str(args.table),
        "catalog": str(args.catalog),
        "forced_image": args.forced_image,
        "provenance": "derived (residual photometry); model_prediction (positions, pred_flux); "
        "simulated (injection recovery)",
        "assumptions": {
            k: v for k, v in globals().items() if k.isupper() and isinstance(v, (int, float, tuple))
        },
    }
    args.out.write_text(json.dumps({"meta": meta, "rows": res}, indent=1) + "\n")
    for r in res:
        print(json.dumps(r))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
