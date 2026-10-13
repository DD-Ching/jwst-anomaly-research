"""Halo-subtracted photometry of compact multiple images on cluster-galaxy (BCG/ICL) halos.

D-032's ``fluxratios`` uses DJA ``mag_auto``, which adds halo light to images that sit on a cluster
galaxy (Abell 2744 600.1b, 602.1c). Here each image is measured on an ``_i2d`` stamp after removing
the smooth halo, by two routes:

- ``medfilt``: the stamp minus a median filter of ``filter_arcsec`` (larger than a compact image,
  smaller than a halo's scale), then a ``r_ap_arcsec`` aperture at the brightest pixel within
  ``recentre_arcsec``;
- ``annulus``: the same aperture minus the median of an ``annulus_arcsec`` ring (a local plane).

Only compact images take part: a resolved arc's aperture flux follows surface brightness, which
lensing conserves (D-032). Compact means concentrated (f(0.2")/f(0.4")) *and* round: a knot on a
thin arc passes the concentration gate (light along a line grows ∝ r), so :func:`axis_ratio` must
also be star-like (``MIN_AXIS_RATIO``, calibrated on Gaia stars per band). For two unresolved
images in one band the aperture loss cancels in the ratio, so
``s = -2.5 log10(flux) + 2.5 log10 |mu|`` is a source magnitude up to a constant.
Each stamp also gets ring checks on the real halo (:func:`ring_checks`): blank apertures measure the
halo light a route leaves behind (and give the background error), and injected PSF-width Gaussians
measure how much of a compact image the route removes. A band's ``flux_class`` is ``flux_outlier``
only when both routes exceed ``OUTLIER_MAG`` with the same sign; ``one_route`` (one route passed its
gates; reported, never flagged) and ``untested`` say that no two-route test ran.

Thresholds and filter scales are ASSUMPTIONs; outputs are derived (injections: simulated).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

FILTER_ARCSEC = (
    1.0  # median-filter size: > 3x a compact image's aperture, < a halo scale (ASSUMPTION)
)
R_AP_ARCSEC = 0.2  # aperture radius (ASSUMPTION)
RECENTRE_ARCSEC = 0.15  # peak search radius around the catalogued position (ASSUMPTION)
ANNULUS_ARCSEC = (0.5, 0.7)  # local-plane ring of the second route, outside 2 r_ap (ASSUMPTION)
CHECK_RING_ARCSEC = 1.0  # blank and injection apertures sit this far out, on 8 position angles
# Gaussian stand-ins for a compact image: about the NIRCam PSF FWHM per band (ASSUMPTION)
PSF_FWHM_ARCSEC = {
    "f115w": 0.04,
    "f150w": 0.05,
    "f200w": 0.066,
    "f277w": 0.092,
    "f356w": 0.116,
    "f410m": 0.137,
    "f444w": 0.145,
}
MIN_SNR = 10.0  # on flux / (ERR ⊕ blank-aperture scatter), as D-032 (ASSUMPTION)
MIN_CONCENTRATION = 0.6  # f(0.2")/f(0.4") below this: resolved, as exotic_screens (ASSUMPTION)
MAX_ABS_MU = 20.0  # as D-032's flux test (ASSUMPTION)
OUTLIER_MAG = 0.75  # as D-032's flux_outlier (ASSUMPTION)
ROUTES = ("medfilt", "annulus")
MIN_RECOVERY = 0.5  # a route keeping less of an injected compact image is dropped (ASSUMPTION)
# Axis-ratio gate per band (ASSUMPTION): between the known thin arcs (El Gordo 1a-c: q 0.57-0.68)
# and unsaturated Gaia DR3 stars (5th percentile 0.91 F150W / 0.94 F444W, 55 / 58 stars on the three
# fields' i2d mosaics, `halo_phot_stars.py`), so mildly stretched compact images still take part.
# Bands without a star calibration are refused.
MIN_AXIS_RATIO = {"f150w": 0.80, "f444w": 0.80}
WINDOW_SIGMA_ARCSEC = 0.15  # axis-ratio moments' Gaussian window, ~1-2.5 PSF FWHM (ASSUMPTION)
INJECT_ARCSEC = 0.45  # --inject-mag scales the filtered light within this of the image (ASSUMPTION)


def _fill_nan(img: np.ndarray) -> np.ndarray:
    out = np.array(img, float)
    bad = ~np.isfinite(out)
    if bad.any():
        out[bad] = np.nanmedian(out) if np.isfinite(out).any() else 0.0
    return out


def halo_filtered(
    sci: np.ndarray, pix_arcsec: float, filter_arcsec: float = FILTER_ARCSEC
) -> np.ndarray:
    """The stamp minus its ``filter_arcsec`` median filter (NaNs filled first)."""
    img = _fill_nan(sci)
    size = int(round(filter_arcsec / pix_arcsec)) | 1
    return img - ndimage.median_filter(img, size=size, mode="reflect")


def measure(
    sci: np.ndarray,
    err: np.ndarray,
    xx: np.ndarray,
    yy: np.ndarray,
    route: str = "medfilt",
    cx: float = 0.0,
    cy: float = 0.0,
    pix_arcsec: float | None = None,
    filter_arcsec: float = FILTER_ARCSEC,
    r_ap: float = R_AP_ARCSEC,
    recentre: float = RECENTRE_ARCSEC,
    annulus: tuple[float, float] = ANNULUS_ARCSEC,
    filtered: np.ndarray | None = None,
) -> tuple[float, float, float]:
    """Halo-subtracted aperture flux, its error and the concentration f(r_ap) / f(2 r_ap) near
    ``(cx, cy)`` (arcsec offsets, as ``xx``/``yy``). ``recentre = 0`` measures at ``(cx, cy)``;
    otherwise the aperture moves to the peak of the median-filtered image within ``recentre`` in
    both routes, so a halo gradient cannot pull it. The error is the i2d ERR only;
    :func:`measure_image` adds the ring-check terms and the recovery. ``filtered`` (from
    :func:`halo_filtered`) skips recomputing the median filter for a stamp measured again.

    Returns NaNs when either aperture holds a non-finite pixel."""
    if pix_arcsec is None:
        pix_arcsec = _pix(xx, yy)
    img = _fill_nan(sci)
    if route not in ROUTES:
        raise ValueError(f"unknown route {route!r}")
    if filtered is None:
        filtered = halo_filtered(img, pix_arcsec, filter_arcsec)
    res = filtered if route == "medfilt" else img
    px, py = _peak(filtered, xx, yy, cx, cy, recentre)
    r = np.hypot(xx - px, yy - py)
    ap, big = r <= r_ap, r <= 2 * r_ap
    if not np.isfinite(sci[big]).all():
        return float("nan"), float("nan"), float("nan")
    ring = (r >= annulus[0]) & (r <= annulus[1]) & np.isfinite(sci)
    level = float(np.median(img[ring]))
    flux = float(np.sum(res[ap]) - (level * ap.sum() if route == "annulus" else 0.0))
    # concentration on the unfiltered light in both routes: the median filter strips a
    # resolved source's extended light, which would let it pass as compact
    f_small = float(np.sum(img[ap]) - level * ap.sum())
    f_big = float(np.sum(img[big]) - level * big.sum())
    conc = f_small / f_big if f_big > 0 else float("nan")
    e = lc.ERR_SCALE * err[ap]  # ERR underestimates the noise (D-027 amendment)
    return flux, float(np.sqrt(np.sum(np.where(np.isfinite(e), e, 0.0) ** 2))), conc


def _peak(filtered, xx, yy, cx=0.0, cy=0.0, recentre=RECENTRE_ARCSEC) -> tuple[float, float]:
    """The brightest point of the (3-pixel smoothed) filtered stamp within ``recentre`` of
    ``(cx, cy)``; ``(cx, cy)`` itself when ``recentre = 0``."""
    if recentre <= 0:
        return cx, cy
    d = np.hypot(xx - cx, yy - cy)
    near = np.where(d <= recentre, ndimage.uniform_filter(filtered, size=3), -np.inf)
    k = int(np.argmax(near))
    return float(xx.flat[k]), float(yy.flat[k])


def axis_ratio(filtered, xx, yy, cx=0.0, cy=0.0, sigma_w=WINDOW_SIGMA_ARCSEC, n_iter=5) -> float:
    """Minor/major axis ratio from the Gaussian-windowed second moments of the halo-filtered light
    (negative pixels clipped) near ``(cx, cy)``. The window starts at the peak and moves to the
    windowed centroid ``n_iter`` times, so a point source off a pixel centre is not made
    elongated by the pixel grid. A point source gives ≈ 1; a thin arc that a bright knot lifts
    over the concentration gate gives ≪ 1 (its light extends along the arc, inside the window).
    Clipping makes residual noise a round positive weight, so faint sources read rounder."""
    px, py = _peak(filtered, xx, yy, cx, cy)
    light = np.clip(np.nan_to_num(filtered), 0, None)
    for _ in range(n_iter + 1):
        dx, dy = xx - px, yy - py
        w = np.exp(-0.5 * (dx**2 + dy**2) / sigma_w**2) * light
        if w.sum() <= 0:
            return float("nan")
        mx, my = float(np.sum(w * dx) / w.sum()), float(np.sum(w * dy) / w.sum())
        px, py = px + mx, py + my
    dx, dy = dx - mx, dy - my  # moments about the final centroid
    mxx, myy, mxy = (float(np.sum(w * a) / w.sum()) for a in (dx * dx, dy * dy, dx * dy))
    lam = np.linalg.eigvalsh([[mxx, mxy], [mxy, myy]])
    return float(np.sqrt(max(lam[0], 0.0) / lam[1])) if lam[1] > 0 else float("nan")


def gaussian(xx, yy, cx, cy, flux, fwhm):
    """A unit-normalised (on the pixel grid) Gaussian of total ``flux`` at ``(cx, cy)``."""
    s = fwhm / 2.3548
    g = np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / s**2)
    return flux * g / g.sum()


def _pix(xx: np.ndarray, yy: np.ndarray) -> float:
    """Pixel size (arcsec) of a stamp's offset grid."""
    return float(np.hypot(xx[0, 1] - xx[0, 0], yy[0, 1] - yy[0, 0]))


def ring_checks(sci, err, xx, yy, fluxes, fwhm, n_angles=8, ring=CHECK_RING_ARCSEC, filtered=None):
    """Two checks on apertures placed on a ring around the image, on the real halo, for each
    route in ``fluxes`` (route -> the image's measured flux):

    - ``blank``: the halo-subtracted flux of the empty apertures (no injection). Halo light the
      route leaves behind shows up here; its median is the residual bias per aperture and its
      spread (half the 16–84 % range) the background error. Neighbours can fall on the ring, so
      only robust statistics are used. The ring samples the halo 1″ out, not at the image, so
      :func:`measure_image` also counts the whole correction as an error.
    - ``recovery``: (flux with a Gaussian injected − ``blank``) / injected flux: how much of a
      compact image the route itself removes (the filter or ring eating its wings). One Gaussian
      (the routes' median flux) per angle serves both routes, so each injected stamp is filtered
      once.

    ``filtered`` is the base stamp's :func:`halo_filtered`, if already computed.
    Returns ``{route: (blank, recovery)}`` arrays."""
    pix = _pix(xx, yy)
    base = halo_filtered(sci, pix) if filtered is None else filtered
    flux = float(np.median(list(fluxes.values())))
    out = {route: ([], []) for route in fluxes}
    for t in np.arange(n_angles) * 2 * np.pi / n_angles:
        cx, cy = ring * np.cos(t), ring * np.sin(t)
        inj = sci + gaussian(xx, yy, cx, cy, flux, fwhm)
        inj_filt = halo_filtered(inj, pix)
        for route in fluxes:
            # same peak search as the real measurement, so its positive bias shows up here too
            b = measure(sci, err, xx, yy, route, cx, cy, filtered=base)[0]
            f = measure(inj, err, xx, yy, route, cx, cy, filtered=inj_filt)[0]
            out[route][0].append(b)
            out[route][1].append((f - b) / flux)
    return {k: (np.asarray(b, float), np.asarray(r, float)) for k, (b, r) in out.items()}


def measure_image(sci, err, xx, yy, mu: float, band: str) -> dict:
    """Both routes and their ring checks for one image in one band (a process-pool task).

    The flux is corrected for the leftover halo (blank median) and divided by the recovery; the
    error adds, in quadrature, the blank scatter and the whole halo correction (ASSUMPTION: the
    ring 1″ out may misjudge the leftover halo at the image, where the cusp is steeper). A route
    that keeps less than ``MIN_RECOVERY`` of an injected image is dropped (NaN flux). Ring checks
    are skipped for images the gates drop anyway (concentration, axis ratio, |μ|)."""
    filt = halo_filtered(sci, _pix(xx, yy))
    meas = {route: measure(sci, err, xx, yy, route, filtered=filt) for route in ROUTES}
    conc = meas["medfilt"][2]  # the same unfiltered measure in both routes
    q = axis_ratio(filt, xx, yy)
    out = {f"{band}_conc": conc, f"{band}_axis_ratio": q}
    gated = _compact(conc, q, band) and abs(mu) <= MAX_ABS_MU
    todo = {r: m[0] for r, m in meas.items() if gated and np.isfinite(m[0]) and m[0] > 0}
    checks = (
        ring_checks(sci, err, xx, yy, todo, PSF_FWHM_ARCSEC[band], filtered=filt) if todo else {}
    )
    for route, (f, e, _) in meas.items():
        blank_med, blank_sd, rec = np.nan, 0.0, np.nan
        if route in checks:
            blank, inj = checks[route]
            blank = blank[np.isfinite(blank)]
            if blank.size:
                blank_med = float(np.median(blank))
                blank_sd = float(np.diff(np.percentile(blank, [16, 84]))[0] / 2)
            rec = float(np.nanmedian(inj))
        corr = f - blank_med if np.isfinite(blank_med) else f
        scale = rec if np.isfinite(rec) else 1.0  # no check ran: the gates drop this image
        bias = abs(blank_med) if np.isfinite(blank_med) else 0.0
        if np.isfinite(rec) and rec < MIN_RECOVERY:
            corr = np.nan  # the route removes too much of a compact image to be corrected
        out[f"{band}_{route}_flux"] = corr / scale
        out[f"{band}_{route}_err"] = float(np.sqrt(e**2 + blank_sd**2 + bias**2)) / scale
        out[f"{band}_{route}_blank_frac"] = blank_med / f if f > 0 else np.nan
        out[f"{band}_{route}_recovery"] = rec
    return out


def _compact(conc, q, band: str):
    """The compactness gates: concentration and axis ratio (NaN fails both)."""
    with np.errstate(invalid="ignore"):
        return (np.asarray(conc) >= MIN_CONCENTRATION) & (np.asarray(q) >= MIN_AXIS_RATIO[band])


def _init_worker() -> None:
    """One BLAS/OpenMP thread per pool worker (scripts/CLAUDE.md), also on fork."""
    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:  # optional; median_filter itself is single-threaded
        return
    threadpool_limits(1)


def _scaled_image(stamp, inject_mag: float):
    """``stamp`` with the image's own halo-subtracted light (the filtered stamp within
    ``INJECT_ARCSEC``) scaled so its flux changes by ``inject_mag`` mag. The real profile and the
    halo stay; the pixel noise and any neighbour inside the radius are scaled too."""
    sci, err, xx, yy = stamp
    near = np.hypot(xx, yy) <= INJECT_ARCSEC
    image = np.where(near, halo_filtered(sci, _pix(xx, yy)), 0.0)
    return sci + (10 ** (-0.4 * inject_mag) - 1) * image, err, xx, yy


def classify(rows: list[dict], bands) -> None:
    """Gates, leave-one-out residuals and ``flux_class`` per band, in place (see module doc)."""
    systems_arr = np.array([r["system"] for r in rows])
    mu = np.array([r["mu"] for r in rows], float)
    for band in bands:
        for route in ROUTES:
            f = np.array([r[f"{band}_{route}_flux"] for r in rows], float)
            e = np.array([r[f"{band}_{route}_err"] for r in rows], float)
            c = np.array([r[f"{band}_conc"] for r in rows], float)
            q = np.array([r[f"{band}_axis_ratio"] for r in rows], float)
            ok = np.isfinite(f) & (f > 0) & (f / np.where(e > 0, e, np.inf) >= MIN_SNR)
            ok &= _compact(c, q, band) & (np.abs(mu) <= MAX_ABS_MU)
            src = np.where(
                ok, -2.5 * np.log10(np.where(ok, f, 1.0)) + 2.5 * np.log10(np.abs(mu)), np.nan
            )
            serr = np.where(ok, 1.0857 * e / np.where(ok, f, 1.0), np.nan)
            res, rerr, worst = lc._loo_residual(src, serr, ok, systems_arr)
            for row, a, b, w in zip(rows, res, rerr, worst, strict=True):
                row[f"{band}_{route}_resid"], row[f"{band}_{route}_resid_err"] = float(a), float(b)
                row[f"{band}_{route}_worst"] = bool(w)
    for row in rows:  # an outlier in a band needs both routes beyond the threshold, same sign
        for band in bands:
            a, b = row[f"{band}_medfilt_resid"], row[f"{band}_annulus_resid"]
            if np.isfinite(a) and np.isfinite(b):
                # only a system's worst image (both of a pair), as D-032: one bad image shifts
                # its siblings' references, which must not flag them too
                worst = row[f"{band}_medfilt_worst"] and row[f"{band}_annulus_worst"]
                big = abs(a) > OUTLIER_MAG and abs(b) > OUTLIER_MAG and np.sign(a) == np.sign(b)
                row[f"{band}_flux_class"] = "flux_outlier" if big and worst else "consistent"
            elif np.isfinite(a) or np.isfinite(b):  # one route only: reported, never flagged
                row[f"{band}_flux_class"] = "one_route"
            else:
                row[f"{band}_flux_class"] = "untested"


def run(
    model_name: str,
    systems: list[str],
    bands: dict[str, str],
    half_arcsec: float = 3.0,
    inject_mags: tuple[float, ...] = (),
    workers: int | None = None,
):
    """Measure and classify ``systems``' images. Each of ``inject_mags`` is an end-to-end
    injection: each image in turn is changed by that many mag (negative: brighter) with
    :func:`_scaled_image`, and the whole chain (both routes, ring checks, gates, residual, class)
    re-runs for it; the others keep their base measurements. Each record holds the target's class
    before and after, per band (``flux_outlier`` after = recovered; a target untested before was
    never testable and does not count as a miss). Returns ``(rows, injections)``."""
    if not lc.has_image_list(model_name) or lc.is_map_model(model_name):
        raise SystemExit(f"error: {model_name}: needs a Lenstool model with a usable image list")
    known = sorted(set(PSF_FWHM_ARCSEC) & set(MIN_AXIS_RATIO))
    missing = [b for b in bands if b not in known]
    if missing:
        raise SystemExit(
            f"error: no PSF width or star-calibrated axis-ratio gate for {missing}; known: {known}"
        )
    model, files, par = lc.load_model(model_name)
    images, _ = lc.image_list(model_name, files, par)
    lc.apply_frame_offset(model_name, model, images)  # into the JWST frame (D-034)
    known = set(np.asarray(images["system"], str))
    unknown = [s for s in systems if s not in known]
    if unknown:
        raise SystemExit(f"error: {model_name}: no systems {unknown}")
    images = images[np.isin(np.asarray(images["system"], str), systems)]
    z = lensmodel.image_redshifts(images, par["z_m_limit"])
    mu = np.asarray(model.evaluate(images["ra"], images["dec"], z)["magnification"], float)
    rows = [
        {"image_id": str(r["image_id"]), "system": str(r["system"]), "mu": float(m)}
        for r, m in zip(images, mu, strict=True)
    ]
    stamps = {}
    for band, uri in bands.items():
        stamp, close = lc.image_stamper(uri, half_arcsec)
        try:
            stamps[band] = [stamp(float(r["ra"]), float(r["dec"])) for r in images]
        finally:
            close()
    # CPU-bound: one process per core, one BLAS thread each (scripts/CLAUDE.md)
    with ProcessPoolExecutor(
        max_workers=workers or os.cpu_count(), initializer=_init_worker
    ) as pool:
        base_f = {
            (band, k): pool.submit(measure_image, *st, mu[k], band)
            for band, sts in stamps.items()
            for k, st in enumerate(sts)
        }
        # every injection is independent: submit them all before waiting on any
        inj_f = {
            (dm, k, band): pool.submit(
                measure_image, *_scaled_image(stamps[band][k], dm), mu[k], band
            )
            for dm in inject_mags
            for k in range(len(rows))
            for band in bands
        }
        for (_, k), fut in base_f.items():
            rows[k].update(fut.result())
        classify(rows, bands)
        injections = []
        for dm in inject_mags:
            for k, target in enumerate(rows):
                out = [dict(r) for r in rows]
                for band in bands:
                    out[k].update(inj_f[(dm, k, band)].result())
                classify(out, bands)
                injections.append(
                    {
                        "image_id": target["image_id"],
                        "inject_mag": dm,
                        "before": {b: target[f"{b}_flux_class"] for b in bands},
                        "after": {b: out[k][f"{b}_flux_class"] for b in bands},
                    }
                )
    return rows, injections


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", required=True)
    ap.add_argument("--systems", required=True, help="comma-separated system ids")
    ap.add_argument("--band", action="append", required=True, help="BAND=i2d-uri (repeatable)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--inject-mag",
        type=float,
        action="append",
        default=[],
        help="end-to-end injection per image, mag (repeatable; negative = brighter)",
    )
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args(argv)
    bands = {k.strip().lower(): v for k, v in (b.split("=", 1) for b in args.band)}
    systems = [s.strip() for s in args.systems.split(",") if s.strip()]
    rows, injections = run(
        args.model, systems, bands, inject_mags=tuple(args.inject_mag), workers=args.workers
    )
    clean = [{k: (lc._r(v, 3) if isinstance(v, float) else v) for k, v in r.items()} for r in rows]
    summary = {
        "model": args.model,
        "provenance": "derived (injections: simulated)",
        "bands": bands,
        "assumptions": {
            "filter_arcsec": FILTER_ARCSEC,
            "r_ap_arcsec": R_AP_ARCSEC,
            "recentre_arcsec": RECENTRE_ARCSEC,
            "annulus_arcsec": list(ANNULUS_ARCSEC),
            "check_ring_arcsec": CHECK_RING_ARCSEC,
            "psf_fwhm_arcsec": PSF_FWHM_ARCSEC,
            "outlier_mag": OUTLIER_MAG,
            "min_concentration": MIN_CONCENTRATION,
            "min_axis_ratio": MIN_AXIS_RATIO,
            "window_sigma_arcsec": WINDOW_SIGMA_ARCSEC,
            "max_abs_mu": MAX_ABS_MU,
            "min_snr": MIN_SNR,
            "min_recovery": MIN_RECOVERY,
            "inject_arcsec": INJECT_ARCSEC,
        },
        "images": clean,
        "injections": injections,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "halo_phot.json").write_text(json.dumps(summary, indent=1) + "\n")
    for r in clean:
        print(
            r["image_id"],
            {
                k: v
                for k, v in r.items()
                if k.endswith(("resid", "blank_frac", "recovery", "class"))
            },
        )
    for inj in injections:
        print("inject", inj)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
