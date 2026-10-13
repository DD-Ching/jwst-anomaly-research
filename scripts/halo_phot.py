"""Halo-subtracted photometry of compact multiple images on cluster-galaxy (BCG/ICL) halos.

D-032's ``fluxratios`` uses DJA ``mag_auto``, which adds halo light to images that sit on a cluster
galaxy (Abell 2744 600.1b, 602.1c). Here each image is measured on an ``_i2d`` stamp after removing
the smooth halo, by two routes:

- ``medfilt``: the stamp minus a median filter of ``filter_arcsec`` (larger than a compact image,
  smaller than a halo's scale), then a ``r_ap_arcsec`` aperture at the brightest pixel within
  ``recentre_arcsec``;
- ``annulus``: the same aperture minus the median of an ``annulus_arcsec`` ring (a local plane).

Only compact images take part: a resolved arc's aperture flux follows surface brightness, which
lensing conserves (D-032). For two unresolved images in one band the aperture loss cancels in the
ratio, so ``s = -2.5 log10(flux) + 2.5 log10 |mu|`` is a source magnitude up to a constant.
Each stamp also gets ring checks on the real halo (:func:`ring_checks`): blank apertures measure the
halo light a route leaves behind (and give the background error), and injected PSF-width Gaussians
measure how much of a compact image the route removes. An image is a ``flux_outlier`` in a band only
when both routes exceed ``OUTLIER_MAG`` with the same sign.

Thresholds and filter scales are ASSUMPTIONs; outputs are derived (injections: simulated).
"""

from __future__ import annotations

import argparse
import json
import sys
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


def _fill_nan(img: np.ndarray) -> np.ndarray:
    out = np.array(img, float)
    bad = ~np.isfinite(out)
    if bad.any():
        out[bad] = np.nanmedian(out) if np.isfinite(out).any() else 0.0
    return out


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
) -> tuple[float, float, float]:
    """Halo-subtracted aperture flux, its error and the concentration f(r_ap) / f(2 r_ap) near
    ``(cx, cy)`` (arcsec offsets, as ``xx``/``yy``). ``recentre = 0`` measures at ``(cx, cy)``;
    otherwise the aperture moves to the peak of the median-filtered image within ``recentre`` in
    both routes, so a halo gradient cannot pull it. The error is the i2d ERR only; :func:`run`
    adds the blank-aperture scatter.

    Returns NaNs when either aperture holds a non-finite pixel."""
    if pix_arcsec is None:
        pix_arcsec = float(np.hypot(xx[0, 1] - xx[0, 0], yy[0, 1] - yy[0, 0]))
    img = _fill_nan(sci)
    if route not in ("medfilt", "annulus"):
        raise ValueError(f"unknown route {route!r}")
    size = int(round(filter_arcsec / pix_arcsec)) | 1
    filtered = img - ndimage.median_filter(img, size=size, mode="reflect")
    res = filtered if route == "medfilt" else img
    px, py = cx, cy
    if recentre > 0:
        d = np.hypot(xx - cx, yy - cy)
        near = np.where(d <= recentre, ndimage.uniform_filter(filtered, size=3), -np.inf)
        k = int(np.argmax(near))
        px, py = float(xx.flat[k]), float(yy.flat[k])
    r = np.hypot(xx - px, yy - py)
    ap, big = r <= r_ap, r <= 2 * r_ap
    if not np.isfinite(sci[big]).all():
        return float("nan"), float("nan"), float("nan")
    level = 0.0
    if route == "annulus":
        ring = (r >= annulus[0]) & (r <= annulus[1]) & np.isfinite(sci)
        level = float(np.median(res[ring]))
    flux = float(np.sum(res[ap]) - level * ap.sum())
    f_big = float(np.sum(res[big]) - level * big.sum())
    conc = flux / f_big if f_big > 0 else float("nan")
    e = err[ap]
    return flux, float(np.sqrt(np.sum(np.where(np.isfinite(e), e, 0.0) ** 2))), conc


def gaussian(xx, yy, cx, cy, flux, fwhm):
    """A unit-normalised (on the pixel grid) Gaussian of total ``flux`` at ``(cx, cy)``."""
    s = fwhm / 2.3548
    g = np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / s**2)
    return flux * g / g.sum()


def ring_checks(sci, err, xx, yy, flux, route, fwhm, n_angles=8, ring=CHECK_RING_ARCSEC):
    """Two checks on apertures placed on a ring around the image, on the real halo:

    - ``blank``: the halo-subtracted flux of the empty apertures (no injection). Halo light the
      route leaves behind shows up here; its median is the residual bias per aperture and its
      spread (half the 16–84 % range) the background error. Neighbours can fall on the ring, so
      only robust statistics are used.
    - ``recovery``: (flux with a Gaussian of ``flux`` injected − ``blank``) / ``flux``: how much
      of a compact image the route itself removes (the filter or ring eating its wings).

    Returns ``(blank, recovery)`` arrays."""
    blank, rec = [], []
    for t in np.arange(n_angles) * 2 * np.pi / n_angles:
        cx, cy = ring * np.cos(t), ring * np.sin(t)
        b = measure(sci, err, xx, yy, route, cx, cy, recentre=0.0)[0]
        f = measure(
            sci + gaussian(xx, yy, cx, cy, flux, fwhm), err, xx, yy, route, cx, cy, recentre=0.0
        )[0]
        blank.append(b)
        rec.append((f - b) / flux)
    return np.asarray(blank, float), np.asarray(rec, float)


def run(model_name: str, systems: list[str], bands: dict[str, str], half_arcsec: float = 3.0):
    files = lc.model_files(model_name)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    images = images[np.isin(np.asarray(images["system"], str), systems)]
    z = lensmodel.image_redshifts(images, par["z_m_limit"])
    mu = np.asarray(model.evaluate(images["ra"], images["dec"], z)["magnification"], float)
    rows = [
        {"image_id": str(r["image_id"]), "system": str(r["system"]), "mu": float(m)}
        for r, m in zip(images, mu, strict=True)
    ]
    for band, uri in bands.items():
        stamp, close = lc.image_stamper(uri, half_arcsec)
        try:
            for row, r in zip(rows, images, strict=True):
                sci, err, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
                for route in ("medfilt", "annulus"):
                    f, e, c = measure(sci, err, xx, yy, route)
                    row[f"{band}_{route}_conc"] = c
                    blank_med, blank_sd, rec = np.nan, 0.0, np.nan
                    if np.isfinite(f) and f > 0:
                        blank, inj = ring_checks(sci, err, xx, yy, f, route, PSF_FWHM_ARCSEC[band])
                        blank = blank[np.isfinite(blank)]
                        if blank.size:
                            blank_med = float(np.median(blank))
                            blank_sd = float(np.diff(np.percentile(blank, [16, 84]))[0] / 2)
                        rec = float(np.nanmedian(inj))
                    row[f"{band}_{route}_flux"] = f
                    row[f"{band}_{route}_err"] = float(np.hypot(e, blank_sd))
                    row[f"{band}_{route}_blank_frac"] = blank_med / f if f > 0 else np.nan
                    row[f"{band}_{route}_recovery"] = rec
        finally:
            close()
    systems_arr = np.array([r["system"] for r in rows])
    for band in bands:
        for route in ("medfilt", "annulus"):
            f = np.array([r[f"{band}_{route}_flux"] for r in rows])
            e = np.array([r[f"{band}_{route}_err"] for r in rows])
            c = np.array([r[f"{band}_{route}_conc"] for r in rows])
            ok = np.isfinite(f) & (f > 0) & (f / np.where(e > 0, e, np.inf) >= MIN_SNR)
            ok &= (c >= MIN_CONCENTRATION) & (np.abs(mu) <= MAX_ABS_MU)
            src = np.where(
                ok, -2.5 * np.log10(np.where(ok, f, 1.0)) + 2.5 * np.log10(np.abs(mu)), np.nan
            )
            serr = np.where(ok, 1.0857 * e / np.where(ok, f, 1.0), np.nan)
            res, rerr, _ = lc._loo_residual(src, serr, ok, systems_arr)
            for row, a, b in zip(rows, res, rerr, strict=True):
                row[f"{band}_{route}_resid"], row[f"{band}_{route}_resid_err"] = float(a), float(b)
    for row in rows:  # an outlier in a band needs both routes beyond the threshold, same sign
        for band in bands:
            a, b = row[f"{band}_medfilt_resid"], row[f"{band}_annulus_resid"]
            row[f"{band}_flux_outlier"] = bool(
                abs(a) > OUTLIER_MAG and abs(b) > OUTLIER_MAG and np.sign(a) == np.sign(b)
            )
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", required=True)
    ap.add_argument("--systems", required=True, help="comma-separated system ids")
    ap.add_argument("--band", action="append", required=True, help="BAND=i2d-uri (repeatable)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    bands = dict(b.split("=", 1) for b in args.band)
    rows = run(args.model, args.systems.split(","), bands)
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
            "max_abs_mu": MAX_ABS_MU,
            "min_snr": MIN_SNR,
        },
        "images": clean,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "halo_phot.json").write_text(json.dumps(summary, indent=1) + "\n")
    for r in clean:
        print(
            r["image_id"],
            {
                k: v
                for k, v in r.items()
                if k.endswith(("resid", "blank_frac", "recovery", "outlier"))
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
