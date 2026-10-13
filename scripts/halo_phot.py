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
Each stamp also gets an injection check: a Gaussian of the measured flux placed on a ring around the
image and measured by the same route, so the halo-subtraction bias is known on the real background.

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
ANNULUS_ARCSEC = (0.4, 0.6)  # local-plane ring of the second route (ASSUMPTION)
INJECT_RING_ARCSEC = 0.7  # injections sit this far from the image, on 8 position angles
INJECT_FWHM_ARCSEC = 0.06  # Gaussian stand-in for a compact image (ASSUMPTION)
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
) -> tuple[float, float]:
    """Halo-subtracted aperture flux, its error and the concentration f(r_ap) / f(2 r_ap) near
    ``(cx, cy)`` (arcsec offsets, as ``xx``/``yy``). ``recentre = 0`` measures at ``(cx, cy)``.

    Returns NaNs when either aperture holds a non-finite pixel."""
    if pix_arcsec is None:
        pix_arcsec = float(np.hypot(xx[0, 1] - xx[0, 0], yy[0, 1] - yy[0, 0]))
    img = _fill_nan(sci)
    if route == "medfilt":
        size = int(round(filter_arcsec / pix_arcsec)) | 1
        res = img - ndimage.median_filter(img, size=size, mode="reflect")
    elif route == "annulus":
        res = img
    else:
        raise ValueError(f"unknown route {route!r}")
    px, py = cx, cy
    if recentre > 0:
        d = np.hypot(xx - cx, yy - cy)
        near = np.where(d <= recentre, ndimage.uniform_filter(res, size=3), -np.inf)
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


def gaussian(xx, yy, cx, cy, flux, fwhm=INJECT_FWHM_ARCSEC, pix_arcsec=None):
    """A unit-normalised (on the pixel grid) Gaussian of total ``flux`` at ``(cx, cy)``."""
    s = fwhm / 2.3548
    g = np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / s**2)
    return flux * g / g.sum()


def inject_check(sci, err, xx, yy, flux, route, n_angles=8, ring=INJECT_RING_ARCSEC):
    """Recovered / injected flux for Gaussians on a ring around the image (measured minus the
    same aperture without the injection, so the local scene cancels). Returns the ratios."""
    out = []
    for t in np.arange(n_angles) * 2 * np.pi / n_angles:
        cx, cy = ring * np.cos(t), ring * np.sin(t)
        base = measure(sci, err, xx, yy, route, cx, cy, recentre=0.0)[0]
        f = measure(sci + gaussian(xx, yy, cx, cy, flux), err, xx, yy, route, cx, cy, recentre=0.0)[
            0
        ]
        out.append((f - base) / flux)
    return np.asarray(out, float)


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
                    row[f"{band}_{route}_flux"], row[f"{band}_{route}_err"] = f, e
                    row[f"{band}_{route}_conc"] = c
                    if np.isfinite(f) and f > 0:
                        inj = inject_check(sci, err, xx, yy, f, route)
                        row[f"{band}_{route}_inj_median"] = float(np.median(inj))
                        row[f"{band}_{route}_inj_spread"] = float(np.ptp(inj))
        finally:
            close()
    systems_arr = np.array([r["system"] for r in rows])
    for band in bands:
        for route in ("medfilt", "annulus"):
            f = np.array([r[f"{band}_{route}_flux"] for r in rows])
            e = np.array([r[f"{band}_{route}_err"] for r in rows])
            c = np.array([r[f"{band}_{route}_conc"] for r in rows])
            ok = np.isfinite(f) & (f > 0) & (f / np.where(e > 0, e, np.inf) >= 10)
            ok &= (c >= MIN_CONCENTRATION) & (np.abs(mu) <= MAX_ABS_MU)
            src = np.where(
                ok, -2.5 * np.log10(np.where(ok, f, 1.0)) + 2.5 * np.log10(np.abs(mu)), np.nan
            )
            serr = np.where(ok, 1.0857 * e / np.where(ok, f, 1.0), np.nan)
            res, rerr, _ = lc._loo_residual(src, serr, ok, systems_arr)
            for row, a, b in zip(rows, res, rerr, strict=True):
                row[f"{band}_{route}_resid"], row[f"{band}_{route}_resid_err"] = float(a), float(b)
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
            "inject_ring_arcsec": INJECT_RING_ARCSEC,
            "inject_fwhm_arcsec": INJECT_FWHM_ARCSEC,
            "outlier_mag": OUTLIER_MAG,
            "min_concentration": MIN_CONCENTRATION,
            "max_abs_mu": MAX_ABS_MU,
            "min_snr": 10,
        },
        "images": clean,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "halo_phot.json").write_text(json.dumps(summary, indent=1) + "\n")
    for r in clean:
        print(
            r["image_id"],
            {k: v for k, v in r.items() if k.endswith(("resid", "inj_median", "_conc"))},
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
