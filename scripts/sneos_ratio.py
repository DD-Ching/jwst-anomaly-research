"""SN Eos 101.2 / 101.1 difference-image flux ratio vs the VENUS AstroLensPy v2.1 model.

Frozen plan: ``results/transients/sneos_ratio/PLAN.md``. This module holds the plan's rules as
functions (resampling, registration check, aperture ratio, injection, decision); the run driver
(empirical-PSF star, S3 stamps of the 9 mosaics, contact sheet) is the next step. R_MODEL is a
model_prediction from ``results/transients/sneos/summary.json``.
"""

from __future__ import annotations

import numpy as np
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve

R_MODEL = 1.162  # |μ_101.2 / μ_101.1|, sneos/summary.json (PLAN.md)
APER_ARCSEC = 0.20  # ASSUMPTION (PLAN.md rule 2)
ANNULUS_ARCSEC = (0.35, 0.60)  # PLAN.md rule 2
MAX_SHIFT_ARCSEC = 0.03  # ASSUMPTION (PLAN.md rule 1)
MIN_SNR = 5.0  # ASSUMPTION (PLAN.md rule 3)
INJ_TOL = 0.10  # ASSUMPTION (PLAN.md rule 4)
INJ_RATIOS = (1.16, 0.70)  # PLAN.md rule 4
FLOOR_LN = 0.30  # ASSUMPTION (PLAN.md rule 6)
SN_MASK_ARCSEC = 0.4  # PLAN.md rule 1


def resample(src: np.ndarray, src_wcs, dst_wcs, shape: tuple[int, int]) -> np.ndarray:
    """``src`` (surface brightness) bilinearly resampled onto the ``dst_wcs`` pixel grid."""
    jj, ii = np.indices(shape, dtype=float)
    ra, dec = dst_wcs.pixel_to_world_values(ii, jj)
    si, sj = src_wcs.world_to_pixel_values(ra, dec)
    return map_coordinates(src, [sj, si], order=1, cval=np.nan)


def register_offset(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """Sub-pixel (dx, dy) of ``b`` relative to ``a`` from the cross-correlation peak (parabolic
    refinement); pixels in ``mask`` and non-finite pixels are set to zero after mean removal."""

    def prep(x):
        x = np.where(mask | ~np.isfinite(x), np.nan, x)
        x = x - np.nanmean(x)
        return np.nan_to_num(x)

    pa, pb = prep(a), prep(b)
    cc = fftconvolve(pb, pa[::-1, ::-1], mode="same")
    j, i = np.unravel_index(np.argmax(cc), cc.shape)

    def refine(m, c, p):
        den = m - 2 * c + p
        return 0.0 if den == 0 else 0.5 * (m - p) / den

    di = refine(cc[j, i - 1], cc[j, i], cc[j, i + 1]) if 0 < i < cc.shape[1] - 1 else 0.0
    dj = refine(cc[j - 1, i], cc[j, i], cc[j + 1, i]) if 0 < j < cc.shape[0] - 1 else 0.0
    cy, cx = cc.shape[0] // 2, cc.shape[1] // 2  # zero lag of a "same" correlation
    return float(i + di - cx), float(j + dj - cy)


def aperture_pair(
    diff: np.ndarray, err: np.ndarray, p1: tuple[float, float], p2: tuple[float, float], pix: float
) -> dict:
    """PLAN.md rule 2: background-subtracted aperture fluxes at ``p1`` (101.1) and ``p2``
    (101.2), pixel positions (x, y); the annulus excludes the other image's aperture."""
    r = APER_ARCSEC / pix
    rin, rout = (a / pix for a in ANNULUS_ARCSEC)
    yy, xx = np.indices(diff.shape)
    out = {}
    for name, p, q in (("f1", p1, p2), ("f2", p2, p1)):
        other = np.hypot(xx - q[0], yy - q[1]) <= r
        ann = CircularAnnulus(p, rin, rout).to_mask(method="center").to_image(diff.shape) > 0
        sel = ann & ~other & np.isfinite(diff)
        bkg = float(np.median(diff[sel]))
        ap = CircularAperture(p, r)
        bad = ~np.isfinite(diff) | ~np.isfinite(err)
        t = aperture_photometry(diff - bkg, ap, error=np.where(bad, 0, err), mask=bad)
        out[name] = float(t["aperture_sum"][0])
        out["s" + name] = float(t["aperture_sum_err"][0])
    f1, f2, s1, s2 = out["f1"], out["f2"], out["sf1"], out["sf2"]
    ratio = f2 / f1
    out["snr1"], out["snr2"] = f1 / s1 if s1 else np.inf, f2 / s2 if s2 else np.inf
    out["ratio"] = ratio
    out["sratio"] = abs(ratio) * float(np.hypot(s1 / f1, s2 / f2))
    out["detected"] = bool(out["snr1"] >= MIN_SNR and out["snr2"] >= MIN_SNR)
    return out


def inject_pair(
    img: np.ndarray, psf: np.ndarray, p1: tuple[float, float], sep: tuple[float, float], f1, ratio
) -> np.ndarray:
    """``img`` plus a unit-sum ``psf`` scaled to ``f1`` at ``p1`` and ``ratio * f1`` at
    ``p1 + sep`` (pixels; integer-pixel placement after a sub-pixel bilinear shift)."""
    out = np.array(img, float, copy=True)
    k = psf / psf.sum()
    h0, h1 = k.shape[0] // 2, k.shape[1] // 2
    for (x, y), f in (
        (p1, f1),
        ((p1[0] + sep[0], p1[1] + sep[1]), ratio * f1),
    ):
        ix, iy = int(np.floor(x)), int(np.floor(y))
        fx, fy = x - ix, y - iy
        jj, ii = np.indices(k.shape, dtype=float)
        ks = map_coordinates(k, [jj - fy, ii - fx], order=1, cval=0.0)
        out[iy - h0 : iy - h0 + k.shape[0], ix - h1 : ix - h1 + k.shape[1]] += f * ks
    return out


def recovered(measured: float, injected: float) -> bool:
    """PLAN.md rule 4: ratio recovered within INJ_TOL (relative)."""
    return bool(abs(measured / injected - 1.0) <= INJ_TOL)


def decide(ln_r: float, sig_ln: float, r_model: float = R_MODEL) -> str:
    """PLAN.md rule 6: 'differs' when |Δ| > 3σ and > FLOOR_LN; otherwise 'agrees'."""
    if not np.isfinite(ln_r) or not np.isfinite(sig_ln):
        return "undetermined"
    d = abs(ln_r - np.log(r_model))
    return "differs" if (d > 3.0 * sig_ln and d > FLOOR_LN) else "agrees"


def combine(ln_r: np.ndarray, sig: np.ndarray) -> tuple[float, float]:
    """PLAN.md rule 5: inverse-variance mean of ln R; σ = max(propagated, scatter / √n)."""
    ln_r, sig = np.asarray(ln_r, float), np.asarray(sig, float)
    if ln_r.size == 0:
        return float("nan"), float("nan")
    w = 1.0 / sig**2
    mean = float(np.sum(w * ln_r) / np.sum(w))
    prop = float(1.0 / np.sqrt(np.sum(w)))
    scat = float(np.std(ln_r, ddof=1) / np.sqrt(ln_r.size)) if ln_r.size > 1 else 0.0
    return mean, max(prop, scat)
