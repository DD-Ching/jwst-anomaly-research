"""Two-epoch difference imaging of a small sky box from level-3 mosaics (arcs; D-027 addendum).

The catalog grid (``transient_grid.py``) measures only catalogued positions in 0.15″ apertures, so
transients on bright lensed arcs (caustic-crossing stars, SNe on their host) fall outside it. This
script reads a box around a sky position from both epochs' ``_i2d`` mosaics by byte range (full
rows, one range read per extension), resamples epoch 2 onto epoch 1's pixel grid, registers the
two by FFT cross-correlation, matches the PSF width with a one-parameter Gaussian blur of the
sharper epoch (the blur that minimises the robust rms of the difference), removes the smooth
residual of the difference with a median filter, and forms a matched-filter significance map.
``--kernel-px N`` replaces the Gaussian with a least-squares pixel kernel fitted on the box
(``kernel_fit``); ``--kernel-stamps`` adds a second kernel for compact-source stamps, fitted
leaving each stamp out (``kernel_match``). Hits within ``--core-arcsec`` of a pixel above the
``--core-pct`` percentile of epoch 1 are flagged ``near_bright_core`` (PSF-rotation rosettes).

Noise: the matched-filter significance uses the pipeline ERR of both epochs; its robust std over
the box is the noise scale (never below 1), as in ``transient_grid.py``. A **detection** is a local
maximum of |S| ≥ ``--min-sigma`` in one band; a **hit** is a group of detections within
``--link-arcsec`` in ≥ ``--min-bands`` bands (either sign). ``--inject`` adds Gaussian point
sources of given AB magnitudes in epoch 2 on pixels brighter than the box's 90th percentile in
epoch 1 (the arcs and bright galaxies) and reports the recovered fraction through the same chain.
``--inject-grid N`` instead adds N stars (``--inject-grid-mags`` in turn) at common sky positions
on those pixels, ``--inject-grid-sep`` apart, to epoch 2 of every band before detection; hits, maps
and the per-star match (``injected.csv``) then carry them, so later stages (``arc_diff_morph.py
--known injected.csv``) are calibrated through the whole chain. Its outputs are injection runs,
never a screen. Every number is ``derived``; thresholds are ASSUMPTIONs.

Validation targets (``--known``): a CSV with name, ra, dec, shifted by ``--frame-offset`` (published
frame → pipeline frame) and matched to hits within ``--match-arcsec``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import pickle
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy import ndimage

# PSF FWHM (arcsec) per band, approximate NIRCam values (JWST user documentation); ASSUMPTION.
PSF_FWHM = {
    "F090W": 0.034,
    "F115W": 0.040,
    "F150W": 0.050,
    "F200W": 0.066,
    "F277W": 0.092,
    "F356W": 0.116,
    "F410M": 0.137,
    "F444W": 0.145,
}


def _wcs(header: fits.Header) -> WCS:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(header)


def read_box_raw(uri: str, ra: float, dec: float, half_arcsec: float) -> dict:
    """SCI, ERR, WCS (of the box), pixel scale and PIXAR_SR of a box of a level-3 mosaic.

    Remote URIs (``s3://``, ``https://``) are opened with fsspec; the rows of the box are read in
    one range per extension and the columns sliced afterwards.
    """
    if "://" in uri:
        import fsspec

        opts = {"anon": True} if uri.startswith("s3://") else {}
        fs, path = fsspec.core.url_to_fs(uri, **opts)
        fobj = fs.open(path, "rb", block_size=2**22, cache_type="readahead")
    else:
        fobj = open(uri, "rb")  # noqa: SIM115 (closed below)
    try:
        with fits.open(fobj, lazy_load_hdus=True) as h:
            hdr = h["SCI"].header
            wcs = _wcs(hdr)
            scale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
            ny, nx = int(hdr["NAXIS2"]), int(hdr["NAXIS1"])
            x, y = (float(v) for v in wcs.all_world2pix(ra, dec, 0))
            half = int(np.ceil(half_arcsec / scale))
            x0, x1 = max(int(round(x)) - half, 0), min(int(round(x)) + half + 1, nx)
            y0, y1 = max(int(round(y)) - half, 0), min(int(round(y)) + half + 1, ny)
            if x1 - x0 < 16 or y1 - y0 < 16:
                raise ValueError(f"{uri}: box outside the mosaic")
            out = {}
            for ext in ("SCI", "ERR"):
                rows = np.asarray(h[ext].section[y0:y1, :], dtype=float)
                out[ext] = rows[:, x0:x1]
            pixar = float(hdr.get("PIXAR_SR", (scale / 206264.806) ** 2))
    finally:
        fobj.close()
    return {
        "sci": out["SCI"],
        "err": out["ERR"],
        "header": hdr.tostring(),
        "x0": x0,
        "y0": y0,
        "scale": scale,
        "pixar_sr": pixar,
    }


def box_from_raw(raw: dict) -> dict:
    """``read_box`` output from ``read_box_raw`` output (the WCS is rebuilt from the header text,
    so a cached raw box gives exactly the same WCS as a fresh read)."""
    sub = _wcs(fits.Header.fromstring(raw["header"]))
    sub.wcs.crpix = sub.wcs.crpix - np.array([raw["x0"], raw["y0"]])
    return {
        "sci": raw["sci"],
        "err": raw["err"],
        "wcs": sub,
        "scale": raw["scale"],
        "pixar_sr": raw["pixar_sr"],
    }


def read_box(uri: str, ra: float, dec: float, half_arcsec: float) -> dict:
    """SCI, ERR, WCS (of the box), pixel scale and PIXAR_SR of a box of a level-3 mosaic."""
    return box_from_raw(read_box_raw(uri, ra, dec, half_arcsec))


CACHE_VERSION = "2026-10-17.1"  # bump when a cached function's output changes


def cached(cache_dir: str | None, fn, *args):
    """``fn(*args)``, stored as a pickle under ``cache_dir`` keyed by ``fn``'s name and ``args``.

    For raw-data products only (archive boxes, empirical PSFs, the injection-free registration):
    nothing that depends on injected sources may go through this cache. Every argument of ``fn``
    is part of the key, so a changed setting never serves a stale entry. Each call returns a fresh
    copy; a WCS is never pickled (``read_box_raw`` keeps the header text instead). Code changes
    to a cached function need a ``CACHE_VERSION`` bump.
    """
    if not cache_dir:
        return fn(*args)
    h = hashlib.sha1(repr((CACHE_VERSION, fn.__name__, *args)).encode()).hexdigest()[:20]
    path = Path(cache_dir) / f"{fn.__name__}_{h}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    out = fn(*args)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f".{os.getpid()}.tmp")
    tmp.write_bytes(pickle.dumps(out))
    tmp.replace(path)
    return out


def resample(img: np.ndarray, wcs_from: WCS, wcs_to: WCS, shape: tuple[int, int]) -> np.ndarray:
    """``img`` (on ``wcs_from``) resampled onto ``wcs_to``'s pixel grid (cubic, NaN outside)."""
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    ra, dec = wcs_to.all_pix2world(xx.ravel(), yy.ravel(), 0)
    x, y = wcs_from.all_world2pix(ra, dec, 0)
    good = np.isfinite(img)
    filled = np.where(good, img, 0.0)
    val = ndimage.map_coordinates(filled, [y, x], order=3, mode="constant", cval=np.nan)
    ok = ndimage.map_coordinates(good.astype(float), [y, x], order=1, mode="constant", cval=0.0)
    val[ok < 0.999] = np.nan
    return val.reshape(shape)


def register(a: np.ndarray, b: np.ndarray, max_shift: float = 3.0) -> tuple[float, float]:
    """Sub-pixel (dy, dx) such that ``b`` shifted by it matches ``a`` (FFT cross-correlation).

    Both images are high-passed (minus a 9-px median) so the correlation is set by compact
    structure; the peak is refined with a parabola in each axis.
    """

    def prep(im: np.ndarray) -> np.ndarray:
        im = np.nan_to_num(im - np.nanmedian(im))
        im = im - ndimage.median_filter(im, size=9)
        lim = np.percentile(np.abs(im), 99.5)
        return np.clip(im, -lim, lim)

    pa, pb = prep(a), prep(b)
    cc = np.fft.irfft2(np.fft.rfft2(pa) * np.conj(np.fft.rfft2(pb)), s=pa.shape)
    cc = np.fft.fftshift(cc)
    cy, cx = np.array(cc.shape) // 2
    r = int(np.ceil(max_shift))
    win = cc[cy - r : cy + r + 1, cx - r : cx + r + 1]
    iy, ix = np.unravel_index(np.argmax(win), win.shape)

    def parab(m1: float, c: float, p1: float) -> float:
        d = m1 - 2 * c + p1
        return 0.0 if d == 0 else 0.5 * (m1 - p1) / d

    dy = iy - r
    dx = ix - r
    if 0 < iy < win.shape[0] - 1:
        dy += parab(win[iy - 1, ix], win[iy, ix], win[iy + 1, ix])
    if 0 < ix < win.shape[1] - 1:
        dx += parab(win[iy, ix - 1], win[iy, ix], win[iy, ix + 1])
    return float(dy), float(dx)


def register_median(
    a: np.ndarray, b: np.ndarray, tile: int, max_shift: float = 3.0
) -> tuple[float, float, int]:
    """Global (dy, dx) as the median of ``register`` over non-overlapping ``tile``-px tiles.

    The whole-box correlation can be pulled to its window edge by a few tiles holding a bright
    saturated star and its spikes (El Gordo F150W, 2026-10-17: (1.71, -3.00) px against a tile
    median of (0.64, -0.52) px); the median follows the majority. Falls back to the whole-box
    value when fewer than 4 usable tiles remain (tiles < 90 % finite in either image, or pinned at
    the window edge, are dropped). Returns dy, dx and the number of tiles used.
    """
    ny, nx = a.shape
    if tile < 8:
        return (*register(a, b, max_shift), 0)
    offs = []
    for y in range(0, ny - tile + 1, tile):
        for x in range(0, nx - tile + 1, tile):
            ta, tb = a[y : y + tile, x : x + tile], b[y : y + tile, x : x + tile]
            if np.isfinite(ta).mean() < 0.9 or np.isfinite(tb).mean() < 0.9:
                continue  # chip gaps, footprint edges: an empty correlation peaks at a corner
            d = register(ta, tb, max_shift)
            if max(abs(d[0]), abs(d[1])) < np.ceil(max_shift) - 1e-6:
                offs.append(d)  # tiles pinned at the window edge carry no position
    if len(offs) < 4:
        return (*register(a, b, max_shift), 0)
    dy, dx = np.median(np.array(offs), axis=0)
    return float(dy), float(dx), len(offs)


def shift(img: np.ndarray, dy: float, dx: float) -> np.ndarray:
    good = np.isfinite(img)
    out = ndimage.shift(np.where(good, img, 0.0), (dy, dx), order=3, mode="constant", cval=np.nan)
    ok = ndimage.shift(good.astype(float), (dy, dx), order=1, mode="constant", cval=0.0)
    out[ok < 0.999] = np.nan
    return out


def local_shift(a: np.ndarray, b: np.ndarray, tile: int, *others: np.ndarray) -> list[np.ndarray]:
    """``b`` (and ``others``, e.g. its variance) warped onto ``a`` by a smooth tiled shift field.

    A single global shift leaves a positive/negative ridge along every thin arc (sub-pixel
    distortion differences between the epochs' mosaics; first MACS0416 run, 2026-10-15). Shifts
    are measured on half-overlapping tiles, median-filtered over 3 × 3 tiles and interpolated.
    """
    ny, nx = a.shape
    step = max(tile // 2, 8)
    cy = np.arange(tile // 2, ny - tile // 2 + 1, step)
    cx = np.arange(tile // 2, nx - tile // 2 + 1, step)
    if cy.size < 2 or cx.size < 2:
        return [b, *others]
    dyg = np.zeros((cy.size, cx.size))
    dxg = np.zeros_like(dyg)
    for i, y in enumerate(cy):
        for j, x in enumerate(cx):
            sl = (slice(y - tile // 2, y + tile // 2), slice(x - tile // 2, x + tile // 2))
            dyg[i, j], dxg[i, j] = register(a[sl], b[sl], max_shift=2.0)
    dyg = ndimage.median_filter(dyg, size=3, mode="nearest")
    dxg = ndimage.median_filter(dxg, size=3, mode="nearest")
    yy, xx = np.mgrid[0:ny, 0:nx]
    gy = np.interp(np.arange(ny), cy, np.arange(cy.size))
    gx = np.interp(np.arange(nx), cx, np.arange(cx.size))
    gyy, gxx = np.meshgrid(gy, gx, indexing="ij")
    fy = ndimage.map_coordinates(dyg, [gyy, gxx], order=1, mode="nearest")
    fx = ndimage.map_coordinates(dxg, [gyy, gxx], order=1, mode="nearest")
    out = []
    for img in (b, *others):
        good = np.isfinite(img)
        w = ndimage.map_coordinates(
            np.where(good, img, 0.0), [yy - fy, xx - fx], order=3, mode="constant", cval=np.nan
        )
        ok = ndimage.map_coordinates(
            good.astype(float), [yy - fy, xx - fx], order=1, mode="constant", cval=0.0
        )
        w[ok < 0.999] = np.nan
        out.append(w)
    return out


def _blur(img: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return img
    good = np.isfinite(img)
    num = ndimage.gaussian_filter(np.where(good, img, 0.0), sigma)
    den = ndimage.gaussian_filter(good.astype(float), sigma)
    out = num / np.where(den > 0, den, np.nan)
    out[~good] = np.nan
    return out


def _robust_std(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def psf_match(
    e1: np.ndarray, e2: np.ndarray, grid: np.ndarray
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Blur the sharper epoch (and rescale epoch 2) to minimise ``e2 - e1`` on bright pixels.

    ``grid`` holds signed blur sigmas in pixels: > 0 blurs epoch 1, < 0 blurs epoch 2. Returns the
    chosen signed sigma, the flux scale k (epoch 2 is divided by it) and the two matched images.
    A Gaussian kernel only matches the PSF core width; rotated spikes of bright stars stay in the
    difference (vetted by eye).
    """
    # Score on bright pixels only: over the whole box, blurring an epoch also lowers its noise,
    # so the global rms prefers over-blurring and leaves a ridge along every arc (2026-10-15).
    ok = np.isfinite(e1)
    bright = ok & (np.where(ok, e1, -np.inf) > np.percentile(e1[ok], 90.0))
    best = (np.inf, 0.0, 1.0)
    for s in grid:
        a = _blur(e1, s) if s > 0 else e1
        b = _blur(e2, -s) if s < 0 else e2
        k = flux_scale(a, b)
        rms = _robust_std((b / k - a)[bright])
        if rms < best[0]:
            best = (rms, float(s), k)
    _, s, k = best
    return s, k, (_blur(e1, s) if s > 0 else e1), (_blur(e2, -s) if s < 0 else e2) / k


def kernel_fit(
    ref: np.ndarray,
    tgt: np.ndarray,
    var: np.ndarray,
    half: int,
    clip: float = 4.0,
    n_iter: int = 3,
    max_pix: int = 200_000,
    seed: int = 0,
    mask: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    """Least-squares pixel kernel K (plus a constant) with ``K * ref ≈ tgt`` (Bramich 2008 style).

    Each kernel pixel is a free parameter, so K absorbs the PSF difference (core and wings),
    sub-pixel shifts and the flux scale at once. Fitted on the box itself with weights 1 / var;
    pixels with |residual| > ``clip`` σ are dropped and the fit repeated (transients, saturated
    cores, spikes). At most ``max_pix`` pixels enter the fit (a random subset), which keeps the
    ~(2·half+1)² parameters many thousand times over-determined, so a point transient barely
    moves K. ``mask`` (optional) restricts the fit to those pixels (e.g. compact-source stamps).
    Returns K of shape (2·half+1, 2·half+1) and the clipped χ² per fitted pixel; the constant is
    discarded (the median filter of the difference removes it).
    """
    if half < 1:
        raise ValueError("kernel half-width must be >= 1")
    n = 2 * half + 1
    good = np.isfinite(ref) & np.isfinite(tgt) & np.isfinite(var) & (var > 0)
    good[:half, :] = good[-half:, :] = good[:, :half] = good[:, -half:] = False
    if mask is not None:
        good &= mask
    r0 = np.where(np.isfinite(ref), ref, 0.0)
    # every pixel of the kernel footprint must be finite in ref
    fin = ndimage.minimum_filter(np.isfinite(ref).astype(np.uint8), size=n) > 0
    ys, xs = np.nonzero(good & fin)
    if ys.size > max_pix:
        pick = np.random.default_rng(seed).choice(ys.size, max_pix, replace=False)
        ys, xs = ys[pick], xs[pick]
    offs = [(dy, dx) for dy in range(-half, half + 1) for dx in range(-half, half + 1)]
    # convolution: (K * ref)[y, x] = sum K[i, j] ref[y - (i - half), x - (j - half)]
    x_mat = np.empty((ys.size, len(offs) + 1))
    for c, (dy, dx) in enumerate(offs):
        x_mat[:, c] = r0[ys - dy, xs - dx]
    x_mat[:, -1] = 1.0
    t = tgt[ys, xs]
    w = 1.0 / np.sqrt(var[ys, xs])
    if ys.size < 10 * (len(offs) + 1):
        raise ValueError(f"too few pixels ({ys.size}) for a {n}x{n} kernel")
    keep = np.ones(ys.size, bool)
    for it in range(max(n_iter, 1)):
        coef, *_ = np.linalg.lstsq(x_mat[keep] * w[keep, None], t[keep] * w[keep], rcond=None)
        z = (x_mat @ coef - t) * w
        if it == n_iter - 1:  # the returned kernel is the fit on the final mask
            break
        sd = _robust_std(z[keep])
        new = np.abs(z) < clip * max(sd, 1e-12)
        if np.array_equal(new, keep):
            break
        keep = new
    k = np.zeros((n, n))
    for c, (dy, dx) in enumerate(offs):
        k[dy + half, dx + half] = coef[c]
    return k, float(np.mean(z[keep] ** 2))


def compact_stamps(img: np.ndarray, half: int, pct: float = 97.0, radius: int = 6) -> np.ndarray:
    """Pixels within ``radius`` of local maxima of the background-subtracted image above ``pct``.

    Compact sources (stars, knots, galaxy cores) carry the PSF core that a whole-box fit
    under-weights; saturated or masked cores (non-finite within the stamp) are left out.
    """
    fin = np.isfinite(img)
    c = np.where(fin, img, 0.0)
    c = c - ndimage.median_filter(c, size=8 * half + 1)
    mx = (c == ndimage.maximum_filter(c, size=2 * half + 1)) & (c > np.percentile(c[fin], pct))
    disk = np.hypot(*np.mgrid[-radius : radius + 1, -radius : radius + 1]) <= radius
    stamps = ndimage.binary_dilation(mx, structure=disk)
    bad = ndimage.binary_dilation(~fin, structure=disk)
    return stamps & ~bad


def kernel_match(
    e1: np.ndarray,
    e2: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    blur: float,
    half: int,
    stamps: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, int]:
    """Convolve the sharper epoch (sign of the Gaussian ``blur``) with a fitted pixel kernel.

    The kernel is fitted on the whole box. With ``stamps``, pixels around compact sources are
    convolved with a second kernel fitted on the other compact-source stamps. Variances follow
    the stamp kernels inside their stamps. Returns matched epoch 1, matched epoch 2, the
    difference variance, the whole-box kernel, the flux scale (epoch 2 / epoch 1) and the
    convolved epoch (1 or 2).
    """
    var = v1 + v2
    conv1 = blur >= 0
    ref, tgt = (e1, e2) if conv1 else (e2, e1)
    vref, vtgt = (v1, v2) if conv1 else (v2, v1)
    k, _ = kernel_fit(ref, tgt, var, half)
    filled = np.where(np.isfinite(ref), ref, 0.0)
    vfill = np.where(np.isfinite(vref), vref, 0.0)
    m = ndimage.convolve(filled, k, mode="constant")
    vm = ndimage.convolve(vfill, k**2, mode="constant")
    if stamps:
        # No single kernel fits both: the whole-box kernel leaves +10-30 % peak residuals on
        # compact sources in SW, a stamp-only kernel leaves the arc light (2026-10-15). Compact
        # stamps therefore get their own kernel, fitted leaving the stamp itself out, because a
        # transient on a compact source would otherwise train the kernel that subtracts it.
        stamp = compact_stamps(e1, half) & compact_stamps(e2, half)
        lab, nlab = ndimage.label(stamp)
        # a stamp's own pixels always belong to it; the grown rim goes to the highest label
        owner = np.where(lab > 0, lab, ndimage.grey_dilation(lab, size=2 * half + 1))
        for i in range(1, nlab + 1):
            sel = owner == i
            if not sel.any():
                continue
            try:
                ki, _ = kernel_fit(ref, tgt, var, half, mask=stamp & (lab != i))
            except ValueError:  # too few other stamps: keep the whole-box kernel here
                continue
            ys, xs = np.nonzero(sel)
            y0, y1 = max(ys.min() - half, 0), ys.max() + half + 1
            x0, x1 = max(xs.min() - half, 0), xs.max() + half + 1
            win = sel[y0:y1, x0:x1]
            m[y0:y1, x0:x1][win] = ndimage.convolve(filled[y0:y1, x0:x1], ki, mode="constant")[win]
            vm[y0:y1, x0:x1][win] = ndimage.convolve(vfill[y0:y1, x0:x1], ki**2, mode="constant")[
                win
            ]
    m[~np.isfinite(ref)] = np.nan
    var = vm + vtgt
    ksum = float(k.sum())
    if conv1:
        return m, e2, var, k, ksum, 1
    return e1, m, var, k, 1.0 / ksum, 2


def flux_scale(a: np.ndarray, b: np.ndarray) -> float:
    """Robust ratio b / a on the brightest 10 % of pixels (calibration or pipeline differences).

    Without it a 1 % zero-point difference leaves every bright arc and galaxy core in the
    difference as a same-sign residual in every band (first MACS0416 run, 2026-10-15).
    """
    good = np.isfinite(a) & np.isfinite(b)
    lo, hi = np.percentile(a[good], [90.0, 99.9])
    sel = good & (a > lo) & (a < hi) & (lo > 0)
    if sel.sum() < 50:
        return 1.0
    return float(np.median(b[sel] / a[sel]))


def significance(
    diff: np.ndarray, var: np.ndarray, sigma_px: float, bg_px: int
) -> tuple[np.ndarray, np.ndarray]:
    """Matched-filter significance (uncalibrated) and the background-subtracted difference."""
    good = np.isfinite(diff) & np.isfinite(var)
    d = np.where(good, diff, 0.0)
    bg = ndimage.median_filter(d, size=bg_px)
    d = np.where(good, d - bg, 0.0)
    r = int(np.ceil(3 * sigma_px))
    yy, xx = np.mgrid[-r : r + 1, -r : r + 1]
    k = np.exp(-0.5 * (xx**2 + yy**2) / sigma_px**2)
    k /= k.sum()
    num = ndimage.convolve(d, k, mode="constant")
    den = np.sqrt(ndimage.convolve(np.where(good, var, 0.0), k**2, mode="constant"))
    cover = ndimage.convolve(good.astype(float), np.ones_like(k), mode="constant") / k.size
    s = num / np.where(den > 0, den, np.nan)
    s[cover < 0.999] = np.nan
    d[~good] = np.nan
    return s, d


def peaks(s: np.ndarray, min_sigma: float, min_sep_px: int) -> list[tuple[int, int, float]]:
    """Local maxima of |s| ≥ ``min_sigma`` (pixel y, x and signed s)."""
    a = np.nan_to_num(np.abs(s))
    mx = ndimage.maximum_filter(a, size=2 * min_sep_px + 1)
    ys, xs = np.nonzero((a == mx) & (a >= min_sigma))
    return [(int(y), int(x), float(s[y, x])) for y, x in zip(ys, xs, strict=True)]


def gaussian_star(shape: tuple[int, int], y: float, x: float, sigma: float, flux: float):
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / sigma**2)
    return flux * g / g.sum()


def mag_to_mjysr(mag: float, pixar_sr: float) -> float:
    """Total flux of an AB magnitude in MJy/sr × pixel (sum of pixel values)."""
    jy = 3631.0 * 10 ** (-0.4 * mag)
    return jy / 1e6 / pixar_sr


def add_star(
    img: np.ndarray, y: float, x: float, sigma: float, flux: float, kern=None, pad: int = 0
) -> None:
    """Add a Gaussian point source (optionally convolved with ``kern``) to ``img`` in place."""
    w = int(np.ceil(4 * sigma)) + pad + (kern.shape[0] // 2 if kern is not None else 0) + 2
    yi, xi = int(round(y)), int(round(x))
    sl = (slice(max(yi - w, 0), yi + w + 1), slice(max(xi - w, 0), xi + w + 1))
    star = gaussian_star(img[sl].shape, y - sl[0].start, x - sl[1].start, sigma, flux)
    if kern is not None:
        star = ndimage.convolve(star, kern, mode="constant")
    img[sl] += star


def empirical_psfs(
    b1: dict, b2: dict, fwhm_arcsec: float, half_arcsec: float, n: int = 5, min_snr: float = 50.0
) -> tuple[np.ndarray, np.ndarray, list[tuple[float, float, float]]]:
    """Per-epoch empirical point-source stamps (sum 1) from the same point-like sources.

    Sources are chosen on epoch 1 only: background-subtracted local maxima with peak / ERR ≥
    ``min_snr``, no brighter maximum within the stamp, finite stamps in both epochs, ranked by
    core concentration (peak / sum within ±r, r = max(3, 2 FWHM) px). Concentrations above that of
    the band-FWHM Gaussian (hot pixels, cosmic rays) are dropped; the ``n`` most concentrated are
    kept. Each epoch's stamps are cut at the same sky positions on that epoch's grid, recentred on
    their own centroid, normalised and median-stacked. ASSUMPTION: these sources are stars or
    unresolved knots and carry each epoch's PSF; wings beyond the stamp are truncated.
    """

    def bgsub(b):
        fin = np.isfinite(b["sci"])
        c = np.where(fin, b["sci"], 0.0)
        k = int(round(4 * fwhm_arcsec / b["scale"])) * 2 + 1
        return c - ndimage.median_filter(c, size=max(k, 9)), fin

    c1, f1 = bgsub(b1)
    c2, f2 = bgsub(b2)
    h1 = int(round(half_arcsec / b1["scale"]))
    r = max(3, int(round(2 * fwhm_arcsec / b1["scale"])))
    g = gaussian_star((2 * r + 1, 2 * r + 1), r, r, fwhm_arcsec / 2.3548 / b1["scale"], 1.0)
    conc_max = float(g.max())
    mx = c1 == ndimage.maximum_filter(c1, size=2 * h1 + 1)
    snr = c1 / np.where(b1["err"] > 0, b1["err"], np.inf)
    cand = []
    for y, x in zip(*np.nonzero(mx & (snr >= min_snr) & f1), strict=True):
        if not (h1 <= y < c1.shape[0] - h1 and h1 <= x < c1.shape[1] - h1):
            continue
        if not f1[y - h1 : y + h1 + 1, x - h1 : x + h1 + 1].all():
            continue
        conc = float(c1[y, x] / c1[y - r : y + r + 1, x - r : x + r + 1].sum())
        if conc <= conc_max:
            ra, dec = b1["wcs"].all_pix2world(x, y, 0)
            cand.append((-conc, float(ra), float(dec)))
    cand.sort()

    def stamp(c, fin, b, ra, dec):
        x, y = (float(v) for v in b["wcs"].all_world2pix(ra, dec, 0))
        h = int(round(half_arcsec / b["scale"]))
        yi, xi = int(round(y)), int(round(x))
        if not (h <= yi < c.shape[0] - h and h <= xi < c.shape[1] - h):
            return None
        sl = (slice(yi - h, yi + h + 1), slice(xi - h, xi + h + 1))
        if not fin[sl].all():
            return None
        st = c[sl]
        cy, cx = ndimage.center_of_mass(np.clip(st, 0.0, None))
        st = ndimage.shift(st, (h - cy, h - cx), order=3, mode="constant")
        return st / st.sum()

    s1, s2, used = [], [], []
    for nc, ra, dec in cand:
        a, b = stamp(c1, f1, b1, ra, dec), stamp(c2, f2, b2, ra, dec)
        if a is not None and b is not None:
            s1.append(a)
            s2.append(b)
            used.append((ra, dec, round(-nc, 4)))
        if len(used) == n:
            break
    if not used:
        raise ValueError("no compact sources for the empirical PSF")
    p1, p2 = np.median(s1, axis=0), np.median(s2, axis=0)
    return p1 / p1.sum(), p2 / p2.sum(), used


def add_stamp(img: np.ndarray, y: float, x: float, psf: np.ndarray, flux: float) -> None:
    """Add ``flux`` × ``psf`` (sum 1, centred) to ``img`` at sub-pixel (y, x), in place."""
    h = psf.shape[0] // 2
    yi, xi = int(round(y)), int(round(x))
    st = flux * ndimage.shift(psf, (y - yi, x - xi), order=3, mode="constant")
    y0, x0 = yi - h, xi - h
    ys, xs = max(y0, 0), max(x0, 0)
    ye, xe = min(y0 + psf.shape[0], img.shape[0]), min(x0 + psf.shape[1], img.shape[1])
    if ye > ys and xe > xs:
        img[ys:ye, xs:xe] += st[ys - y0 : ye - y0, xs - x0 : xe - x0]


def grid_positions(
    img: np.ndarray, w: WCS, scale: float, n: int, min_sep: float, margin: float, pct: float, seed
) -> list[tuple[float, float]]:
    """Up to ``n`` sky positions on pixels above the ``pct`` percentile of ``img`` (arcs and
    galaxies), at least ``min_sep`` arcsec apart and ``margin`` arcsec from the box border."""
    rng = np.random.default_rng(seed)
    m = max(int(np.ceil(margin / scale)), 1)  # ok[-0:] would clear the whole mask
    ok = np.nan_to_num(img) > np.nanpercentile(img, pct)
    ok[:m, :] = ok[-m:, :] = ok[:, :m] = ok[:, -m:] = False
    ys, xs = np.nonzero(ok)
    chosen: list[tuple[float, float]] = []
    for i in rng.permutation(ys.size):
        y, x = ys[i] + rng.uniform(-0.5, 0.5), xs[i] + rng.uniform(-0.5, 0.5)
        if all(np.hypot(y - cy, x - cx) * scale >= min_sep for cy, cx in chosen):
            chosen.append((y, x))
            if len(chosen) == n:
                break
    return [tuple(float(v) for v in w.all_pix2world(x, y, 0)) for y, x in chosen]


def grid_job(uri, ra, dec, half, n, min_sep, margin, seed) -> list[tuple[float, float]]:
    """``grid_positions`` on the epoch-1 box of ``uri`` (bright = above the 90th percentile)."""
    b = read_box(uri, ra, dec, half)
    return grid_positions(b["sci"], b["wcs"], b["scale"], n, min_sep, margin, 90.0, seed)


PSF_N = 5  # ASSUMPTION: empirical PSF stack size (PLAN addendum 3)
PSF_MIN_SNR = 50.0


def global_shift(a: np.ndarray, b: np.ndarray, scale: float, reg_tile_arcsec: float):
    """The chain's global registration: tile median (``--reg-tile-arcsec`` > 0) or whole box."""
    if reg_tile_arcsec > 0:
        return register_median(a, b, int(round(reg_tile_arcsec / scale)))
    return (*register(a, b), 0)


def pair_shift(
    e1: str, e2: str, ra: float, dec: float, half: float, reg_tile_arcsec: float, cache_dir=None
):
    """Global registration (dy, dx) of the injection-free boxes, as the chain measures it.

    Known limitation (final review 2026-10-17): the chain re-measures the registration with the
    injected sources present, which moves it by 0.02-0.15 px (F150W) from this value
    (``results/transients/fade_calib/inject_shift_residual.json``), so injected sources keep that
    small misregistration."""
    b1 = box_from_raw(cached(cache_dir, read_box_raw, e1, ra, dec, half))
    b2 = box_from_raw(cached(cache_dir, read_box_raw, e2, ra, dec, half + 1.0))
    s2 = resample(b2["sci"], b2["wcs"], b1["wcs"], b1["sci"].shape)
    dy, dx, _ = global_shift(b1["sci"], s2, b1["scale"], reg_tile_arcsec)
    return dy, dx


def pair_psfs(e1, e2, ra, dec, box, half, fwhm, n, min_snr):
    """Per-epoch empirical PSFs from ``box``-half boxes of the two mosaics (``empirical_psfs``)."""
    q1 = read_box(e1, ra, dec, box)
    q2 = read_box(e2, ra, dec, box + 1.0)
    return empirical_psfs(q1, q2, fwhm, half, n=n, min_snr=min_snr)


def epoch2_pixel(w1: WCS, w2: WCS, ra: float, dec: float, dy: float, dx: float):
    """(y, x) in epoch 2 of a sky-fixed source at (ra, dec), given the chain's registration.

    Real sources sit at epoch-1 pixel p - (dy, dx) in WCS-resampled epoch 2 (``shift`` then moves
    them by +(dy, dx) onto p), so the source is placed at the sky position of that pixel.
    """
    x1, y1 = (float(v) for v in w1.all_world2pix(ra, dec, 0))
    ra2, dec2 = w1.all_pix2world(x1 - dx, y1 - dy, 0)
    x2, y2 = (float(v) for v in w2.all_world2pix(ra2, dec2, 0))
    return y2, x2


def band_diff(job: dict) -> dict:
    """Read, align, match and difference one band; detections, maps and injections."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    t0 = time.perf_counter()
    cd = job.get("cache_dir")
    b1 = box_from_raw(cached(cd, read_box_raw, job["e1"], job["ra"], job["dec"], job["half"]))
    b2 = box_from_raw(cached(cd, read_box_raw, job["e2"], job["ra"], job["dec"], job["half"] + 1.0))
    t_read = time.perf_counter() - t0
    # Pre-matching injections (``--inject-pre``): sources added to the raw boxes of each epoch,
    # on that epoch's own grid and flux scale, before resampling, registration and matching.
    # With ``inj_pre_psf == "empirical"`` each epoch gets its own empirical PSF (chosen before
    # any injection), so an unchanged source has each epoch's real shape.
    psfs, psf_info, pre_shift = None, None, None
    if job.get("inj_pre"):
        # Where a real sky-fixed source lands in epoch 2: the chain shifts resampled epoch 2 by
        # the registration (dy, dx), so a source at epoch-1 pixel p sits at p - (dy, dx) there.
        # Placing it at the WCS position instead would leave every injection misregistered by
        # (dy, dx) after the shift. The shift is measured on the injection-free boxes.
        pre_shift = cached(
            cd,
            pair_shift,
            job["e1"],
            job["e2"],
            job["ra"],
            job["dec"],
            job["half"],
            job.get("reg_tile_arcsec", 0.0),
            cd,
        )
    if job.get("inj_pre") and job.get("inj_pre_psf") == "empirical":
        fw = PSF_FWHM[job["band"]]
        box = job["inj_pre_psf_box"]
        p1, p2, used = cached(
            cd,
            pair_psfs,
            job["e1"],
            job["e2"],
            job["ra"],
            job["dec"],
            box,
            job["inj_pre_psf_half"],
            fw,
            PSF_N,
            PSF_MIN_SNR,
        )
        psfs = {1: p1, 2: p2}
        psf_info = {
            "sources_ra_dec_conc": used,
            "half_arcsec": job["inj_pre_psf_half"],
            "box_half_arcsec": box,
            "peak_frac": [round(float(p1.max()), 4), round(float(p2.max()), 4)],
        }
    for ra, dec, mag1, mag2 in job.get("inj_pre") or []:
        for ep, b, mag in ((1, b1, mag1), (2, b2, mag2)):
            if mag is None:
                continue
            if ep == 1:
                x, y = (float(v) for v in b1["wcs"].all_world2pix(ra, dec, 0))
            else:
                y, x = epoch2_pixel(b1["wcs"], b2["wcs"], ra, dec, *pre_shift)
            if 0 <= y < b["sci"].shape[0] and 0 <= x < b["sci"].shape[1]:
                flux = mag_to_mjysr(mag, b["pixar_sr"])
                if psfs is not None:
                    add_stamp(b["sci"], y, x, psfs[ep], flux)
                else:
                    sig = PSF_FWHM[job["band"]] / 2.3548 / b["scale"]
                    add_star(b["sci"], y, x, sig, flux)
    t1 = time.perf_counter()
    t_pre = t1 - t0 - t_read
    shape = b1["sci"].shape
    s2 = resample(b2["sci"], b2["wcs"], b1["wcs"], shape)
    v2 = resample(b2["err"] ** 2, b2["wcs"], b1["wcs"], shape)
    dy, dx, n_reg_tiles = global_shift(b1["sci"], s2, b1["scale"], job.get("reg_tile_arcsec", 0.0))
    s2 = shift(s2, dy, dx)
    v2 = shift(v2, dy, dx)
    if job["tile_arcsec"] > 0:
        tile = int(round(job["tile_arcsec"] / b1["scale"]))
        s2, v2 = local_shift(b1["sci"], s2, tile, v2)
    sig_psf = PSF_FWHM[job["band"]] / 2.3548 / b1["scale"]
    grid = np.round(np.linspace(-3.0, 3.0, 25) * sig_psf, 3)
    blur, kscale, m1, m2 = psf_match(b1["sci"], s2, grid)
    var = b1["err"] ** 2 + np.clip(v2, 0.0, None) / kscale**2
    kern = None
    conv2 = blur < 0
    if job["kernel_px"] > 0:
        m1, m2, var, kern, kscale, conv_epoch = kernel_match(
            b1["sci"],
            s2,
            b1["err"] ** 2,
            np.clip(v2, 0.0, None),
            blur,
            job["kernel_px"],
            job["kernel_stamps"],
        )
        conv2 = conv_epoch == 2
    bg_px = int(np.ceil(job["bg_arcsec"] / b1["scale"])) | 1
    sep = max(int(round(2 * sig_psf)), 1)
    t_match = time.perf_counter() - t1

    # Whole-chain injection grid: stars at common sky positions in every band, added to epoch 2
    # before detection, so hits, maps and any later stage (arc_diff_morph.py) see them.
    if job.get("inj_sky"):
        m2 = m2.copy()
        sig_m2 = sig_psf if (kern is not None or blur >= 0) else float(np.hypot(sig_psf, -blur))
        for ra, dec, mag in job["inj_sky"]:
            x, y = (float(v) for v in b1["wcs"].all_world2pix(ra, dec, 0))
            if 0 <= y < shape[0] and 0 <= x < shape[1]:
                add_star(
                    m2,
                    y,
                    x,
                    sig_m2,
                    mag_to_mjysr(mag, b1["pixar_sr"]),
                    kern if (kern is not None and conv2) else None,
                )

    def detect(e2img: np.ndarray, scale: float | None = None):
        sraw, d = significance(e2img - m1, var, sig_psf, bg_px)
        sc = scale if scale is not None else max(_robust_std(sraw), 1.0)
        return sraw / sc, d, sc

    t2 = time.perf_counter()
    s, d, noise_scale = detect(m2)
    det = []
    for y, x, v in peaks(s, job["min_sigma"], sep):
        ra, dec = b1["wcs"].all_pix2world(x, y, 0)
        det.append(
            {"band": job["band"], "ra": float(ra), "dec": float(dec), "x": x, "y": y, "sigma": v}
        )
    out = {
        "band": job["band"],
        "e1": job["e1"].rsplit("/", 1)[-1],
        "e2": job["e2"].rsplit("/", 1)[-1],
        "scale_arcsec": b1["scale"],
        "shift_px": [dy, dx],
        "blur_sigma_px": blur,
        "kernel_sum": None if kern is None else float(kern.sum()),
        "convolved_epoch": 2 if conv2 else 1,
        "flux_scale": kscale,
        "noise_scale": noise_scale,
        "detections": det,
        "s": s.astype(np.float32),
        "d": d.astype(np.float32),
        "e1img": b1["sci"].astype(np.float32),
        "wcs": b1["wcs"].to_header_string(relax=True),
    }
    t_detect = time.perf_counter() - t2
    t3 = time.perf_counter()
    # Injections: epoch 2 only, on bright pixels (arcs); same chain, same noise scale.
    rng = np.random.default_rng(job["seed"])
    bright = np.nan_to_num(b1["sci"]) > np.nanpercentile(b1["sci"], 90)
    edge = int(np.ceil(bg_px / 2)) + 3
    bright[:edge, :] = bright[-edge:, :] = bright[:, :edge] = bright[:, -edge:] = False
    ys, xs = np.nonzero(bright)
    inj = {}
    inj_unflagged = {}
    core_thr = np.nanpercentile(b1["sci"], job["core_pct"])
    core_r = int(np.ceil(job["core_arcsec"] / b1["scale"]))
    if job["inject"] and ys.size:
        sig_m2 = sig_psf if blur >= 0 else float(np.hypot(sig_psf, -blur))
        for mag in job["inject_mags"]:
            found = 0
            found_unflagged = 0
            for _ in range(job["inject"]):
                i = rng.integers(ys.size)
                y, x = ys[i] + rng.uniform(-0.5, 0.5), xs[i] + rng.uniform(-0.5, 0.5)
                # Local window: the star changes S only within the filter + median footprint.
                yi, xi = int(round(y)), int(round(x))
                w = bg_px + int(np.ceil(3 * sig_psf)) + 4
                sl = (slice(max(yi - w, 0), yi + w + 1), slice(max(xi - w, 0), xi + w + 1))
                star = gaussian_star(
                    m2[sl].shape,
                    y - sl[0].start,
                    x - sl[1].start,
                    sig_psf if kern is not None else sig_m2,
                    mag_to_mjysr(mag, b1["pixar_sr"]),
                )
                if kern is not None and conv2:  # epoch 2 is the convolved side
                    star = ndimage.convolve(star, kern, mode="constant")
                si, _ = significance(m2[sl] + star - m1[sl], var[sl], sig_psf, bg_px)
                si = si / noise_scale
                r = sep + 1
                cy, cx = yi - sl[0].start, xi - sl[1].start
                win = si[max(cy - r, 0) : cy + r + 1, max(cx - r, 0) : cx + r + 1]
                if np.nanmax(win, initial=-np.inf) >= job["min_sigma"]:
                    found += 1
                    core = b1["sci"][
                        max(yi - core_r, 0) : yi + core_r + 1, max(xi - core_r, 0) : xi + core_r + 1
                    ]
                    if not np.nanmax(core, initial=-np.inf) >= core_thr:
                        found_unflagged += 1
            inj[f"{mag:.2f}"] = found / job["inject"]
            inj_unflagged[f"{mag:.2f}"] = found_unflagged / job["inject"]
    out["injection"] = inj
    out["injection_unflagged"] = inj_unflagged
    out["inject_pre_psf"] = psf_info
    out["inject_pre_shift_px"] = list(pre_shift) if pre_shift is not None else None
    out["reg_tiles_used"] = n_reg_tiles
    out["timing_s"] = {
        "read": round(t_read, 2),
        "pre_injections": round(t_pre, 2),
        "resample_register_match": round(t_match, 2),
        "detect": round(t_detect, 2),
        "post_injections": round(time.perf_counter() - t3, 2),
    }
    return out


def group(dets: list[dict], link_arcsec: float) -> list[list[dict]]:
    """Friends-of-friends groups of detections across bands."""
    n = len(dets)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    ra = np.array([d["ra"] for d in dets])
    dec = np.array([d["dec"] for d in dets])
    for i in range(n):
        dra = (ra[i + 1 :] - ra[i]) * np.cos(np.radians(dec[i])) * 3600
        ddec = (dec[i + 1 :] - dec[i]) * 3600
        for j in np.nonzero(np.hypot(dra, ddec) <= link_arcsec)[0]:
            parent[find(i)] = find(i + 1 + int(j))
    groups: dict[int, list[dict]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(dets[i])
    return list(groups.values())


def match_targets(targets: list[dict], res: list[dict], hits: list[dict], r_arcsec: float):
    """Per target: the signed S of largest |S| within ``r_arcsec`` in each band, and the nearest
    hit within ``r_arcsec`` (its id, or None)."""
    out = []
    for k in targets:
        entry = {kk: k[kk] for kk in k if kk not in ("ra", "dec")}
        for r in res:
            w = _wcs(fits.Header.fromstring(r["wcs"]))
            x, y = w.all_world2pix(k["ra"], k["dec"], 0)
            xi, yi = int(round(float(x))), int(round(float(y)))
            rr = int(np.ceil(r_arcsec / r["scale_arcsec"]))
            s = r["s"]
            win = s[max(yi - rr, 0) : yi + rr + 1, max(xi - rr, 0) : xi + rr + 1]
            v = win[np.isfinite(win)]
            entry[r["band"]] = float(v[np.argmax(np.abs(v))]) if v.size else None
        best = None
        for h in hits:
            sep = (
                np.hypot((h["ra"] - k["ra"]) * np.cos(np.radians(k["dec"])), h["dec"] - k["dec"])
                * 3600
            )
            if sep <= r_arcsec and (best is None or sep < best[1]):
                best = (h["hit_id"], sep)
        entry["hit"] = best[0] if best else None
        out.append(entry)
    return out


def load_known(path: Path, offset: tuple[float, float]) -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r["name"].startswith("#"):
                continue
            dec = float(r["dec"]) + offset[1] / 3600
            ra = float(r["ra"]) + offset[0] / 3600 / np.cos(np.radians(dec))
            rows.append({**r, "ra": ra, "dec": dec})
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--band", nargs=3, action="append", required=True, metavar=("BAND", "E1", "E2"))
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--half-arcsec", type=float, default=12.0)
    ap.add_argument("--bg-arcsec", type=float, default=0.6, help="median filter of the difference")
    ap.add_argument("--tile-arcsec", type=float, default=3.0, help="local registration; 0 = off")
    ap.add_argument(
        "--reg-tile-arcsec",
        type=float,
        default=0.0,
        help="global shift = median over tiles of this size (0 = whole-box correlation)",
    )
    ap.add_argument(
        "--kernel-px", type=int, default=0, help="half-width of a fitted pixel kernel; 0 = Gaussian"
    )
    ap.add_argument(
        "--kernel-stamps", action="store_true", help="fit the kernel on compact-source stamps"
    )
    ap.add_argument("--min-sigma", type=float, default=5.0)
    ap.add_argument("--core-arcsec", type=float, default=0.8, help="bright-core flag radius")
    ap.add_argument("--core-pct", type=float, default=99.9, help="bright-core flag level")
    ap.add_argument("--min-bands", type=int, default=2)
    ap.add_argument("--link-arcsec", type=float, default=0.1)
    ap.add_argument("--known", type=Path)
    ap.add_argument("--frame-offset", type=float, nargs=2, default=[0.0, 0.0])
    ap.add_argument("--match-arcsec", type=float, default=0.15)
    ap.add_argument("--inject", type=int, default=0)
    ap.add_argument("--inject-mags", type=float, nargs="+", default=[27.5, 28.0, 28.5])
    ap.add_argument(
        "--inject-grid",
        type=int,
        default=0,
        help="N stars at common sky positions in every band, through the whole chain",
    )
    ap.add_argument("--inject-grid-mags", type=float, nargs="+", default=[27.5, 28.0])
    ap.add_argument("--inject-grid-sep", type=float, default=1.5, help="arcsec between stars")
    ap.add_argument(
        "--inject-pre",
        type=Path,
        help="CSV name,ra,dec,mag1,mag2: sources added to both raw epochs before matching",
    )
    ap.add_argument(
        "--inject-pre-psf",
        choices=("gauss", "empirical"),
        default="gauss",
        help="shape of --inject-pre sources: band-FWHM Gaussian, or per-epoch empirical PSF",
    )
    ap.add_argument(
        "--inject-pre-psf-half", type=float, default=0.6, help="stamp half-size (arcsec)"
    )
    ap.add_argument(
        "--inject-pre-psf-box", type=float, default=40.0, help="PSF-star box half-size (arcsec)"
    )
    ap.add_argument(
        "--cache-dir",
        type=Path,
        help="local cache of raw boxes and empirical PSFs (untracked; never injection-dependent)",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--procs", type=int, default=min(os.cpu_count() or 1, 4))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--maps", type=Path, help="S, difference and epoch-1 maps (.npz; untracked)")
    a = ap.parse_args(argv)
    t_start = time.perf_counter()
    pre_rows = []
    if a.inject_pre:
        with open(a.inject_pre, newline="") as f:
            pre_rows = list(csv.DictReader(f))

    def _mag(v):
        return float(v) if v not in (None, "") else None

    inj_pre = [
        (float(r["ra"]), float(r["dec"]), _mag(r["mag1"]), _mag(r["mag2"])) for r in pre_rows
    ]

    jobs = [
        {
            "band": b.upper(),
            "e1": e1,
            "e2": e2,
            "ra": a.ra,
            "dec": a.dec,
            "half": a.half_arcsec,
            "bg_arcsec": a.bg_arcsec,
            "tile_arcsec": a.tile_arcsec,
            "kernel_px": a.kernel_px,
            "kernel_stamps": a.kernel_stamps,
            "core_arcsec": a.core_arcsec,
            "core_pct": a.core_pct,
            "min_sigma": a.min_sigma,
            "inject": a.inject,
            "inject_mags": a.inject_mags,
            "seed": a.seed + i,
            "inj_sky": [],
            "inj_pre": inj_pre,
            "inj_pre_psf": a.inject_pre_psf,
            "inj_pre_psf_half": a.inject_pre_psf_half,
            "inj_pre_psf_box": a.inject_pre_psf_box,
            "reg_tile_arcsec": a.reg_tile_arcsec,
            "cache_dir": str(a.cache_dir) if a.cache_dir else None,
        }
        for i, (b, e1, e2) in enumerate(a.band)
    ]
    inj_sky: list[tuple[float, float, float]] = []
    with ProcessPoolExecutor(max_workers=a.procs) as ex:
        if a.inject_grid:
            # Read in a worker: a remote read in the parent leaves an event loop that forked
            # workers inherit.
            pos = ex.submit(
                grid_job,
                a.band[0][1],
                a.ra,
                a.dec,
                a.half_arcsec,
                a.inject_grid,
                a.inject_grid_sep,
                a.bg_arcsec + 0.5,
                a.seed,
            ).result()
            mags = a.inject_grid_mags
            inj_sky = [(ra, dec, mags[i % len(mags)]) for i, (ra, dec) in enumerate(pos)]
            for j in jobs:
                j["inj_sky"] = inj_sky
        res = list(ex.map(band_diff, jobs))
    a.out.mkdir(parents=True, exist_ok=True)
    dets = [d for r in res for d in r["detections"]]
    groups = group(dets, a.link_arcsec)
    hits = []
    for g in groups:
        bands = sorted({d["band"] for d in g})
        if len(bands) < a.min_bands:
            continue
        signs = {np.sign(d["sigma"]) for d in g}
        best = max(g, key=lambda d: abs(d["sigma"]))
        hits.append(
            {
                "ra": float(np.mean([d["ra"] for d in g])),
                "dec": float(np.mean([d["dec"] for d in g])),
                "n_bands": len(bands),
                "bands": ",".join(bands),
                "max_abs_sigma": abs(best["sigma"]),
                "sign": "mixed" if len(signs) > 1 else ("+" if signs == {1.0} else "-"),
            }
        )
    # Bright-core flag (ASSUMPTION: radius and level): PSF rotation between epochs leaves rosettes
    # around stars and compact galaxy cores that no single kernel removes (2026-10-15).
    core_ref = [
        (
            _wcs(fits.Header.fromstring(r["wcs"])),
            r["e1img"],
            np.nanpercentile(r["e1img"], a.core_pct),
            int(np.ceil(a.core_arcsec / r["scale_arcsec"])),
        )
        for r in res
    ]
    for h in hits:
        h["near_bright_core"] = False
        for w, img, thr, rr in core_ref:
            x, y = w.all_world2pix(h["ra"], h["dec"], 0)
            xi, yi = int(round(float(x))), int(round(float(y)))
            win = img[max(yi - rr, 0) : yi + rr + 1, max(xi - rr, 0) : xi + rr + 1]
            if np.nanmax(win, initial=-np.inf) >= thr:
                h["near_bright_core"] = True
                break
    hits.sort(key=lambda h: -h["max_abs_sigma"])
    for i, h in enumerate(hits):
        h["hit_id"] = f"h{i:03d}"
    known_out = []
    if a.known:
        known_out = match_targets(
            [
                {"name": k["name"], "pair": k.get("pair", ""), "ra": k["ra"], "dec": k["dec"]}
                for k in load_known(a.known, tuple(a.frame_offset))
            ],
            res,
            hits,
            a.match_arcsec,
        )
    injected = []
    if inj_sky:
        flagged = {h["hit_id"] for h in hits if h["near_bright_core"]}
        injected = match_targets(
            [
                {"name": f"inj{i:03d}", "ra": ra, "dec": dec, "mag": mag}
                for i, (ra, dec, mag) in enumerate(inj_sky)
            ],
            res,
            hits,
            a.match_arcsec,
        )
        for e, (ra, dec, _) in zip(injected, inj_sky, strict=True):
            e["near_bright_core"] = e["hit"] in flagged
            e["ra"], e["dec"] = ra, dec
        with open(a.out / "injected.csv", "w", newline="") as f:
            wr = csv.DictWriter(
                f, fieldnames=["name", "ra", "dec", "mag", "hit", "near_bright_core"]
            )
            wr.writeheader()
            for e in injected:
                wr.writerow({k: e[k] for k in wr.fieldnames})
    if pre_rows:
        flagged = {h["hit_id"] for h in hits if h["near_bright_core"]}
        hit_sign = {h["hit_id"]: h["sign"] for h in hits}
        pre = match_targets(
            [{"name": r["name"], "ra": float(r["ra"]), "dec": float(r["dec"])} for r in pre_rows],
            res,
            hits,
            a.match_arcsec,
        )
        fields = list(pre_rows[0].keys()) + ["hit", "hit_sign", "near_bright_core"]
        fields += [r["band"] for r in res]
        with open(a.out / "injected_pre.csv", "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=fields)
            wr.writeheader()
            for row, e in zip(pre_rows, pre, strict=True):
                out = dict(row)
                out["hit"] = e["hit"] or ""
                out["hit_sign"] = hit_sign.get(e["hit"], "")
                out["near_bright_core"] = e["hit"] in flagged
                for r in res:
                    out[r["band"]] = e[r["band"]]
                wr.writerow(out)
    tab = Table(rows=hits, names=list(hits[0].keys())) if hits else Table()
    tab.meta = {"provenance": "derived", "thresholds": "ASSUMPTION: see summary.json"}
    tab.write(a.out / "hits.ecsv", overwrite=True)
    if a.maps:
        np.savez_compressed(
            a.maps,
            **{f"{r['band']}_{k}": r[k] for r in res for k in ("s", "d", "e1img")},
            **{f"{r['band']}_wcs": np.array(r["wcs"]) for r in res},
        )
    summary = {
        "command": "arc_diff.py " + " ".join(argv if argv is not None else os.sys.argv[1:]),
        "box": {"ra": a.ra, "dec": a.dec, "half_arcsec": a.half_arcsec},
        "thresholds_assumption": {
            "min_sigma": a.min_sigma,
            "min_bands": a.min_bands,
            "link_arcsec": a.link_arcsec,
            "bg_arcsec": a.bg_arcsec,
            "tile_arcsec": a.tile_arcsec,
            "reg_tile_arcsec": a.reg_tile_arcsec,
            "inject_pre_psf": a.inject_pre_psf,
            "kernel_px": a.kernel_px,
            "kernel_stamps": a.kernel_stamps,
            "match_arcsec": a.match_arcsec,
            "core_arcsec": a.core_arcsec,
            "core_pct": a.core_pct,
            "psf_fwhm_arcsec": {r["band"]: PSF_FWHM[r["band"]] for r in res},
        },
        "bands": {
            r["band"]: {
                k: r[k]
                for k in (
                    "e1",
                    "e2",
                    "scale_arcsec",
                    "shift_px",
                    "blur_sigma_px",
                    "kernel_sum",
                    "convolved_epoch",
                    "flux_scale",
                    "noise_scale",
                )
            }
            | {
                "n_detections": len(r["detections"]),
                "injection_recovered": r["injection"],
                "injection_recovered_not_near_core_this_band": r["injection_unflagged"],
                "timing_s": r["timing_s"],
                "inject_pre_psf": r.get("inject_pre_psf"),
                "inject_pre_shift_px": r.get("inject_pre_shift_px"),
                "reg_tiles_used": r.get("reg_tiles_used"),
            }
            for r in res
        },
        "n_hits": len(hits),
        "n_hits_not_near_bright_core": sum(not h["near_bright_core"] for h in hits),
        "known": known_out,
        "inject_pre": {"n": len(pre_rows), "csv": str(a.inject_pre) if a.inject_pre else None},
        "wall_s": round(time.perf_counter() - t_start, 1),
        "inject_grid": {
            "n": len(injected),
            "n_hit": sum(e["hit"] is not None for e in injected),
            "n_hit_not_near_bright_core": sum(
                e["hit"] is not None and not e["near_bright_core"] for e in injected
            ),
        },
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in ("n_hits", "n_hits_not_near_bright_core", "known", "inject_grid")
            },
            indent=1,
        )
    )
    print(json.dumps(summary["bands"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
