"""W1 blank-field responsivity R by metacalibration of real F200W images.

Plan and rules (frozen first): ``results/w1_r_realmetacal/PLAN.md``. ``validate`` runs the
simulated gate (V1 metacal vs direct shearing, V2 sign); ``real`` streams 1024-px tiles of a JADES
GOODS-N F200W i2d from MAST S3 and measures R with the same chain. Shapes follow
``w1_r_sim.measure`` (kernel 2 px, 3σ smoothed-noise threshold, 5 px, no deblending) on whole
scenes; ε from ``exotic_screens.intrinsic_ellipticity``.

    python scripts/w1_r_realmetacal.py validate --out results/w1_r_realmetacal/validate.json
    python scripts/w1_r_realmetacal.py real --out results/w1_r_realmetacal/r_real.json
    python scripts/w1_r_realmetacal.py posthoc --r-real results/w1_r_realmetacal/r_real.json \
        --out results/w1_r_realmetacal/r_real_posthoc_excl_tile.json
    python scripts/w1_r_realmetacal.py sheet --out results/w1_r_realmetacal/sheet.png
    python scripts/w1_r_realmetacal.py validate2 --out results/w1_r_realmetacal/v2/validate.json
    python scripts/w1_r_realmetacal.py real2 --out results/w1_r_realmetacal/v2/r_real.json
    python scripts/w1_r_realmetacal.py epsf3 --out results/w1_r_realmetacal/v3/epsf.json
    python scripts/w1_r_realmetacal.py validate3 --epsf results/w1_r_realmetacal/v3/epsf.json \
        --out results/w1_r_realmetacal/v3/validate.json
    python scripts/w1_r_realmetacal.py real3 --epsf results/w1_r_realmetacal/v3/epsf.json \
        --out results/w1_r_realmetacal/v3/r_real.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exotic_screens as es  # noqa: E402
import w1_r_sim  # noqa: E402
from w1_psf_smooth import KERNEL_SIGMA_PX  # noqa: E402
from w1_r_sim import PIX_ARCSEC, PSF_FWHM_ARCSEC  # noqa: E402

I2D = "s3://stpubdata/jwst/public/jw01181/L3/t/o001/jw01181-o001_t001_nircam_clear-f200w_i2d.fits"
TILE = 1024
N_TILES = 6
MIN_GOOD = 0.98
G = 0.01
PSF_R_FWHM = 0.10  # arcsec, round Gaussian reconvolution PSF (ASSUMPTION)
MATCH_PX = 1.5
EDGE_PX = 8
MIN_SNR = 10.0
# v2 (PLAN_v2.md; ASSUMPTIONs)
PSF_SIGMA_FLOOR_PX = 0.071 / PIX_ARCSEC / 2.3548  # measured PSF_d FWHM (D-101) as σ, ≈ 0.97 px
MASK_DILATE_PX = 2
MASK_PAD_PX = 16
SHEARS = {"noshear": (0.0, 0.0), "1p": (G, 0.0), "1m": (-G, 0.0), "2p": (0.0, G), "2m": (0.0, -G)}


def _gsp():
    import galsim

    return galsim.GSParams(maximum_fft_size=32768)


def psf_d(fwhm=PSF_FWHM_ARCSEC):
    """Deconvolution PSF: ``w1_r_sim.psf`` (Airy ⊗ Gaussian of FWHM ``fwhm`` / 2; measured FWHM
    0.071″ at the default 0.066″, D-101, not stored); ASSUMPTION for the drizzled PSF."""
    return w1_r_sim.psf(fwhm)


def psf_r():
    import galsim

    return galsim.Gaussian(fwhm=PSF_R_FWHM)


def psf_d_effective(spec):
    """Pixel-convolved PSF_d: a float is the analytic ``psf_d`` FWHM (⊗ pixel); an array is an
    empirical drizzled PSF image (PLAN_v3), which already includes the pixel."""
    import galsim

    if isinstance(spec, np.ndarray):
        return galsim.InterpolatedImage(
            galsim.Image(np.ascontiguousarray(spec, dtype=float), scale=PIX_ARCSEC),
            x_interpolant=galsim.Lanczos(15),
            gsparams=_gsp(),
        )
    return galsim.Convolve(psf_d(spec), galsim.Pixel(PIX_ARCSEC))


def metacal(arr: np.ndarray, g1: float, g2: float, psf_d_fwhm=PSF_FWHM_ARCSEC) -> np.ndarray:
    """``arr`` deconvolved by PSF_d ⊗ pixel (``psf_d_fwhm``: see :func:`psf_d_effective`), sheared
    by (g1, g2) about the image centre, reconvolved by PSF_r ⊗ pixel, on the same pixel grid."""
    import galsim

    gsp = _gsp()
    ny, nx = arr.shape
    ii = galsim.InterpolatedImage(
        galsim.Image(np.ascontiguousarray(arr, dtype=float), scale=PIX_ARCSEC),
        x_interpolant=galsim.Lanczos(15),
        gsparams=gsp,
    )
    dec = galsim.Deconvolve(psf_d_effective(psf_d_fwhm), gsparams=gsp)
    obj = galsim.Convolve(ii, dec, gsparams=gsp)
    if g1 or g2:
        obj = obj.shear(g1=g1, g2=g2)
    obj = galsim.Convolve(obj, psf_r(), gsparams=gsp)
    return obj.drawImage(nx=nx, ny=ny, scale=PIX_ARCSEC).array.astype(float)


def metacal_set(arr: np.ndarray, noise_sigma: float | None, seed: int, psf_d_fwhm=PSF_FWHM_ARCSEC):
    """The five metacal images; with ``noise_sigma`` a 90°-rotated white-noise image put through
    the same operation is added to each (fixnoise)."""
    noise = None
    if noise_sigma:
        noise = np.rot90(np.random.default_rng(seed).normal(0.0, noise_sigma, arr.shape))
    out = {}
    for k, (g1, g2) in SHEARS.items():
        im = metacal(arr, g1, g2, psf_d_fwhm)
        if noise is not None:
            im = im + np.rot90(metacal(noise, g1, g2, psf_d_fwhm), k=-1)
        out[k] = im
    return out


def clipped_sigma(arr: np.ndarray) -> float:
    from astropy.stats import sigma_clipped_stats

    return float(sigma_clipped_stats(arr[np.isfinite(arr)], sigma=3.0, maxiters=10)[2])


def catalogue(im: np.ndarray, kernel_sigma=KERNEL_SIGMA_PX, nsig=3.0, npixels=5) -> dict:
    """Scene version of ``w1_r_sim.measure`` (same threshold and settings): arrays x, y, a, b,
    pa, snr."""
    from photutils.segmentation import SourceCatalog, detect_sources
    from scipy import ndimage

    sig = clipped_sigma(im)
    conv = ndimage.gaussian_filter(im, kernel_sigma, mode="constant")
    thr = nsig * sig / (2 * np.sqrt(np.pi) * kernel_sigma)
    segm = detect_sources(conv, thr, n_pixels=npixels)
    if segm is None:
        return {k: np.zeros(0) for k in ("x", "y", "a", "b", "pa", "snr")}
    cat = SourceCatalog(im, segm, convolved_data=conv)
    area = np.asarray(cat.area.value, float)
    return {
        "x": np.asarray(cat.x_centroid, float),
        "y": np.asarray(cat.y_centroid, float),
        "a": np.asarray(cat.semimajor_axis.value, float),
        "b": np.asarray(cat.semiminor_axis.value, float),
        "pa": np.asarray(cat.orientation.value, float),
        "snr": np.asarray(cat.segment_flux, float) / (sig * np.sqrt(area)),
    }


def _match(
    ref: dict, cat: dict, g1: float, g2: float, centre: tuple[float, float], unique: bool = False
):
    """Index in ``cat`` of each ``ref`` source at its first-order sheared position (−1 if none
    within MATCH_PX; with ``unique``, also −1 where several ``ref`` sources share a match)."""
    from scipy.spatial import cKDTree

    if not len(cat["x"]):
        return np.full(len(ref["x"]), -1)
    dx, dy = ref["x"] - centre[0], ref["y"] - centre[1]
    xp = centre[0] + (1 + g1) * dx + g2 * dy
    yp = centre[1] + g2 * dx + (1 - g1) * dy
    d, i = cKDTree(np.c_[cat["x"], cat["y"]]).query(np.c_[xp, yp])
    idx = np.where(d <= MATCH_PX, i, -1)
    if unique:  # PLAN_v2 change 4: a target shared by several references is dropped for all
        tgt, cnt = np.unique(idx[idx >= 0], return_counts=True)
        idx[np.isin(idx, tgt[cnt > 1])] = -1
    return idx


def responses(
    images: dict,
    move_positions: bool = True,
    psf_sigma: float | None = None,
    v2: bool = False,
    near_mask: np.ndarray | None = None,
) -> dict:
    """Per-source R_11, R_22 from the five images (PLAN.md cuts); ``move_positions`` False when each
    galaxy was sheared about its own centre (V1 direct). ``v2`` adds the PLAN_v2 PSF σ floor and
    unique matching; ``near_mask`` (bool image) drops sources centred on True pixels."""
    cats = {k: catalogue(v) for k, v in images.items()}
    ref = cats["noshear"]
    ny, nx = images["noshear"].shape
    centre = ((nx - 1) / 2.0, (ny - 1) / 2.0)
    shapes = {
        "snr": ref["snr"],
        "semimajor_px": ref["a"],
        "ellipticity": 1.0 - ref["b"] / np.where(ref["a"] > 0, ref["a"], np.nan),
    }
    psig_raw = es.psf_sigma_px(shapes) if psf_sigma is None else psf_sigma
    psig = max(psig_raw, PSF_SIGMA_FLOOR_PX) if v2 and psf_sigma is None else psig_raw

    def eps(c, idx):
        out = np.full(len(idx), np.nan + 0j)
        ok = idx >= 0
        out[ok] = es.intrinsic_ellipticity(c["a"][idx[ok]], c["b"][idx[ok]], c["pa"][idx[ok]], psig)
        return out

    e = {"noshear": eps(ref, np.arange(len(ref["x"])))}
    for k, (g1, g2) in SHEARS.items():
        if k == "noshear":
            continue
        idx = (
            _match(ref, cats[k], g1, g2, centre, unique=v2)
            if move_positions
            else _match(ref, cats[k], 0, 0, centre, unique=v2)
        )
        e[k] = eps(cats[k], idx)
    inside = (
        (ref["x"] > EDGE_PX)
        & (ref["x"] < nx - 1 - EDGE_PX)
        & (ref["y"] > EDGE_PX)
        & (ref["y"] < ny - 1 - EDGE_PX)
    )
    ok = inside & (ref["snr"] >= MIN_SNR) & np.all([np.isfinite(v) for v in e.values()], axis=0)
    if near_mask is not None and len(ok):
        iy = np.clip(np.round(ref["y"]).astype(int), 0, ny - 1)
        ix = np.clip(np.round(ref["x"]).astype(int), 0, nx - 1)
        ok &= ~near_mask[iy, ix]
    r11 = (e["1p"].real - e["1m"].real) / (2 * G)
    r22 = (e["2p"].imag - e["2m"].imag) / (2 * G)
    return {
        "r11": r11[ok],
        "r22": r22[ok],
        "snr": ref["snr"][ok],
        "a": ref["a"][ok],
        "psf_sigma_px": psig,
        "psf_sigma_raw_px": psig_raw,
        "n_detected": int(len(ref["x"])),
        "mean_eps_noshear": complex(np.nanmean(e["noshear"][ok])) if ok.any() else complex("nan"),
    }


def summary(r11, r22, n_boot=1000, seed=1) -> dict:
    r = 0.5 * (np.asarray(r11) + np.asarray(r22))
    if not len(r):
        return {"n": 0, "R": None, "err": None}
    rng = np.random.default_rng(seed)
    boot = [r[rng.integers(0, len(r), len(r))].mean() for _ in range(n_boot)]
    return {
        "n": int(len(r)),
        "R": float(r.mean()),
        "err": float(np.std(boot)),
        "R_median": float(np.median(r)),
        "R11": float(np.mean(r11)),
        "R22": float(np.mean(r22)),
    }


def rule_outcome(r: float | None, err: float | None, r_borrowed: float = 0.414) -> str:
    if r is None or err is None:
        return "undetermined"
    return "conservative" if r >= r_borrowed - 2 * err else "optimistic: re-run lower"


# ---- validation (simulated) ----


def scene(
    g1: float, g2: float, the_psf, noise_sigma: float, seed: int, spacing=48, size=TILE,
    no_pixel=False,
):  # fmt: skip
    """Grid of Sérsic n = 1 galaxies (hlr cycles {2.5, 4, 6} px, S/N {20, 50, 200}), each sheared
    about its own centre by (g1, g2); the same ``seed`` gives the same galaxies and noise.
    ``no_pixel`` for a PSF that already includes the pixel (an empirical PSF image)."""
    import galsim

    rng = np.random.default_rng(seed)
    img = galsim.Image(size, size, scale=PIX_ARCSEC)
    hlrs, snrs = (2.5, 4.0, 6.0), (20.0, 50.0, 200.0)
    k = 0
    for yc in range(spacing // 2, size - spacing // 2 + 1, spacing):
        for xc in range(spacing // 2, size - spacing // 2 + 1, spacing):
            hlr, snr = hlrs[k % 3], snrs[(k // 3) % 3]
            k += 1
            e = rng.normal(0, 0.25, 2)
            while np.hypot(*e) >= 0.8:
                e = rng.normal(0, 0.25, 2)
            flux = snr * noise_sigma * np.sqrt(np.pi) * 2 * hlr  # rough S/N mapping (ASSUMPTION)
            gal = galsim.Sersic(1.0, half_light_radius=hlr * PIX_ARCSEC, flux=flux).shear(
                e1=e[0], e2=e[1]
            )
            if g1 or g2:
                gal = gal.shear(g1=g1, g2=g2)
            stamp = galsim.Convolve(gal, the_psf).drawImage(
                nx=spacing,
                ny=spacing,
                scale=PIX_ARCSEC,
                center=galsim.PositionI(xc + 1, yc + 1),
                method="no_pixel" if no_pixel else "auto",
            )
            b = stamp.bounds & img.bounds
            img[b] += stamp[b]
    arr = img.array.astype(float)
    return arr + np.random.default_rng(seed + 10_000).normal(0, noise_sigma, arr.shape)


def validate_pooled(seeds=(7, 8, 9, 10), noise_sigma=1.0) -> dict:
    """Post hoc (not the frozen gate): V1 on several scenes, per-source responses pooled."""
    per = [validate(noise_sigma, s, keep_arrays=True) for s in seeds]
    pooled = {}
    for k in ("metacal", "direct"):
        r11 = np.concatenate([p.pop(f"_{k}_r11") for p in per])
        r22 = np.concatenate([p.pop(f"_{k}_r22") for p in per])
        pooled[k] = summary(r11, r22)
    return {"label": "simulated", "seeds": list(seeds), "per_seed": per, "pooled": pooled}


def validate(noise_sigma=1.0, seed=7, keep_arrays=False, v2=False, epsf=None) -> dict:
    """PLAN.md gate (``v2``: PLAN_v2 changes); with ``epsf`` (PLAN_v3) the scene PSF and PSF_d are
    both that empirical PSF image."""
    t0 = time.time()
    if epsf is None:
        the_psf, spec, nopix = psf_d(), PSF_FWHM_ARCSEC, False
    else:
        the_psf, spec, nopix = psf_d_effective(epsf), epsf, True
    mc_in = scene(0, 0, the_psf, noise_sigma, seed, no_pixel=nopix)
    mc_imgs = metacal_set(mc_in, noise_sigma, seed + 1, spec)
    mc = responses(mc_imgs, v2=v2)
    out_sigma = clipped_sigma(mc_imgs["noshear"])
    # V1 direct: the same galaxies drawn with PSF_r at ±g, white noise at the metacal output level
    direct_imgs = {k: scene(g1, g2, psf_r(), out_sigma, seed) for k, (g1, g2) in SHEARS.items()}
    di = responses(direct_imgs, move_positions=False, psf_sigma=mc["psf_sigma_px"], v2=v2)
    s_mc, s_di = summary(mc["r11"], mc["r22"]), summary(di["r11"], di["r22"])
    if s_mc["R"] is None or s_di["R"] is None:
        tol, v1 = float("nan"), False
    else:
        tol = 0.10 * s_di["R"] + 2 * np.hypot(s_mc["err"], s_di["err"])
        v1 = abs(s_mc["R"] - s_di["R"]) <= tol
    # V2 sign: true g1 = +0.03 through the noshear metacal image vs g = 0
    plus = responses(
        metacal_set(
            scene(0.03, 0, the_psf, noise_sigma, seed, no_pixel=nopix), noise_sigma, seed + 1, spec
        ),
        v2=v2,
    )
    shift = plus["mean_eps_noshear"].real - mc["mean_eps_noshear"].real
    extra = {}
    if keep_arrays:
        extra = {"_metacal_r11": mc["r11"], "_metacal_r22": mc["r22"], "_direct_r11": di["r11"],
                 "_direct_r22": di["r22"]}  # fmt: skip
    return {
        **extra,
        "label": "simulated",
        "seed": seed,
        "v2": v2,
        "psf": "analytic" if epsf is None else "epsf",
        "metacal": s_mc,
        "direct": s_di,
        "psf_sigma_px": mc["psf_sigma_px"],
        "output_noise_sigma": out_sigma,
        "v1_tolerance": float(tol),
        "v1_pass": bool(v1),
        "v2_mean_re_eps_shift": float(shift),
        "v2_pass": bool(shift > 0),
        "gate_pass": bool(v1 and shift > 0),
        "wall_s": round(time.time() - t0, 1),
    }


# ---- real data ----


def good_tiles(wht_rows, shape, n=N_TILES, tile=TILE, min_good=MIN_GOOD):
    """First ``n`` tile origins (row-major from the lower left) whose finite, WHT > 0 fraction
    is ≥ ``min_good``; ``wht_rows(y0, x0)`` returns (sci, wht) for the tile."""
    out = []
    ny, nx = shape
    for y0 in range(0, ny - tile + 1, tile):
        for x0 in range(0, nx - tile + 1, tile):
            sci, wht = wht_rows(y0, x0)
            good = np.isfinite(sci) & (wht > 0)
            if good.mean() >= min_good:
                out.append((y0, x0))
                if len(out) == n:
                    return out
    return out


def _real_tile(job) -> dict:
    arr, y0, x0, seed, psf_fwhm, fixnoise = job
    good = np.isfinite(arr)
    from astropy.stats import sigma_clipped_stats

    _, sky, sig = (float(v) for v in sigma_clipped_stats(arr[good], sigma=3.0, maxiters=10))
    a = np.where(good, arr - sky, 0.0)
    try:
        res = responses(metacal_set(a, sig if fixnoise else None, seed, psf_fwhm))
    except ValueError as exc:  # e.g. no S/N > 50 source for the PSF σ: keep the other tiles
        empty = np.zeros(0)
        res = {"r11": empty, "r22": empty, "snr": empty, "a": empty, "psf_sigma_px": None,
               "n_detected": None, "error": str(exc)}  # fmt: skip
    return {"y0": y0, "x0": x0, "sky": sky, "sigma": sig, **res}


def _init_worker():
    os.environ["OMP_NUM_THREADS"] = "1"  # scripts/CLAUDE.md parallel topology


def real(out: Path, uri=I2D, workers=None) -> dict:
    from astropy.io import fits

    t0 = time.time()
    with fits.open(uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as h:
        sci, wht = h["SCI"], h["WHT"]

        def rows(y0, x0):
            return (
                np.asarray(sci.section[y0 : y0 + TILE, x0 : x0 + TILE], float),
                np.asarray(wht.section[y0 : y0 + TILE, x0 : x0 + TILE], float),
            )

        cache = {}

        def rows_cached(y0, x0):
            sci_t, wht_t = rows(y0, x0)
            if (np.isfinite(sci_t) & (wht_t > 0)).mean() >= MIN_GOOD:
                cache[(y0, x0)] = (sci_t, wht_t)  # keep accepted tiles only
            return sci_t, wht_t

        tiles = good_tiles(rows_cached, sci.shape)
        arrs = {t: np.where(cache[t][1] > 0, cache[t][0], np.nan) for t in tiles}
        mosaic_shape = list(sci.shape)
    t_read = time.time() - t0
    variants = {"plan": (PSF_FWHM_ARCSEC, True), "no_fixnoise": (PSF_FWHM_ARCSEC, False),
                "psf_d_x1.15": (1.15 * PSF_FWHM_ARCSEC, True)}  # fmt: skip
    jobs = [
        (arrs[t], t[0], t[1], 100 + i, fw, fx, name)
        for name, (fw, fx) in variants.items()
        for i, t in enumerate(tiles)
    ]
    with ProcessPoolExecutor(max_workers=workers or os.cpu_count(), initializer=_init_worker) as ex:
        res = list(ex.map(_real_tile, [j[:6] for j in jobs]))
    result = {
        "label": "derived (metacal of observed MAST JWST 1181 o001 F200W i2d tiles)",
        "uri": uri,
        "mosaic_shape": mosaic_shape,
        "tiles": [list(t) for t in tiles],
        "settings": {
            "g": G,
            "psf_d_fwhm_arcsec": PSF_FWHM_ARCSEC,
            "psf_r_fwhm_arcsec": PSF_R_FWHM,
            "match_px": MATCH_PX,
            "edge_px": EDGE_PX,
            "min_snr": MIN_SNR,
        },  # fmt: skip
        "variants": {},
        "read_s": round(t_read, 1),
    }
    for name in variants:
        rr = [r for r, j in zip(res, jobs, strict=True) if j[6] == name]
        r11 = np.concatenate([r["r11"] for r in rr])
        r22 = np.concatenate([r["r22"] for r in rr])
        s = summary(r11, r22)
        s["rule"] = rule_outcome(s["R"], s["err"])
        s["per_tile"] = [
            {
                "y0": r["y0"],
                "x0": r["x0"],
                "psf_sigma_px": r["psf_sigma_px"],
                "n_detected": r["n_detected"],
                "sigma": r["sigma"],
                **summary(r["r11"], r["r22"], n_boot=200),
            }  # fmt: skip
            for r in rr
        ]
        snr = np.concatenate([r["snr"] for r in rr])
        s["by_snr"] = {
            f"{lo:g}-{hi:g}": summary(
                r11[(snr >= lo) & (snr < hi)], r22[(snr >= lo) & (snr < hi)], n_boot=200
            )
            for lo, hi in ((10, 20), (20, 50), (50, 200), (200, np.inf))
        }
        result["variants"][name] = s
    result["wall_s"] = round(time.time() - t0, 1)
    out.write_text(json.dumps(result, indent=1) + "\n")
    return result


def sheet(out: Path, tiles=((0, 1024), (1024, 1024)), uri=I2D, v2=False, seeds=None) -> None:
    """Contact sheet: original, metacal unsheared and 1p − 1m (plan settings) for ``tiles``;
    ``v2`` prepares tiles with ``prep_tile_v2`` (noise-filled mask) and ``seeds`` (the run's
    per-tile seeds, ``real2``: 100 + tile index) reproduces the measured realisations."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats

    fig, ax = plt.subplots(len(tiles), 3, figsize=(12, 4 * len(tiles)))
    with fits.open(uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as h:
        for i, (y0, x0) in enumerate(tiles):
            a = np.asarray(h["SCI"].section[y0 : y0 + TILE, x0 : x0 + TILE], float)
            w = np.asarray(h["WHT"].section[y0 : y0 + TILE, x0 : x0 + TILE], float)
            a = np.where(w > 0, a, np.nan)
            seed = seeds[i] if seeds else 1
            if v2:
                b, _, _, sd, _ = prep_tile_v2(a, seed + 50_000 if seeds else 1)
            else:
                good = np.isfinite(a)
                _, med, sd = sigma_clipped_stats(a[good], sigma=3.0, maxiters=10)
                b = np.where(good, a - med, 0.0)
            ims = metacal_set(b, sd, seed)
            panels = [(b, "original"), (ims["noshear"], "metacal noshear"),
                      (ims["1p"] - ims["1m"], "1p - 1m")]  # fmt: skip
            for j, (im, title) in enumerate(panels):
                s = sigma_clipped_stats(im)[2]
                ax[i, j].imshow(im, origin="lower", cmap="gray", vmin=-3 * s,
                                vmax=(10 if j < 2 else 3) * s)  # fmt: skip
                ax[i, j].set_title(f"tile y0={y0} x0={x0}: {title}", fontsize=8)
                ax[i, j].axis("off")
    plt.tight_layout()
    plt.savefig(out, dpi=60)


def posthoc(r_real: Path, out: Path, exclude=((0, 1024),), n_boot=2000, seed=1) -> dict:
    """Post hoc (deviation) summaries from ``r_real.json``'s per-tile numbers: tiles in ``exclude``
    dropped (n-weighted mean of the per-tile means; error from the per-tile bootstrap errors in
    quadrature, approximate), and for every tile set a bootstrap over tiles (tiles as the
    independent unit)."""
    res = json.loads(r_real.read_text())
    excl = {tuple(t) for t in exclude}
    rng = np.random.default_rng(seed)
    result = {
        "label": "post hoc (deviation), derived from r_real.json per-tile summaries",
        "excluded_tiles": [list(t) for t in exclude],
        "variants": {},
    }
    for name, v in res["variants"].items():
        entry = {}
        for key, tiles in (
            ("all_tiles", [t for t in v["per_tile"] if t["n"]]),
            ("excluded", [t for t in v["per_tile"] if t["n"] and (t["y0"], t["x0"]) not in excl]),
        ):
            n = np.array([t["n"] for t in tiles], float)
            r = np.array([t["R"] for t in tiles], float)
            e = np.array([t["err"] for t in tiles], float)
            mean = float((n * r).sum() / n.sum())
            idx = rng.integers(0, len(r), (n_boot, len(r)))
            boot = (n[idx] * r[idx]).sum(axis=1) / n[idx].sum(axis=1)
            entry[key] = {
                "n_tiles": len(tiles),
                "n": int(n.sum()),
                "R": mean,
                "err_sources": float(np.sqrt(((n * e) ** 2).sum()) / n.sum()),
                "err_tile_bootstrap": float(np.std(boot)),
            }
            err = max(entry[key]["err_sources"], entry[key]["err_tile_bootstrap"])
            entry[key]["rule_larger_error"] = rule_outcome(mean, err)
        result["variants"][name] = entry
    out.write_text(json.dumps(result, indent=1) + "\n")
    return result


# ---- v2 (PLAN_v2.md) ----


def _disk(r: int) -> np.ndarray:
    y, x = np.mgrid[-r : r + 1, -r : r + 1]
    return x * x + y * y <= r * r


def prep_tile_v2(arr: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Sky-subtracted tile with masked pixels (non-finite, dilated by a MASK_DILATE_PX disk) filled
    with white noise of the clipped σ; returns (image, near-mask bool image within MASK_PAD_PX
    (Euclidean) of a masked pixel, sky, σ, masked-pixel fraction)."""
    from astropy.stats import sigma_clipped_stats
    from scipy import ndimage

    bad = ~np.isfinite(arr)
    _, sky, sig = (float(v) for v in sigma_clipped_stats(arr[~bad], sigma=3.0, maxiters=10))
    if bad.any():
        fill = ndimage.binary_dilation(bad, structure=_disk(MASK_DILATE_PX))
        near = ndimage.distance_transform_edt(~bad) <= MASK_PAD_PX
    else:
        fill = near = bad
    noise = np.random.default_rng(seed).normal(0.0, sig, arr.shape)
    return np.where(fill, noise, arr - sky), near, sky, sig, float(bad.mean())


def _real_tile_v2(job) -> dict:
    """One tile, prepared once, through every variant ``{name: (psf_fwhm, fixnoise)}``."""
    arr, y0, x0, seed, variants = job
    a, near, sky, sig, bad_frac = prep_tile_v2(arr, seed + 50_000)
    out = {}
    for name, (psf_fwhm, fixnoise) in variants.items():
        try:
            res = responses(metacal_set(a, sig if fixnoise else None, seed, psf_fwhm), v2=True,
                            near_mask=near)  # fmt: skip
        except ValueError as exc:
            empty = np.zeros(0)
            res = {"r11": empty, "r22": empty, "snr": empty, "a": empty, "psf_sigma_px": None,
                   "psf_sigma_raw_px": None, "n_detected": None, "error": str(exc)}  # fmt: skip
        out[name] = {"y0": y0, "x0": x0, "sky": sky, "sigma": sig, "masked_frac": bad_frac,
                     "pad_frac": float(near.mean()), **res}  # fmt: skip
    return out


def tile_gate(per_tile: list[dict]) -> list[bool]:
    """PLAN_v2 change 5: keep a tile unless its raw PSF σ is below the floor or its detection
    count is below half the median over tiles with a count (deviation: the plan's "accepted
    tiles" is circular; the median over all coverage-passing tiles is used, ``gate_median``)."""
    med = gate_median(per_tile)
    return [
        bool(t.get("n_detected"))
        and t["n_detected"] >= 0.5 * med
        and t.get("psf_sigma_raw_px") is not None
        and t["psf_sigma_raw_px"] >= PSF_SIGMA_FLOOR_PX
        for t in per_tile
    ]


def gate_median(per_tile: list[dict]) -> float:
    counts = [t["n_detected"] for t in per_tile if t.get("n_detected")]
    return float(np.median(counts)) if counts else 0.0


def tile_bootstrap(per_tile_r: list[np.ndarray], n_boot=2000, seed=1) -> float | None:
    """1σ of the pooled mean R over bootstrap resamples of tiles (sources pooled per draw)."""
    tiles = [r for r in per_tile_r if len(r)]
    if len(tiles) < 2:
        return None
    rng = np.random.default_rng(seed)
    sums = np.array([r.sum() for r in tiles])
    ns = np.array([len(r) for r in tiles])
    draws = rng.integers(0, len(tiles), (n_boot, len(tiles)))
    return float(np.std(sums[draws].sum(1) / ns[draws].sum(1)))


def rule_outcome_v2(r: float | None, err: float | None, adopt=0.80, r_borrowed=0.414) -> str:
    if r is None or err is None:
        return "undetermined"
    if r - 2 * err >= adopt:
        return "adoptable"
    return "conservative" if r >= r_borrowed - 2 * err else "optimistic: re-run lower"


def real2(
    out: Path, uri=I2D, workers=None, variants=None, plan="results/w1_r_realmetacal/PLAN_v2.md",
    label="derived (metacal v2 of observed MAST JWST 1181 o001 F200W i2d tiles)", extra=None,
) -> dict:  # fmt: skip
    from astropy.io import fits

    t0 = time.time()
    with fits.open(uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as h:
        sci, wht = h["SCI"], h["WHT"]
        cache = {}

        def rows_cached(y0, x0):
            sci_t = np.asarray(sci.section[y0 : y0 + TILE, x0 : x0 + TILE], float)
            wht_t = np.asarray(wht.section[y0 : y0 + TILE, x0 : x0 + TILE], float)
            if (np.isfinite(sci_t) & (wht_t > 0)).mean() >= MIN_GOOD:
                cache[(y0, x0)] = np.where(wht_t > 0, sci_t, np.nan)
            return sci_t, wht_t

        tiles = good_tiles(rows_cached, sci.shape, n=10**6)
        mosaic_shape = list(sci.shape)
    t_read = time.time() - t0
    variants = variants or {"plan": (PSF_FWHM_ARCSEC, True),
                            "no_fixnoise": (PSF_FWHM_ARCSEC, False),
                            "psf_d_x1.15": (1.15 * PSF_FWHM_ARCSEC, True)}  # fmt: skip
    jobs = [(cache[t], t[0], t[1], 100 + i, variants) for i, t in enumerate(tiles)]
    with ProcessPoolExecutor(max_workers=workers or os.cpu_count(), initializer=_init_worker) as ex:
        res = list(ex.map(_real_tile_v2, jobs))
    keep = tile_gate([r["plan"] for r in res])  # decided once, on the plan variant
    result = {
        "label": label,
        "plan": plan,
        **(extra or {}),
        "uri": uri,
        "mosaic_shape": mosaic_shape,
        "tiles": [list(t) for t in tiles],
        "settings": {
            "g": G,
            "psf_d_fwhm_arcsec": PSF_FWHM_ARCSEC,
            "psf_r_fwhm_arcsec": PSF_R_FWHM,
            "match_px": MATCH_PX,
            "edge_px": EDGE_PX,
            "min_snr": MIN_SNR,
            "psf_sigma_floor_px": PSF_SIGMA_FLOOR_PX,
            "mask_dilate_px": MASK_DILATE_PX,
            "mask_pad_px": MASK_PAD_PX,
        },  # fmt: skip
        "variants": {},
        "gate": {"variant": "plan", "median_n_detected": gate_median([r["plan"] for r in res])},
        "read_s": round(t_read, 1),
    }
    for name in variants:
        rr = [r[name] for r in res]
        kept = [r for r, k in zip(rr, keep, strict=True) if k]
        r11 = np.concatenate([r["r11"] for r in kept]) if kept else np.zeros(0)
        r22 = np.concatenate([r["r22"] for r in kept]) if kept else np.zeros(0)
        s = summary(r11, r22)
        s["err_source_bootstrap"] = s.pop("err")
        s["err"] = tile_bootstrap([0.5 * (r["r11"] + r["r22"]) for r in kept])
        s["n_tiles"] = len(kept)
        s["rule"] = rule_outcome_v2(s["R"], s["err"])
        s["excluded_tiles"] = [[r["y0"], r["x0"]] for r, k in zip(rr, keep, strict=True) if not k]
        s["per_tile"] = [
            {
                "y0": r["y0"],
                "x0": r["x0"],
                "kept": k,
                "psf_sigma_raw_px": r["psf_sigma_raw_px"],
                "n_detected": r["n_detected"],
                "sigma": r["sigma"],
                "masked_frac": r["masked_frac"],
                "pad_frac": r["pad_frac"],
                **summary(r["r11"], r["r22"], n_boot=200),
            }  # fmt: skip
            for r, k in zip(rr, keep, strict=True)
        ]
        snr = np.concatenate([r["snr"] for r in kept]) if kept else np.zeros(0)
        s["by_snr"] = {
            f"{lo:g}-{hi:g}": summary(
                r11[(snr >= lo) & (snr < hi)], r22[(snr >= lo) & (snr < hi)], n_boot=200
            )
            for lo, hi in ((10, 20), (20, 50), (50, 200), (200, np.inf))
        }
        result["variants"][name] = s
    result["wall_s"] = round(time.time() - t0, 1)
    out.write_text(json.dumps(result, indent=1) + "\n")
    return result


# ---- v3: empirical drizzled PSF_d (PLAN_v3.md; ASSUMPTIONs) ----
EPSF_HALF = 16
EPSF_G = (17.0, 21.0)
EPSF_RUWE = 1.4
EPSF_EDGE_PX = 20
EPSF_PEAK_PX = 3
EPSF_FLAT_CORE = 0.85
EPSF_SIZE_TOL = 0.20
EPSF_MIN_STARS = 5
EPSF_FWHM_RANGE = (0.060, 0.090)  # arcsec, gate E


def gaia_in_footprint(wcs, shape) -> list[dict]:
    """Gaia DR3 sources inside the mosaic (cone around its centre, then pixel bounds), PLAN_v3 1."""
    from astroquery.gaia import Gaia

    ny, nx = shape
    c = wcs.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    corners = ((0, 0), (nx - 1, ny - 1), (0, ny - 1), (nx - 1, 0))
    r = max(c.separation(wcs.pixel_to_world(x, y)).deg for x, y in corners)
    q = (
        "SELECT source_id, ra, dec, phot_g_mean_mag, ruwe, astrometric_params_solved FROM "
        "gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),"
        f"CIRCLE({c.ra.deg},{c.dec.deg},{r + 0.002}))"
    )
    t = Gaia.launch_job(q).get_results()
    x, y = wcs.world_to_pixel_values(np.asarray(t["ra"], float), np.asarray(t["dec"], float))
    return [
        {"source_id": int(t["source_id"][i]), "g": float(t["phot_g_mean_mag"][i]),
         "ruwe": None if np.ma.is_masked(t["ruwe"][i]) else float(t["ruwe"][i]),
         "params": int(t["astrometric_params_solved"][i]), "x": float(x[i]), "y": float(y[i])}
        for i in range(len(t))
        if EPSF_EDGE_PX <= x[i] < nx - EPSF_EDGE_PX and EPSF_EDGE_PX <= y[i] < ny - EPSF_EDGE_PX
    ]  # fmt: skip


def star_ok(star: dict) -> bool:
    return (
        EPSF_G[0] <= star["g"] <= EPSF_G[1]
        and star["ruwe"] is not None
        and star["ruwe"] < EPSF_RUWE
        and star["params"] > 3
    )


def prep_star(st: np.ndarray, wh: np.ndarray) -> tuple[np.ndarray | None, str]:
    """PLAN_v3 2: one star stamp → (recentred unit-sum stamp, "ok") or (None, reason)."""
    from scipy import ndimage

    h = st.shape[0] // 2
    if not (np.isfinite(st).all() and (wh > 0).all()):
        return None, "masked pixel"
    st = st - np.median(np.r_[st[0], st[-1], st[:, 0], st[:, -1]])
    py, px = np.unravel_index(np.argmax(st), st.shape)
    if max(abs(py - h), abs(px - h)) > EPSF_PEAK_PX:
        return None, "peak off centre"
    core = st[h - 1 : h + 2, h - 1 : h + 2]  # stamp centre: a plateau's argmax is its corner
    if core.max() <= 0 or core.min() > EPSF_FLAT_CORE * core.max():
        return None, "flat core"
    cy, cx = ndimage.center_of_mass(np.clip(st, 0, None) * (st > 0.1 * st.max()))
    st = ndimage.shift(st, (h - cy, h - cx), order=3, mode="nearest")
    if st.sum() <= 0:
        return None, "no flux"
    return st / st.sum(), "ok"


def moment_sigma(st: np.ndarray, r=5.0) -> float:
    """Flux-weighted second-moment σ (px) within ``r`` px of the stamp centre (positive pixels)."""
    h = st.shape[0] // 2
    y, x = np.mgrid[-h : h + 1, -h : h + 1]
    w = np.clip(st, 0, None) * (x * x + y * y <= r * r)
    return float(np.sqrt((w * (x * x + y * y)).sum() / w.sum() / 2))


def hm_fwhm_arcsec(psf: np.ndarray) -> float:
    """Gate E's FWHM: √(4A/π) pixels of the half-maximum area A, in arcsec."""
    return float(np.sqrt(4 * (psf >= 0.5 * psf.max()).sum() / np.pi) * PIX_ARCSEC)


def gate_e(n_stars: int, fwhm_arcsec: float) -> bool:
    return n_stars >= EPSF_MIN_STARS and EPSF_FWHM_RANGE[0] <= fwhm_arcsec <= EPSF_FWHM_RANGE[1]


def stack_epsf(stamps: list[np.ndarray]) -> tuple[np.ndarray | None, list[bool]]:
    """PLAN_v3 3–4: drop size outliers, median-stack, clip at 0, unit sum."""
    if not stamps:
        return None, []
    sig = np.array([moment_sigma(s) for s in stamps])
    keep = list(np.abs(sig / np.median(sig) - 1) <= EPSF_SIZE_TOL)
    kept = [s for s, k in zip(stamps, keep, strict=True) if k]
    psf = np.clip(np.median(kept, axis=0), 0, None)
    return psf / psf.sum(), keep


def epsf_sheet(stamps: list[np.ndarray], psf: np.ndarray | None, path: Path) -> None:
    """Star stamps, the ePSF and the analytic PSF_d ⊗ pixel on the same grid (log stretch)."""
    import galsim
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = 2 * EPSF_HALF + 1
    ana = (
        galsim.Convolve(psf_d(), galsim.Pixel(PIX_ARCSEC))
        .drawImage(nx=n, ny=n, scale=PIX_ARCSEC, method="no_pixel")
        .array
    )
    panels = [(f"star {i}", s) for i, s in enumerate(stamps)]
    panels += ([("ePSF (median)", psf)] if psf is not None else []) + [("analytic PSF_d", ana)]
    fig, axs = plt.subplots(1, len(panels), figsize=(2.6 * len(panels), 2.8))
    for ax, (t, im) in zip(np.atleast_1d(axs), panels, strict=True):
        ax.imshow(np.log10(np.clip(im / im.max(), 1e-4, None)), origin="lower", vmin=-4, vmax=0)
        ax.set_title(t, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def epsf3(out: Path, uri=I2D) -> dict:
    """Build the ePSF from the mosaic (PLAN_v3), write ``out`` (JSON with the image) and gate E."""
    import galsim
    from astropy.io import fits
    from astropy.wcs import WCS

    t0 = time.time()
    h = EPSF_HALF
    with fits.open(uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as f:
        sci, wht = f["SCI"], f["WHT"]
        wcs = WCS(sci.header)
        stars = gaia_in_footprint(wcs, sci.shape)
        stamps, rows = [], []
        for st in stars:
            row = {**st, "selected": star_ok(st)}
            if row["selected"]:
                xi, yi = int(round(st["x"])), int(round(st["y"]))
                a = np.asarray(sci.section[yi - h : yi + h + 1, xi - h : xi + h + 1], float)
                w = np.asarray(wht.section[yi - h : yi + h + 1, xi - h : xi + h + 1], float)
                s, why = prep_star(a, w)
                row["stamp"] = why
                if s is not None:
                    row["stamp_index"] = len(stamps)
                    stamps.append(s)
            rows.append(row)
    psf, keep = stack_epsf(stamps)
    epsf_sheet(stamps, psf, out.with_name("epsf_sheet.png"))
    for r in rows:
        if "stamp_index" in r:
            r["size_ok"] = bool(keep[r["stamp_index"]])
    n = int(sum(keep))
    fwhm = hm_fwhm_arcsec(psf) if psf is not None else float("nan")
    gs_fwhm = (
        float(psf_d_effective(psf).calculateFWHM()) if psf is not None else float("nan")
    )  # reported
    res = {
        "label": "derived (empirical PSF, observed JWST 1181 o001 F200W i2d, Gaia DR3 stars)",
        "plan": "results/w1_r_realmetacal/PLAN_v3.md",
        "uri": uri,
        "n_gaia_in_footprint": len(stars),
        "n_selected": int(sum(r["selected"] for r in rows)),
        "n_stamps": len(stamps),
        "n_stars": n,
        "fwhm_halfmax_area_arcsec": fwhm,
        "fwhm_galsim_arcsec_reported": gs_fwhm,
        "moment_sigma_px": moment_sigma(psf) if psf is not None else None,
        "analytic_psf_d_fwhm_galsim_arcsec": float(
            galsim.Convolve(psf_d(), galsim.Pixel(PIX_ARCSEC)).calculateFWHM()
        ),
        "gate_e_pass": gate_e(n, fwhm),
        "stars": rows,
        "epsf": None if psf is None else np.round(psf, 8).tolist(),
        "wall_s": round(time.time() - t0, 1),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1) + "\n")
    return res


def load_epsf(path: Path, posthoc: bool = False) -> np.ndarray:
    """The ePSF image; refused when gate E failed unless ``posthoc`` (diagnostic, not the rule)."""
    d = json.loads(path.read_text())
    if not d["gate_e_pass"] and not posthoc:
        raise SystemExit("gate E failed: PLAN_v3 stops before any R")
    return np.asarray(d["epsf"], float)


def real3(out: Path, epsf: np.ndarray, uri=I2D, workers=None, posthoc=False) -> dict:
    return real2(
        out, uri, workers,
        variants={"plan": (epsf, True), "epsf_no_fixnoise": (epsf, False),
                  "analytic_v2": (PSF_FWHM_ARCSEC, True)},
        plan="results/w1_r_realmetacal/PLAN_v3.md",
        label="derived (metacal v3, empirical PSF_d, observed MAST JWST 1181 o001 F200W i2d tiles)",
        extra={"psf_d": "plan, epsf_no_fixnoise: v3/epsf.json; analytic_v2: settings",
               "posthoc_gate_e_failed": posthoc},
    )  # fmt: skip


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--out", type=Path, required=True)
    v.add_argument("--pooled", action="store_true", help="post hoc: four scenes pooled")
    v2 = sub.add_parser("validate2", help="PLAN_v2 gate: V1/V2 with the floor and unique matching")
    v2.add_argument("--out", type=Path, required=True)
    r2 = sub.add_parser("real2", help="PLAN_v2 real run: all tiles, mask, floor, unique, tile gate")
    r2.add_argument("--out", type=Path, required=True)
    r2.add_argument("--workers", type=int, default=None)
    e3 = sub.add_parser("epsf3", help="PLAN_v3: build the empirical PSF and gate E")
    e3.add_argument("--out", type=Path, required=True)
    v3 = sub.add_parser("validate3", help="PLAN_v3 simulated gate with the ePSF")
    v3.add_argument("--epsf", type=Path, required=True)
    v3.add_argument("--out", type=Path, required=True)
    v3.add_argument("--posthoc", action="store_true", help="run although gate E failed")
    r3 = sub.add_parser("real3", help="PLAN_v3 real run with the ePSF as PSF_d")
    r3.add_argument("--epsf", type=Path, required=True)
    r3.add_argument("--out", type=Path, required=True)
    r3.add_argument("--workers", type=int, default=None)
    r3.add_argument("--posthoc", action="store_true", help="run although gate E failed")
    r = sub.add_parser("real")
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--workers", type=int, default=None)
    ph = sub.add_parser("posthoc")
    ph.add_argument("--r-real", type=Path, required=True)
    ph.add_argument("--out", type=Path, required=True)
    ph.add_argument("--exclude", type=int, nargs="+", default=None, help="y0 x0 pairs to drop")
    sh = sub.add_parser("sheet")
    sh.add_argument("--out", type=Path, required=True)
    sh.add_argument("--v2", action="store_true", help="PLAN_v2 tile preparation")
    sh.add_argument("--tiles", type=int, nargs="+", default=None, help="y0 x0 pairs")
    sh.add_argument("--r-real", type=Path, default=None, help="real2 JSON: use the run's seeds")
    a = p.parse_args(argv)
    if a.cmd == "posthoc":
        ex = list(zip(a.exclude[::2], a.exclude[1::2], strict=True)) if a.exclude else None
        print(json.dumps(posthoc(a.r_real, a.out, **({"exclude": ex} if ex else {})), indent=1))
        return 0
    if a.cmd == "sheet":
        t = list(zip(a.tiles[::2], a.tiles[1::2], strict=True)) if a.tiles else None
        if a.r_real and not t:
            t = [(0, 1024), (1024, 1024)]  # the default tiles, with the run's seeds
        seeds = None
        if a.r_real and t:
            order = [tuple(x) for x in json.loads(a.r_real.read_text())["tiles"]]
            seeds = [100 + order.index(x) for x in t]
        sheet(a.out, **({"tiles": t} if t else {}), v2=a.v2, seeds=seeds)
        return 0
    if a.cmd == "epsf3":
        res = epsf3(a.out)
        print(json.dumps({k: v for k, v in res.items() if k not in ("stars", "epsf")}, indent=1))
        return 0
    if a.cmd == "validate3":
        res = validate(v2=True, epsf=load_epsf(a.epsf, a.posthoc))
        res["posthoc_gate_e_failed"] = a.posthoc
        a.out.write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res, indent=1))
        return 0
    if a.cmd == "real3":
        res = real3(a.out, load_epsf(a.epsf, a.posthoc), workers=a.workers, posthoc=a.posthoc)
        print(json.dumps({k: (v["R"], v["err"], v["n_tiles"], v["rule"])
                          for k, v in res["variants"].items()}, indent=1))  # fmt: skip
        return 0
    if a.cmd == "validate2":
        res = validate(v2=True)
        a.out.write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res, indent=1))
        return 0
    if a.cmd == "real2":
        res = real2(a.out, workers=a.workers)
        print(json.dumps({k: (v["R"], v["err"], v["n_tiles"], v["rule"])
                          for k, v in res["variants"].items()}, indent=1))  # fmt: skip
        return 0
    if a.cmd == "validate":
        res = validate_pooled() if a.pooled else validate()
        a.out.write_text(json.dumps(res, indent=1) + "\n")
    else:
        res = real(a.out, workers=a.workers)
    print(
        json.dumps(
            res
            if a.cmd == "validate"
            else {k: (v["R"], v["err"], v["rule"]) for k, v in res["variants"].items()},
            indent=1,
        )
    )  # noqa: E501
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
