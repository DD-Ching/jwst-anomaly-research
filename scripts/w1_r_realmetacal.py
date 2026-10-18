"""W1 blank-field responsivity R by metacalibration of real F200W images.

Plan and rules (frozen first): ``results/w1_r_realmetacal/PLAN.md``. ``validate`` runs the
simulated gate (V1 metacal vs direct shearing, V2 sign); ``real`` streams 1024-px tiles of a JADES
GOODS-N F200W i2d from MAST S3 and measures R with the same chain. Shapes follow
``w1_r_sim.measure`` (kernel 2 px, 3σ smoothed-noise threshold, 5 px, no deblending) on whole
scenes; ε from ``exotic_screens.intrinsic_ellipticity``.

    python scripts/w1_r_realmetacal.py validate --out results/w1_r_realmetacal/validate.json
    python scripts/w1_r_realmetacal.py real --out results/w1_r_realmetacal/r_real.json
    python scripts/w1_r_realmetacal.py sheet --out results/w1_r_realmetacal/sheet.png
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
SHEARS = {"noshear": (0.0, 0.0), "1p": (G, 0.0), "1m": (-G, 0.0), "2p": (0.0, G), "2m": (0.0, -G)}


def _gsp():
    import galsim

    return galsim.GSParams(maximum_fft_size=32768)


def psf_d(fwhm=PSF_FWHM_ARCSEC):
    """Deconvolution PSF: ``w1_r_sim.psf`` (Airy ⊗ Gaussian); ASSUMPTION for the drizzled PSF."""
    import galsim

    airy = galsim.Airy(lam=2000.0, diam=6.5, scale_unit=galsim.arcsec)
    return galsim.Convolve(airy, galsim.Gaussian(fwhm=0.5 * fwhm))


def psf_r():
    import galsim

    return galsim.Gaussian(fwhm=PSF_R_FWHM)


def metacal(arr: np.ndarray, g1: float, g2: float, psf_d_fwhm=PSF_FWHM_ARCSEC) -> np.ndarray:
    """``arr`` deconvolved by PSF_d ⊗ pixel, sheared by (g1, g2) about the image centre,
    reconvolved by PSF_r ⊗ pixel, on the same pixel grid."""
    import galsim

    gsp = _gsp()
    ny, nx = arr.shape
    ii = galsim.InterpolatedImage(
        galsim.Image(np.ascontiguousarray(arr, dtype=float), scale=PIX_ARCSEC),
        x_interpolant=galsim.Lanczos(15),
        gsparams=gsp,
    )
    dec = galsim.Deconvolve(
        galsim.Convolve(psf_d(psf_d_fwhm), galsim.Pixel(PIX_ARCSEC)), gsparams=gsp
    )
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


def _match(ref: dict, cat: dict, g1: float, g2: float, centre: tuple[float, float]):
    """Index in ``cat`` of each ``ref`` source at its first-order sheared position (−1 if none
    within MATCH_PX)."""
    from scipy.spatial import cKDTree

    if not len(cat["x"]):
        return np.full(len(ref["x"]), -1)
    dx, dy = ref["x"] - centre[0], ref["y"] - centre[1]
    xp = centre[0] + (1 + g1) * dx + g2 * dy
    yp = centre[1] + g2 * dx + (1 - g1) * dy
    d, i = cKDTree(np.c_[cat["x"], cat["y"]]).query(np.c_[xp, yp])
    return np.where(d <= MATCH_PX, i, -1)


def responses(images: dict, move_positions: bool = True, psf_sigma: float | None = None) -> dict:
    """Per-source R_11, R_22 from the five images (PLAN.md cuts); ``move_positions`` False when each
    galaxy was sheared about its own centre (V1 direct)."""
    cats = {k: catalogue(v) for k, v in images.items()}
    ref = cats["noshear"]
    ny, nx = images["noshear"].shape
    centre = ((nx - 1) / 2.0, (ny - 1) / 2.0)
    shapes = {
        "snr": ref["snr"],
        "semimajor_px": ref["a"],
        "ellipticity": 1.0 - ref["b"] / np.where(ref["a"] > 0, ref["a"], np.nan),
    }
    psig = es.psf_sigma_px(shapes) if psf_sigma is None else psf_sigma

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
            _match(ref, cats[k], g1, g2, centre)
            if move_positions
            else _match(ref, cats[k], 0, 0, centre)
        )
        e[k] = eps(cats[k], idx)
    inside = (
        (ref["x"] > EDGE_PX)
        & (ref["x"] < nx - 1 - EDGE_PX)
        & (ref["y"] > EDGE_PX)
        & (ref["y"] < ny - 1 - EDGE_PX)
    )
    ok = inside & (ref["snr"] >= MIN_SNR) & np.all([np.isfinite(v) for v in e.values()], axis=0)
    r11 = (e["1p"].real - e["1m"].real) / (2 * G)
    r22 = (e["2p"].imag - e["2m"].imag) / (2 * G)
    return {
        "r11": r11[ok],
        "r22": r22[ok],
        "snr": ref["snr"][ok],
        "a": ref["a"][ok],
        "psf_sigma_px": psig,
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


def scene(g1: float, g2: float, the_psf, noise_sigma: float, seed: int, spacing=48, size=TILE):
    """Grid of Sérsic n = 1 galaxies (hlr cycles {2.5, 4, 6} px, S/N {20, 50, 200}), each sheared
    about its own centre by (g1, g2); the same ``seed`` gives the same galaxies and noise."""
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
                nx=spacing, ny=spacing, scale=PIX_ARCSEC, center=galsim.PositionI(xc + 1, yc + 1)
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


def validate(noise_sigma=1.0, seed=7, keep_arrays=False) -> dict:
    t0 = time.time()
    mc_in = scene(0, 0, psf_d(), noise_sigma, seed)
    mc_imgs = metacal_set(mc_in, noise_sigma, seed + 1)
    mc = responses(mc_imgs)
    out_sigma = clipped_sigma(mc_imgs["noshear"])
    # V1 direct: the same galaxies drawn with PSF_r at ±g, white noise at the metacal output level
    direct_imgs = {k: scene(g1, g2, psf_r(), out_sigma, seed) for k, (g1, g2) in SHEARS.items()}
    di = responses(direct_imgs, move_positions=False, psf_sigma=mc["psf_sigma_px"])
    s_mc, s_di = summary(mc["r11"], mc["r22"]), summary(di["r11"], di["r22"])
    tol = 0.10 * s_di["R"] + 2 * np.hypot(s_mc["err"], s_di["err"])
    v1 = abs(s_mc["R"] - s_di["R"]) <= tol
    # V2 sign: true g1 = +0.03 through the noshear metacal image vs g = 0
    plus = responses(metacal_set(scene(0.03, 0, psf_d(), noise_sigma, seed), noise_sigma, seed + 1))
    shift = plus["mean_eps_noshear"].real - mc["mean_eps_noshear"].real
    extra = {}
    if keep_arrays:
        extra = {"_metacal_r11": mc["r11"], "_metacal_r22": mc["r22"], "_direct_r11": di["r11"],
                 "_direct_r22": di["r22"]}  # fmt: skip
    return {
        **extra,
        "label": "simulated",
        "seed": seed,
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
    sig = clipped_sigma(arr)
    from astropy.stats import sigma_clipped_stats

    sky = float(sigma_clipped_stats(arr[good], sigma=3.0, maxiters=10)[1])
    a = np.where(good, arr - sky, 0.0)
    res = responses(metacal_set(a, sig if fixnoise else None, seed, psf_fwhm))
    return {"y0": y0, "x0": x0, "sky": sky, "sigma": sig, **res}


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
            cache[(y0, x0)] = rows(y0, x0)
            return cache[(y0, x0)]

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
    with ProcessPoolExecutor(max_workers=workers or os.cpu_count()) as ex:
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


def sheet(out: Path, tiles=((0, 1024), (1024, 1024)), uri=I2D) -> None:
    """Contact sheet: original, metacal unsheared and 1p − 1m (plan settings) for ``tiles``."""
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
            good = np.isfinite(a)
            _, med, sd = sigma_clipped_stats(a[good], sigma=3.0, maxiters=10)
            b = np.where(good, a - med, 0.0)
            ims = metacal_set(b, sd, 1)
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


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--out", type=Path, required=True)
    v.add_argument("--pooled", action="store_true", help="post hoc: four scenes pooled")
    r = sub.add_parser("real")
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--workers", type=int, default=None)
    sh = sub.add_parser("sheet")
    sh.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    if a.cmd == "sheet":
        sheet(a.out)
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
