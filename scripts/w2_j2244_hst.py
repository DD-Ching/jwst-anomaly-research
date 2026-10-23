"""J2244-0550 public HST test (results/w2_j2244_hst/PLAN.md): two point sources +/- extended light.

Runs in an isolated environment with lenstronomy (see PLAN.md); the frozen rules live in
w2_j2244_hst_rules.py. Data are read as byte-range cutouts of the public HST drizzled products on
AWS S3 (stpubdata), so no full ACS frame is downloaded.

  stars   Gaia DR3 stars in the programme fields of a band -> star stamps (PSF library input)
  fit     one system stamp: M0 (two point sources + uniform) and M1 (M0 + Sersic), PSF iteration
          in both, delta chi^2 on the common noise model, status by the frozen rules
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w2_j2244_hst_rules as rules  # noqa: E402

S3 = "https://stpubdata.s3.amazonaws.com/hst/public"
# fsspec's default 5 MB read-ahead block made each ACS star cutout cost ~12 MB of traffic in this
# round (RESULT "Cost"); small blocks keep a byte-range cutout close to its own size.
FSSPEC = {"block_size": 2**16}
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
BAND = {
    # stamp and star-stamp half-sizes (px), core mask radius (arcsec), pixel (arcsec), PSF size (px)
    "F160W": {"half": 24, "star_half": 24, "core": 0.30, "pix": 0.12825, "psf": 21},
    "F814W": {"half": 50, "star_half": 50, "core": 0.15, "pix": 0.05, "psf": 41},
}


def url_of(obs_id: str) -> str:
    """Pipeline association drizzled product (DRZ for WFC3/IR, DRC for ACS/WFC) on S3."""
    ext = "drz" if obs_id.startswith("i") else "drc"
    return f"{S3}/{obs_id[:4]}/{obs_id}/{obs_id}_{ext}.fits"


def cutout(obs_id: str, ra: float, dec: float, half: int) -> dict:
    """SCI (e-/s) and WHT cutout centred on (ra, dec) by byte-range reads."""
    from astropy.io import fits
    from astropy.wcs import WCS

    with fits.open(url_of(obs_id), use_fsspec=True, fsspec_kwargs=FSSPEC) as h:
        hdr = h["SCI"].header
        w = WCS(hdr)
        x, y = (float(v) for v in w.all_world2pix(ra, dec, 0))
        xi, yi = int(round(x)), int(round(y))
        ny, nx = hdr["NAXIS2"], hdr["NAXIS1"]
        if not (half <= xi < nx - half and half <= yi < ny - half):
            return {"ok": False, "why": "edge"}
        sl = (slice(yi - half, yi + half + 1), slice(xi - half, xi + half + 1))
        sci = np.asarray(h["SCI"].section[sl], dtype=float)
        wht = np.asarray(h["WHT"].section[sl], dtype=float)
        p0 = h[0].header
        return {
            "ok": True,
            "sci": sci,
            "wht": wht,
            "x0": xi - half,
            "y0": yi - half,
            "dx": x - xi,
            "dy": y - yi,
            "exptime": float(p0["EXPTIME"]),
            "ncombine": int(hdr.get("NCOMBINE", p0.get("NCOMBINE", 1))),
            "photflam": float(hdr.get("PHOTFLAM", p0.get("PHOTFLAM", np.nan))),
            "photplam": float(hdr.get("PHOTPLAM", p0.get("PHOTPLAM", np.nan))),
            "obs_id": obs_id,
        }


def zp_ab(photflam: float, photplam: float) -> float:
    return -2.5 * np.log10(photflam) - 5 * np.log10(photplam) - 2.408


def sky_stats(img: np.ndarray, mask: np.ndarray | None = None) -> tuple[float, float]:
    from astropy.stats import sigma_clipped_stats

    good = np.isfinite(img) if mask is None else (np.isfinite(img) & mask)
    _, med, std = sigma_clipped_stats(img[good], sigma=3.0, maxiters=10)
    return float(med), float(std)


def gaia_stars(ra: float, dec: float, radius_deg: float) -> list[dict]:
    q = (
        "SELECT source_id, ra, dec, phot_g_mean_mag, parallax_over_error, pmra, pmdec, pmra_error, "
        "pmdec_error, ruwe FROM gaiadr3.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))=1"
    )
    for k in range(3):
        r = requests.post(
            GAIA_TAP,
            data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": q},
            timeout=180,
        )
        if r.status_code == 200:
            j = r.json()
            cols = [c["name"] for c in j["metadata"]]
            return [dict(zip(cols, row, strict=True)) for row in j["data"]]
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"Gaia TAP {r.status_code}")


def is_star(s: dict) -> bool:
    """Gaia DR3 star: significant parallax or proper motion, and a clean single-star fit."""
    if s.get("ruwe") is None or s["ruwe"] > 1.4:
        return False
    par = (s.get("parallax_over_error") or 0) >= 3
    pm = 0.0
    if s.get("pmra") is not None and s.get("pmra_error"):
        pm = float(np.hypot(s["pmra"] / s["pmra_error"], s["pmdec"] / s["pmdec_error"]))
    return bool(par or pm >= 5)


def model_psf_stamp(stamps: list[np.ndarray], size: int) -> np.ndarray:
    """Empirical PSF from star stamps with photutils EPSFBuilder (oversampling 1, drizzled grid)."""
    from astropy.nddata import NDData
    from astropy.table import Table
    from photutils.psf import EPSFBuilder, extract_stars

    n = stamps[0].shape[0]
    big = []
    for s in stamps:
        big.append(NDData(s))
    tabs = [Table({"x": [(n - 1) / 2], "y": [(n - 1) / 2]}) for _ in stamps]
    stars = extract_stars(big, tabs, size=size)
    epsf, _ = EPSFBuilder(oversampling=1, maxiters=10, progress_bar=False).build_epsf(stars)
    k = np.clip(np.asarray(epsf.data, dtype=float), 0, None)
    return k / k.sum()


def peaks(
    img: np.ndarray, pix: float, n: int = 2, min_sep: float = 0.3
) -> list[tuple[float, float]]:
    """Brightest local maxima (pixel x, y), at least min_sep apart."""
    from scipy.ndimage import gaussian_filter, maximum_filter

    sm = gaussian_filter(np.nan_to_num(img), 1.0)
    loc = (sm == maximum_filter(sm, size=3)) & np.isfinite(img)
    ys, xs = np.nonzero(loc)
    order = np.argsort(sm[ys, xs])[::-1]
    out: list[tuple[float, float]] = []
    for i in order:
        x, y = float(xs[i]), float(ys[i])
        if all(np.hypot(x - a, y - b) * pix >= min_sep for a, b in out):
            out.append((x, y))
        if len(out) == n:
            break
    return out


def fit(
    img: np.ndarray,
    wht: np.ndarray,
    psf: np.ndarray,
    pix: float,
    bg_rms: float,
    core: float,
    with_sersic: bool,
    seed: int,
    pso: tuple[int, int] = (60, 120),
) -> dict:
    """lenstronomy fit of two unlensed point sources + uniform background (+ one Sersic)."""
    from lenstronomy.Data.imaging_data import ImageData
    from lenstronomy.Data.psf import PSF
    from lenstronomy.ImSim.image_linear_solve import ImageLinearFit
    from lenstronomy.LightModel.light_model import LightModel
    from lenstronomy.PointSource.point_source import PointSource
    from lenstronomy.Workflow.fitting_sequence import FittingSequence

    np.random.seed(seed)  # noqa: NPY002 - lenstronomy's PSO draws from the global generator
    n = img.shape[0]
    c = (n - 1) / 2
    good = np.isfinite(img) & (wht > 0.5 * np.median(wht[wht > 0]))
    data = np.where(good, img, 0.0)
    expmap = np.where(good, wht, np.median(wht[wht > 0]))
    kwargs_data = {
        "image_data": data,
        "background_rms": bg_rms,
        "exposure_time": expmap,
        "ra_at_xy_0": -c * pix,
        "dec_at_xy_0": -c * pix,
        "transform_pix2angle": np.array([[pix, 0.0], [0.0, pix]]),
    }
    kwargs_psf = {"psf_type": "PIXEL", "kernel_point_source": psf}
    kwargs_numerics = {"supersampling_factor": 1}
    pk = peaks(np.where(good, img, np.nan), pix)
    if len(pk) < 2:
        return {"ok": False, "why": "fewer than 2 peaks"}
    xs = [(x - c) * pix for x, _ in pk]
    ys = [(y - c) * pix for _, y in pk]
    light = ["UNIFORM"] + (["SERSIC_ELLIPSE"] if with_sersic else [])
    kwargs_model = {"lens_light_model_list": light, "point_source_model_list": ["UNLENSED"]}
    mid = (np.mean(xs), np.mean(ys))
    ll_init = [{"amp": 0.0}]
    ll_sig = [{"amp": 1.0}]
    ll_fix = [{}]
    ll_lo = [{"amp": -1e3}]
    ll_up = [{"amp": 1e3}]
    if with_sersic:
        ll_init.append(
            {
                "R_sersic": 0.3,
                "n_sersic": 2.0,
                "e1": 0.0,
                "e2": 0.0,
                "center_x": mid[0],
                "center_y": mid[1],
            }
        )
        ll_sig.append(
            {
                "R_sersic": 0.2,
                "n_sersic": 1.0,
                "e1": 0.1,
                "e2": 0.1,
                "center_x": 0.1,
                "center_y": 0.1,
            }
        )
        ll_fix.append({})
        ll_lo.append(
            {
                "R_sersic": 0.04,
                "n_sersic": 0.5,
                "e1": -0.5,
                "e2": -0.5,
                "center_x": mid[0] - 0.8,
                "center_y": mid[1] - 0.8,
            }
        )
        ll_up.append(
            {
                "R_sersic": 2.0,
                "n_sersic": 6.0,
                "e1": 0.5,
                "e2": 0.5,
                "center_x": mid[0] + 0.8,
                "center_y": mid[1] + 0.8,
            }
        )
    ps_init = [{"ra_image": np.array(xs), "dec_image": np.array(ys)}]
    ps_sig = [{"ra_image": np.full(2, 0.02), "dec_image": np.full(2, 0.02)}]
    ps_lo = [{"ra_image": np.array(xs) - 0.15, "dec_image": np.array(ys) - 0.15}]
    ps_up = [{"ra_image": np.array(xs) + 0.15, "dec_image": np.array(ys) + 0.15}]
    kwargs_params = {
        "lens_light_model": [ll_init, ll_sig, ll_fix, ll_lo, ll_up],
        "point_source_model": [ps_init, ps_sig, [{}], ps_lo, ps_up],
    }
    kwargs_likelihood = {
        "check_bounds": True,
        "check_positive_flux": True,  # no negative point-source or Sersic amplitudes (v2)
        "image_likelihood_mask_list": [good.astype(int)],
    }
    fs = FittingSequence(
        {
            "multi_band_list": [[kwargs_data, kwargs_psf, kwargs_numerics]],
            "multi_band_type": "single-band",
        },
        kwargs_model,
        {"num_point_source_list": [2]},
        kwargs_likelihood,
        kwargs_params,
    )
    seq = [
        ["PSO", {"sigma_scale": 1.0, "n_particles": pso[0], "n_iterations": pso[1]}],
        [
            "psf_iteration",
            {
                "num_iter": 30,
                "psf_iter_factor": 0.5,
                "stacking_method": "median",
                "psf_symmetry": 1,
                "block_center_neighbour": 0.2,
                # v3 (review): False keeps the variance map the iteration estimates; True put back
                # the (absent, i.e. zero) initial map, so v2 ran without a PSF error term
                "keep_psf_variance_map": False,
                "error_map_radius": 0.6,
            },
        ],
        ["PSO", {"sigma_scale": 0.5, "n_particles": pso[0], "n_iterations": pso[1]}],
    ]
    t0 = time.time()
    fs.fit_sequence(seq)
    res = fs.best_fit()
    kpsf = fs.multi_band_list[0][1]
    # model image with the iterated PSF, linear amplitudes solved
    im = ImageLinearFit(
        ImageData(**kwargs_data),
        PSF(
            psf_type="PIXEL",
            kernel_point_source=kpsf["kernel_point_source"],
            psf_variance_map=kpsf.get("psf_variance_map"),
        ),
        lens_light_model_class=LightModel(light),
        point_source_class=PointSource(["UNLENSED"]),
        kwargs_numerics=kwargs_numerics,
        likelihood_mask=good.astype(int),
    )
    model, errmap, _, _ = im.image_linear_solve(
        kwargs_lens_light=res["kwargs_lens_light"], kwargs_ps=res["kwargs_ps"], inv_bool=False
    )
    ps = res["kwargs_ps"][0]
    out = {
        "ok": True,
        "kwargs_result": _jsonable(res),
        "model": model,
        "errmap": np.asarray(errmap, dtype=float),
        "psf_final": np.asarray(kpsf["kernel_point_source"]),
        "ps_x": [float(v) for v in ps["ra_image"]],
        "ps_y": [float(v) for v in ps["dec_image"]],
        "init_x": xs,
        "init_y": ys,
        "good": good,
        "seconds": round(time.time() - t0, 1),
    }
    if with_sersic:
        sl = res["kwargs_lens_light"][1]
        lm = LightModel(["SERSIC_ELLIPSE"])
        # total Sersic flux (e-/s): integrate the unconvolved profile on a 5x finer, 3x larger grid
        f = 5
        g = (np.arange(3 * n * f) - (3 * n * f - 1) / 2) * pix / f
        gx, gy = np.meshgrid(g, g)
        # lenstronomy amplitudes are surface brightness per arcsec^2: times the sub-pixel area.
        # (Fixed after the validation and injection runs; earlier outputs are stored as *_units_v1.)
        out["sersic_flux"] = float(lm.surface_brightness(gx, gy, [sl]).sum() * (pix / f) ** 2)
        out["sersic"] = {k: float(v) for k, v in sl.items()}
    return out


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [float(v) for v in o.ravel()]
    if isinstance(o, np.floating | np.integer):
        return float(o)
    return o


def measure(img, wht, psf, pix, core, seed, t_band=None, band=None):
    """M0 and M1 on one stamp; delta chi^2 on the M1 cores; status if a threshold is given."""
    med_w = np.median(wht[wht > 0])
    expmap = np.where(np.isfinite(img) & (wht > 0.5 * med_w), wht, med_w)
    _, bg = sky_stats(img, mask=_outer(img.shape[0]))
    r0 = fit(img, wht, psf, pix, bg, core, False, seed)
    r1 = fit(img, wht, psf, pix, bg, core, True, seed)
    if not (r0.get("ok") and r1.get("ok")):
        return {"ok": False, "why": r0.get("why") or r1.get("why")}
    cores = list(zip(r1["ps_x"], r1["ps_y"], strict=True))
    good = r0["good"] & r1["good"]

    # one noise model for both: M1's model and PSF error map. M0 and M1 are fitted separately
    # (own PSF iteration), so they are not nested and delta chi^2 can be negative.
    def _chi(model):
        n = img.shape[0]
        c_ = (n - 1) / 2
        yy, xx = np.mgrid[0:n, 0:n]
        m = good.copy()
        for x, y in cores:
            m &= np.hypot((xx - c_) * pix - x, (yy - c_) * pix - y) > core
        var = bg**2 + np.clip(r1["model"], 0, None) / expmap + r1["errmap"]
        return float(((img - model) ** 2 / var)[m].sum()), int(m.sum())

    c0, npx = _chi(r0["model"])
    c1, _ = _chi(r1["model"])
    nfree = 4 + 2 + 1 + 7  # positions, point amplitudes, uniform, Sersic (6 non-linear + amp)
    redchi = c1 / max(npx - nfree, 1)
    mid = (np.mean(r1["ps_x"]), np.mean(r1["ps_y"]))
    sx = r1["sersic"]["center_x"] - mid[0]
    sy = r1["sersic"]["center_y"] - mid[1]
    off = float(np.hypot(sx, sy))
    pos_ok = all(
        np.hypot(a - b, c_ - d) <= rules.P.pos_tol
        for a, b, c_, d in zip(r1["ps_x"], r1["init_x"], r1["ps_y"], r1["init_y"], strict=True)
    )
    sep = float(np.hypot(r1["ps_x"][0] - r1["ps_x"][1], r1["ps_y"][0] - r1["ps_y"][1]))
    out = {
        "ok": True,
        "bg_rms": bg,
        "chi2_m0": c0,
        "chi2_m1": c1,
        "n_pix": npx,
        "dchi2": c0 - c1,
        "redchi_m1": redchi,
        "sersic_offset": off,
        "sersic": r1["sersic"],
        "sersic_flux": r1["sersic_flux"],
        "pos_ok": bool(pos_ok),
        "sep": sep,
        "ps_x": r1["ps_x"],
        "ps_y": r1["ps_y"],
        "seconds": r0["seconds"] + r1["seconds"],
        "_r0": r0,
        "_r1": r1,
    }
    if t_band is not None:
        out["status"] = rules.status(
            out["dchi2"],
            t_band,
            off,
            redchi,
            pos_ok,
            sersic_r=r1["sersic"]["R_sersic"],
            sersic_flux=r1["sersic_flux"],
            band=band,
        )
    return out


def _outer(n: int, frac: float = 0.3) -> np.ndarray:
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2
    return np.hypot(xx - c, yy - c) > (0.5 - frac / 2) * n * 0.9


STAR_RULES = {
    # Gaia G range, field radius (deg), isolation (arcsec), saturation (e-/exposure), sat_g_max
    "F160W": {"g": (16.0, 21.5), "radius": 35 / 3600, "iso": 2.5, "sat": 1e12, "sat_g_max": 0.0},
    "F814W": {
        "g": (16.5, 21.5),
        "radius": 100 / 3600,
        "iso": 2.5,
        "sat": 50000.0,
        "sat_g_max": 19.0,
    },
}


def star_stamps(band: str, fields: list[dict], exclude: list[tuple[float, float]]) -> list[dict]:
    """Background-subtracted star stamps from the given programme fields."""
    from astropy.io import fits
    from astropy.wcs import WCS
    from scipy.ndimage import gaussian_filter

    sr = STAR_RULES[band]
    half = BAND[band]["star_half"]
    out = []
    for f in fields:
        try:
            g = gaia_stars(f["ra"], f["dec"], sr["radius"])
        except RuntimeError:
            continue
        cand = [
            s
            for s in g
            if s["phot_g_mean_mag"]
            and sr["g"][0] <= s["phot_g_mean_mag"] <= sr["g"][1]
            and is_star(s)
        ]
        if not cand:
            continue
        try:
            h = fits.open(url_of(f["obs_id"]), use_fsspec=True, fsspec_kwargs=FSSPEC)
        except (FileNotFoundError, OSError):
            continue
        with h:
            hdr = h["SCI"].header
            w = WCS(hdr)
            ny, nx = hdr["NAXIS2"], hdr["NAXIS1"]
            # v3 (review): NCOMBINE is in the SCI header of these products, not the primary
            ncomb = int(hdr.get("NCOMBINE", h[0].header.get("NCOMBINE", 1)))
            per_exp = float(h[0].header["EXPTIME"]) / max(ncomb, 1)
            for s in cand:
                if any(
                    np.hypot((s["ra"] - a) * np.cos(np.radians(s["dec"])), s["dec"] - b) * 3600
                    < 3.0
                    for a, b in exclude
                ):
                    continue
                near = [
                    o
                    for o in g
                    if o["source_id"] != s["source_id"]
                    and np.hypot(
                        (o["ra"] - s["ra"]) * np.cos(np.radians(s["dec"])), o["dec"] - s["dec"]
                    )
                    * 3600
                    < sr["iso"]
                ]
                if near:
                    continue
                x, y = (float(v) for v in w.all_world2pix(s["ra"], s["dec"], 0))
                xi, yi = int(round(x)), int(round(y))
                m = 6  # re-centre on the brightest pixel within +/- 6 px (Gaia epoch, WCS)
                if not (
                    half + m + 2 <= xi < nx - half - m - 2
                    and half + m + 2 <= yi < ny - half - m - 2
                ):
                    continue
                box = np.asarray(
                    h["SCI"].section[yi - m : yi + m + 1, xi - m : xi + m + 1], dtype=float
                )
                if not np.all(np.isfinite(box)):
                    continue
                # smoothed peak, so a saturated (zero-weight) core still centres on the star
                by, bx = np.unravel_index(np.argmax(gaussian_filter(box, 1.5)), box.shape)
                yi, xi = yi - m + int(by), xi - m + int(bx)
                sl = (slice(yi - half, yi + half + 1), slice(xi - half, xi + half + 1))
                sci = np.asarray(h["SCI"].section[sl], dtype=float)
                wht = np.asarray(h["WHT"].section[sl], dtype=float)
                medw = np.median(wht)
                core = wht[half - 3 : half + 4, half - 3 : half + 4]
                sat = bool(np.any(core <= 0.5 * medw))
                if not np.all(np.isfinite(sci)) or medw <= 0 or np.mean(wht <= 0.5 * medw) > 0.05:
                    continue
                med, std = sky_stats(sci, mask=_outer(sci.shape[0]))
                st = sci - med
                pk = float(gaussian_filter(st, 1.5)[half, half])
                if not sat:
                    pk = st[half, half]
                    sat = bool(pk * per_exp > sr["sat"])
                if sat and s["phot_g_mean_mag"] > sr["sat_g_max"]:
                    continue  # a saturated core is accepted only for stars brighter than sat_g_max
                if pk < 20 * std or (not sat and pk < st.max()):
                    continue
                # a second local maximum above 5 % of the peak means a neighbour or a galaxy
                p2 = peaks(st, BAND[band]["pix"], n=2, min_sep=0.25)
                if len(p2) > 1 and st[int(p2[1][1]), int(p2[1][0])] > 0.05 * pk:
                    continue
                out.append(
                    {
                        "obs_id": f["obs_id"],
                        "source_id": str(s["source_id"]),
                        "ra": float(s["ra"]),
                        "dec": float(s["dec"]),
                        "g": s["phot_g_mean_mag"],
                        "stamp": st,
                        "wht": wht,
                        "sky_std": std,
                        "peak": float(pk),
                        "sat": sat,
                    }
                )
    return out


def make_pair(s1: dict, s2: dict, sep_arcsec: float, pix: float, pa_deg: float) -> dict:
    """Negative control: two real star stamps, the second shifted by an integer-pixel vector of
    length close to sep_arcsec (no resampling), summed with their real noise; WHT = the smaller."""
    dx = int(round(sep_arcsec / pix * np.cos(np.radians(pa_deg)) / 2))
    dy = int(round(sep_arcsec / pix * np.sin(np.radians(pa_deg)) / 2))
    a = np.roll(np.roll(s1["stamp"], -dy, axis=0), -dx, axis=1)
    b = np.roll(np.roll(s2["stamp"], dy, axis=0), dx, axis=1)
    img = a + b
    wa = np.roll(np.roll(s1["wht"], -dy, axis=0), -dx, axis=1)
    wb = np.roll(np.roll(s2["wht"], dy, axis=0), dx, axis=1)
    edge = max(abs(dx), abs(dy)) + 1
    img[:edge, :] = np.nan
    img[-edge:, :] = np.nan
    img[:, :edge] = np.nan
    img[:, -edge:] = np.nan
    return {
        "sci": img,
        "wht": np.minimum(wa, wb),
        "sep_true": float(2 * np.hypot(dx, dy) * pix),
        "flux_ratio": float(s1["stamp"].sum() / s2["stamp"].sum()),
        "dx": dx,
        "dy": dy,
    }


def inject(
    img: np.ndarray,
    psf: np.ndarray,
    pix: float,
    x: float,
    y: float,
    mag: float,
    zp: float,
    r_eff: float,
    n_s: float,
    q: float,
    seed: int,
) -> np.ndarray:
    """Add a PSF-convolved Sersic of total AB magnitude mag (lenstronomy profile, 5x subsampled)."""
    from lenstronomy.LightModel.light_model import LightModel
    from lenstronomy.Util import param_util
    from scipy.signal import fftconvolve

    rng = np.random.default_rng(seed)
    e1, e2 = param_util.phi_q2_ellipticity(float(rng.uniform(0, np.pi)), q)
    n = img.shape[0]
    f = 5
    g = (np.arange(n * f) - (n * f - 1) / 2) * pix / f
    gx, gy = np.meshgrid(g, g)
    lm = LightModel(["SERSIC_ELLIPSE"])
    kw = [
        {
            "amp": 1.0,
            "R_sersic": r_eff,
            "n_sersic": n_s,
            "e1": e1,
            "e2": e2,
            "center_x": x,
            "center_y": y,
        }
    ]
    sub = lm.surface_brightness(gx, gy, kw)
    im = sub.reshape(n, f, n, f).sum(axis=(1, 3))
    # normalise to the total flux of the profile (integral on a 3x larger grid)
    g3 = (np.arange(3 * n * f) - (3 * n * f - 1) / 2) * pix / f
    gx3, gy3 = np.meshgrid(g3, g3)
    tot = lm.surface_brightness(gx3, gy3, kw).sum()
    flux = 10 ** (-0.4 * (mag - zp))
    im = im / tot * flux
    conv = fftconvolve(im, psf, mode="same")
    return img + conv  # source Poisson noise neglected (faint, sky-dominated)


OUTDIR = Path(__file__).resolve().parents[1] / "results" / "w2_j2244_hst"
FIELD_RADIUS = {"17455": 35 / 3600, "17199": 60 / 3600, "17308": 100 / 3600}


def build_library(
    band: str, prog: str, fields_json: str, out_pkl: str, max_fields: int = 80
) -> None:
    import pickle

    fields = json.load(open(fields_json))[f"{prog}_{band}"]
    seen, chosen = set(), []
    for f in fields:  # first association per target, positive exposure time, obs_id order
        if f["t_exptime"] < 100 or f["target"] in seen:
            continue
        seen.add(f["target"])
        chosen.append(f)
    chosen = chosen[:max_fields]
    STAR_RULES[band] = {**STAR_RULES[band], "radius": FIELD_RADIUS[prog]}
    st = star_stamps(band, chosen, exclude=[])
    st.sort(key=lambda s: -s["peak"])
    pickle.dump(st, open(out_pkl, "wb"))
    summary = [
        {
            "obs_id": s["obs_id"],
            "gaia_dr3": s["source_id"],
            "g": s["g"],
            "peak_e_s": s["peak"],
            "sat": s["sat"],
        }
        for s in st
    ]
    (OUTDIR / f"library_{prog}_{band}.json").write_text(
        json.dumps(
            {
                "label": "observed (Gaia DR3 stars in programme fields; star rules in PLAN.md)",
                "fields_used": [f["obs_id"] for f in chosen],
                "stars": summary,
            },
            indent=1,
        )
        + "\n"
    )
    print(band, prog, "fields", len(chosen), "stars", len(st))


def psf_from(lib: list[dict], band: str, skip: set) -> np.ndarray:
    stamps = [s["stamp"] for s in lib if s["source_id"] not in skip and not s.get("sat")]
    return model_psf_stamp(stamps, BAND[band]["psf"])


def png(path: Path, img, r0, r1, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    good = r1["good"]
    v = np.nanpercentile(np.where(good, img, np.nan), [1, 99.5])
    res0 = np.where(good, img - r0["model"], np.nan)
    res1 = np.where(good, img - r1["model"], np.nan)
    s = np.nanstd(res1) * 5
    fig, ax = plt.subplots(1, 4, figsize=(9.6, 2.6))
    for a, d, t, lim in (
        (ax[0], np.arcsinh(img / max(v[1] * 0.01, 1e-6)), "data (asinh)", None),
        (ax[1], np.arcsinh(r1["model"] / max(v[1] * 0.01, 1e-6)), "M1 model", None),
        (ax[2], res0, "data - M0", s),
        (ax[3], res1, "data - M1", s),
    ):
        if lim is None:
            a.imshow(d, origin="lower", cmap="gray")
        else:
            a.imshow(d, origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim)
        a.set_title(t, fontsize=8)
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle(title, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=70)
    plt.close(fig)


def run_job(job: dict) -> dict:
    """One job: a positive control, a negative pair, an injection or the target (same chain)."""
    import pickle

    band = job["band"]
    pix = BAND[band]["pix"]
    lib = pickle.load(open(job["library"], "rb"))
    out = {k: v for k, v in job.items() if k not in ("library",)}
    if job["kind"] in ("positive", "target", "target_injection"):
        c = cutout(job["obs_id"], job["ra"], job["dec"], BAND[band]["half"])
        if not c["ok"]:
            return {**out, "status": "no data", "why": c["why"]}
        img, wht = c["sci"], c["wht"]
        out["zp_ab"] = float(zp_ab(c["photflam"], c["photplam"]))
        near = {
            st["source_id"]
            for st in lib
            if np.hypot(
                (st["ra"] - job["ra"]) * np.cos(np.radians(job["dec"])), st["dec"] - job["dec"]
            )
            * 3600
            < 3.0
        }
        psf = psf_from(lib, band, near)
    else:  # negative pair from library stars i, j
        s1, s2 = lib[job["i"]], lib[job["j"]]
        pr = make_pair(s1, s2, job["sep"], pix, job["pa"])
        img, wht = pr["sci"], pr["wht"]
        out.update({k: pr[k] for k in ("sep_true", "flux_ratio", "dx", "dy")})
        out["stars"] = [s1["source_id"], s2["source_id"]]
        psf = psf_from(lib, band, {s1["source_id"], s2["source_id"]})
    if job.get("inject"):
        ij = job["inject"]
        # image positions from a peak search; injection 0.4 of the way from the fainter one
        pk = peaks(np.where(np.isfinite(img), img, np.nan), pix)
        n = img.shape[0]
        cc = (n - 1) / 2
        (xa, ya), (xb, yb) = pk[0], pk[1]  # brighter first
        x = ((xb + ij["frac"] * (xa - xb)) - cc) * pix
        y = ((yb + ij["frac"] * (ya - yb)) - cc) * pix
        img = inject(
            img, psf, pix, x, y, ij["mag"], job["zp_ab"], ij["r_eff"], 4.0, 0.8, ij["seed"]
        )
        out["inj_xy"] = [float(x), float(y)]
        out["inj_flux"] = float(10 ** (-0.4 * (ij["mag"] - job["zp_ab"])))
    r = measure(
        img,
        wht,
        psf,
        pix,
        BAND[band]["core"],
        seed=job.get("seed", 1),
        t_band=job.get("T"),
        band=band,
    )
    if not r.get("ok"):
        return {**out, "status": "fit unreliable", "why": r.get("why")}
    out.update({k: v for k, v in r.items() if not k.startswith("_")})
    if job.get("png"):
        png(
            OUTDIR / job["png"],
            img,
            r["_r0"],
            r["_r1"],
            f"{job['name']} {band} dchi2 {r['dchi2']:.0f}",
        )
    return out


def run_jobs(jobs_json: str, out_json: str, procs: int) -> None:
    """Run jobs in a process pool; each result is appended to <out_json>.jsonl as it finishes
    (a restart loses only running jobs), then all are written to out_json in job order."""
    import multiprocessing as mp

    jobs = json.load(open(jobs_json))
    part = Path(out_json + "l")
    done = {}
    if part.exists():
        for ln in part.read_text().splitlines():
            r = json.loads(ln)
            done[r["_k"]] = r
    todo = [dict(j, _k=k) for k, j in enumerate(jobs) if k not in done]
    with mp.get_context("fork").Pool(procs) as pool, part.open("a") as fh:
        for r in pool.imap_unordered(_run_keyed, todo, chunksize=1):
            fh.write(json.dumps(_jsonable(r)) + "\n")
            fh.flush()
            done[r["_k"]] = r
            print(
                r.get("name"),
                r.get("band"),
                r.get("status"),
                round(r.get("dchi2", float("nan")), 1),
                round(r.get("redchi_m1", float("nan")), 2),
                r.get("sersic_offset"),
                flush=True,
            )
    res = [done[k] for k in range(len(jobs))]
    Path(out_json).write_text(json.dumps(_jsonable(res), indent=1) + "\n")


def _run_keyed(job: dict) -> dict:
    k = job.pop("_k")
    r = run_job(job)
    r["_k"] = k
    return r


def gate_report(val_json: str, inj_json: str | None) -> dict:
    """Per-band gate on validation controls (PLAN section 4); positives must meet the pool geometry
    rule (fitted separation within rules.P.sep_range); injection limits and absorption (sec. 5)."""
    val = json.load(open(val_json))
    inj = json.load(open(inj_json)) if inj_json and Path(inj_json).exists() else []
    out = {"label": "derived (validation controls and injections; frozen v2 rules)", "bands": {}}
    for band in ("F160W", "F814W"):
        pos = [r for r in val if r["band"] == band and r["kind"] == "positive"]
        neg = [r for r in val if r["band"] == band and r["kind"] == "negative"]
        elig = [r for r in pos if rules.geometry_ok(r.get("sep"))]
        excl = [
            {"name": r["name"], "sep": r.get("sep"), "status": r.get("status")}
            for r in pos
            if r not in elig
        ]
        base = {r["name"]: r for r in neg}
        rows = []
        for r in (x for x in inj if x["band"] == band):
            b = base.get(r["name"].split("_")[0], {})
            rows.append(
                {
                    "pair": r["name"].split("_")[0],
                    "r_eff": r["inject"]["r_eff"],
                    "mag": r["inject"]["mag"],
                    "status": r.get("status"),
                    "dchi2": r.get("dchi2"),
                    "recovered": rules.detected(r.get("status", "")),
                    "flux_ratio_raw": (r.get("sersic_flux") or 0) / r["inj_flux"],
                    "flux_ratio_minus_baseline": (
                        (r.get("sersic_flux") or 0) - (b.get("sersic_flux") or 0)
                    )
                    / r["inj_flux"],
                }
            )
        limits = {}
        absorb_ok = bool(rows)
        for pair in sorted({x["pair"] for x in rows}):
            for reff in rules.P.inj_reff:
                sel = sorted(
                    (x for x in rows if x["pair"] == pair and x["r_eff"] == reff),
                    key=lambda x: x["mag"],
                )
                limits[f"{pair}_r{reff}"] = rules.contiguous_limit(
                    [x["mag"] for x in sel], [x["recovered"] for x in sel]
                )
                if sel and sel[0]["flux_ratio_raw"] < rules.P.absorb_min:
                    absorb_ok = False
        g = rules.gate(
            [r.get("status", "fit unreliable") for r in neg],
            [r.get("status", "fit unreliable") for r in elig],
            absorb_ok,
        )
        out["bands"][band] = {
            "T": next((r["T"] for r in val if r["band"] == band), None),
            "gate": g,
            "positives_eligible": [
                {"name": r["name"], "sep": r["sep"], "dchi2": r["dchi2"], "status": r["status"]}
                for r in elig
            ],
            "positives_excluded_geometry": excl,
            "negatives": [
                {
                    "name": r["name"],
                    "g": r.get("g"),
                    "dchi2": r.get("dchi2"),
                    "redchi": r.get("redchi_m1"),
                    "status": r.get("status"),
                }
                for r in neg
            ],
            "injections": rows,
            "injection_limits": limits,
        }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("library")
    p.add_argument("band")
    p.add_argument("prog")
    p.add_argument("fields_json")
    p.add_argument("out_pkl")
    p.add_argument("--max-fields", type=int, default=80)
    p = sub.add_parser("gate")
    p.add_argument("val_json")
    p.add_argument("out_json")
    p.add_argument("--inj", default=None)
    p = sub.add_parser("run")
    p.add_argument("jobs_json")
    p.add_argument("out_json")
    p.add_argument("--procs", type=int, default=4)
    a = ap.parse_args(argv)
    if a.cmd == "gate":
        rep = gate_report(a.val_json, a.inj)
        Path(a.out_json).write_text(json.dumps(_jsonable(rep), indent=1) + "\n")
        for b, v in rep["bands"].items():
            print(b, v["gate"], v["injection_limits"])
        return
    if a.cmd == "library":
        build_library(a.band, a.prog, a.fields_json, a.out_pkl, a.max_fields)
    else:
        run_jobs(a.jobs_json, a.out_json, a.procs)


if __name__ == "__main__":
    main()
