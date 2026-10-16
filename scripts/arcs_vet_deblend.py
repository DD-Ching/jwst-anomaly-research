"""Vetting step 2 for Abell S1063 ``arcs`` anti rows: deblend 1316, arc-fragment test for 2184.

Frozen plan: ``results/arcs/abells1063-cats/vetting/PLAN_deblend.md`` (written before the run).

1316: one elliptical Gaussian (M1) against two round Gaussians (M2) on the 1″-median-filtered
F200W cutout; ΔBIC = BIC(M1) − BIC(M2); rule ``deblend_class``, calibrated by injections of both
kinds at empty positions (gate). If two knots: separation PA against the z_s = 2 prediction and
the knot flux ratio across the SW bands (rule ``ratio_consistency``).

2184: 0.2″ aperture SEDs of 2184 and of two points on the nearby giant arc (rule ``sed_match``;
the two arc points are the control).

    python scripts/arcs_vet_deblend.py

Output: ``results/arcs/abells1063-cats/vetting/deblend.json`` and ``deblend.png`` (derived;
pixels observed, injections simulations, predictions model_prediction).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import stats
from scipy.ndimage import gaussian_filter, median_filter
from scipy.optimize import least_squares

from jwst_anomaly import cutouts, paths

_spec = importlib.util.spec_from_file_location(
    "arcs_bgsub_shapes", Path(__file__).resolve().parent / "arcs_bgsub_shapes.py"
)
abs_ = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(abs_)

URI = "s3://stpubdata/jwst/public/jw03293/L3/t/o001/jw03293-o001_t001_nircam_clear-{b}_i2d.fits"
SW = ["f090w", "f115w", "f150w", "f200w"]
BANDS = SW + ["f277w", "f356w", "f410m", "f444w", "f480m"]


# ---- frozen rules ---------------------------------------------------------------------------
def deblend_class(dbic: float, snr1: float, snr2: float, sep_arcsec: float) -> str:
    """``two-knot`` / ``single`` / ``undecided`` (PLAN_deblend.md)."""
    if not np.isfinite(dbic):
        return "undecided"
    if dbic >= 10 and min(snr1, snr2) >= 3 and sep_arcsec >= 0.10:
        return "two-knot"
    if dbic <= -10:
        return "single"
    return "undecided"


def gate_ok(frac_single_ok: float, frac_two_ok: float, min_frac: float = 0.8) -> bool:
    return bool(frac_single_ok >= min_frac and frac_two_ok >= min_frac)


def ratio_consistency(f1, e1, f2, e2, snr_min: float = 2.0) -> tuple[float, int]:
    """χ² p-value that f1/f2 is one constant across bands (bands with both S/N ≥ ``snr_min``)."""
    f1, e1, f2, e2 = (np.asarray(v, float) for v in (f1, e1, f2, e2))
    ok = (f1 / e1 >= snr_min) & (f2 / e2 >= snr_min)
    return sed_match(f1[ok], e1[ok], f2[ok], e2[ok], snr_min=-np.inf)


def sed_match(fa, ea, fb, eb, snr_min: float = 2.0, floor: float = 0.0) -> tuple[float, int]:
    """χ² p-value that SED a = A × SED b for one free A; bands with both S/N ≥ ``snr_min``.

    ``floor`` adds a fractional error in quadrature (0 in the frozen rule). Returns (p, number
    of bands used); p is NaN with fewer than 2 bands.
    """
    fa, ea, fb, eb = (np.asarray(v, float) for v in (fa, ea, fb, eb))
    ea, eb = np.hypot(ea, floor * np.abs(fa)), np.hypot(eb, floor * np.abs(fb))
    ok = (fa / ea >= snr_min) & (fb / eb >= snr_min) & np.isfinite(fa + fb + ea + eb)
    fa, ea, fb, eb = fa[ok], ea[ok], fb[ok], eb[ok]
    n = int(fa.size)
    if n < 2:
        return float("nan"), n

    def chi2(a):
        return float(np.sum((fa - a * fb) ** 2 / (ea**2 + a**2 * eb**2)))

    a0 = float(np.sum(fa * fb / ea**2) / np.sum(fb**2 / ea**2))
    grid = a0 * np.exp(np.linspace(-1.5, 1.5, 601))
    c = min(chi2(a) for a in grid)
    return float(stats.chi2.sf(c, n - 1)), n


def sed_class(p: float) -> str:
    if not np.isfinite(p):
        return "undetermined"
    return "SED-consistent" if p >= 0.003 else "SED-different"


# ---- models ---------------------------------------------------------------------------------
def _gauss(xx, yy, x0, y0, flux, sig, q=1.0, th=0.0):
    c, s = np.cos(th), np.sin(th)
    dx, dy = xx - x0, yy - y0
    uu, vv = c * dx + s * dy, -s * dx + c * dy
    sb = sig * q
    return flux / (2 * np.pi * sig * sb) * np.exp(-0.5 * ((uu / sig) ** 2 + (vv / sb) ** 2))


def flux_scale(D) -> tuple[float, float]:
    """Start flux and upper bound for the fits (PLAN_fixes.md fix 3).

    The start is the filtered aperture sum when it is at least 10 % of the sum of the positive
    pixels, else the positive-pixel sum; the bound is 10 × the larger of the two. (Before the fix
    both came from the sum floored at 1e-6, which pinned faint sources whose filtered sum is
    ≤ 0 or near 0.)
    """
    D = np.asarray(D, float)
    s, s_pos = float(np.sum(D)), float(np.sum(np.clip(D, 0.0, None)))
    start = s if s >= 0.1 * s_pos and s > 0 else max(s_pos, 1e-6)
    return start, 10.0 * max(s, s_pos, 1e-6)


def fit_models(img, err, x0, y0, r_px):
    """Fit M1 and M2 within ``r_px`` of (x0, y0); returns dict with ΔBIC, knots, S/N."""
    yy, xx = np.indices(img.shape, dtype=float)
    m = ((xx - x0) ** 2 + (yy - y0) ** 2 <= r_px**2) & np.isfinite(img) & (err > 0)
    X, Y, D, E = xx[m], yy[m], img[m], err[m]
    n = D.size
    f0, f_up = flux_scale(D)

    def r1(p):
        return (_gauss(X, Y, p[0], p[1], p[2], p[3], p[4], p[5]) - D) / E

    best1 = None
    for th in np.linspace(0, np.pi, 4, endpoint=False):
        s = least_squares(
            r1,
            [x0, y0, f0, 2.0, 0.6, th],
            bounds=([x0 - 4, y0 - 4, 0, 0.5, 0.1, -10], [x0 + 4, y0 + 4, f_up, 8, 1, 10]),
        )
        if best1 is None or s.cost < best1.cost:
            best1 = s

    def r2(p):
        return (_gauss(X, Y, p[0], p[1], p[2], p[3]) + _gauss(X, Y, p[4], p[5], p[6], p[7]) - D) / E

    best2 = None
    th1 = best1.x[5]
    for sep in (2.0, 3.0, 4.0, 6.0):
        for dth in (0.0, np.pi / 4, -np.pi / 4, np.pi / 2):
            dx, dy = 0.5 * sep * np.cos(th1 + dth), 0.5 * sep * np.sin(th1 + dth)
            s = least_squares(
                r2,
                [x0 + dx, y0 + dy, f0 / 2, 1.2, x0 - dx, y0 - dy, f0 / 2, 1.2],
                bounds=(
                    [x0 - 6, y0 - 6, 0, 0.5, x0 - 6, y0 - 6, 0, 0.5],
                    [x0 + 6, y0 + 6, f_up, 5, x0 + 6, y0 + 6, f_up, 5],
                ),
            )
            if best2 is None or s.cost < best2.cost:
                best2 = s
    chi1, chi2_ = 2 * best1.cost, 2 * best2.cost
    bic1, bic2 = chi1 + 6 * np.log(n), chi2_ + 8 * np.log(n)
    p = best2.x
    try:
        cov = np.linalg.inv(best2.jac.T @ best2.jac)
        sn1, sn2 = p[2] / np.sqrt(cov[2, 2]), p[6] / np.sqrt(cov[6, 6])
    except np.linalg.LinAlgError:
        sn1 = sn2 = float("nan")
    return {
        "dbic": float(bic1 - bic2),
        "n_pix": int(n),
        "chi2_m1": float(chi1),
        "chi2_m2": float(chi2_),
        "m1": [float(v) for v in best1.x],
        "m2": [float(v) for v in p],
        "snr1": float(sn1),
        "snr2": float(sn2),
        "sep_px": float(np.hypot(p[0] - p[4], p[1] - p[5])),
    }


def forced_fluxes(img, err, knots, r_px, x0, y0):
    """Linear fluxes (and σ) of fixed round Gaussians plus a constant."""
    yy, xx = np.indices(img.shape, dtype=float)
    m = ((xx - x0) ** 2 + (yy - y0) ** 2 <= r_px**2) & np.isfinite(img) & (err > 0)
    A = np.stack(
        [_gauss(xx[m], yy[m], kx, ky, 1.0, ks) for kx, ky, ks in knots] + [np.ones(m.sum())], 1
    )
    w = 1.0 / err[m]
    Aw, bw = A * w[:, None], img[m] * w
    sol, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    cov = np.linalg.inv(Aw.T @ Aw)
    return sol[:-1], np.sqrt(np.diag(cov))[:-1]


# ---- helpers --------------------------------------------------------------------------------
def load(path, filt_arcsec=1.0):
    with fits.open(path) as h:
        sci = np.asarray(h["SCI"].data, float)
        err = np.asarray(h["ERR"].data, float)
        wcs = WCS(h["SCI"].header)
    px = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))) * 3600.0
    sub = sci - median_filter(np.nan_to_num(sci), size=int(round(filt_arcsec / px)) | 1)
    return sci, sub, err, wcs, px


def sky_pa(wcs, x0, y0, x1, y1) -> float:
    return float(wcs.pixel_to_world(x0, y0).position_angle(wcs.pixel_to_world(x1, y1)).deg % 180)


def to_ujy(px):
    """MJy/sr summed over pixels → µJy."""
    return 1e6 * (px / 206264.806) ** 2 * 1e6


def aperture(sub, err, wcs, ra, dec, r_arcsec, px):
    x, y = wcs.world_to_pixel_values(ra, dec)
    yy, xx = np.indices(sub.shape, dtype=float)
    m = (xx - x) ** 2 + (yy - y) ** 2 <= (r_arcsec / px) ** 2
    k = to_ujy(px)
    return float(np.nansum(sub[m]) * k), float(np.sqrt(np.nansum(err[m] ** 2)) * k)


def gate(sci, err, xt, yt, px, r_px, sizes, n_inject, rng, sep_arcsec=0.16) -> dict:
    """Injections at empty positions 1.5–3.5″ away; ``sizes`` gives (flux, σ, q) for ``single``
    and (flux1, σ1, flux2, σ2) for ``two`` (knots ``sep_arcsec`` apart, random axis), image
    units / px."""
    filt = int(round(1.0 / px)) | 1
    half = filt // 2 + int(np.ceil(r_px)) + 6
    res = {"single": [], "two": []}
    drawn = 0
    while min(len(res["single"]), len(res["two"])) < n_inject and drawn < 40 * n_inject:
        drawn += 1
        rr = rng.uniform(1.5, 3.5) / px
        ph = rng.uniform(0, 2 * np.pi)
        x, y = int(round(xt + rr * np.cos(ph))), int(round(yt + rr * np.sin(ph)))
        if x - half < 0 or y - half < 0 or x + half >= sci.shape[1] or y + half >= sci.shape[0]:
            continue
        box = sci[y - half : y + half + 1, x - half : x + half + 1]
        ebox = err[y - half : y + half + 1, x - half : x + half + 1]
        if not np.all(np.isfinite(box)):
            continue
        s0 = box - median_filter(box, size=filt)
        sig = 1.4826 * np.median(np.abs(s0 - np.median(s0)))
        k = 4
        if np.max(np.abs(s0[half - k : half + k + 1, half - k : half + k + 1])) >= 3 * sig:
            continue
        yy, xx = np.indices(box.shape, dtype=float)
        th = rng.uniform(0, np.pi)
        kind = "single" if len(res["single"]) <= len(res["two"]) else "two"
        if kind == "single":
            stamp = _gauss(xx, yy, half, half, *sizes["single"], th)
        else:
            d = 0.5 * sep_arcsec / px
            f1, s1, f2, s2 = sizes["two"]
            c, sn = np.cos(th), np.sin(th)
            stamp = _gauss(xx, yy, half + d * c, half + d * sn, f1, s1) + _gauss(
                xx, yy, half - d * c, half - d * sn, f2, s2
            )
        b = box + stamp
        bs = b - median_filter(b, size=filt)
        w = gaussian_filter(bs, 1.0)[half - 3 : half + 4, half - 3 : half + 4]
        jj, ii = np.unravel_index(np.argmax(w), w.shape)
        f = fit_models(bs, ebox, half - 3 + ii, half - 3 + jj, r_px)
        res[kind].append(deblend_class(f["dbic"], f["snr1"], f["snr2"], f["sep_px"] * px))
    fs = float(np.mean([c != "two-knot" for c in res["single"]])) if res["single"] else float("nan")
    ft = float(np.mean([c == "two-knot" for c in res["two"]])) if res["two"] else float("nan")
    return {
        "n_single": len(res["single"]),
        "n_two": len(res["two"]),
        "frac_single_not_two_knot": round(fs, 3),
        "frac_two_called_two_knot": round(ft, 3),
        "single_classes": {c: res["single"].count(c) for c in sorted(set(res["single"]))},
        "two_classes": {c: res["two"].count(c) for c in sorted(set(res["two"]))},
        "passed": gate_ok(fs, ft),
    }


# ---- 1316 -----------------------------------------------------------------------------------
def run_1316(paths_by_band, row, a, rng, gate_sep_observed=False) -> dict:
    """Deblend one target (1316 in step 2; 1557 in step 3 with ``gate_sep_observed``: the fitted
    gate uses the observed knot separation and is the one that counts, PLAN_1557.md)."""
    sci, sub, err, wcs, px = load(paths_by_band["f200w"])
    xt, yt = wcs.world_to_pixel_values(row["ra"], row["dec"])
    r_px = a.fit_arcsec / px
    yy, xx = np.indices(sub.shape, dtype=float)
    near = np.hypot(xx - xt, yy - yt) <= 3
    j, i = np.unravel_index(
        np.argmax(np.where(near, gaussian_filter(np.nan_to_num(sub), 1.0), -np.inf)), sub.shape
    )
    obs = fit_models(sub, err, float(i), float(j), r_px)
    p = obs["m2"]
    sep_arcsec = obs["sep_px"] * px
    cls = deblend_class(obs["dbic"], obs["snr1"], obs["snr2"], sep_arcsec)
    out = {
        "pixel_arcsec": round(px, 4),
        "observed": obs,
        "sep_arcsec": round(sep_arcsec, 3),
        "class": cls,
    }
    pa_sep = sky_pa(wcs, p[0], p[1], p[4], p[5])
    out["pa_sep"] = round(pa_sep, 1)
    out["offset_sep_pred_z2"] = round(float(abs_.axis_offset(pa_sep, row["pa_pred_z2"])), 1)
    m1 = obs["m1"]
    pa_m1 = sky_pa(wcs, m1[0], m1[1], m1[0] + np.cos(m1[5]), m1[1] + np.sin(m1[5]))
    out["pa_m1"] = round(pa_m1, 1)
    out["offset_m1_pred_z2"] = round(float(abs_.axis_offset(pa_m1, row["pa_pred_z2"])), 1)
    out["m1_q"] = round(m1[4], 3)

    # gate (frozen): catalogue flux, semimajor and ellipticity; knots 0.16″ apart, σ knot_sigma
    tot = 3631.0 * 10 ** (-0.4 * row["mag"]) / to_ujy(px) * 1e6
    sizes = {
        "single": (tot, row["semimajor_px"], 1 - row["ellipticity"]),
        "two": (tot / 2, a.knot_sigma_arcsec / px, tot / 2, a.knot_sigma_arcsec / px),
    }
    out["gate"] = gate(sci, err, xt, yt, px, r_px, sizes, a.n_inject, rng)
    # post hoc (deviation): the fitted total is below the catalogue flux, so repeat the gate with
    # the fitted M1 shape and the fitted M2 knot fluxes and σ (separation still 0.16″)
    out["catalogue_flux_image_units"] = tot
    rr = np.hypot(xx - xt, yy - yt) * px
    bkg = float(np.nanmedian(sci[(rr >= 0.6) & (rr <= 1.0)]))
    out["raw_aperture_0p3_image_units"] = {
        "annulus_bkg": float(np.nansum(sci[rr <= 0.3] - bkg)),
        "filtered": float(np.nansum(sub[rr <= 0.3])),
    }
    out["fitted_flux_image_units"] = {"m1": m1[2], "m2": p[2] + p[6]}
    sizes_fit = {"single": (m1[2], m1[3], m1[4]), "two": (p[2], p[3], p[6], p[7])}
    if gate_sep_observed:
        sep_g = round(sep_arcsec, 2)
        out["gate_fitted"] = gate(
            sci, err, xt, yt, px, r_px, sizes_fit, a.n_inject, rng, sep_arcsec=sep_g
        )
        out["gate_fitted"]["sep_arcsec"] = sep_g
        valid = out["gate_fitted"]["passed"]
    else:
        out["gate_fitted_post_hoc"] = gate(sci, err, xt, yt, px, r_px, sizes_fit, a.n_inject, rng)
        valid = out["gate"]["passed"]

    # knot SEDs (fixed positions and σ from F200W) in the SW bands
    knots_sky = [wcs.pixel_to_world(p[0], p[1]), wcs.pixel_to_world(p[4], p[5])]
    sig_arcsec = [p[3] * px, p[7] * px]
    sed = {}
    for bnd in SW:
        _, s_b, e_b, w_b, px_b = load(paths_by_band[bnd])
        kn = []
        for c, sa in zip(knots_sky, sig_arcsec, strict=True):
            kx, ky = w_b.world_to_pixel(c)
            kn.append((float(kx), float(ky), sa / px_b))
        cx, cy = 0.5 * (kn[0][0] + kn[1][0]), 0.5 * (kn[0][1] + kn[1][1])
        f, e = forced_fluxes(s_b, e_b, kn, a.fit_arcsec / px_b, cx, cy)
        k = to_ujy(px_b)
        sed[bnd] = {"f1_ujy": f[0] * k, "e1_ujy": e[0] * k, "f2_ujy": f[1] * k, "e2_ujy": e[1] * k}
    pr, nb = ratio_consistency(
        *[[sed[b][key] for b in SW] for key in ("f1_ujy", "e1_ujy", "f2_ujy", "e2_ujy")]
    )
    out["knot_sed"] = {b: {k: round(v, 4) for k, v in d.items()} for b, d in sed.items()}
    out["ratio_p"] = pr
    out["ratio_bands"] = nb
    out["ratio_class"] = (
        "undetermined"
        if not np.isfinite(pr)
        else ("consistent SED" if pr >= 0.003 else "different SED")
    )
    if not valid:
        outcome = "test invalid (gate failed)"
    elif cls != "two-knot":
        outcome = "blend not supported; stays open"
    elif out["offset_sep_pred_z2"] >= 60 and out["ratio_class"] == "different SED":
        outcome = "ordinary: blend of two sources"
    elif out["ratio_class"] == "consistent SED":
        outcome = "clump pair; stays open"
    else:
        outcome = "two knots, SED test undetermined or separation not anti; stays open"
    out["outcome"] = outcome
    out["_knots_pix_f200w"] = [(p[0], p[1]), (p[4], p[5])]
    return out


# ---- 2184 -----------------------------------------------------------------------------------
def run_2184(paths_by_band, row, a) -> dict:
    sci, sub, err, wcs, px = load(paths_by_band["f200w"])
    xt, yt = wcs.world_to_pixel_values(row["ra"], row["dec"])
    sm = gaussian_filter(np.nan_to_num(sub), 2.0)
    yy, xx = np.indices(sm.shape, dtype=float)
    rr = np.hypot(xx - xt, yy - yt) * px
    cand = np.where((rr <= 3.0) & (rr >= 0.5), sm, -np.inf)
    pts = []
    for _ in range(2):
        j, i = np.unravel_index(np.argmax(cand), cand.shape)
        pts.append((float(i), float(j)))
        cand[np.hypot(xx - i, yy - j) * px < 0.5] = -np.inf
    pts_sky = [wcs.pixel_to_world(x, y) for x, y in pts]
    bkg = float(np.nanmedian(sci[(rr >= 0.6) & (rr <= 1.0)]))
    cat_flux = 3631.0 * 10 ** (-0.4 * row["mag"]) / to_ujy(px) * 1e6
    flux_check = {
        "catalogue_flux_image_units": cat_flux,
        "raw_aperture_0p3_image_units": float(np.nansum(sci[rr <= 0.3] - bkg)),
    }
    seds = {"2184": [], "arc_a": [], "arc_b": []}
    for bnd in BANDS:
        _, s_b, e_b, w_b, px_b = load(paths_by_band[bnd])
        seds["2184"].append(aperture(s_b, e_b, w_b, row["ra"], row["dec"], a.ap_arcsec, px_b))
        for key, c in zip(("arc_a", "arc_b"), pts_sky, strict=True):
            seds[key].append(aperture(s_b, e_b, w_b, c.ra.deg, c.dec.deg, a.ap_arcsec, px_b))

    def cmp(k1, k2, floor=0.0):
        f1, e1 = zip(*seds[k1], strict=True)
        f2, e2 = zip(*seds[k2], strict=True)
        p, n = sed_match(f1, e1, f2, e2, floor=floor)
        return {"p": p, "n_bands": n, "class": sed_class(p)}

    # post hoc (deviation): a 10 % per-band floor, the plan's own 0.1 mag PSF ASSUMPTION
    post = {
        "floor": 0.10,
        "control_arc_a_vs_b": cmp("arc_a", "arc_b", 0.10),
        "2184_vs_arc_a": cmp("2184", "arc_a", 0.10),
        "2184_vs_arc_b": cmp("2184", "arc_b", 0.10),
    }

    ctrl = cmp("arc_a", "arc_b")
    ta, tb = cmp("2184", "arc_a"), cmp("2184", "arc_b")
    if ctrl["class"] != "SED-consistent":
        outcome = "test invalid (control failed)"
    elif ta["class"] == "SED-different" and tb["class"] == "SED-different":
        outcome = "fragment excluded"
    else:
        outcome = "fragment possible; stays open"
    return {
        "arc_points_radec": [[round(c.ra.deg, 6), round(c.dec.deg, 6)] for c in pts_sky],
        "arc_points_sep_arcsec": [
            round(float(c.separation(wcs.pixel_to_world(xt, yt)).arcsec), 2) for c in pts_sky
        ],
        "sed_ujy": {
            k: {b: [round(f, 4), round(e, 4)] for b, (f, e) in zip(BANDS, v, strict=True)}
            for k, v in seds.items()
        },
        "control_arc_a_vs_b": ctrl,
        "2184_vs_arc_a": ta,
        "2184_vs_arc_b": tb,
        "outcome": outcome,
        "post_hoc_floor": post,
        "flux_check": flux_check,
        "_pts_pix_f200w": pts,
    }


def figure(cut, o1316, o2184, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.5))
    for k, (uid, r) in enumerate((("1316", 0.6), ("2184", 3.5))):
        _, sub, _, wcs, px = load(cut[uid]["f200w"])
        c = sub.shape[0] // 2
        h = int(r / px)
        im = sub[c - h : c + h + 1, c - h : c + h + 1]
        lo, hi = np.nanpercentile(im, [5, 99.5])
        ax[k].imshow(im, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        ax[k].set_title(f"{uid} F200W, 1″-filtered, {2 * r:.1f}″ box")
        if uid == "1316":
            for x, y in o1316["_knots_pix_f200w"]:
                ax[k].plot(x - c + h, y - c + h, "r+", ms=12)
        else:
            for x, y in o2184["_pts_pix_f200w"]:
                ax[k].plot(x - c + h, y - c + h, "co", mfc="none", ms=14)
            ax[k].plot(h, h, "y+", ms=14)
    seds = o2184["sed_ujy"]
    wl = [0.9, 1.15, 1.5, 2.0, 2.77, 3.56, 4.1, 4.44, 4.8]
    for key, col in (("2184", "y"), ("arc_a", "c"), ("arc_b", "b")):
        f = np.array([seds[key][b][0] for b in BANDS])
        e = np.array([seds[key][b][1] for b in BANDS])
        nrm = np.nansum(f * (f / e) ** 2) / np.nansum((f / e) ** 2)
        ax[2].errorbar(wl, f / nrm, e / abs(nrm), fmt="o-", color=col, label=key)
    ax[2].set_xlabel("λ (µm)")
    ax[2].set_ylabel("flux / weighted mean")
    ax[2].legend()
    ax[2].set_title("0.2″ aperture SEDs")
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def _single_target(a, d, cut, row, rng) -> int:
    uid = str(row["label"])
    o = run_1316(cut[uid], row, a, rng, gate_sep_observed=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, len(SW), figsize=(3.2 * len(SW), 3.4))
    for k, bnd in enumerate(SW):
        _, sub, _, wcs, px = load(cut[uid][bnd])
        x, y = wcs.world_to_pixel_values(row["ra"], row["dec"])
        h = int(0.6 / px)
        xi, yi = int(round(float(x))), int(round(float(y)))
        im = sub[yi - h : yi + h + 1, xi - h : xi + h + 1]
        lo, hi = np.nanpercentile(im, [5, 99.5])
        ax[k].imshow(im, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        if bnd == "f200w":
            for kx, ky in o["_knots_pix_f200w"]:
                ax[k].plot(kx - xi + h, ky - yi + h, "r+", ms=12)
        ax[k].set_title(f"{uid} {bnd.upper()}, 1.2″, 1″-filtered")
    fig.tight_layout()
    fig.savefig(d / "vetting" / f"deblend_{uid}.png", dpi=110)
    out = {
        "plan": f"PLAN_{uid}.md",
        "label": "ordinary-explanation test (derived); pixels observed; injections simulations",
        "assumptions": {
            k: getattr(a, k) for k in ("fit_arcsec", "knot_sigma_arcsec", "n_inject", "seed")
        }
        | {"filter_arcsec": 1.0, "images": URI.format(b="<band>")},
        uid: {k: v for k, v in o.items() if not k.startswith("_")},
    }
    (d / "vetting" / f"deblend_{uid}.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n"
    )
    print(json.dumps(out, indent=1, default=float))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", default="abells1063-cats")
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--fit-arcsec", type=float, default=0.4, help="ASSUMPTION")
    p.add_argument("--knot-sigma-arcsec", type=float, default=0.04, help="ASSUMPTION")
    p.add_argument("--ap-arcsec", type=float, default=0.2, help="ASSUMPTION")
    p.add_argument("--n-inject", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261016)
    p.add_argument("--cut-dir", type=Path, default=None, help="reuse cutouts (else a temp dir)")
    p.add_argument(
        "--target", default=None, help="deblend only this label (step 3: 1557, PLAN_1557.md)"
    )
    a = p.parse_args(argv)
    d = a.arcs_dir / a.model
    t = Table.read(d / "arcs_anti.ecsv")
    labels = (a.target,) if a.target else ("1316", "2184")
    rows = {str(r["label"]): r for r in t if str(r["label"]) in labels}
    tg = Table(
        {
            "source_uid": list(rows),
            "ra": [rows[k]["ra"] for k in rows],
            "dec": [rows[k]["dec"] for k in rows],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        base = a.cut_dir or Path(tmp)
        cut = {k: {} for k in rows}
        for b in BANDS:
            c = cutouts.make_cutouts(
                URI.format(b=b), tg, size_arcsec=8.0, out_dir=base / b, extensions=["ERR", "WHT"]
            )
            for uid, pth in zip(c["source_uid"], c["path"], strict=True):
                cut[str(uid)][b] = str(pth)
        rng = np.random.default_rng(a.seed)
        if a.target:
            return _single_target(a, d, cut, rows[a.target], rng)
        o1 = run_1316(cut["1316"], rows["1316"], a, rng)
        o2 = run_2184(cut["2184"], rows["2184"], a)
        figure(cut, o1, o2, d / "vetting" / "deblend.png")
    out = {
        "plan": "PLAN_deblend.md",
        "label": "ordinary-explanation tests (derived); pixels observed; injections simulations",
        "assumptions": {
            k: getattr(a, k)
            for k in ("fit_arcsec", "knot_sigma_arcsec", "ap_arcsec", "n_inject", "seed")
        }
        | {"filter_arcsec": 1.0, "images": URI.format(b="<band>")},
        "1316": {k: v for k, v in o1.items() if not k.startswith("_")},
        "2184": {k: v for k, v in o2.items() if not k.startswith("_")},
    }
    (d / "vetting" / "deblend.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
