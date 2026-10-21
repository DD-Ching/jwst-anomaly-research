"""Pure functions of the line-C HST test (frozen in results/w2_hst_blend/PLAN.md before any target).

N point sources (2-4 quasar images) with one shared Moffat profile plus a constant background,
generalising `w12_niq_hst.py` (D-064 addendum) from a pair to a quad; residual flux inside the
smallest circle about the images' centroid that holds every image (D-056 search circle), outside
0.3" image cores, after the per-image halo correction. The selection rule, the status rule and the
control gate are functions here so that they, not a reading by eye, decide the outcome.
Thresholds are ASSUMPTIONs (``Params``); fits and residuals are ``derived``; injections
``simulated``.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.modeling import fitting, models
from astropy.table import Table
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_niq_hst as hst  # noqa: E402


@dataclass(frozen=True)
class Params:
    # --- selection (computable; PLAN.md "Selection") ---
    min_t814_s: float = 1000.0  # summed F814W ACS/WFC + WFC3/UVIS exposure containing the position
    instruments: tuple = ("ACS/WFC", "WFC3/UVIS")  # 0.04"/px HAP-MVM drizzle, F814W
    max_targets: int = 10
    max_controls: int = 10
    # --- chain ---
    size_px: int = 250  # hapcut cutout, 10" at 0.04"/px
    search_arcsec: float = 3.0  # quasar images within D-056 image_radius of the position
    peak_sigma: float = 20.0  # an image peak is >= 20 sigma_pix above the background
    peak_frac: float = 0.01  # and >= 1 % of the brightest peak
    peak_min_sep_arcsec: float = 0.15  # local maxima closer than this are one source
    conc_r_in_px: float = 2.0  # concentration = flux(r <= 2 px) / flux(r <= 5 px)
    conc_r_out_px: float = 5.0
    conc_frac: float = 0.7  # point-like: concentration >= 0.7 x that of the brightest peak
    max_images: int = 4  # galaxy-scale lensed quasars: 2 or 4 images (brightest kept)
    image_mask_arcsec: float = 0.3  # D-064: image cores excluded from the residual
    min_config_arcsec: float = 0.8  # largest image separation below this: no aperture left
    detect_sigma: float = 5.0  # D-064
    sep_match_arcsec: float = 0.5  # D-064 / D-056 PAIR_SEP_TOL for an expected pair separation
    fit_margin_px: int = 25  # fit stamp: images' extent + 1"
    # --- control gate (D-060) ---
    gate_eff: float = 0.9
    gate_min_decided: int = 5
    # --- injections ---
    inj_mags: tuple = (19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0)
    inj_reff_arcsec: tuple = (0.3, 0.6)
    inj_n: float = 4.0
    inj_ellip: float = 0.2
    # --- ordinary-lens prediction (D-056 Faber-Jackson in LS z; D-064 F814W - z = 0-0.6) ---
    fj: dict = field(
        default_factory=lambda: {"a": 20.41, "k": 0.67, "slope": -10.0, "rms": 0.88, "n": 605}
    )
    f814w_minus_z: tuple = (0.0, 0.6)
    z_l_default: float = 0.5
    z_s_default: float = 2.0


DECIDED = ("light between images", "not detected between images")


# ----------------------------------------------------------------------------- selection
_SVM = re.compile(r"^hst_\d+_[0-9a-z]{2}_.+_([0-9a-z]{6})$")  # HAP single-visit mosaic: one visit
_ASN = re.compile(r"^([0-9a-z]{6})\d{2}0$")  # pipeline association product, e.g. jf5626010
_EXP = re.compile(r"^([0-9a-z]{6})[0-9a-z]{2}[a-z]$")  # individual exposure, e.g. jf5626ryq


def visit_exposure(obs_ids, instruments, exptimes) -> float:
    """F814W depth of one system with each exposure counted once (fix of the v1 rule, which
    summed every level-3 CAOM row and so counted the same exposures 1.9-7 times).

    HST level-3 products overlap: HAP multi-visit mosaics (``hst_skycell-...``: all visits, and
    the same visits again in every overlapping skycell), HAP single-visit mosaics
    (``hst_<prop>_<visit>_..._<ipppss>``: one visit), pipeline association products
    (``<ipppss>nn0``) and individual exposures (``<ipppss>nn<letter>``). Per distinct visit
    (ipppss, instrument) the depth is, in this order of preference: the SVM exposure; else the sum
    of the visit's association products; else the sum of its individual exposures. Multi-visit
    mosaics are never counted (their visits appear as SVM or pipeline rows)."""
    svm: dict[tuple, float] = {}
    asn: dict[tuple, float] = {}
    exp: dict[tuple, float] = {}
    for o, ins, t in zip(obs_ids, instruments, exptimes, strict=True):
        o, t = str(o).lower(), float(t) if np.isfinite(float(t)) else 0.0
        if o.startswith("hst_skycell"):
            continue
        m = _SVM.match(o)
        if m:
            k = (m.group(1), str(ins))
            svm[k] = max(svm.get(k, 0.0), t)
            continue
        for rx, acc in ((_ASN, asn), (_EXP, exp)):
            m = rx.match(o)
            if m:
                k = (m.group(1), str(ins))
                acc[k] = acc.get(k, 0.0) + t
                break
    visits = set(svm) | set(asn) | set(exp)
    return float(sum(svm.get(k, asn.get(k, exp.get(k, 0.0))) for k in visits))


def select_systems(classes: Table, coverage: Table, p: Params) -> Table:
    """Frozen selection rule (only imaging coverage, depth and catalogued model information).

    Pool: D-056 quasar "blended"/"too close" systems whose position lies inside (``contains``) a
    calibrated F814W image from ``p.instruments``, with F814W exposure >= ``min_t814_s``, each
    exposure counted once (``visit_exposure``; v1 summed overlapping products, see RESULT.md).
    Targets: no catalogued lens redshift (``lens_z_known`` False); ordered by summed F814W
    exposure (deepest first), ties by ``system_id``; the first ``max_targets``.
    Controls: a catalogued lens redshift (a published deflector), F814W from a programme that
    also supplied a selected target's F814W (same imaging type); same order; first
    ``max_controls``. Returns one row per selected system with ``role``."""
    cov = Table(coverage, masked=False)
    ok = (
        np.asarray(cov["contains"], bool)
        & (np.asarray(cov["filters"], str) == "F814W")
        & np.isin(np.asarray(cov["instrument"], str), p.instruments)
    )
    cov = cov[ok]
    sid_cov = np.asarray(cov["system_id"], str)
    rows = []
    for r in classes:
        if str(r["status"]) not in ("blended", "too close") or str(r["selection"]) != "quasar":
            continue
        m = sid_cov == str(r["system_id"])
        texp = visit_exposure(
            np.asarray(cov["obs_id"], str)[m],
            np.asarray(cov["instrument"], str)[m],
            np.asarray(cov["t_exptime"], float)[m],
        )
        if texp < p.min_t814_s:
            continue
        progs = sorted({x for x in np.asarray(cov["proposal_id"], str)[m] if x not in ("", "0")})
        rows.append(
            {
                "system_id": str(r["system_id"]),
                "name": str(r["name"]),
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "status_ls": str(r["status"]),
                "ls_sep": float(r["ls_sep"]),
                "theta_e_cat": float(r["theta_e"]),
                "z_lens": float(r["z_lens"]),
                "z_source": float(r["z_source"]),
                "lens_z_known": bool(r["lens_z_known"]),
                "t814_s": texp,
                "programmes": ";".join(progs),
            }
        )
    if not rows:
        return Table(names=("system_id", "role"), dtype=(str, str))
    pool = Table(rows=rows)
    order = np.lexsort((np.asarray(pool["system_id"], str), -np.asarray(pool["t814_s"])))
    pool = pool[order]
    tgt = pool[~np.asarray(pool["lens_z_known"], bool)][: p.max_targets]
    tprog = {x for s in tgt["programmes"] for x in str(s).split(";") if x}
    ctl_pool = pool[np.asarray(pool["lens_z_known"], bool)]
    same = np.array([bool(tprog & set(str(s).split(";"))) for s in ctl_pool["programmes"]], bool)
    ctl = ctl_pool[same][: p.max_controls]
    tgt["role"] = "target"
    ctl["role"] = "control"
    from astropy.table import vstack

    out = vstack([ctl, tgt])
    out.meta.update(provenance="derived", source="w2_hst_blend_chain.select_systems")
    return out


# ----------------------------------------------------------------------------- chain
def concentration(img: np.ndarray, x: float, y: float, bg: float, p: Params) -> float:
    yy, xx = np.indices(img.shape)
    r = np.hypot(xx - x, yy - y)
    a = img - bg
    inner = float(a[r <= p.conc_r_in_px].sum())
    outer = float(a[r <= p.conc_r_out_px].sum())
    return inner / outer if outer > 0 else np.nan


def find_point_sources(
    img: np.ndarray, centre: tuple[float, float], pix: float, sigma_pix: float, p: Params
) -> list[tuple[float, float]]:
    """Quasar images: local maxima (1-px Gaussian smoothing, maxima within the minimum
    separation) inside ``search_arcsec`` of ``centre``, >= ``peak_sigma`` sigma and >=
    ``peak_frac`` of the brightest, and point-like (concentration >= ``conc_frac`` x that of the
    brightest peak, ASSUMPTION: the brightest peak of a lensed quasar is an image). Extended peaks
    (a lens galaxy) are never images. At most ``max_images``, brightest first."""
    bg = float(np.median(img))
    sm = ndimage.gaussian_filter(img, 1.0)
    w = max(3, int(round(2 * p.peak_min_sep_arcsec / pix)) | 1)
    mx = ndimage.maximum_filter(sm, size=w)
    yy, xx = np.indices(img.shape)
    inside = np.hypot(xx - centre[0], yy - centre[1]) <= p.search_arcsec / pix
    cand = (sm == mx) & inside & (sm - bg >= p.peak_sigma * sigma_pix)
    ys, xs = np.nonzero(cand)
    if not len(xs):
        return []
    val = sm[ys, xs] - bg
    o = np.argsort(-val)
    xs, ys, val = xs[o], ys[o], val[o]
    keep = val >= p.peak_frac * val[0]
    xs, ys, val = xs[keep], ys[keep], val[keep]
    c = np.array([concentration(img, x, y, bg, p) for x, y in zip(xs, ys, strict=True)])
    point = np.isfinite(c) & (c >= p.conc_frac * c[0])
    pts = [(float(x), float(y)) for x, y in zip(xs[point], ys[point], strict=True)]
    return pts[: p.max_images]


def n_psf_fit(img: np.ndarray, peaks: list[tuple[float, float]]):
    """N Moffat2D point sources with one shared shape (gamma, alpha) and a constant background
    (``w12_niq_hst.two_psf_fit`` for N = 2). Returns (positions, amplitudes, fitted model,
    residual)."""
    yy, xx = np.indices(img.shape)
    bg = float(np.median(img))
    m = models.Const2D(bg)
    for x0, y0 in peaks:
        amp = float(img[int(round(y0)), int(round(x0))] - bg)
        m = m + models.Moffat2D(
            max(amp, 1e-6),
            x0,
            y0,
            gamma=2.0,
            alpha=2.5,
            bounds={"gamma": (0.5, 10.0), "alpha": (1.2, 8.0), "amplitude": (0, None)},
        )
    for k in range(2, len(peaks) + 1):
        getattr(m, f"gamma_{k}").tied = lambda mm: mm.gamma_1
        getattr(m, f"alpha_{k}").tied = lambda mm: mm.alpha_1
    fit = fitting.LevMarLSQFitter()(m, xx, yy, img, maxiter=1000)
    pos = [
        (float(getattr(fit, f"x_0_{k}").value), float(getattr(fit, f"y_0_{k}").value))
        for k in range(1, len(peaks) + 1)
    ]
    amps = [float(getattr(fit, f"amplitude_{k}").value) for k in range(1, len(peaks) + 1)]
    return pos, amps, fit, img - fit(xx, yy)


def remove_halos_n(res: np.ndarray, pos: list[tuple[float, float]], nbin: float = 1.0):
    """Per image, subtract the radial residual profile measured on the side facing away from the
    other images' centroid, out to the nearest other image (``w12_niq_hst.remove_halos`` for
    N = 2: PSF mismatch is radially symmetric about an image, light between images is not)."""
    out = res.copy()
    yy, xx = np.indices(res.shape)
    for k, (x, y) in enumerate(pos):
        oth = [q for j, q in enumerate(pos) if j != k]
        ox, oy = np.mean([q[0] for q in oth]), np.mean([q[1] for q in oth])
        dn = min(np.hypot(q[0] - x, q[1] - y) for q in oth)
        r = np.hypot(xx - x, yy - y)
        dd = max(np.hypot(ox - x, oy - y), 1e-9)
        cosang = ((xx - x) * (ox - x) + (yy - y) * (oy - y)) / np.maximum(r * dd, 1e-9)
        away = cosang < 0
        near = r < dn
        if not near.any():
            continue
        rbin = (r / nbin).astype(int)
        prof = np.zeros(rbin[near].max() + 2)
        for b in np.unique(rbin[near]):
            sel = near & away & (rbin == b)
            prof[b] = np.median(res[sel]) if sel.any() else 0.0
        out[near] -= prof[rbin[near]]
    return out


def aperture_mask(shape, pos, pix: float, p: Params, valid=None) -> np.ndarray:
    """Smallest circle about the images' centroid holding every image (the pair-diameter circle
    for N = 2), minus ``image_mask_arcsec`` cores."""
    yy, xx = np.indices(shape)
    xc, yc = np.mean([q[0] for q in pos]), np.mean([q[1] for q in pos])
    rad = max(np.hypot(q[0] - xc, q[1] - yc) for q in pos)
    m = np.hypot(xx - xc, yy - yc) <= rad
    core = p.image_mask_arcsec / pix
    for x, y in pos:
        m &= np.hypot(xx - x, yy - y) > core
    if valid is not None:
        m &= valid
    return m


def max_separation(pos) -> float:
    return float(
        max(np.hypot(a[0] - b[0], a[1] - b[1]) for i, a in enumerate(pos) for b in pos[i + 1 :])
    )


def classify(
    n_img: int, max_sep_arcsec: float, snr: float, sep_expected: float, noise: float, p: Params
) -> str:
    """Status rule (frozen). Decided: "light between images" (S/N >= detect_sigma) or "none
    between images" (|S/N| < detect_sigma). Everything else is undecided."""
    if n_img < 2:
        return "unresolved (< 2 point sources)"
    if not np.isfinite(max_sep_arcsec) or max_sep_arcsec < p.min_config_arcsec:
        return "too compact for the core mask"
    if (
        n_img == 2
        and np.isfinite(sep_expected)
        and abs(max_sep_arcsec - sep_expected) > p.sep_match_arcsec
    ):
        return "pair mismatch"
    if not np.isfinite(noise) or noise <= 0:
        return "noise not estimable"
    if not np.isfinite(snr) or snr <= -p.detect_sigma:
        return "fit inconclusive (over-subtraction)"
    # owner decision 2026-10-21: a sub-threshold residual reads "not detected" (v1: "none")
    return "light between images" if snr >= p.detect_sigma else "not detected between images"


def method_validation(gate: dict, psf_summary: dict | None, p: Params) -> dict:
    """Owner decision 2026-10-21, applied after the plan was frozen (stricter than PLAN.md's
    "PSF mismatch not excluded" wording): the method is validated only if the control gate passed
    AND at least one negative control (no-lens empirical-PSF stamp) exists AND every negative
    control stays below ``detect_sigma``. Otherwise target residuals are not interpreted."""
    fp = (psf_summary or {}).get("epsf_no_lens_max_snr")
    if not gate.get("passed"):
        verdict = "failed (control gate)"
    elif fp is None or not np.isfinite(fp):
        verdict = "incomplete (no negative control with an empirical PSF)"
    elif fp >= p.detect_sigma:
        verdict = "failed (PSF-mismatch false positive on a negative control)"
    else:
        verdict = "passed"
    return {
        "verdict": verdict,
        "validated": verdict == "passed",
        "gate_efficiency": gate.get("efficiency"),
        "negative_control_max_snr": fp,
        "detect_sigma": p.detect_sigma,
        "target_interpretation": "interpreted under PLAN.md"
        if verdict == "passed"
        else "not interpreted (method validation failed)",
        "rule": "owner decision 2026-10-21, applied after the frozen plan",
    }


def control_gate(statuses, p: Params) -> dict:
    """D-060 gate: the fraction of decided controls showing light between the images must be
    >= ``gate_eff`` with >= ``gate_min_decided`` decided controls."""
    s = np.asarray(list(statuses), dtype=object)
    n_light = int(np.sum(s == DECIDED[0]))
    n_dec = int(np.isin(s, DECIDED).sum())
    eff = n_light / n_dec if n_dec else np.nan
    return {
        "n_controls": len(s),
        "n_decided": n_dec,
        "n_light": n_light,
        "efficiency": eff,
        "passed": bool(n_dec >= p.gate_min_decided and eff >= p.gate_eff),
    }


def measure(
    img, valid, centre, pix, p: Params, sigma_pix, noise_fn, sep_expected=np.nan, peaks=None
) -> dict:
    """Whole chain on one stamp: point sources, N-Moffat fit, halo correction, residual S/N in
    the aperture, status. ``noise_fn(npix)`` returns the per-pixel aperture noise for an aperture
    of ``npix`` pixels (``w12_niq_hst`` rule: the larger of MAD x correlation and the empirical
    sky-aperture scatter)."""
    pk = find_point_sources(img, centre, pix, sigma_pix, p) if peaks is None else peaks
    out = {"n_img": len(pk), "max_sep_arcsec": np.nan, "resid_snr": np.nan, "resid_flux": np.nan}
    out.update(resid_err=np.nan, n_pix=0, resid_snr_raw=np.nan, pos=pk, amps=[])
    if len(pk) < 2:
        out["status"] = classify(len(pk), np.nan, np.nan, sep_expected, np.nan, p)
        return out
    # fit on a stamp around the images (the cutout holds unrelated sources)
    xs, ys = [q[0] for q in pk], [q[1] for q in pk]
    x0 = max(int(min(xs)) - p.fit_margin_px, 0)
    y0 = max(int(min(ys)) - p.fit_margin_px, 0)
    x1 = min(int(max(xs)) + p.fit_margin_px + 1, img.shape[1])
    y1 = min(int(max(ys)) + p.fit_margin_px + 1, img.shape[0])
    sub, vs = img[y0:y1, x0:x1], valid[y0:y1, x0:x1]
    pos, amps, _fit, res = n_psf_fit(sub, [(x - x0, y - y0) for x, y in pk])
    msep = max_separation(pos) * pix
    ap = aperture_mask(sub.shape, pos, pix, p, vs)
    npix = int(ap.sum())
    noise = noise_fn(npix) if npix else np.nan
    clean = remove_halos_n(res, pos)
    flux = float(clean[ap].sum())
    err = float(noise * np.sqrt(npix)) if npix else np.nan
    snr = flux / err if err and np.isfinite(err) and err > 0 else np.nan
    raw = float(res[ap].sum()) / err if err and np.isfinite(err) and err > 0 else np.nan
    out.update(
        psf_gamma=float(_fit.gamma_1.value),
        psf_alpha=float(_fit.alpha_1.value),
        max_sep_arcsec=msep,
        resid_snr=snr,
        resid_flux=flux,
        resid_err=err,
        n_pix=npix,
        resid_snr_raw=raw,
        pos=[(x + x0, y + y0) for x, y in pos],
        amps=amps,
        stamp=(y0, y1, x0, x1),
        _sub=sub,
        _res=res,
        _clean=clean,
    )
    out["status"] = classify(len(pk), msep, snr, sep_expected, noise, p)
    return out


def lens_position(pos, amps) -> tuple[float, float]:
    """Injection position: for a double the SIS position on the line between the images at
    sep * f / (1 + f) from the fainter one (f = faint/bright, as D-064); else the centroid."""
    if len(pos) == 2:
        (xb, yb), (xf, yf) = (pos[0], pos[1]) if amps[0] >= amps[1] else (pos[1], pos[0])
        q = min(amps) / max(amps)
        return xf + (xb - xf) * q / (1 + q), yf + (yb - yf) * q / (1 + q)
    return float(np.mean([q[0] for q in pos])), float(np.mean([q[1] for q in pos]))


def moffat_kernel(gamma: float, alpha: float, half: int = 25) -> np.ndarray:
    yy, xx = np.indices((2 * half + 1, 2 * half + 1))
    k = models.Moffat2D(1.0, half, half, gamma=gamma, alpha=alpha)(xx, yy)
    return k / k.sum()


def sersic_unit(shape, x, y, reff_px, n, ellip, theta=0.5, psf=None) -> np.ndarray:
    """Unit-flux Sersic lens, 5x5 subsampled per pixel (the n = 4 cusp is not a point source)
    and convolved with ``psf`` (the chain's fitted Moffat) when given: a real lens is seen
    through the PSF. (D-064 injected an unconvolved, pixel-centre Sersic.)"""
    ss = 5
    yy, xx = np.indices((shape[0] * ss, shape[1] * ss))
    xs, ys = (xx + 0.5) / ss - 0.5, (yy + 0.5) / ss - 0.5
    g = models.Sersic2D(1, reff_px, n, x_0=x, y_0=y, ellip=ellip, theta=theta)(xs, ys)
    g = g.reshape(shape[0], ss, shape[1], ss).mean(axis=(1, 3))
    if psf is not None:
        from scipy.signal import fftconvolve

        g = fftconvolve(g, psf, mode="same")
    return g / g.sum()


def recovered(base: dict, inj: dict, p: Params) -> bool:
    """An injected lens is recovered when the chain ends in "light between images", finds the
    same number of images, and the S/N rises >= detect_sigma above the no-injection baseline."""
    b = base["resid_snr"] if np.isfinite(base["resid_snr"]) else 0.0
    return bool(
        inj["status"] == DECIDED[0]
        and inj["n_img"] == base["n_img"]
        and inj["resid_snr"] - b >= p.detect_sigma
    )


def contiguous_limit(passed) -> float:
    return hst.contiguous_limit(passed)


def predicted_f814w(theta_e, z_l, z_s, p: Params) -> dict:
    """Ordinary lens F814W for an SIS of Einstein radius theta_E (D-056 FJ in LS z,
    ``model_prediction``), typical and 2-rms under-luminous, with F814W - z = 0.6 (the fainter end
    of D-064's 0-0.6 assumption)."""
    from jwst_anomaly import lenscats

    cal = lenscats.FJCalibration(band="ls_z", **p.fj)
    zl = z_l if np.isfinite(z_l) else p.z_l_default
    zs = z_s if np.isfinite(z_s) else p.z_s_default
    if not zl < zs - 0.1:
        return {"m_typ": np.nan, "m_faint": np.nan}
    sig = lenscats.sis_sigma(theta_e, zl, zs)
    m = float(cal.mag(sig, zl)) + p.f814w_minus_z[1]
    return {"m_typ": m, "m_faint": m + 2 * cal.rms}
