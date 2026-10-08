"""W3 in the OGLE-IV Mróz et al. microlensing samples: ordinary vs exotic fits and limits (D-057).

Subcommands (outputs under ``$JWST_ANOMALY_DATA/derived/w3_ogle/``, never in git):

- ``fit``: every event of a sample (``jwst_anomaly.ogle``) is fitted with ordinary models — PSPL;
  FSPL (MulensModel, uniform disk) when the PSPL u0 < ``Params.fspl_u0_max``; PSPL with annual
  parallax (MulensModel) when t_E ≥ ``Params.parallax_te_min`` — and with exotic models from
  ``exotic_sim``: ``N1neg`` (n = 1, ε < 0, finite source, ρ free), ``E2pos`` (Ellis, n = 2) and
  ``E2neg`` (n = 2, ε < 0, ρ free). Every model shares the source trajectory (straight line; the
  same as ``MulensModel.Model.get_trajectory`` without parallax, asserted in the tests) and gets
  its source and blend fluxes from the same weighted linear least squares (F_b ≥ −F_min,
  F_s ≥ 0, as in the published selection). One ``derived`` row per event: χ², dof, BIC and
  parameters per model, and ΔBIC = BIC(exotic) − min BIC(ordinary).
- ``vet``: flags (ΔBIC < ``Params.flag_dbic``, an ASSUMPTION) through ordinary tests, cheapest
  first.
- ``inject``: W3 events (``exotic_sim.inject_light_curve``) injected into real light curves of the
  sample, classified by this fitter and by an emulation of the published selection.
- ``limit``: 95 % upper limit on the W3 rate from the published efficiencies and the injections.

Exotic physics is a hypothesis. A better exotic fit is an anomaly to vet, never a discovery.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import minimize

from jwst_anomaly import exotic_sim as es
from jwst_anomaly import ogle, paths, schema

LN10 = math.log(10.0)


@dataclass(frozen=True)
class Params:
    """Every threshold is an ASSUMPTION unless it is a published selection cut (marked)."""

    flag_dbic: float = -10.0  # ΔBIC below which an exotic fit flags an event
    fspl_u0_max: float = 0.1  # fit FSPL when PSPL u0 is below this (ρ constrained only then)
    parallax_te_min: float = 20.0  # days; annual parallax fitted for longer events
    pie_max: float = 5.0  # ASSUMPTION: |π_E| bound; unbounded fits reached 30–1,400 (D-058)
    fs_near: float = 10.0  # exact finite-source integral within fs_near·ρ of a singular radius
    beta_far: float = 50.0  # A = 1 beyond this impact parameter (|A − 1| < 1e-5)
    log_rho_bounds: tuple = (-3.5, 0.0)
    te_bounds: tuple = (0.1, 5000.0)
    maxfev: int = 800
    # published selection (Mróz et al. 2019, Table 2): the cuts below are theirs
    sel_chi2_out: float = 2.0
    sel_n3: int = 3
    sel_chi3: float = 32.0
    sel_amp_mag: float = 0.1
    sel_chi2_fit: float = 2.0
    sel_u0_max: float = 1.0
    sel_te_max: float = 300.0
    sel_is_max: float = 21.0
    sel_fs_min: float = 0.01


P = Params()
ORDINARY = ("PSPL", "FSPL", "PAR")
EXOTIC = {"N1neg": (1.0, -1), "E2pos": (2.0, 1), "E2neg": (2.0, -1)}
N_FLUX = 2  # fs, fb from the linear solve
N_NONLIN = {"PSPL": 3, "FSPL": 4, "PAR": 5, "N1neg": 4, "E2pos": 4, "E2neg": 4}


def out_dir() -> Path:
    d = paths.data_root() / "derived" / "w3_ogle"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ----------------------------------------------------------------------------- data


def mag_to_flux(mag, err):
    """OGLE flux (1 = I of 18 mag, the published convention) and its error."""
    f = 10.0 ** (-0.4 * (np.asarray(mag, float) - ogle.OGLE_FLUX_ZP))
    return f, f * 0.4 * LN10 * np.asarray(err, float)


def flux_to_mag(f):
    return ogle.OGLE_FLUX_ZP - 2.5 * np.log10(f)


# ----------------------------------------------------------------------------- magnification


def pspl(beta):
    """Paczyński point-source point-lens magnification (MulensModel's ``point_source``)."""
    b = np.maximum(np.asarray(beta, float), 1e-12)
    b2 = b * b
    return (b2 + 2.0) / (b * np.sqrt(b2 + 4.0))


@lru_cache(maxsize=8)
def _n2_table(sign: int):
    """Point-source total magnification for n = 2 on a dense log grid (exotic_sim, cached)."""
    if sign == 1:
        lb = np.linspace(-6.0, 3.0, 20001)
        a = es.total_magnification(10.0**lb, 2.0, 1)
        return lb, np.log10(a)
    bc = es.caustic_beta(2.0)
    ld = np.linspace(-12.0, 3.0, 20001)  # log10(beta - beta_c)
    a = es.total_magnification(bc + 10.0**ld, 2.0, -1)
    return ld, np.log10(a)


def point_magnification(beta, n: float, sign: int) -> np.ndarray:
    """Fast point-source total magnification equal to ``exotic_sim.total_magnification``."""
    b = np.asarray(beta, float)
    if n == 1 and sign == 1:
        return pspl(b)
    if n == 1 and sign == -1:
        out = np.zeros_like(b)
        ok = b > 2.0
        bo = b[ok]
        out[ok] = (bo * bo - 2.0) / (bo * np.sqrt(bo * bo - 4.0))
        return out
    if n == 2 and sign == 1:
        lb, la = _n2_table(1)
        x = np.log10(np.maximum(b, 1e-300))
        out = 10.0 ** np.interp(x, lb, la)
        lo = x < lb[0]
        out[lo] = 10.0 ** la[0] * 10.0 ** lb[0] / np.maximum(b[lo], 1e-300)  # A ∝ 1/β
        out[x > lb[-1]] = 1.0
        return out
    if n == 2 and sign == -1:
        ld, la = _n2_table(-1)
        d = b - es.caustic_beta(2.0)
        out = np.zeros_like(b)
        ok = d > 0
        x = np.log10(d[ok])
        v = 10.0 ** np.interp(x, ld, la)
        v[x < ld[0]] = 10.0 ** la[0] * np.sqrt(10.0 ** ld[0] / d[ok][x < ld[0]])
        v[x > ld[-1]] = 1.0
        out[ok] = v
        return out
    return es.total_magnification(b, n, sign)


def exotic_magnification(beta, n: float, sign: int, rho: float, near: float | None = None):
    """Finite-source total magnification: exact ``exotic_sim`` integral near singular radii
    (the caustic of a repulsive lens, β ≈ 0 of an attractive one), point source elsewhere."""
    near = P.fs_near if near is None else near
    b = np.asarray(beta, float)
    a = np.ones_like(b)
    close = b < P.beta_far
    a[close] = point_magnification(b[close], n, sign)
    if rho > 0:
        if sign == -1:
            sel = np.abs(b - es.caustic_beta(n)) < near * rho
        else:
            sel = b < near * rho
        if sel.any():
            a[sel] = es.finite_source_magnification(
                b[sel],
                rho,
                n,
                sign,
                n_nodes=32,
                point_magnification=lambda x: point_magnification(x, n, sign),
            )
    return a


def straight_beta(t, t0, te, u0):
    tau = (np.asarray(t, float) - t0) / te
    return np.sqrt(u0 * u0 + tau * tau)


# ----------------------------------------------------------------------------- linear fluxes


def linear_fluxes(a, f, w):
    """Weighted least-squares F = fs·A + fb with fb ≥ −F_min and fs ≥ 0 (published convention).

    ``a`` has shape (N,) or (K, N); returns fs, fb, χ² with shape () or (K,).
    """
    a = np.atleast_2d(a)
    sw = w.sum()
    sa = a @ w
    saa = (a * a) @ w
    sf = f @ w
    saf = a @ (w * f)
    det = saa * sw - sa * sa
    with np.errstate(divide="ignore", invalid="ignore"):
        fs = np.where(det > 0, (saf * sw - sa * sf) / det, 0.0)
        fb = (sf - fs * sa) / sw
        low = fb < -ogle.F_MIN
        fb = np.where(low, -ogle.F_MIN, fb)
        fs = np.where(low, (saf + ogle.F_MIN * sa) / np.where(saa > 0, saa, 1.0), fs)
    neg = ~(fs >= 0)
    fs = np.where(neg, 0.0, fs)
    fb = np.where(neg, sf / sw, fb)
    r = f - fs[:, None] * a - fb[:, None]
    chi2 = (r * r) @ w
    if chi2.size == 1:
        return float(fs[0]), float(fb[0]), float(chi2[0])
    return fs, fb, chi2


# ----------------------------------------------------------------------------- models


class LightCurve:
    """Flux light curve of one event, with the published RA/Dec for parallax."""

    def __init__(self, t, f, sf, ra=None, dec=None, event_id=""):
        self.t = np.asarray(t, float)
        self.f = np.asarray(f, float)
        self.sf = np.asarray(sf, float)
        self.w = 1.0 / self.sf**2
        self.ra, self.dec, self.event_id = ra, dec, event_id
        self.seasons = None  # one-hot season design (with_season_offsets), else one blend flux
        self.n_extra = 0  # linear parameters beyond fs, fb (BIC)

    def with_season_offsets(self, gap_days: float = 60.0, trend: bool = False) -> LightCurve:
        """Copy whose blend flux is free per observing season (gaps > ``gap_days`` split them),
        plus a linear drift per season when ``trend``: the ordinary-systematics model for
        season-to-season zero points and slow baseline (blend or source) variability."""
        out = LightCurve(self.t, self.f, self.sf, self.ra, self.dec, self.event_id)
        sid = np.concatenate([[0], np.cumsum(np.diff(self.t) > gap_days)])
        onehot = (sid[None, :] == np.arange(sid.max() + 1)[:, None]).astype(float)
        cols = [onehot]
        if trend:
            mid = np.array([self.t[sid == k].mean() for k in range(sid.max() + 1)])
            cols.append(onehot * (self.t - mid[sid])[None, :] / 100.0)
        out.seasons = np.vstack(cols)
        out.n_extra = out.seasons.shape[0] - 1
        return out

    @classmethod
    def from_mag(cls, t, mag, err, **kw):
        f, sf = mag_to_flux(mag, err)
        return cls(t, f, sf, **kw)

    def subset(self, keep):
        out = LightCurve(
            self.t[keep], self.f[keep], self.sf[keep], self.ra, self.dec, self.event_id
        )
        if self.seasons is not None:
            out.seasons = self.seasons[:, keep]
            out.seasons = out.seasons[out.seasons.sum(axis=1) > 0]
            out.n_extra = out.seasons.shape[0] - 1
        return out


def _mm():
    import MulensModel as mm  # optional extra `mulens` (D-054)

    return mm


def have_mm() -> bool:
    try:
        _mm()
    except ImportError:
        return False
    return True


def mm_crosscheck(model: str, lc: LightCurve, r: dict) -> float:
    """Max |A_here − A_MulensModel| at the best fit (PSPL, PAR: trajectory; FSPL: WittMao94)."""
    p = {k: r[k] for k in ("t0", "tE", "u0", "rho", "pi_E_N", "pi_E_E", "t0_par") if k in r}
    a = magnification(model, lc, p)
    ref = np.asarray(mm_model(lc, p, rho=(model == "FSPL")).get_magnification(lc.t), float)
    return float(np.max(np.abs(a - ref) / ref))


def mm_model(lc: LightCurve, p: dict, rho: bool = False):
    """The MulensModel model with these parameters (used for trajectories and cross-checks)."""
    from astropy.coordinates import SkyCoord

    mm = _mm()
    par = {"t_0": p["t0"], "u_0": p["u0"], "t_E": p["tE"]}
    kw = {}
    if "pi_E_N" in p:
        par.update(pi_E_N=p["pi_E_N"], pi_E_E=p["pi_E_E"], t_0_par=p["t0_par"])
        kw["coords"] = SkyCoord(lc.ra, lc.dec, unit="deg")
    if rho:
        par["rho"] = p["rho"]
    mod = mm.Model(par, **kw)
    if rho:
        half = p["tE"] * (abs(p["u0"]) + 10.0 * p["rho"] + 0.05)
        mod.set_magnification_methods(
            [p["t0"] - half, "finite_source_uniform_WittMao94", p["t0"] + half]
        )
    return mod


def parallax_basis(lc: LightCurve, t0_par: float) -> np.ndarray:
    """d(x, y)/d(pi_E_N, pi_E_E) of the MulensModel annual-parallax trajectory, shape (4, N).

    MulensModel's trajectory is affine in (pi_E_N, pi_E_E) at fixed t_0_par and coordinates, so
    three trajectories give it exactly; ``trajectory`` then reproduces ``get_trajectory``.
    """
    key = round(t0_par, 4)
    cache = lc.__dict__.setdefault("_pbasis", {})
    if key not in cache:
        base = {"t0": t0_par, "tE": 100.0, "u0": 0.1, "t0_par": t0_par}
        tr = [
            mm_model(lc, {**base, "pi_E_N": a, "pi_E_E": b}).get_trajectory(lc.t)
            for a, b in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        ]
        cache[key] = np.array(
            [tr[1].x - tr[0].x, tr[1].y - tr[0].y, tr[2].x - tr[0].x, tr[2].y - tr[0].y]
        )
    return cache[key]


def trajectory(lc: LightCurve, p: dict) -> tuple[np.ndarray, np.ndarray]:
    """Source position (x, y) in Einstein radii: MulensModel's convention (y = u0 at t0)."""
    x = (lc.t - p["t0"]) / p["tE"]
    y = np.full_like(x, p["u0"])
    if "pi_E_N" in p:
        b = parallax_basis(lc, p["t0_par"])
        x = x + b[0] * p["pi_E_N"] + b[2] * p["pi_E_E"]
        y = y + b[1] * p["pi_E_N"] + b[3] * p["pi_E_E"]
    return x, y


def magnification(model: str, lc: LightCurve, p: dict) -> np.ndarray:
    """Magnification of ``model`` at the epochs of ``lc`` on the shared trajectory.

    PSPL/PAR: Paczyński; FSPL: uniform disk (``exotic_sim`` n = 1, ε > 0, cross-checked against
    MulensModel's WittMao94 at the best fit); exotic: ``exotic_magnification``.
    """
    x, y = trajectory(lc, p)
    beta = np.hypot(x, y)
    if model in ("PSPL", "PAR"):
        return pspl(beta)
    if model == "FSPL":
        return exotic_magnification(beta, 1.0, 1, p["rho"])
    n, sign = EXOTIC[model]
    return exotic_magnification(beta, n, sign, p["rho"])


def _unpack(model: str, x, t0_par=None) -> dict:
    p = {"t0": x[0], "tE": 10.0 ** x[1], "u0": x[2]}
    if model != "PAR":
        p["u0"] = abs(x[2])
    if model in ("FSPL", *EXOTIC):
        p["rho"] = 10.0 ** np.clip(x[3], *P.log_rho_bounds)
    if model == "PAR":
        p.update(pi_E_N=x[3], pi_E_E=x[4], t0_par=t0_par)
    return p


def chi2_of(model: str, lc: LightCurve, p: dict):
    a = magnification(model, lc, p)
    if not np.all(np.isfinite(a)):
        return np.inf, 0.0, 0.0
    if lc.seasons is None:
        fs, fb, chi2 = linear_fluxes(a, lc.f, lc.w)
        return chi2, fs, fb
    m = np.vstack([a, lc.seasons]).T
    sw = np.sqrt(lc.w)
    coef, *_ = np.linalg.lstsq(m * sw[:, None], lc.f * sw, rcond=None)
    if coef[0] < 0:
        return np.inf, 0.0, 0.0
    r = lc.f - m @ coef
    return float(np.sum(r * r * lc.w)), float(coef[0]), float(np.mean(coef[1:]))


def _objective(model, lc, t0_par):
    lo_te, hi_te = np.log10(P.te_bounds[0]), np.log10(P.te_bounds[1])
    tmin, tmax = lc.t.min() - 1000.0, lc.t.max() + 1000.0

    def fun(x):
        if not (lo_te <= x[1] <= hi_te and tmin <= x[0] <= tmax):
            return 1e30
        if (
            len(x) > 3
            and model != "PAR"
            and not (P.log_rho_bounds[0] - 0.5 <= x[3] <= P.log_rho_bounds[1])
        ):
            return 1e30
        if model == "PAR" and math.hypot(x[3], x[4]) > P.pie_max:
            return 1e30
        try:
            c = chi2_of(model, lc, _unpack(model, x, t0_par))[0]
        except (ValueError, ZeroDivisionError, FloatingPointError):
            return 1e30
        return c if np.isfinite(c) else 1e30

    return fun


def _simplex_steps(model, x0):
    te = 10.0 ** x0[1]
    steps = [0.2 * te, 0.15, max(0.05, 0.3 * abs(x0[2]))]
    if model in ("FSPL", *EXOTIC):
        steps.append(0.4)
    if model == "PAR":
        steps += [0.1, 0.1]
    return np.array(steps)


def optimise(model, lc, starts, t0_par=None, n_best=2):
    """Nelder-Mead from the ``n_best`` best starting points, then one restart from the best."""
    fun = _objective(model, lc, t0_par)
    vals = np.array([fun(np.asarray(s, float)) for s in starts])
    order = np.argsort(vals)[:n_best]
    best_x, best_v = None, np.inf
    for i in order:
        x0 = np.asarray(starts[i], float)
        sim = np.vstack([x0, x0 + np.diag(_simplex_steps(model, x0))])
        r = minimize(
            fun,
            x0,
            method="Nelder-Mead",
            options={
                "initial_simplex": sim,
                "maxfev": P.maxfev,
                "xatol": 1e-5,
                "fatol": 1e-3,
            },
        )
        if r.fun < best_v:
            best_x, best_v = r.x, r.fun
    if best_x is not None:  # restart: Nelder-Mead often stalls on a fresh simplex
        sim = np.vstack([best_x, best_x + np.diag(0.3 * _simplex_steps(model, best_x))])
        r = minimize(
            fun,
            best_x,
            method="Nelder-Mead",
            options={"initial_simplex": sim, "maxfev": P.maxfev, "xatol": 1e-6, "fatol": 1e-4},
        )
        if r.fun < best_v:
            best_x, best_v = r.x, r.fun
    p = _unpack(model, best_x, t0_par)
    chi2, fs, fb = chi2_of(model, lc, p)
    k = N_NONLIN[model] + N_FLUX + lc.n_extra
    n = lc.t.size
    return {
        "model": model,
        "chi2": chi2,
        "dof": n - k,
        "k": k,
        "bic": chi2 + k * math.log(n),
        "fs": fs,
        "fb": fb,
        **{key: float(v) for key, v in p.items() if v is not None},
    }


def pspl_starts(t0, te, u0):
    lt = np.log10(max(te, 0.2))
    out = [(t0, lt, u0)]
    for du in (0.02, 0.1, 0.3, 0.7, 1.0):
        for dl in (-0.3, 0.0, 0.3):
            out.append((t0, lt + dl, du))
    return out


def exotic_starts(model, t0, te, u0_pspl):
    """Starting grid for an exotic model around a PSPL-like bump (t0, t_E from the PSPL fit)."""
    n, sign = EXOTIC[model]
    lt = np.log10(max(te, 0.2))
    if sign == -1:
        bc = es.caustic_beta(n)
        u0s = bc * np.array([0.2, 0.6, 0.9, 0.98, 1.01, 1.04, 1.1, 1.25, 1.5])
        # one caustic spike can play the bump: the spike sits t_E sqrt(bc² - u0²) from t0
        offs = (0.0, 1.0, -1.0)
    else:
        u0s = np.array([0.5, 1.0, 2.0]) * max(u0_pspl, 0.01)
        u0s = np.concatenate([u0s, [0.05, 0.3, 0.8]])
        offs = (0.0,)
    out = []
    for u0 in u0s:
        for dl in (-0.5, 0.0, 0.5):
            for lr in (-2.5, -1.2):
                for o in offs:
                    tt = te * 10**dl
                    shift = o * tt * math.sqrt(max(bc**2 - u0**2, 0.0)) if sign == -1 else 0.0
                    out.append((t0 + shift, lt + dl, u0, lr))
    if sign == -1:  # time scales independent of the PSPL fit (it can be far off for W3 shapes)
        for te_abs in (1.0, 3.0, 10.0, 30.0, 100.0, 300.0):
            for u0 in bc * np.array([0.5, 0.95, 0.99, 1.02]):
                for lr in (-2.5, -1.2):
                    for o in (0.0, 1.0, -1.0):  # t0 guess on either caustic spike
                        shift = o * te_abs * math.sqrt(max(bc**2 - u0**2, 0.0))
                        out.append((t0 + shift, math.log10(te_abs), u0, lr))
    return out


def spike_pair_starts(lc: LightCurve, bc: float, n_peaks: int = 5) -> list:
    """Starts that put the two caustic spikes of a repulsive lens on pairs of light-curve maxima."""
    k = np.convolve(lc.f, np.ones(3) / 3, mode="same")
    peaks = []
    for i in np.argsort(k)[::-1]:
        if all(abs(i - j) > 3 for j in peaks):
            peaks.append(int(i))
        if len(peaks) == n_peaks:
            break
    out = []
    for a in range(len(peaks)):
        for b in range(a + 1, len(peaks)):
            ta, tb = sorted((lc.t[peaks[a]], lc.t[peaks[b]]))
            half = 0.5 * (tb - ta)
            if half <= 0:
                continue
            for u0 in bc * np.array([0.3, 0.7, 0.95, 0.99]):
                te = half / math.sqrt(bc**2 - u0**2)
                if P.te_bounds[0] < te < P.te_bounds[1]:
                    for lr in (-2.5, -1.2):
                        out.append((0.5 * (ta + tb), math.log10(te), u0, lr))
    return out


def fit_event(
    lc: LightCurve,
    t0_guess: float,
    te_guess: float,
    u0_guess: float,
    models=None,
    te_grid=(),
):
    """All models for one light curve. Returns {model: result dict}.

    ``te_grid``: extra PSPL starting time scales (injections have no published fit to start from).
    """
    models = models or ("PSPL", "FSPL", "PAR", *EXOTIC)
    res = {}
    starts = pspl_starts(t0_guess, te_guess, u0_guess)
    for te in te_grid:
        starts += pspl_starts(t0_guess, te, u0_guess)
    ps = optimise("PSPL", lc, starts, n_best=3)
    res["PSPL"] = ps
    if "FSPL" in models and ps["u0"] < P.fspl_u0_max:
        lt = math.log10(ps["tE"])
        starts = [(ps["t0"], lt, ps["u0"], math.log10(r)) for r in (0.002, 0.01, 0.03, 0.1)]
        res["FSPL"] = optimise("FSPL", lc, starts)
    if "PAR" in models and ps["tE"] >= P.parallax_te_min and lc.ra is not None and have_mm():
        lt = math.log10(ps["tE"])
        starts = [
            (ps["t0"], lt, s * ps["u0"], pn, pe)
            for s in (1, -1)
            for pn, pe in ((0.0, 0.0), (0.2, 0.0), (-0.2, 0.0), (0.0, 0.2), (0.0, -0.2))
        ]
        res["PAR"] = optimise("PAR", lc, starts, t0_par=round(ps["t0"], 1), n_best=2)
    for m in EXOTIC:
        if m in models:
            starts = exotic_starts(m, ps["t0"], ps["tE"], ps["u0"])
            if EXOTIC[m][1] == -1:
                starts += spike_pair_starts(lc, es.caustic_beta(EXOTIC[m][0]))
            res[m] = optimise(m, lc, starts, n_best=3)
    if have_mm():
        for m in ("PSPL", "FSPL", "PAR"):
            if m in res:
                res[m]["mm_check"] = mm_crosscheck(m, lc, res[m])
    return res


def summarise(event_id: str, res: dict, n: int) -> dict:
    """One flat row: χ², dof, BIC and parameters per model; ΔBIC per exotic model."""
    row = {"event_id": event_id, "n_points": n}
    best_ord = min(res[m]["bic"] for m in ORDINARY if m in res)
    row["best_ordinary"] = min((m for m in ORDINARY if m in res), key=lambda m: res[m]["bic"])
    for m, r in res.items():
        for key in (
            "chi2",
            "dof",
            "bic",
            "t0",
            "tE",
            "u0",
            "rho",
            "pi_E_N",
            "pi_E_E",
            "fs",
            "fb",
            "mm_check",
        ):
            row[f"{m}_{key}"] = r.get(key, np.nan)
    for m in EXOTIC:
        row[f"dbic_{m}"] = res[m]["bic"] - best_ord if m in res else np.nan
    row["dbic_min"] = np.nanmin([row[f"dbic_{m}"] for m in EXOTIC])
    return row


# ----------------------------------------------------------------------------- fit stage


def _fit_worker(job):
    ev, t, mag, err = job
    t1 = time.time()
    lc = LightCurve.from_mag(t, mag, err, ra=ev["ra"], dec=ev["dec"], event_id=ev["event_id"])
    try:
        res = fit_event(lc, ev["t0_pub"], ev["tE_pub"], ev["u0_pub"])
        row = summarise(ev["event_id"], res, lc.t.size)
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001 — one bad light curve must not stop the run
        row = {"event_id": ev["event_id"], "n_points": lc.t.size, "error": repr(exc)[:200]}
    row["seconds"] = time.time() - t1
    return row


def _jobs(sample: ogle.OgleMrozSample, limit: int | None, ids=None, skip=(), chunk=None):
    ev = sample.events()
    if chunk is not None:  # (k, n): every n-th event from k, so each chunk spans all fields
        ev = ev[chunk[0] :: chunk[1]]
    if ids is not None:
        ev = ev[np.isin(ev["event_id"], list(ids))]
    if limit:
        ev = ev[:limit]
    if skip:
        ev = ev[~np.isin(ev["event_id"], list(skip))]
    for r in ev:
        try:
            lc = sample.light_curve(r["event_id"])
        except KeyError:  # listed in the table but no photometry file: report, don't stop the run
            print(f"skipped {r['event_id']}: no light curve", flush=True)
            continue
        d = {k: (r[k].item() if hasattr(r[k], "item") else r[k]) for k in ev.colnames}
        yield d, np.asarray(lc["time"]), np.asarray(lc["mag"]), np.asarray(lc["mag_err"])


def load_checkpoint(path: Path) -> list[dict]:
    """Rows of an interrupted ``fit`` (one JSON object per line); a torn last line is dropped."""
    rows = []
    if path.exists():
        for line in path.read_text().splitlines():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                break
    return rows


def params_tag() -> str:
    """Short hash of the fitting ``Params``; checkpoint rows from other Params are not reused."""
    return hashlib.sha256(json.dumps(asdict(P), sort_keys=True).encode()).hexdigest()[:12]


def parse_chunk(text: str | None) -> tuple[int, int] | None:
    """``"K/N"`` (1-based K) -> ``(K - 1, N)``; ``None`` passes through."""
    if text is None:
        return None
    k, n = (int(x) for x in text.split("/"))
    if not 1 <= k <= n:
        raise ValueError(f"chunk {text!r}: need 1 <= K <= N")
    return k - 1, n


def run_fit(
    sample_key: str,
    limit: int | None,
    procs: int,
    fresh: bool = False,
    chunk: tuple[int, int] | None = None,
) -> Path:
    """Fit the sample. Each row is appended to ``fits_<key>.partial.jsonl`` as it finishes, so an
    interrupted run resumes where it stopped (a cloud session ends after ~40 min; the bulge sample
    takes ~4 h on 4 cores, ~10 s per event); ``fresh`` starts over. ``chunk=(k, n)`` fits only
    events k, k+n, ... (0-based), so separate sessions fit disjoint, field-balanced parts."""
    sample = ogle.OgleMrozSample(sample_key)
    t1 = time.time()
    ckpt = out_dir() / f"fits_{sample_key}.partial.jsonl"
    if fresh:
        ckpt.unlink(missing_ok=True)
    tag = params_tag()
    rows = load_checkpoint(ckpt)
    stale = [r for r in rows if r.get("params_tag") != tag]
    if stale:  # fitted under other Params (e.g. before the D-058 π_E bound): refit them
        print(f"dropping {len(stale)} checkpointed events fitted with other Params", flush=True)
        rows = [r for r in rows if r.get("params_tag") == tag]
    ckpt.write_text("".join(json.dumps(r, default=float) + "\n" for r in rows))
    if chunk is not None:  # other chunks' rows stay in the checkpoint but not in this chunk's table
        ids = set(sample.events()["event_id"][chunk[0] :: chunk[1]].tolist())
        rows = [r for r in rows if r["event_id"] in ids]
    done = {r["event_id"] for r in rows}
    if done:
        print(f"resuming: {len(done)} events already fitted", flush=True)
    with Pool(procs) as pool, ckpt.open("a") as fh:
        jobs = _jobs(sample, limit, skip=done, chunk=chunk)
        for i, row in enumerate(pool.imap_unordered(_fit_worker, jobs, 4)):
            row["params_tag"] = tag
            rows.append(row)
            fh.write(json.dumps(row, default=float) + "\n")
            fh.flush()
            if (i + 1) % 250 == 0:
                print(f"{i + 1} events, {time.time() - t1:.0f} s", flush=True)
    keys = sorted({k for r in rows for k in r} - {"params_tag"}, key=lambda k: (k != "event_id", k))
    tab = Table(
        {
            k: [r.get(k, np.nan if k not in ("error", "best_ordinary") else "") for r in rows]
            for k in keys
        }
    )
    ev = sample.events()
    tab = _join_pub(tab, ev)
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/w3_microlensing.py fit on {sample.name} ({sample.spec.reference})",
        params=json.dumps(asdict(P)),
        chunk="" if chunk is None else f"{chunk[0] + 1}/{chunk[1]}",
        wall_time_s=round(time.time() - t1, 1),  # this invocation only
        cpu_time_s=round(float(sum(r.get("seconds", 0.0) for r in rows)), 1),
        procs=procs,
    )
    path = out_dir() / f"fits_{sample_key}.ecsv"  # what `vet` / `sheet` / `summary` read
    tab.write(path, overwrite=True)
    if chunk is not None:  # keep each chunk's table; the next chunk replaces `path`
        tab.write(
            out_dir() / f"fits_{sample_key}_chunk{chunk[0] + 1}of{chunk[1]}.ecsv", overwrite=True
        )
    print(f"wrote {path}: {len(tab)} events, {time.time() - t1:.0f} s wall")
    return path


def no_error(tab: Table) -> np.ndarray:
    """Rows without a fit error (ECSV reads empty strings back as masked)."""
    col = tab["error"]
    vals = col.filled("") if hasattr(col, "filled") else col
    return np.asarray(vals == "", bool)


def _join_pub(tab: Table, ev: Table) -> Table:
    idx = {e: i for i, e in enumerate(ev["event_id"])}
    j = np.array([idx[e] for e in tab["event_id"]])
    for c in ("field", "ra", "dec", "t0_pub", "tE_pub", "u0_pub", "Is_pub", "fs_pub", "alt_id"):
        tab[c] = ev[c][j]
    return tab


# ----------------------------------------------------------------------------- published selection


def published_selection(
    lc: LightCurve, fit: dict | None = None, sample_key: str = "bulge2019"
) -> dict:
    """Emulation of the Mróz et al. 2019 low-cadence selection (their Table 2) on one light curve.

    Emulated: χ²_out/dof ≤ 2 outside a 720 d window (360 d if the light curve spans < 6 yr) centred
    on the brightest run, ≥ 3 consecutive points ≥ 3σ_base above the 5σ-clipped baseline, χ_3+ ≥ 32,
    amplitude ≥ 0.1 mag, one bump, and the PSPL fit cuts (χ²/dof ≤ 2 overall and within t_E of t0,
    t0 window, u0 ≤ 1, t_E ≤ 300 d, I_s ≤ 21, F_b > −F_min, f_s > 0.01). Not emulated: n_DIA ≥ 3
    (centroid on the difference images), the s < 0.4 artifact statistic, and the removal of
    neighbours brightened together. How bumps are counted is not published: here a second bump is
    another run of ≥ 3 consecutive ≥ 3σ points separated from the first by ≥ 3 consecutive points
    below F_base + 3σ_base whose minimum is below F_base (the flux returns to the baseline).
    """
    t, f, sf = lc.t, lc.f, lc.sf
    out = {}
    span = t.max() - t.min()
    half = 360.0 if span >= 6 * 365.25 else 180.0
    # brightest 3-point running mean defines the window centre
    k = np.convolve(f, np.ones(3) / 3, mode="same")
    tc = t[int(np.argmax(k))]
    outside = np.abs(t - tc) > half
    if outside.sum() < 10:
        outside = np.abs(t - tc) > 0.5 * half
    fo = f[outside]
    for _ in range(5):
        med, sd = np.mean(fo), np.std(fo)
        keep = np.abs(fo - med) < 5 * sd
        if keep.all():
            break
        fo = fo[keep]
    fbase, sbase = float(np.mean(fo)), float(np.std(fo))
    so = sf[outside]
    out["chi2_out"] = float(np.sum((f[outside] - fbase) ** 2 / so**2) / max(outside.sum() - 1, 1))
    hi = f >= fbase + 3.0 * sbase
    runs = _runs(hi)
    good = [(a, b) for a, b in runs if b - a >= P.sel_n3]
    out["chi3"] = float(sum(np.sum((f[a:b] - fbase) / sf[a:b]) for a, b in good))
    out["n_bump"] = _count_bumps(good, f, fbase, sbase)
    peak = f.max()
    out["amp_mag"] = float(2.5 * np.log10(peak / fbase)) if fbase > 0 and peak > 0 else 0.0
    passed = {
        "chi2_out": out["chi2_out"] <= P.sel_chi2_out,
        "n3": len(good) > 0,
        "chi3": out["chi3"] >= P.sel_chi3,
        "amp": out["amp_mag"] >= P.sel_amp_mag,
        "one_bump": out["n_bump"] == 1,
    }
    if fit is not None:
        n = t.size
        a = pspl(straight_beta(t, fit["t0"], fit["tE"], fit["u0"]))
        model = fit["fs"] * a + fit["fb"]
        near = np.abs(t - fit["t0"]) < fit["tE"]
        chi_all = np.sum((f - model) ** 2 / sf**2) / max(n - 5, 1)
        chi_te = np.sum((f[near] - model[near]) ** 2 / sf[near] ** 2) / max(near.sum() - 5, 1)
        i_s = flux_to_mag(fit["fs"]) if fit["fs"] > 0 else np.inf
        fs_frac = fit["fs"] / (fit["fs"] + fit["fb"]) if fit["fs"] + fit["fb"] > 0 else 0.0
        passed.update(
            chi2_fit=chi_all <= P.sel_chi2_fit,
            chi2_fit_te=(near.sum() <= 5) or chi_te <= P.sel_chi2_fit,
            t0=ogle.SAMPLES[sample_key].t_min <= fit["t0"] <= ogle.SAMPLES[sample_key].t_max,
            u0=fit["u0"] <= P.sel_u0_max,
            te=fit["tE"] <= P.sel_te_max,
            i_s=i_s <= P.sel_is_max,
            fb=_blend_ok(lc, fit),
            fs=fs_frac > P.sel_fs_min,
        )
        out.update(chi2_fit=float(chi_all), chi2_fit_te=float(chi_te), Is=float(i_s))
    out["passed"] = passed
    out["selected"] = all(passed.values())
    return out


def _blend_ok(lc: LightCurve, fit: dict) -> bool:
    """F_b > −F_min, or, at the bound, the F_b = 0 (four-parameter) model is worse by Δχ² < 9:
    Mróz et al. 2019 then keep the four-parameter fit."""
    if fit["fb"] > -ogle.F_MIN + 1e-9:
        return True
    a = pspl(straight_beta(lc.t, fit["t0"], fit["tE"], fit["u0"]))
    fs4 = float(np.sum(lc.w * a * lc.f) / np.sum(lc.w * a * a))
    chi4 = float(np.sum(lc.w * (lc.f - fs4 * a) ** 2))
    chi5 = float(np.sum(lc.w * (lc.f - fit["fs"] * a - fit["fb"]) ** 2))
    return chi4 - chi5 < 9.0


def _runs(mask):
    """(start, stop) of runs of True."""
    m = np.concatenate([[False], mask, [False]]).astype(int)
    d = np.diff(m)
    return list(zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0], strict=True))


def _count_bumps(good, f, fbase, sbase):
    if not good:
        return 0
    n = 1
    for (_, b1), (a2, _) in zip(good[:-1], good[1:], strict=True):
        gap = f[b1:a2]
        low = gap < fbase + 3.0 * sbase
        if any(b - a >= 3 for a, b in _runs(low)) and gap.min() < fbase:
            n += 1
    return n


# ----------------------------------------------------------------------------- injection


def baseline_from_fit(lc: LightCurve, fit: dict) -> tuple[LightCurve, float, float]:
    """Real cadence and residual noise with the fitted PSPL event removed (``derived``).

    Returns the constant-baseline light curve, the baseline flux and the published-model source
    flux. Errors of magnified points are scaled back to the baseline flux (Poisson-like;
    ASSUMPTION).
    """
    a = pspl(straight_beta(lc.t, fit["t0"], fit["tE"], fit["u0"]))
    model = fit["fs"] * a + fit["fb"]
    base = fit["fs"] + fit["fb"]
    scale = np.sqrt(np.clip(base / np.maximum(model, 1e-6), 0.0, 1.0))
    resid = (lc.f - model) * scale
    sf = lc.sf * scale
    return LightCurve(lc.t, base + resid, sf, lc.ra, lc.dec, lc.event_id), base, fit["fs"]


def inject_w3(
    base_lc: LightCurve,
    base_flux: float,
    source_flux: float,
    t0,
    te,
    u0,
    rho,
    n=1.0,
    sign=-1,
):
    """Exotic event on a constant-baseline real light curve (``simulated``).

    ``exotic_sim.inject_light_curve`` gives the noiseless flux; the residual noise of each epoch is
    rescaled by sqrt(F_new / F_base) with a floor of 0.3 (sky-limited when the source vanishes;
    ASSUMPTION), errors likewise.
    """
    blend = min(max(source_flux / base_flux, 0.0), 1.0)
    inj = es.inject_light_curve(base_lc.t, base_flux, t0, te, u0, n, sign, rho, blend)
    fnew = np.asarray(inj["flux"], float)
    scale = np.sqrt(np.clip(fnew / base_flux, 0.09, None))
    noise = (base_lc.f - base_flux) * scale
    lc = LightCurve(
        base_lc.t, fnew + noise, base_lc.sf * scale, base_lc.ra, base_lc.dec, base_lc.event_id
    )
    return lc


INJ_TE = (3.0, 10.0, 30.0, 100.0, 300.0)  # days (ASSUMPTION: grid)
INJ_RHO = (0.01, 0.1)
INJ_PER_CELL = 150


def _pick_bases(fits: Table, n: int, seed: int) -> list[str]:
    """Events whose PSPL fit is good (χ²/dof ≤ 1.5) and whose source is bright enough to matter."""
    ok = (
        no_error(fits)
        & (fits["PSPL_chi2"] / fits["PSPL_dof"] <= 1.5)
        & (fits["PSPL_fs"] > 0)
        & (fits["n_points"] >= 100)
    )
    ids = np.asarray(fits["event_id"][ok])
    rng = np.random.default_rng(seed)
    return list(rng.choice(ids, size=min(n, ids.size), replace=False))


def _inject_worker(job):
    kind, base_lc, base_flux, fs, prm = job
    if kind == "W3":
        lc = inject_w3(base_lc, base_flux, fs, prm["t0"], prm["tE"], prm["u0"], prm["rho"])
    else:  # PSPL control: same machinery, n = 1, ε > 0, point source
        lc = inject_w3(base_lc, base_flux, fs, prm["t0"], prm["tE"], prm["u0"], 0.0, 1.0, 1)
    k = np.convolve(lc.f, np.ones(3) / 3, mode="same")
    t0g = float(lc.t[int(np.argmax(k))])
    row = {"kind": kind, "event_id": base_lc.event_id, **prm}
    try:
        models = ("PSPL", "FSPL", "PAR", "N1neg") if kind == "W3" else ("PSPL",)
        res = fit_event(lc, t0g, 10.0, 0.3, models=models, te_grid=(3.0, 30.0, 100.0))
        sel = published_selection(lc, res["PSPL"])
        row.update(selected=sel["selected"], **{f"sel_{k}": v for k, v in sel["passed"].items()})
        if kind == "W3":
            best = min(res[m]["bic"] for m in ORDINARY if m in res)
            row["dbic_N1neg"] = res["N1neg"]["bic"] - best
            row["flagged"] = row["dbic_N1neg"] < P.flag_dbic
            row["flag_vetted"] = row["flagged"] and _survives_cheap_vetting(lc, res)
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row["error"] = repr(exc)[:200]
    return row


def par_starts(ps: dict) -> list:
    """Parallax starts around a PSPL fit (both u0 signs, small π_E offsets)."""
    lt = math.log10(ps["tE"])
    return [
        (ps["t0"], lt, s * ps["u0"], a, b)
        for s in (1, -1)
        for a, b in ((0, 0), (0.3, 0), (-0.3, 0), (0, 0.3), (0, -0.3))
    ]


def _survives_cheap_vetting(lc: LightCurve, res: dict) -> bool:
    """The vetting tests that remove most real-data flags, applied to an injected W3 flag:
    errors rescaled to χ²/dof = 1 of the best ordinary model, and per-season offset + drift."""
    best = min((m for m in ORDINARY if m in res), key=lambda m: res[m]["bic"])
    scale = max(res[best]["chi2"] / res[best]["dof"], 1.0)
    d = (res["N1neg"]["chi2"] - res[best]["chi2"]) / scale + (
        res["N1neg"]["k"] - res[best]["k"]
    ) * math.log(lc.t.size)
    if not d < P.flag_dbic:
        return False
    lct = lc.with_season_offsets(trend=True)
    ps, rx = res["PSPL"], res["N1neg"]
    best_t = optimise("PSPL", lct, [(ps["t0"], math.log10(ps["tE"]), ps["u0"])])["bic"]
    if "PAR" in res:  # as in _vet_one: the season-trend refit always includes parallax
        r = res["PAR"]
        start = (r["t0"], math.log10(r["tE"]), r["u0"], r["pi_E_N"], r["pi_E_E"])
        best_t = min(best_t, optimise("PAR", lct, [start], t0_par=r["t0_par"])["bic"])
    elif have_mm() and lc.ra is not None:  # short t_E: fit_event skipped PAR, _vet_one does not
        best_t = min(best_t, optimise("PAR", lct, par_starts(ps), t0_par=round(ps["t0"], 1))["bic"])
    e = optimise("N1neg", lct, [(rx["t0"], math.log10(rx["tE"]), rx["u0"], math.log10(rx["rho"]))])
    return bool(e["bic"] - best_t < P.flag_dbic)


def run_inject(procs: int, per_cell: int, seed: int = 55) -> Path:
    """Injection-recovery on real bulge light curves (``simulated`` events, ``derived`` result)."""
    sample = ogle.OgleMrozSample("bulge2019")
    fits = Table.read(out_dir() / "fits_bulge2019.ecsv")
    ev = sample.events()
    pub = {r["event_id"]: r for r in ev}
    fit_by_id = {r["event_id"]: r for r in fits}
    rng = np.random.default_rng(seed)
    bases = _pick_bases(fits, 400, seed)
    base_cache = {}
    spec = sample.spec

    def base(eid):
        if eid not in base_cache:
            lc0 = sample.light_curve(eid)
            lc = LightCurve.from_mag(
                lc0["time"],
                lc0["mag"],
                lc0["mag_err"],
                ra=pub[eid]["ra"],
                dec=pub[eid]["dec"],
                event_id=eid,
            )
            r = fit_by_id[eid]
            fit = {k: float(r[f"PSPL_{k}"]) for k in ("t0", "tE", "u0", "fs", "fb")}
            base_cache[eid] = baseline_from_fit(lc, fit)
        return base_cache[eid]

    jobs = []
    for te in INJ_TE:
        for kind, rhos in (("W3", INJ_RHO), ("PSPL", (0.0,))):
            for rho in rhos:
                for _ in range(per_cell):
                    eid = bases[int(rng.integers(len(bases)))]
                    blc, bflux, fs = base(eid)
                    umax = 2.0 if kind == "W3" else 1.0
                    prm = {
                        "tE": te,
                        "rho": rho,
                        "u0": float(rng.uniform(0.0, umax)),
                        "t0": float(rng.uniform(spec.t_min, spec.t_max)),
                        "Is_base": float(flux_to_mag(bflux)),
                        "fs_frac": float(min(fs / bflux, 1.0)),
                    }
                    jobs.append((kind, blc, bflux, fs, prm))
    t1 = time.time()
    rows = []
    with Pool(procs) as pool:
        for i, row in enumerate(pool.imap_unordered(_inject_worker, jobs, 2)):
            rows.append(row)
            if (i + 1) % 200 == 0:
                print(f"{i + 1}/{len(jobs)} injections, {time.time() - t1:.0f} s", flush=True)
    keys = sorted({k for r in rows for k in r})
    fill = {"kind": "", "event_id": "", "error": ""}
    tab = Table({k: [r.get(k, fill.get(k, np.nan)) for r in rows] for k in keys})
    for c in tab.colnames:
        if tab[c].dtype == object:
            tab[c] = [np.nan if v is None else v for v in tab[c]]
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=(
            "scripts/w3_microlensing.py inject: exotic_sim.inject_light_curve (simulated) on "
            f"PSPL-subtracted real light curves of {sample.name}"
        ),
        params=json.dumps(asdict(P)),
        wall_time_s=round(time.time() - t1, 1),
        seed=seed,
    )
    path = out_dir() / "injections_bulge2019.ecsv"
    tab.write(path, overwrite=True)
    print(f"wrote {path}: {len(tab)} rows, {time.time() - t1:.0f} s")
    return path


def selection_audit(procs: int) -> dict:
    """Pass rate of the emulated selection on the real sample (all passed the real one)."""
    sample = ogle.OgleMrozSample("bulge2019")
    fits = Table.read(out_dir() / "fits_bulge2019.ecsv")
    ev = sample.events()
    pub = {r["event_id"]: r for r in ev}
    jobs = []
    for r in fits[no_error(fits)]:
        lc0 = sample.light_curve(r["event_id"])
        e = pub[r["event_id"]]
        fit = {k: float(r[f"PSPL_{k}"]) for k in ("t0", "tE", "u0", "fs", "fb")}
        jobs.append(
            (lc0["time"].data, lc0["mag"].data, lc0["mag_err"].data, e["ra"], e["dec"], fit)
        )
    with Pool(procs) as pool:
        out = pool.map(_audit_worker, jobs, 16)
    keys = out[0].keys()
    rates = {k: float(np.mean([o[k] for o in out])) for k in keys}
    rates["n"] = len(out)
    return rates


def _audit_worker(job):
    t, m, e, ra, dec, fit = job
    lc = LightCurve.from_mag(t, m, e, ra=ra, dec=dec)
    sel = published_selection(lc, fit)
    return {"selected": sel["selected"], **sel["passed"]}


# ----------------------------------------------------------------------------- vetting


def bs_magnification(lc: LightCurve, p: dict) -> np.ndarray:
    """Two PSPL sources with one t_E (binary source; also a xallarap proxy): shape (2, N)."""
    a1 = pspl(straight_beta(lc.t, p["t01"], p["tE"], p["u01"]))
    a2 = pspl(straight_beta(lc.t, p["t02"], p["tE"], p["u02"]))
    return np.vstack([a1, a2])


def linear_fluxes_n(a: np.ndarray, f, w):
    """Unconstrained weighted least squares F = Σ_k fs_k A_k + fb; returns (coeffs, χ²)."""
    m = np.vstack([a, np.ones_like(f)]).T
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(m * sw[:, None], f * sw, rcond=None)
    r = f - m @ coef
    return coef, float(np.sum(r * r * w))


def fit_binary_source(lc: LightCurve, ps: dict) -> dict:
    def fun(x):
        p = {"t01": x[0], "u01": abs(x[1]), "t02": x[2], "u02": abs(x[3]), "tE": 10 ** x[4]}
        coef, chi2 = linear_fluxes_n(bs_magnification(lc, p), lc.f, lc.w)
        return chi2 if coef[0] >= 0 and coef[1] >= 0 else chi2 + 1e6

    te = ps["tE"]
    starts = [
        (ps["t0"], ps["u0"], ps["t0"] + d * te, u2, math.log10(te))
        for d in (-2.0, -1.0, -0.3, 0.3, 1.0, 2.0)
        for u2 in (0.01, 0.1, 0.5)
    ]
    vals = [fun(np.array(s)) for s in starts]
    best = None
    for i in np.argsort(vals)[:3]:
        r = minimize(
            fun,
            np.array(starts[i]),
            method="Nelder-Mead",
            options={"maxfev": 2000, "xatol": 1e-5, "fatol": 1e-3},
        )
        if best is None or r.fun < best.fun:
            best = r
    n = lc.t.size
    k = 5 + 3
    return {
        "model": "BS",
        "chi2": float(best.fun),
        "dof": n - k,
        "k": k,
        "bic": float(best.fun) + k * math.log(n),
    }


def fit_binary_lens(lc: LightCurve, ps: dict, maxfev: int = 600) -> dict:
    """Binary lens (MulensModel + VBMicrolensing), grid over s, q, alpha then Nelder-Mead."""
    mm = _mm()
    t0, te, u0 = ps["t0"], ps["tE"], ps["u0"]
    near = (t0 - 3 * te, t0 + 3 * te)

    def chi2(x):
        p = {
            "t_0": x[0],
            "u_0": x[1],
            "t_E": 10 ** x[2],
            "rho": 10 ** np.clip(x[3], -4, -0.5),
            "s": 10 ** x[4],
            "q": 10 ** np.clip(x[5], -5, 0),
            "alpha": x[6],
        }
        if not (P.te_bounds[0] < p["t_E"] < P.te_bounds[1]):
            return 1e30
        try:
            mod = mm.Model(p)
            mod.set_magnification_methods([near[0], "VBBL", near[1]])
            a = np.asarray(mod.get_magnification(lc.t), float)
        except Exception:  # noqa: BLE001 — VBBL can fail far from sensible parameters
            return 1e30
        if not np.all(np.isfinite(a)):
            return 1e30
        return linear_fluxes(a, lc.f, lc.w)[2]

    starts = [
        (t0, u0, math.log10(te), -2.0, ls, lq, al)
        for ls in (-0.3, 0.0, 0.3)
        for lq in (-3.0, -1.5, 0.0)
        for al in (30.0, 90.0, 150.0, 210.0, 270.0, 330.0)  # degrees, as MulensModel's alpha
    ]
    # near-PSPL limits (tiny or distant companion), so the fit is never worse than PSPL
    starts += [(t0, u0, math.log10(te), -2.0, ls, -4.5, 60.0) for ls in (-0.5, 0.0, 0.5)]
    starts += [(t0, u0, math.log10(te), -2.0, 0.9, lq, 60.0) for lq in (-2.0, -1.0)]
    vals = np.array([chi2(np.array(s)) for s in starts])
    best = None
    for i in np.argsort(vals)[:3]:
        r = minimize(
            chi2,
            np.array(starts[i]),
            method="Nelder-Mead",
            options={"maxfev": maxfev, "xatol": 1e-5, "fatol": 1e-3},
        )
        if best is None or r.fun < best.fun:
            best = r
    n = lc.t.size
    k = 7 + 2
    return {
        "model": "BL",
        "chi2": float(best.fun),
        "dof": n - k,
        "k": k,
        "bic": float(best.fun) + k * math.log(n),
        "x": [float(v) for v in best.x],
    }


def isolated_outliers(lc: LightCurve, model_flux: np.ndarray, nsig=4.0, nadj=3.0) -> np.ndarray:
    """Mróz et al.'s rule: > 4σ from the model while both neighbours are within 3σ."""
    z = np.abs(lc.f - model_flux) / lc.sf
    prev = np.concatenate([[0.0], z[:-1]])
    nxt = np.concatenate([z[1:], [0.0]])
    return (z > nsig) & (prev < nadj) & (nxt < nadj)


def model_flux(model: str, lc: LightCurve, r: dict) -> np.ndarray:
    p = {k: r[k] for k in ("t0", "tE", "u0", "rho", "pi_E_N", "pi_E_E", "t0_par") if k in r}
    p = {k: float(v) for k, v in p.items() if v is not None and np.isfinite(v)}
    return r["fs"] * magnification(model, lc, p) + r["fb"]


def _vet_worker(job):
    """``_vet_one`` with errors recorded: a failed test keeps the flag open (``survives``)."""
    try:
        return _vet_one(job)
    except Exception as exc:  # noqa: BLE001 — one bad flag must not lose the others' vetting
        return {
            "event_id": job[0]["event_id"],
            "tests": [("vet_error", True, repr(exc)[:200])],
            "complete": False,
        }


def _vet_one(job):
    ev, t, mag, err, do_bl = job
    lc = LightCurve.from_mag(t, mag, err, ra=ev["ra"], dec=ev["dec"], event_id=ev["event_id"])
    out = {"event_id": ev["event_id"], "tests": []}
    allm = ("PSPL", "FSPL", "PAR", *EXOTIC)
    res = fit_event(lc, ev["t0_pub"], ev["tE_pub"], ev["u0_pub"], models=allm)
    # always fit FSPL and parallax in vetting, whatever u0 and t_E
    ps = res["PSPL"]
    lt = math.log10(ps["tE"])
    if "FSPL" not in res:
        res["FSPL"] = optimise(
            "FSPL", lc, [(ps["t0"], lt, ps["u0"], math.log10(r)) for r in (0.003, 0.03, 0.3)]
        )
    if "PAR" not in res and have_mm():
        res["PAR"] = optimise("PAR", lc, par_starts(ps), t0_par=round(ps["t0"], 1))
    ex = min(EXOTIC, key=lambda m: res[m]["bic"])
    ordinary = {m: res[m] for m in ORDINARY if m in res}
    best_o = min(ordinary, key=lambda m: ordinary[m]["bic"])
    d0 = res[ex]["bic"] - ordinary[best_o]["bic"]
    out.update(exotic=ex, dbic_all=d0, best_ordinary=best_o)
    out["tests"].append(("refit_all_ordinary", d0 < P.flag_dbic, f"ΔBIC {d0:.1f} vs {best_o}"))
    # 1. isolated outliers + errors rescaled to χ²/dof = 1 of the best ordinary model
    bad = isolated_outliers(lc, model_flux(best_o, lc, ordinary[best_o])) | isolated_outliers(
        lc, model_flux(ex, lc, res[ex])
    )
    lc2 = lc.subset(~bad)
    res2 = fit_event(lc2, ps["t0"], ps["tE"], ps["u0"], models=("PSPL", "FSPL", "PAR", ex))
    if "FSPL" not in res2:
        res2["FSPL"] = optimise("FSPL", lc2, [(ps["t0"], lt, ps["u0"], math.log10(0.03))])
    best_o2 = min((m for m in ORDINARY if m in res2), key=lambda m: res2[m]["bic"])
    scale = max(res2[best_o2]["chi2"] / res2[best_o2]["dof"], 1.0)
    d1 = (res2[ex]["chi2"] - res2[best_o2]["chi2"]) / scale + (
        res2[ex]["k"] - res2[best_o2]["k"]
    ) * math.log(lc2.t.size)
    out["tests"].append(
        (
            "robust_errors_outliers",
            d1 < P.flag_dbic,
            f"{int(bad.sum())} isolated outliers removed, errors ×{math.sqrt(scale):.2f}; "
            f"ΔBIC {d1:.1f}",
        )
    )
    # 2. variable baseline: χ²/dof of a constant outside |t − t0| < 2 t_E + 60 d
    far = np.abs(lc.t - ps["t0"]) > 2 * ps["tE"] + 60.0
    if far.sum() > 20:
        fo, wo = lc.f[far], lc.w[far]
        mean = np.sum(fo * wo) / np.sum(wo)
        chi_out = float(np.sum((fo - mean) ** 2 * wo) / (far.sum() - 1))
    else:
        chi_out = np.nan
    out["chi2_out"] = chi_out
    out["tests"].append(
        (
            "variable_baseline" if np.isfinite(chi_out) else "variable_baseline_untestable",
            not (chi_out > 2.0),
            f"baseline χ²/dof {chi_out:.2f} ({far.sum()} pts)",
        )
    )
    # 2b. where the exotic preference comes from: core (|t − t0| < t_E) or wings
    ao = model_flux(best_o, lc, ordinary[best_o])
    ae = model_flux(ex, lc, res[ex])
    dchi = (lc.f - ae) ** 2 * lc.w - (lc.f - ao) ** 2 * lc.w
    core = np.abs(lc.t - ps["t0"]) < ps["tE"]
    out["dchi2_core"] = float(dchi[core].sum())
    out["dchi2_wings"] = float(dchi[~core].sum())
    # 2c. season-to-season baseline offsets (OGLE zero points, slow blends): refit both families
    lcs = lc.with_season_offsets()
    o_s = {"PSPL": optimise("PSPL", lcs, [(ps["t0"], lt, ps["u0"])])}
    if "PAR" in ordinary:
        r = ordinary["PAR"]
        o_s["PAR"] = optimise(
            "PAR",
            lcs,
            [(r["t0"], math.log10(r["tE"]), r["u0"], r["pi_E_N"], r["pi_E_E"])],
            t0_par=r["t0_par"],
        )
    rx = res[ex]
    e_s = optimise(ex, lcs, [(rx["t0"], math.log10(rx["tE"]), rx["u0"], math.log10(rx["rho"]))])
    d2 = e_s["bic"] - min(r["bic"] for r in o_s.values())
    out["tests"].append(
        (
            "season_offsets",
            d2 < P.flag_dbic,
            f"{lcs.n_extra + 1} seasons, free baseline each; ΔBIC {d2:.1f}",
        )
    )
    # 2d. slow baseline drifts: free offset and linear trend per season, both families
    lct = lc.with_season_offsets(trend=True)
    o_t = optimise("PSPL", lct, [(ps["t0"], lt, ps["u0"])])
    best_t = o_t["bic"]
    if "PAR" in ordinary:
        r = ordinary["PAR"]
        best_t = min(
            best_t,
            optimise(
                "PAR",
                lct,
                [(r["t0"], math.log10(r["tE"]), r["u0"], r["pi_E_N"], r["pi_E_E"])],
                t0_par=r["t0_par"],
            )["bic"],
        )
    e_t = optimise(ex, lct, [(rx["t0"], math.log10(rx["tE"]), rx["u0"], math.log10(rx["rho"]))])
    d2t = e_t["bic"] - best_t
    out["tests"].append(
        ("season_trends", d2t < P.flag_dbic, f"offset + linear drift per season; ΔBIC {d2t:.1f}")
    )
    # 3. binary source (xallarap proxy)
    bs = fit_binary_source(lc, ps)
    d3 = res[ex]["bic"] - min(ordinary[best_o]["bic"], bs["bic"])
    out["tests"].append(
        ("binary_source", d3 < P.flag_dbic, f"BS BIC {bs['bic']:.1f}; ΔBIC {d3:.1f}")
    )
    # 4. binary lens (expensive; only when asked)
    if do_bl and have_mm():
        bl = fit_binary_lens(lc, ps)
        d4 = res[ex]["bic"] - min(ordinary[best_o]["bic"], bs["bic"], bl["bic"])
        out["tests"].append(
            ("binary_lens", d4 < P.flag_dbic, f"BL BIC {bl['bic']:.1f}; ΔBIC {d4:.1f}")
        )
        out["bl"] = bl
    out["survives"] = all(ok for _, ok, _ in out["tests"])
    out["complete"] = bool(do_bl and have_mm())  # binary lens and PAR tested
    out["res"] = {m: {k: v for k, v in r.items() if k != "model"} for m, r in res.items()}
    return out


def arxiv_mentions(names: list[str]) -> dict[str, int]:
    """Number of arXiv records mentioning each name (export.arxiv.org API, all fields)."""
    import urllib.parse
    import urllib.request

    out = {}
    for name in names:
        if not name or name in ("X", "-"):
            continue
        q = urllib.parse.quote(f'all:"{name}"')
        url = f"https://export.arxiv.org/api/query?search_query={q}&max_results=5"
        try:
            with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310
                txt = r.read().decode()
            m = txt.split("<opensearch:totalResults")[1].split(">")[1].split("<")[0]
            out[name] = int(m)
        except Exception:  # noqa: BLE001
            out[name] = -1
        time.sleep(3.0)  # arXiv API etiquette
    return out


def variable_catalogue_matches(ra, dec, radius_arcsec: float = 1.0) -> dict[str, list]:
    """Batched CDS XMatch of all flags against VSX and Gaia DR3 variability classifications."""
    import astropy.units as u
    from astroquery.xmatch import XMatch

    pos = Table({"ra": np.asarray(ra, float), "dec": np.asarray(dec, float)})
    pos["idx"] = np.arange(len(pos))
    out = {}
    for cat in ("vizier:B/vsx/vsx", "vizier:I/358/vclassre"):
        try:
            m = XMatch.query(
                cat1=pos,
                cat2=cat,
                max_distance=radius_arcsec * u.arcsec,
                colRA1="ra",
                colDec1="dec",
            )
            out[cat] = [int(i) for i in m["idx"]]
        except Exception as exc:  # noqa: BLE001
            out[cat] = [f"error: {exc!r}"[:120]]
    return out


def run_vet(sample_key: str, procs: int, binary_lens: bool = True) -> Path:
    sample = ogle.OgleMrozSample(sample_key)
    fits = Table.read(out_dir() / f"fits_{sample_key}.ecsv")
    ok = no_error(fits)
    flags = fits[ok & (fits["dbic_min"] < P.flag_dbic)]
    print(f"{len(flags)} flags of {ok.sum()} fitted events")
    jobs = []
    ev = sample.events()
    pub = {r["event_id"]: r for r in ev}
    for r in flags:
        lc = sample.light_curve(r["event_id"])
        e = pub[r["event_id"]]
        d = {k: (e[k].item() if hasattr(e[k], "item") else e[k]) for k in ev.colnames}
        jobs.append((d, lc["time"].data, lc["mag"].data, lc["mag_err"].data, binary_lens))
    t1 = time.time()
    with Pool(procs) as pool:
        out = pool.map(_vet_worker, jobs, 1)
    names = [str(pub[o["event_id"]]["alt_id"]) for o in out] + [o["event_id"] for o in out]
    lit = arxiv_mentions(names) if out else {}
    var = variable_catalogue_matches(flags["ra"], flags["dec"]) if out else {}
    failed = [c for c, idx in var.items() if any(isinstance(j, str) for j in idx)]
    for i, o in enumerate(out):
        alt = str(pub[o["event_id"]]["alt_id"])
        n_lit = max(lit.get(alt, 0), lit.get(o["event_id"], 0))
        o["tests"].append(("literature_arxiv", True, f"{alt}: {n_lit} arXiv records (read them)"))
        hits = [c for c, idx in var.items() if i in idx]
        o["tests"].append(("variable_catalogues", not hits, f"matches: {hits or 'none'} (1″)"))
        if failed:  # a failed query is not a clean match list: the flag stays unvetted
            o["tests"].append(("variable_catalogues_failed", True, f"queries failed: {failed}"))
            o["complete"] = False
        o["survives"] = all(ok for _, ok, _ in o["tests"])
    path = out_dir() / f"vetting_{sample_key}.json"
    path.write_text(
        json.dumps(
            {
                "provenance": "derived",
                "wall_time_s": time.time() - t1,
                "variable_xmatch": var,
                "arxiv": lit,
                "flags": out,
            },
            indent=1,
            default=float,
        )
    )
    print(f"wrote {path}; survivors: {[o['event_id'] for o in out if o['survives']]}")
    return path


def contact_sheet(sample_key: str, path_png: Path) -> None:
    """Light curves of the flags with every model curve (for visual inspection)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vet = json.loads((out_dir() / f"vetting_{sample_key}.json").read_text())["flags"]
    vet = [o for o in vet if "res" in o]  # flags whose vetting raised have no fits to draw
    sample = ogle.OgleMrozSample(sample_key)
    pub = {r["event_id"]: r for r in sample.events()}
    n = len(vet)
    if not n:
        return
    ncol = 3
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.4 * nrow), squeeze=False)
    colors = {
        "PSPL": "k",
        "FSPL": "0.5",
        "PAR": "tab:blue",
        "N1neg": "tab:red",
        "E2pos": "tab:green",
        "E2neg": "tab:orange",
    }
    for ax, o in zip(axes.ravel(), vet, strict=False):
        e = pub[o["event_id"]]
        lc0 = sample.light_curve(o["event_id"])
        lc = LightCurve.from_mag(lc0["time"], lc0["mag"], lc0["mag_err"], ra=e["ra"], dec=e["dec"])
        r0 = o["res"]["PSPL"]
        win = 4 * r0["tE"] + 30
        sel = np.abs(lc.t - r0["t0"]) < win
        ax.errorbar(
            lc.t[sel] - r0["t0"],
            lc.f[sel],
            lc.sf[sel],
            fmt=".",
            ms=3,
            color="0.25",
            lw=0.5,
            zorder=5,
        )
        tt = np.linspace(r0["t0"] - win, r0["t0"] + win, 3000)
        fine = LightCurve(tt, np.ones_like(tt), np.ones_like(tt), e["ra"], e["dec"])
        for m, r in o["res"].items():
            ax.plot(
                tt - r0["t0"],
                model_flux(m, fine, r),
                color=colors.get(m, "m"),
                lw=0.9,
                label=f"{m} {r['bic'] - o['res'][o['best_ordinary']]['bic']:+.0f}",
            )
        ax.set_title(
            f"{o['event_id']} {e['alt_id']}\n{o['exotic']} ΔBIC {o['dbic_all']:.0f}; "
            f"survives={o['survives']}",
            fontsize=8,
        )
        ax.legend(fontsize=6, loc="upper right")
        ax.set_xlabel("t − t0 (d)", fontsize=7)
        ax.set_ylabel("flux (I = 18 → 1)", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path_png, dpi=110)


# ----------------------------------------------------------------------------- limit


def einstein_time_days(mass_msun, d_l_kpc=4.0, d_s_kpc=8.0, mu_mas_yr=5.0) -> np.ndarray:
    """t_E for |M| (n = 1) in the stated bulge geometry (ASSUMPTION), via ``exotic_sim``."""
    kpc = 1e3 * es.PC_M
    d_l, d_s = d_l_kpc * kpc, d_s_kpc * kpc
    th = np.array(
        [
            es.einstein_radius(es.eps_bar_point_mass(m), 1.0, d_l, d_s, d_s - d_l)
            for m in np.atleast_1d(mass_msun)
        ]
    )
    mas = th * 180.0 / math.pi * 3600e3
    return mas / mu_mas_yr * 365.25


def run_limit(per_cell_min: int = 20) -> Path:
    """95 % limits use the zero-event Poisson numerator 3.0: they need a complete null vetting."""
    path = out_dir() / "vetting_bulge2019.json"
    if not path.exists():
        raise SystemExit(f"no zero-event limit: run `vet` first ({path} missing)")
    vet = json.loads(path.read_text())["flags"]
    open_flags = [o["event_id"] for o in vet if o.get("survives") or o.get("complete") is not True]
    if open_flags:
        raise SystemExit(f"no zero-event limit: flags survive or are unvetted: {open_flags}")
    sample = ogle.OgleMrozSample("bulge2019")
    inj = Table.read(out_dir() / "injections_bulge2019.ecsv")
    inj = inj[no_error(inj)]
    rows = []
    for te in INJ_TE:
        ctrl = inj[(inj["kind"] == "PSPL") & (inj["tE"] == te)]
        p_ctrl = float(np.mean(ctrl["selected"])) if len(ctrl) else np.nan
        exposure = sample.exposure_star_years(te)
        eps_pub = sample.efficiency(te)
        for rho in INJ_RHO:
            w = inj[(inj["kind"] == "W3") & (inj["tE"] == te) & (inj["rho"] == rho)]
            sel = np.asarray(w["selected"], bool)
            flg = np.asarray(w["flag_vetted"], bool)
            p_sel = float(sel.mean())
            p_both = float((sel & flg).mean())
            p_flag = float(flg.mean())
            ratio = p_both / p_ctrl if p_ctrl > 0 else np.nan
            u0 = np.asarray(w["u0"])
            admits = {
                f"sel_u0_{a}_{b}": float(sel[(u0 >= a) & (u0 < b)].mean())
                if ((u0 >= a) & (u0 < b)).any()
                else np.nan
                for a, b in ((0.0, 1.0), (1.0, 1.8), (1.8, 2.0))
            }
            eff_w3 = eps_pub * ratio
            lim = 3.0 / (exposure * ratio) if ratio > 0 else np.inf
            mass = (te / einstein_time_days(1.0)[0]) ** 2
            rows.append(
                {
                    "tE_days": te,
                    "rho": rho,
                    "mass_msun_model": mass,
                    "n_inj": len(w),
                    "n_ctrl": len(ctrl),
                    "p_sel_pspl_ctrl": p_ctrl,
                    "p_sel_w3": p_sel,
                    "p_flag_w3": p_flag,
                    "p_sel_and_flag_w3": p_both,
                    "ratio_to_pspl": ratio,
                    "eff_pub_pspl": eps_pub,
                    "eff_w3": eff_w3,
                    "exposure_star_yr_pspl": exposure,
                    "rate95_per_star_yr": lim,
                    **admits,
                }
            )
    tab = Table(rows)
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w3_microlensing.py limit: published efficiencies × injection ratio",
        assumptions=(
            "W3 event = umbra crossing, u0 < 2 in |ε| Einstein radii; Poisson 95 % for 0 events = "
            "3.0; ε_W3 = ε_pub(t_E) × P(selected ∧ flagged | W3) / P(selected | PSPL, u0 < 1); "
            "mass: n = 1, D_L = 4 kpc, D_S = 8 kpc, μ_rel = 5 mas/yr (model_prediction)"
        ),
    )
    path = out_dir() / "limits_bulge2019.ecsv"
    tab.write(path, overwrite=True)
    tab.pprint(max_width=250, max_lines=50)
    return path


def summarise_fits(sample_key: str) -> dict:
    """ΔBIC distribution and flag counts of a fitted sample (``derived``)."""
    fits = Table.read(out_dir() / f"fits_{sample_key}.ecsv")
    ok = no_error(fits)
    out = {
        "n_events": len(fits),
        "n_fitted": int(ok.sum()),
        "wall_time_s": fits.meta.get("wall_time_s"),
    }
    f = fits[ok]
    out["best_ordinary"] = {m: int(np.sum(f["best_ordinary"] == m)) for m in ORDINARY}
    for m in (*EXOTIC, "min"):
        d = np.asarray(f[f"dbic_{m}"], float)
        out[m] = {
            "quantiles_5_25_50_75_95": [
                round(float(q), 1) for q in np.nanpercentile(d, [5, 25, 50, 75, 95])
            ],
            "n_lt_0": int(np.sum(d < 0)),
            "n_lt_minus10": int(np.sum(d < P.flag_dbic)),
            "min": round(float(np.nanmin(d)), 1),
        }
    return out


def write_manifest() -> Path:
    """Pinned OGLE files as a tracked manifest (URL, sha256, size, retrieval date)."""
    tab = Table(
        {
            "url": list(ogle.FILES),
            "sha256": [v[0] for v in ogle.FILES.values()],
            "size_bytes": [v[1] for v in ogle.FILES.values()],
            "retrieved_utc": ["2026-10-08"] * len(ogle.FILES),
        }
    )
    tab.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source="jwst_anomaly.ogle.FILES (Mróz et al. 2019, 2020 data products)",
    )
    path = paths.manifests_dir() / "ogle_mroz.ecsv"
    tab.write(path, overwrite=True)
    return path


# ----------------------------------------------------------------------------- CLI


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--procs", type=int, default=min(4, os.cpu_count() or 1))
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit", help="fit every event of a sample")
    f.add_argument("--sample", default="bulge2019", choices=sorted(ogle.SAMPLES))
    f.add_argument("--limit", type=int, default=None)
    f.add_argument("--fresh", action="store_true", help="ignore an interrupted run's checkpoint")
    f.add_argument("--chunk", default=None, help="K/N: fit only events K-1, K-1+N, ... (1-based K)")
    v = sub.add_parser("vet", help="vet the flags of a fitted sample")
    v.add_argument("--sample", default="bulge2019", choices=sorted(ogle.SAMPLES))
    v.add_argument("--no-binary-lens", action="store_true")
    c = sub.add_parser("sheet", help="contact sheet of the vetted flags")
    c.add_argument("--sample", default="bulge2019", choices=sorted(ogle.SAMPLES))
    c.add_argument("--out", type=Path, required=True)
    i = sub.add_parser("inject", help="injection-recovery on bulge light curves")
    i.add_argument("--per-cell", type=int, default=INJ_PER_CELL)
    sub.add_parser("audit", help="emulated selection on the real bulge sample")
    sub.add_parser("limit", help="95 %% rate limit from the injections")
    sub.add_parser("manifest", help="write data/manifests/ogle_mroz.ecsv")
    m = sub.add_parser("summary", help="ΔBIC distribution of a fitted sample")
    m.add_argument("--sample", default="bulge2019", choices=sorted(ogle.SAMPLES))
    args = ap.parse_args(argv)
    if args.cmd == "fit":
        run_fit(args.sample, args.limit, args.procs, args.fresh, parse_chunk(args.chunk))
    elif args.cmd == "vet":
        run_vet(args.sample, args.procs, binary_lens=not args.no_binary_lens)
    elif args.cmd == "sheet":
        contact_sheet(args.sample, args.out)
    elif args.cmd == "inject":
        run_inject(args.procs, args.per_cell)
    elif args.cmd == "audit":
        print(json.dumps(selection_audit(args.procs), indent=1))
    elif args.cmd == "limit":
        run_limit()
    elif args.cmd == "summary":
        print(json.dumps(summarise_fits(args.sample), indent=1))
    elif args.cmd == "manifest":
        print(write_manifest())
    return 0


if __name__ == "__main__":
    sys.exit(main())
