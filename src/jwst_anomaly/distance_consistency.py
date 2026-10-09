"""D1 distance self-consistency (hypothesis round 1, owner idea 1).

Per-object tests of whether two independent distance measures along one sightline agree with an
ordinary FLRW model, using published posteriors as the data (derived) and astropy.cosmology for the
predictions (model_prediction).

Strong lenses with stellar kinematics: the time-delay distance D_dt and the lens distance D_d
give the H0-free ratio R = D_dt / ((1 + z_d) D_d) = D_s / D_ds. In a flat universe
R = chi_s / (chi_s - chi_d), so it depends only weakly on Omega_m and w. A lens whose R cannot be
reconciled with its redshifts is an
anomaly, not evidence of anything exotic: the mass-sheet degeneracy (lambda_MST), line-of-sight
convergence, kinematic anisotropy, substructure and wrong redshifts are the ordinary explanations.

FRBs: dispersion measure vs host redshift (Macquart relation). The per-burst residual is judged
against
the predictive distribution of DM_MW + DM_cosmic + DM_host / (1 + z), not a Gaussian.

All functions are pure (numpy arrays in, arrays/tables out) so they can be tested offline.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM, Flatw0waCDM
from scipy import stats

C_KMS = const.c.to(u.km / u.s).value


@dataclass(frozen=True)
class Params:
    om_range: tuple[float, float] = (0.1, 0.5)  # ASSUMPTION: broad flat prior on Omega_m
    w_range: tuple[float, float] = (-2.0, -0.5)  # ASSUMPTION: broad flat prior on w (wCDM variant)
    h0_range: tuple[float, float] = (50.0, 100.0)  # ASSUMPTION: flat prior in ln H0 for D_dt tests
    lns_half_width: float = (
        0.7  # ASSUMPTION: flat prior on a common ln scale (lambda_MST-like) in R tests
    )
    n_om: int = 41
    n_w: int = 31
    n_lns: int = 281
    n_pred: int = 200_000  # predictive draws per lens
    n_obs: int = 200_000  # posterior samples kept per lens (random subsample)
    density_bins: int = 300  # histogram density of ln X used as a 1-D likelihood
    flag_sigma_global: float = 5.0  # ASSUMPTION: trials-corrected flag threshold


DEFAULT = Params()


# --------------------------------------------------------------------------------------------------
# Distances (astropy.cosmology; H0 cancels in every dimensionless quantity below)
# --------------------------------------------------------------------------------------------------


_ZGRID = np.linspace(0.0, 3.0, 3001)


@lru_cache(maxsize=4096)
def _chi_table(om: float, w: float) -> np.ndarray:
    """chi(z) / (c/H0) on _ZGRID from astropy's E(z) (cumulative trapezoid; ~1e-7 relative)."""
    cosmo = Flatw0waCDM(H0=100.0, Om0=om, w0=w, wa=0.0, Tcmb0=0.0)
    inv_e = 1.0 / cosmo.efunc(_ZGRID)
    return np.concatenate([[0.0], np.cumsum(0.5 * (inv_e[1:] + inv_e[:-1]) * np.diff(_ZGRID))])


def comoving_dimensionless(z: Sequence[float], om: float, w: float = -1.0) -> np.ndarray:
    """Comoving distance / (c/H0), flat (w)CDM, radiation neglected, z <= 3."""
    z = np.asarray(z, dtype=float)
    if np.any(z > _ZGRID[-1]):
        raise ValueError("z > 3 not tabulated")
    return np.interp(z, _ZGRID, _chi_table(round(float(om), 10), round(float(w), 10)))


def ratio_model(zd, zs, om: float, w: float = -1.0) -> np.ndarray:
    """R = D_s / D_ds = chi_s / (chi_s - chi_d) (flat)."""
    zd, zs = np.atleast_1d(zd), np.atleast_1d(zs)
    chi_d = comoving_dimensionless(zd, om, w)
    chi_s = comoving_dimensionless(zs, om, w)
    return chi_s / (chi_s - chi_d)


def ddt_model_dimensionless(zd, zs, om: float, w: float = -1.0) -> np.ndarray:
    """D_dt in units of c/H0: chi_d chi_s / (chi_s - chi_d) (flat)."""
    zd, zs = np.atleast_1d(zd), np.atleast_1d(zs)
    chi_d = comoving_dimensionless(zd, om, w)
    chi_s = comoving_dimensionless(zs, om, w)
    return chi_d * chi_s / (chi_s - chi_d)


def ratio_from_samples(dd: np.ndarray, ddt: np.ndarray, zd: float) -> np.ndarray:
    """Per-sample R = D_dt / ((1 + z_d) D_d), keeping the joint (D_d, D_dt) correlation."""
    return np.asarray(ddt, float) / ((1.0 + zd) * np.asarray(dd, float))


def grid_ln_model(zd, zs, kind: str, params: Params = DEFAULT) -> tuple[np.ndarray, np.ndarray]:
    """ln model on an Omega_m grid for each lens: shape (n_lens, n_om). kind: 'ratio' or 'ddt'."""
    om = np.linspace(*params.om_range, params.n_om)
    fn = ratio_model if kind == "ratio" else ddt_model_dimensionless
    grid = np.array([np.log(fn(zd, zs, o)) for o in om]).T
    return om, grid


# --------------------------------------------------------------------------------------------------
# Pulls
# --------------------------------------------------------------------------------------------------


def p_to_sigma(p):
    """Two-sided p-value -> Gaussian sigma."""
    return stats.norm.isf(np.clip(np.asarray(p, float), 1e-300, 1.0) / 2.0)


def sigma_to_p(sigma):
    return 2.0 * stats.norm.sf(np.asarray(sigma, float))


def pull(ln_obs: np.ndarray, ln_pred: np.ndarray, rng: np.random.Generator, w_obs=None) -> dict:
    """Compare observed and predicted ln-quantity samples via Delta = obs - pred (random pairs).

    Returns the Gaussian pull z = mean / std of Delta, and the empirical two-sided tail p (floored
    at
    1 / N, so it saturates near 4.4 sigma for 2e5 draws; z is the number to use beyond that)."""
    n = len(ln_pred)
    prob = None if w_obs is None else np.asarray(w_obs, float) / np.sum(w_obs)
    obs = rng.choice(np.asarray(ln_obs, float), size=n, replace=True, p=prob)
    delta = obs - rng.permutation(np.asarray(ln_pred, float))
    mu, sd = float(np.mean(delta)), float(np.std(delta))
    frac_pos = float(np.mean(delta > 0))
    p_emp = max(2.0 * min(frac_pos, 1.0 - frac_pos), 1.0 / n)
    return {
        "delta_ln": mu,
        "sigma_ln": sd,
        "z": mu / sd,
        "p_emp": p_emp,
        "z_emp": float(np.sign(mu) * p_to_sigma(p_emp)),
    }


def prior_predictive_ratio(
    zd: float, zs: float, rng: np.random.Generator, params: Params = DEFAULT, wcdm: bool = False
) -> np.ndarray:
    """ln R draws for one lens under the broad Omega_m (and optionally w) prior."""
    om = np.linspace(*params.om_range, params.n_om)
    if not wcdm:
        ln_r = np.array([np.log(ratio_model(zd, zs, o))[0] for o in om])
        draws = rng.uniform(*params.om_range, params.n_pred)
        return np.interp(draws, om, ln_r)
    w = np.linspace(*params.w_range, params.n_w)
    ln_r = np.array([[np.log(ratio_model(zd, zs, o, ww))[0] for ww in w] for o in om])
    from scipy.interpolate import RegularGridInterpolator

    f = RegularGridInterpolator((om, w), ln_r)
    pts = np.column_stack(
        [rng.uniform(*params.om_range, params.n_pred), rng.uniform(*params.w_range, params.n_pred)]
    )
    return f(pts)


def ln_density(samples: np.ndarray, weights=None, bins: int = 300):
    """Histogram density of a 1-D sample as a callable log-likelihood (floored, interpolated)."""
    samples = np.asarray(samples, float)
    lo, hi = np.percentile(samples, [0.01, 99.99])
    pad = 0.5 * (hi - lo)
    hist, edges = np.histogram(
        samples, bins=bins, range=(lo - pad, hi + pad), weights=weights, density=True
    )
    hist = np.convolve(hist, np.array([0.25, 0.5, 0.25]), mode="same")
    centres = 0.5 * (edges[1:] + edges[:-1])
    floor = 1e-12 * hist.max()

    def f(x):
        return np.log(np.interp(x, centres, hist, left=0.0, right=0.0) + floor)

    return f


def leave_one_out(
    ln_obs: Sequence[np.ndarray],
    ln_model_grid: np.ndarray,
    om_grid: np.ndarray,
    lns_grid: np.ndarray | None,
    rng: np.random.Generator,
    params: Params = DEFAULT,
    weights: Sequence | None = None,
) -> list[dict]:
    """Leave-one-out predictive pull of each object against the fit to all the others.

    Model: ln X_i = ln g_i(Omega_m) + ln s, with s a common scale (c/H0 for D_dt; a lambda_MST-like
    sample-wide factor for R, or fixed at 1 when lns_grid is None). Flat priors on the grids."""
    n = len(ln_obs)
    weights = weights or [None] * n
    dens = [ln_density(x, w, params.density_bins) for x, w in zip(ln_obs, weights, strict=True)]
    lns = np.zeros(1) if lns_grid is None else lns_grid
    # loglike[i, a, b] = ln L_i(om_a, lns_b)
    loglike = np.array(
        [[dens[i](ln_model_grid[i, a] + lns) for a in range(len(om_grid))] for i in range(n)]
    )
    out = []
    for i in range(n):
        post = loglike[np.arange(n) != i].sum(axis=0)
        post = np.exp(post - post.max()).ravel()
        post /= post.sum()
        idx = rng.choice(post.size, size=params.n_pred, p=post)
        a, b = np.unravel_index(idx, (len(om_grid), len(lns)))
        # jitter within grid cells so the predictive is continuous
        dom = om_grid[1] - om_grid[0]
        om_draw = np.clip(
            om_grid[a] + rng.uniform(-0.5, 0.5, idx.size) * dom, om_grid[0], om_grid[-1]
        )
        lns_draw = lns[b] + (
            rng.uniform(-0.5, 0.5, idx.size) * (lns[1] - lns[0]) if lns.size > 1 else 0.0
        )
        pred = np.interp(om_draw, om_grid, ln_model_grid[i]) + lns_draw
        res = pull(ln_obs[i], pred, rng, weights[i])
        res["lns_post_mean"] = float(np.mean(lns_draw))
        out.append(res)
    return out


def local_sigma_threshold(
    n_trials: int, global_sigma: float = 5.0, one_sided: bool = False
) -> float:
    """Local sigma that corresponds to a global two-sided `global_sigma` after n_trials (Sidak).

    one_sided=True converts the local p to a one-sided sigma (for one-tailed tests)."""
    p_glob = sigma_to_p(global_sigma)
    p_loc = -np.expm1(np.log1p(-p_glob) / n_trials)
    if one_sided:
        return float(stats.norm.isf(p_loc))
    return float(p_to_sigma(p_loc))


def global_sigma(local_sigma: float, n_trials: int) -> float:
    p_loc = sigma_to_p(abs(local_sigma))
    p_glob = -np.expm1(n_trials * np.log1p(-p_loc))
    return float(p_to_sigma(p_glob))


# --------------------------------------------------------------------------------------------------
# FRB DM-z (Macquart relation)
# --------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class FRBParams:
    f_d: float = 0.844  # ASSUMPTION: diffuse-baryon fraction (Macquart et al. 2020 fiducial)
    feedback_F: float = 0.32  # ASSUMPTION: sigma_DM = F z^-1/2 (Macquart+2020 form; James+2022 F)
    alpha: float = 3.0  # Macquart+2020 p(Delta) shape (alpha = beta = 3)
    host_mu: float = float(
        np.log(68.2)
    )  # ASSUMPTION: log-normal DM_host, median 68.2 (Macquart+2020)
    host_sigma: float = 0.88  # ASSUMPTION: log-normal DM_host width (Macquart+2020)
    halo_range: tuple[float, float] = (10.0, 80.0)  # ASSUMPTION: Milky Way halo DM, uniform
    ism_frac_err: float = 0.2  # ASSUMPTION: Gaussian 20 % (1 sigma) error on the NE2001 DM_ISM
    z_sigma_floor: float = 0.02  # ASSUMPTION: sigma = F max(z, floor)^-1/2 (z^-1/2 diverges at 0)
    delta_max: float = 20.0  # ASSUMPTION: p(Delta) truncated at 20 <DM_cosmic> (Delta^-3 tail)
    n_grid: int = (
        20_001  # DM grid points of the convolution (finer if the cosmic or ISM term needs it)
    )
    flag_sigma_global: float = 5.0  # ASSUMPTION: trials-corrected flag threshold


FRB_DEFAULT = FRBParams()


def macquart_mean_dm(z, cosmo=None, f_d: float = 0.844, n: int = 400) -> np.ndarray:
    """<DM_cosmic>(z) in pc cm^-3: (3 c H0 Ob f_d / 8 pi G m_p) * int (1+z) chi_e / E(z) dz,
    chi_e = Y_H + Y_He / 2 with Y_He = 0.25 (He fully ionised)."""
    cosmo = cosmo or FlatLambdaCDM(H0=67.66, Om0=0.30966, Ob0=0.04897, Tcmb0=0.0)  # Planck18 values
    chi_e = 0.75 + 0.25 / 2.0
    pref = (
        (3 * const.c * cosmo.H0 * cosmo.Ob0 * f_d / (8 * np.pi * const.G * const.m_p))
        .to(u.pc / u.cm**3)
        .value
    )
    z = np.atleast_1d(np.asarray(z, float))
    out = np.empty_like(z)
    for k, zz in enumerate(z):
        grid = np.linspace(0.0, zz, n)
        integrand = (1 + grid) * chi_e / cosmo.efunc(grid)
        out[k] = pref * np.trapezoid(integrand, grid)
    return out


def _macquart_unnorm(x, c0: float, sigma: float, alpha: float = 3.0):
    """Macquart+2020 p(Delta) ~ Delta^-alpha exp(-(Delta^-alpha - C0)^2 / (2 alpha^2 sigma^2))."""
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    ok = x > 0
    xa = x[ok] ** (-alpha)
    out[ok] = xa * np.exp(-((xa - c0) ** 2) / (2 * alpha**2 * sigma**2))
    return out


def _macquart_c0(sigma: float, alpha: float = 3.0, delta_max: float = 20.0) -> float:
    """C0 such that <Delta> = 1 for the Macquart+2020 p(Delta) truncated at delta_max."""
    from scipy.optimize import brentq

    x = np.linspace(1e-3, delta_max, 40000)

    def mean_minus_one(c0):
        p = _macquart_unnorm(x, c0, sigma, alpha)
        return np.trapezoid(x * p, x) / np.trapezoid(p, x) - 1.0

    return brentq(mean_minus_one, -30, 30)


def macquart_sigma(z: float, p: FRBParams = FRB_DEFAULT) -> float:
    return p.feedback_F / np.sqrt(max(z, p.z_sigma_floor))


def draw_macquart_delta(
    z: float, rng: np.random.Generator, p: FRBParams = FRB_DEFAULT, n: int = 200_000
) -> np.ndarray:
    """Monte Carlo draws of Delta = DM_cosmic / <DM_cosmic> (used only to cross-check the grid)."""
    sigma = macquart_sigma(z, p)
    c0 = _macquart_c0(sigma, p.alpha, p.delta_max)
    x = np.linspace(1e-3, p.delta_max, 40000)
    cdf = np.cumsum(_macquart_unnorm(x, c0, sigma, p.alpha))
    cdf /= cdf[-1]
    return np.interp(rng.uniform(size=n), cdf, x)


def _uniform_gauss_pdf(y, mean: float, s: float, lo: float, hi: float) -> np.ndarray:
    """PDF of N(mean, s) + U(lo, hi), written to avoid cancellation in both tails."""
    u1 = (np.asarray(y, float) - lo - mean) / s
    u2 = (np.asarray(y, float) - hi - mean) / s
    val = np.where(
        u2 > 0, stats.norm.sf(u2) - stats.norm.sf(u1), stats.norm.cdf(u1) - stats.norm.cdf(u2)
    )
    return np.clip(val, 0.0, None) / (hi - lo)


class FRBPredictive:
    """Predictive distribution of DM_obs = DM_ISM + DM_halo + DM_cosmic + DM_host/(1+z) for one FRB.

    ISM (Gaussian) + halo (uniform) is analytic; the cosmic term is convolved on a uniform DM grid
    (direct sum of positive terms, no FFT round-off); the host log-normal enters via its analytic
    CDF/SF. Tail probabilities are sums of positive terms, so they stay accurate far below 1e-10."""

    def __init__(self, dm_ism: float, z: float, p: FRBParams = FRB_DEFAULT, mean_dm=None):
        self.p, self.z = p, z
        self.mean_dm = float(macquart_mean_dm(z, f_d=p.f_d)[0] if mean_dm is None else mean_dm)
        s = max(p.ism_frac_err * dm_ism, 1e-3)
        lo_h, hi_h = p.halo_range
        y_lo, y_hi = dm_ism + lo_h - 14 * s, dm_ism + hi_h + 14 * s
        span = (y_hi - y_lo) + p.delta_max * self.mean_dm
        dx = min(span / (p.n_grid - 1), self.mean_dm / 50.0, s / 5.0)
        self.dx = dx
        y = y_lo + dx * np.arange(int(np.ceil((y_hi - y_lo) / dx)) + 1)
        p_ih = _uniform_gauss_pdf(y, dm_ism, s, lo_h, hi_h)
        sigma = macquart_sigma(z, p)
        c0 = _macquart_c0(sigma, p.alpha, p.delta_max)
        c = dx * np.arange(int(np.ceil(p.delta_max * self.mean_dm / dx)) + 1)
        p_c = _macquart_unnorm(c / self.mean_dm, c0, sigma, p.alpha)
        p_c /= p_c.sum()
        self.p_ih_w = p_ih * dx / np.sum(p_ih * dx)  # weights of ISM + halo on y
        self.y = y
        w = np.convolve(self.p_ih_w, p_c)
        self.a = y_lo + dx * np.arange(w.size)  # grid of ISM + halo + cosmic
        self.w = w / w.sum()

    def _host_cdf_sf(self, t):
        t = np.asarray(t, float)
        pos = t > 0
        arg = np.full(t.shape, -np.inf)
        arg[pos] = (np.log((1 + self.z) * t[pos]) - self.p.host_mu) / self.p.host_sigma
        return stats.norm.cdf(arg), stats.norm.sf(arg)

    def tails(self, dm_obs: float) -> tuple[float, float]:
        """(P(pred <= obs), P(pred >= obs))."""
        cdf, sf = self._host_cdf_sf(dm_obs - self.a)
        return float(np.sum(self.w * cdf)), float(np.sum(self.w * sf))

    def p_below_floor(self, dm_obs: float) -> float:
        """P(DM_ISM + DM_halo >= DM_obs): the observed DM is below the Milky-Way-only floor."""
        return float(np.sum(self.p_ih_w[self.y >= dm_obs]))

    def quantile_dm(self, q: float) -> float:
        from scipy.optimize import brentq

        hi = self.a[-1] + 1e5
        return float(brentq(lambda d: self.tails(d)[0] - q, self.a[0] - 1.0, hi, xtol=1e-3))

    def detect_limits(self, sigma_one_sided: float) -> tuple[float, float]:
        """DM at which the low (high) one-sided tail reaches `sigma_one_sided`; nan if none."""
        from scipy.optimize import brentq

        p_thr = stats.norm.sf(sigma_one_sided)
        med = self.quantile_dm(0.5)
        lo_end = min(self.a[0] - 1.0, 0.0) - 1e3
        low = brentq(
            lambda d: np.log(max(self.tails(d)[0], 1e-300)) - np.log(p_thr), lo_end, med, xtol=1e-3
        )
        hi_end = self.a[-1] + 1e7
        if self.tails(hi_end)[1] > p_thr:
            high = np.nan
        else:
            high = brentq(
                lambda d: np.log(max(self.tails(d)[1], 1e-300)) - np.log(p_thr),
                med,
                hi_end,
                xtol=1e-3,
            )
        return float(low), float(high)


def frb_predictive_tails(
    dm_obs: float, dm_ism: float, z: float, p: FRBParams = FRB_DEFAULT, mean_dm=None
) -> dict:
    """One-sided predictive tails of the observed DM given z (deterministic grid convolution)."""
    pred = FRBPredictive(dm_ism, z, p, mean_dm)
    p_low, p_high = pred.tails(dm_obs)
    return {
        "dm_pred_median": pred.quantile_dm(0.5),
        "p_low": p_low,
        "p_high": p_high,
        "z_low": float(stats.norm.isf(min(max(p_low, 1e-300), 0.5))),  # obs below prediction
        "z_high": float(stats.norm.isf(min(max(p_high, 1e-300), 0.5))),  # obs above prediction
        "p_below_floor": pred.p_below_floor(dm_obs),
        "mean_dm_cosmic": pred.mean_dm,
    }
