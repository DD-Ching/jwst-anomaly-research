"""Closed-form exotic-lens simulations for calibrating the exotic screens (D-047).

Every output is ``simulated`` (``schema.Provenance.SIMULATED``): synthetic lensing of a synthetic
source by a hypothetical lens. Nothing here is evidence of exotic physics.

Lens family (Kitamura, Nakajima & Asada 2013, arXiv:1211.0379; Izumi et al. 2013, arXiv:1305.5037):
the weak-field deflection angle is ``alpha(b) = eps_bar / b**n``.

- ``n = 1, sign = +1``: Schwarzschild point mass;
- ``n = 1, sign = -1``: negative point mass (Cramer et al. 1995; Safonova, Torres & Romero 2002);
- ``n = 2, sign = +1``: Ellis wormhole, ``eps_bar = pi a**2 / 4`` (Abe 2010, arXiv:1009.6084);
- ``n > 1, sign = +1``: attractive but with negative convergence ``kappa = (1 - n) / (2 x**(n+1))``.

In units of the Einstein radius (``x`` = signed image position on the lens-source axis, positive on
the source's side), the lens equation is ``beta = x - sign * sgn(x) / |x|**n``. The eigenvalues of
the lens Jacobian are ``lambda_t = 1 - sign / |x|**(n+1)`` (tangential) and
``lambda_r = 1 + sign * n / |x|**(n+1)`` (radial); the signed magnification is
``1 / (lambda_t * lambda_r)`` (Izumi et al. 2013, eqs. for lambda_+ and lambda_-).

Attractive lenses (``sign = +1``) always give two images, one on each side of the lens.
Repulsive lenses (``sign = -1``) give two images on the source's side when ``beta > beta_c``
(``caustic_beta``), none inside it (the "umbra"), and divergent magnification at ``beta_c``;
for ``n = 1``, ``beta_c = 2`` (Safonova et al. 2002).

Approximations (ASSUMPTIONs): thin lens, weak field, isolated lens (no external shear or
macro-magnification), point or uniform-disk source, geometric optics.
"""

from __future__ import annotations

import math
from functools import lru_cache

import astropy.constants as const
import numpy as np
from astropy.table import Table

from jwst_anomaly import schema

# Constants (SI) for the physical Einstein-radius helpers.
G_SI = const.G.si.value
C_SI = const.c.si.value
MSUN_KG = const.M_sun.si.value
PC_M = const.pc.si.value

_BISECT_STEPS = 64  # bracket width ≤ beta + 1, so 2^-64 of it is below double precision
_BETA_FLOOR = 1e-9  # a point source exactly behind an attractive lens would sit on the ring


def _check(n: float, sign: int) -> None:
    if not n > 0:
        raise ValueError("n must be positive")
    if sign not in (1, -1):
        raise ValueError("sign must be +1 (attractive, eps > 0) or -1 (repulsive, eps < 0)")


def deflection_integral(n: float) -> float:
    """∫_0^{π/2} cos^n ψ dψ: the factor between the metric parameter ε and eps_bar (KNA13)."""
    return 0.5 * math.sqrt(math.pi) * math.gamma((n + 1) / 2) / math.gamma((n + 2) / 2)


def eps_bar_point_mass(mass_msun: float) -> float:
    """eps_bar = 4 G |M| / c² in metres (Schwarzschild weak-field deflection, n = 1)."""
    return 4.0 * G_SI * abs(mass_msun) * MSUN_KG / C_SI**2


def eps_bar_ellis(throat_radius_m: float) -> float:
    """eps_bar = π a² / 4 in m² (Ellis wormhole leading-order deflection, Abe 2010)."""
    return math.pi * throat_radius_m**2 / 4.0


def einstein_radius(eps_bar: float, n: float, d_l: float, d_s: float, d_ls: float) -> float:
    """Angular Einstein radius in radians: (|eps_bar| D_LS / (D_S D_L^n))^(1/(n+1)).

    Kitamura et al. 2014 (arXiv:1307.6637); for ε < 0 it is the scale radius only (no ring
    forms). ``eps_bar`` in m^n and distances (angular-diameter) in m.
    """
    return (abs(eps_bar) * d_ls / (d_s * d_l**n)) ** (1.0 / (n + 1))


def caustic_beta(n: float) -> float:
    """Source radius of the caustic of a repulsive lens, (n+1) n^(-n/(n+1)) Einstein radii.

    Minimum of beta(x) = x + x^-n, reached at x_c = n^(1/(n+1)), where lam_r vanishes.
    """
    return (n + 1.0) * n ** (-n / (n + 1.0))


def critical_x(n: float) -> float:
    """Image radius of the radial critical curve of a repulsive lens, n^(1/(n+1)) Einstein radii."""
    return n ** (1.0 / (n + 1.0))


def demagnification_onset_approx(n: float) -> float:
    """KNA13's leading-order estimate 2/(n+1) of the demagnification onset.

    It comes from expanding A_tot ≈ 2/(beta (n+1)) for beta < 1 and setting it to 1, so it is only
    an approximation for large n (KNA13: 0.182 against their numerical 0.187 for n = 10). Use
    ``demagnification_onset`` for the exact value of this model.
    """
    return 2.0 / (n + 1.0)


def demagnification_onset(n: float, sign: int = 1) -> float:
    """Smallest beta (Einstein radii) at which the point-source total magnification drops below 1.

    Exact for this lens model (numerical): 1.111 for n = 2, 0.643 for n = 3, 0.187 for n = 10.
    Returns inf when the total magnification never drops below 1 (n ≤ 1, or sign = -1, whose
    total magnification is 0 in the umbra and > 1 outside it).
    """
    _check(n, sign)
    if sign == -1 or n <= 1:
        return math.inf
    # widen the bracket until the magnification drops below 1 (the onset -> inf as n -> 1+);
    # the grid starts well below the large-n estimate 2/(n+1), where A_tot > 1
    lo_b, hi_b = 0.01 * demagnification_onset_approx(n), 20.0
    while True:
        grid = np.geomspace(lo_b, hi_b, 4000)
        below = np.nonzero(total_magnification(grid, n, sign) < 1.0)[0]
        if below.size:
            break
        if hi_b > 1e12:
            return math.inf
        hi_b *= 100.0
    i = below[0]
    if i == 0:
        raise ValueError(
            f"demagnification_onset: no bracket for n = {n} (A_tot < 1 at beta {lo_b})"
        )
    lo, hi = np.array([grid[i - 1]]), np.array([grid[i]])
    return float(_bisect(lambda b: total_magnification(b, n, sign) - 1.0, lo, hi)[0])


def _bisect(f, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Vectorised bisection for a root of f inside [lo, hi] (f(lo), f(hi) of opposite sign)."""
    lo = np.array(lo, dtype=float)
    hi = np.array(hi, dtype=float)
    flo = f(lo)
    for _ in range(_BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        same = np.sign(fm) == np.sign(flo)
        lo = np.where(same, mid, lo)
        flo = np.where(same, fm, flo)
        hi = np.where(same, hi, mid)
    return 0.5 * (lo + hi)


def _eigen(x: np.ndarray, n: float, sign: int) -> tuple[np.ndarray, np.ndarray]:
    ax = np.abs(x) ** (n + 1)
    lam_t = 1.0 - sign / ax
    lam_r = 1.0 + sign * n / ax
    return lam_t, lam_r


def _roots_n1(b: np.ndarray, sign: int) -> np.ndarray:
    """Closed-form n = 1 images: x = (beta ± sqrt(beta² + 4 sign)) / 2 (outer, inner)."""
    x = np.full((b.size, 2), np.nan)
    disc = b**2 + 4.0 * sign
    ok = disc >= 0
    root = np.sqrt(np.where(ok, disc, 0.0))
    x[ok, 0] = 0.5 * (b[ok] + root[ok])
    x[ok, 1] = -sign / x[ok, 0]  # product of the roots is -sign; avoids cancellation at large beta
    return x


def _roots_numeric(b: np.ndarray, n: float, sign: int) -> np.ndarray:
    """Images for any n by bisection on the monotonic branches of the lens equation."""
    x = np.full((b.size, 2), np.nan)
    if sign == 1:
        lo = (b + 1.0) ** (-1.0 / n)
        # outer image, x > 0: x - x^-n = beta (increasing in x)
        x[:, 0] = _bisect(lambda t: t - t**-n - b, lo, b + 1.0)
        # inner image, x = -y < 0: y^-n - y = beta (decreasing in y)
        x[:, 1] = -_bisect(lambda t: t**-n - t - b, lo, np.ones_like(b))
    else:
        bc = caustic_beta(n)
        xc = critical_x(n)
        ok = b >= bc
        bo = np.where(ok, b, bc + 1.0)  # dummy values inside the umbra, masked below
        outer = _bisect(lambda t: t + t**-n - bo, np.full_like(bo, xc), bo)
        inner = _bisect(lambda t: t + t**-n - bo, bo ** (-1.0 / n), np.full_like(bo, xc))
        x[ok, 0] = outer[ok]
        x[ok, 1] = inner[ok]
    return x


def solve_images(beta, n: float = 1.0, sign: int = 1) -> dict[str, np.ndarray]:
    """Images of a point source at radius ``beta`` (Einstein radii, ≥ 0).

    Returns arrays of shape ``(len(beta), 2)``: ``x`` (signed image radius along the lens→source
    direction), ``mu`` (signed magnification), ``lam_t`` and ``lam_r`` (tangential and radial
    Jacobian eigenvalues; the image is stretched by 1/|lam| along each direction). Missing images
    (repulsive lens inside the umbra) are NaN. Column 0 is the outer image.
    """
    _check(n, sign)
    b = np.atleast_1d(np.asarray(beta, dtype=float))
    if np.any(b < 0):
        raise ValueError("beta must be non-negative")
    if sign == 1:
        b = np.maximum(b, _BETA_FLOOR)
    x = _roots_n1(b, sign) if n == 1 else _roots_numeric(b, n, sign)
    lam_t, lam_r = _eigen(x, n, sign)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = 1.0 / (lam_t * lam_r)
    return {"x": x, "mu": mu, "lam_t": lam_t, "lam_r": lam_r}


def total_magnification(beta, n: float = 1.0, sign: int = 1) -> np.ndarray:
    """Σ|μ| over a point source's images: 0 inside a repulsive lens's umbra, NaN for NaN beta."""
    b = np.atleast_1d(np.asarray(beta, dtype=float))
    total = np.nansum(np.abs(solve_images(b, n, sign)["mu"]), axis=1)
    total[~np.isfinite(b)] = np.nan  # a missing time must not look like an umbra (W3)
    return total


def image_plane_magnification(x, n: float = 1.0, sign: int = 1) -> np.ndarray:
    """Signed magnification at image radius ``x`` (Einstein radii); no root finding needed."""
    _check(n, sign)
    lam_t, lam_r = _eigen(np.asarray(x, dtype=float), n, sign)
    with np.errstate(divide="ignore"):
        return 1.0 / (lam_t * lam_r)


def convergence_shear(x, n: float = 1.0, sign: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Effective convergence and shear at image radius ``x`` (Izumi et al. 2013).

    kappa = sign (1 - n) / (2 x^(n+1)): negative for attractive n > 1 (Ellis), zero for n = 1.
    gamma = (lam_t - lam_r) / 2; it is negative (tangential stretch) for attractive lenses and
    positive (radial stretch) for repulsive ones.
    """
    lam_t, lam_r = _eigen(np.asarray(x, dtype=float), n, sign)
    return 1.0 - 0.5 * (lam_t + lam_r), 0.5 * (lam_t - lam_r)


def count_ratio(
    x, cumulative_counts, flux_limit: float, n: float = 1.0, sign: int = 1
) -> np.ndarray:
    """Image-plane number density of sources brighter than ``flux_limit``, relative to no lens.

    Magnification bias: N_obs(>S) = N(>S/|mu|) / |mu| for a locally uniform background (each image
    is a one-to-one map of a source-plane patch). ``cumulative_counts(S)`` is N(>S) per unit area,
    an ASSUMPTION the caller supplies (e.g. measured counts in the same band). Inside about one
    Einstein radius |mu| << 1, so steep bright-end counts give a deficit (a "dark hole").
    """
    mu = np.abs(image_plane_magnification(x, n, sign))
    return cumulative_counts(flux_limit / mu) / mu / cumulative_counts(flux_limit)


@lru_cache(maxsize=16)
def _gauss_legendre(n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Cached Gauss-Legendre nodes and weights (recomputing them dominated fitter run time)."""
    nodes, weights = np.polynomial.legendre.leggauss(n_nodes)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


def _check_rho(rho: float) -> None:
    if not rho >= 0:
        raise ValueError("rho must be non-negative")


def finite_source_magnification(
    beta, rho: float, n: float = 1.0, sign: int = 1, n_nodes: int = 48, point_magnification=None
):
    """Total magnification of a uniform disk of radius ``rho`` (Einstein radii) centred at ``beta``.

    The lens is axisymmetric, so the disk average reduces to a 1-D integral over the source radius
    b of A(b) times the length of the circle of radius b inside the disk, divided by pi rho². The
    b range is split at the disk's inner radii and at the caustic (umbra edge) of a repulsive lens;
    each piece uses Gauss-Legendre in t with b = lo + (hi - lo)(1 - cos t)/2, which absorbs the
    1/sqrt singularities at the caustic and the disk edges. Agrees with inverse ray shooting
    (``tests/test_exotic_sim.py``). ``point_magnification(b)`` replaces
    ``total_magnification(b, n, sign)`` inside the integral (e.g. a faster tabulated equivalent).
    """
    _check_rho(rho)
    if rho == 0:
        return total_magnification(beta, n, sign)
    b = np.atleast_1d(np.asarray(beta, dtype=float))
    lo = np.maximum(b - rho, 0.0)
    hi = b + rho
    bc = caustic_beta(n) if sign == -1 else 0.0
    cuts = np.stack([lo, np.abs(b - rho), np.full_like(b, bc), hi], axis=1)
    cuts = np.sort(np.clip(cuts, lo[:, None], hi[:, None]), axis=1)
    a_, b_ = cuts[:, :-1], cuts[:, 1:]  # (N, 3) sub-intervals
    nodes, weights = _gauss_legendre(int(n_nodes))
    t = 0.5 * np.pi * (nodes + 1.0)
    u = 0.5 * (1.0 - np.cos(t))
    du = 0.25 * np.pi * np.sin(t) * weights  # d u = sin(t)/2 dt, dt = pi/2 d(node)
    # Most sub-intervals have zero width (b ≥ ρ, or the caustic outside the disk); their
    # integrand is exactly +0 (w = 0, amp and arc finite), so it is evaluated only on the others
    # and scattered into a zero array of the full (N, 3, n_nodes) shape: the sum below then runs
    # over the same layout and is bit-identical, at ~1/3 to 1/2 of the cost.
    integrand = np.zeros(a_.shape + (nodes.size,))
    nz = b_ > a_
    if nz.any():
        a0 = a_[nz][:, None]
        width = (b_ - a_)[nz][:, None]
        bb = np.broadcast_to(b[:, None], a_.shape)[nz][:, None]
        rr = a0 + width * u  # (M, n_nodes)
        w = width * du
        full = rr <= rho - bb  # whole circle inside the disk (source disk covers the lens)
        with np.errstate(divide="ignore", invalid="ignore"):
            cosphi = (rr**2 + bb**2 - rho**2) / (2.0 * rr * bb)
        # = clip(nan_to_num(cosphi, nan=1), -1, 1): ±inf clip to ±1 either way
        cosphi = np.clip(np.where(np.isnan(cosphi), 1.0, cosphi), -1.0, 1.0)
        phi = np.where(full, np.pi, np.arccos(cosphi))
        arc = 2.0 * rr * phi
        amp_fn = point_magnification or (lambda x: total_magnification(x, n, sign))
        amp = np.asarray(amp_fn(np.maximum(rr, 0.0).ravel()), float).reshape(rr.shape)
        # a node rounding onto the caustic (measure zero) would give inf * 0
        amp = np.where(np.isfinite(amp), amp, 0.0)
        integrand[nz] = amp * arc * w
    out = np.sum(integrand, axis=(1, 2)) / (np.pi * rho**2)
    out[~np.isfinite(b)] = np.nan
    return out


def impact_track(t, t0: float, t_e: float, u0: float) -> np.ndarray:
    """beta(t) = sqrt(u0² + ((t - t0)/t_E)²) for straight-line relative motion (Abe 2010)."""
    if not t_e > 0:
        raise ValueError("t_e must be positive")
    t = np.asarray(t, dtype=float)
    return np.sqrt(u0**2 + ((t - t0) / t_e) ** 2)


def light_curve(
    t, t0: float, t_e: float, u0: float, n: float = 1.0, sign: int = 1, rho: float = 0.0
) -> np.ndarray:
    """Total magnification A(t); for sign = -1, zero in the umbra and spikes at the caustic."""
    _check_rho(rho)
    beta = impact_track(t, t0, t_e, u0)
    return finite_source_magnification(beta, rho, n, sign)


def inject_images(
    source_dx,
    source_dy,
    theta_e: float,
    n: float = 1.0,
    sign: int = 1,
    flux=1.0,
    source_id=None,
) -> Table:
    """Image-plane injection: one row per lensed image of each point source (``simulated``).

    ``source_dx``, ``source_dy``: unlensed source offsets from the lens (any angular unit, e.g.
    arcsec, in the frame the screen uses); ``theta_e`` in the same unit. Output columns:
    ``source_id``, ``image`` (0 outer, 1 inner), ``dx``, ``dy`` (image offsets), ``mu`` (signed),
    ``flux`` (``flux * |mu|``), ``stretch_t`` and ``stretch_r`` (1/|lam_t|, 1/|lam_r|), ``radial``
    (True when the image is longer radially than tangentially) and ``pa_deg`` (position angle of
    the image's long axis, degrees east of north, assuming dx = east, dy = north). Sources inside a
    repulsive lens's umbra produce no rows.
    """
    try:
        sx, sy = np.broadcast_arrays(
            np.atleast_1d(np.asarray(source_dx, dtype=float)),
            np.atleast_1d(np.asarray(source_dy, dtype=float)),
        )
    except ValueError as err:
        raise ValueError("source_dx and source_dy must broadcast to one length") from err
    fl = np.broadcast_to(np.asarray(flux, dtype=float), sx.shape)
    ids = np.arange(sx.size) if source_id is None else np.atleast_1d(np.asarray(source_id))
    if ids.shape != sx.shape:
        raise ValueError("source_id must have one entry per source")
    if not theta_e > 0:
        raise ValueError("theta_e must be positive")
    beta_ang = np.hypot(sx, sy)
    beta = beta_ang / theta_e
    sol = solve_images(beta, n, sign)
    # unit vector lens → source; arbitrary (east) for a source exactly behind the lens
    safe = np.where(beta_ang > 0, beta_ang, 1.0)
    ux = np.where(beta_ang > 0, sx / safe, 1.0)
    uy = np.where(beta_ang > 0, sy / safe, 0.0)
    rows = {
        k: [] for k in ("source_id", "image", "dx", "dy", "mu", "flux", "stretch_t", "stretch_r")
    }
    for k in range(2):
        xk = sol["x"][:, k]
        good = np.isfinite(xk)
        rows["source_id"].append(ids[good])
        rows["image"].append(np.full(good.sum(), k))
        rows["dx"].append(xk[good] * theta_e * ux[good])
        rows["dy"].append(xk[good] * theta_e * uy[good])
        rows["mu"].append(sol["mu"][good, k])
        rows["flux"].append(fl[good] * np.abs(sol["mu"][good, k]))
        rows["stretch_t"].append(1.0 / np.abs(sol["lam_t"][good, k]))
        rows["stretch_r"].append(1.0 / np.abs(sol["lam_r"][good, k]))
    out = Table({k: np.concatenate(v) for k, v in rows.items()})
    out["radial"] = out["stretch_r"] > out["stretch_t"]
    radial_pa = np.degrees(np.arctan2(out["dx"], out["dy"])) % 180.0
    out["pa_deg"] = np.where(out["radial"], radial_pa, (radial_pa + 90.0) % 180.0)
    out.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source=f"exotic_sim.inject_images(n={n}, sign={sign}, theta_e={theta_e})",
    )
    return out


def inject_light_curve(
    times,
    baseline_flux: float,
    t0: float,
    t_e: float,
    u0: float,
    n: float = 1.0,
    sign: int = 1,
    rho: float = 0.0,
    blend: float = 1.0,
) -> Table:
    """Light-curve injection for the transient screen (``simulated``).

    ``flux = baseline_flux * (blend * A(t) + 1 - blend)``, where ``blend`` is the lensed fraction of
    the baseline flux (1 = unblended). Columns: ``time``, ``beta``, ``magnification``, ``flux``.
    """
    if not 0 <= blend <= 1:
        raise ValueError("blend must be in [0, 1]")
    t = np.atleast_1d(np.asarray(times, dtype=float))
    a = light_curve(t, t0, t_e, u0, n, sign, rho)
    out = Table(
        {
            "time": t,
            "beta": impact_track(t, t0, t_e, u0),
            "magnification": a,
            "flux": baseline_flux * (blend * a + 1.0 - blend),
        }
    )
    out.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source=(
            f"exotic_sim.inject_light_curve(n={n}, sign={sign}, t0={t0}, t_e={t_e}, u0={u0}, "
            f"rho={rho}, blend={blend})"
        ),
    )
    return out
