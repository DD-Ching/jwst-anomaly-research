"""Gap G2 of goal round 6: do achronal null geodesics of the Krasnikov tube violate ANEC?

Goal round 7 (docs/goal_round7/README.md). Offline, symbolic plus small numerical checks; no data.
Geometry: the complete Everett & Roman (gr-qc/9702049, "ER") tube, Eqs. (12)-(13), including the
front theta(t - x - rho) laid at light speed and both end caps:

    ds^2 = -dt^2 + (1 - k) dx dt + k dx^2 + drho^2 + rho^2 dphi^2
         = -(dt - dx)(dt + k dx) + drho^2 + rho^2 dphi^2,
    k = 1 - eta theta(rho_max - rho) theta(t - x - rho) [theta(x) - theta(x + eps - D)].

Two step functions are used: a compact C-infinity step (exactly 0 for s <= 0 and 1 for s >= eps,
ER Eq. 7) and ER's own tanh step (Eq. 35), which has tails.

Checked here [derived]:
1. Outgoing lines u = t - x = const, rho, phi = const: null pregeodesics, R_ll = 0 and zero screen
   tidal tensor (shear-free) for ANY k(u, x, rho); achronal for eta < 1 because u never decreases
   along causal curves (u is a null coordinate, not a time function).
2. Where k depends on x alone the spacetime is flat. With the compact step this is the core
   {rho < rho_max - eps, t - x - rho > eps} plus the exterior.
3. Axis regularity. ER's radial factor theta(R - rho) with the tanh step has k_rho(0) != 0, and
   the front factor theta(t - x - rho) has d/d rho != 0 on the axis while the front crosses it
   (0 < t - x < eps). Both are kinks in Cartesian coordinates: the metric is not C^1 there. The
   axis results below therefore need a radial factor smooth at the axis, k = K(x, rho^2) with K
   at least C^2 (C^4 for the curvature to be C^2), and exclude the front-crossing set.
4. Axial incoming ray, eternal tube, k = K(x, rho^2) with the axis minimising k(x, .) for every x:
   w = t + int k(x, 0) dx never decreases along causal curves, so the axial ray is a complete
   achronal null geodesic. With the affine tangent of unit exterior Killing energy,
   l = (2 k0 / (1 + k0), -2 / (1 + k0), 0, 0), R_ab l^a l^b = -4 k2 / (1 + k0) and
   d lambda = -(1 + k0) dx / 2, so ANEC = -2 int k2(x) dx (k2 = d^2 k / d rho^2 on the axis >= 0).
   Checked on two smooth radial factors: exp(-(rho/R)^2) (k2 > 0, ANEC < 0) and exp(-(rho/R)^4)
   (k2 = 0, ANEC = 0, curved off the axis). These are mathematical checks, not source models.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache

import numpy as np
import sympy as sp


@dataclass(frozen=True)
class Params:
    """Numerical choices for the profile check (ASSUMPTIONs; units of eps)."""

    eta: float = 0.999
    eps: float = 1.0
    rho_max: float = 50.0
    d_len: float = 400.0
    n_rays: int = 64
    seed: int = 7


P = Params()


def _christoffel(g, xs):
    gi = sp.simplify(g.inv())
    n = len(xs)
    return [
        [
            [
                sp.simplify(
                    sum(
                        gi[a, d]
                        * (
                            sp.diff(g[d, b], xs[c])
                            + sp.diff(g[d, c], xs[b])
                            - sp.diff(g[b, c], xs[d])
                        )
                        for d in range(n)
                    )
                    / 2
                )
                for c in range(n)
            ]
            for b in range(n)
        ]
        for a in range(n)
    ]


def _ricci(gam, xs, b, c):
    n = len(xs)
    return sp.simplify(
        sum(
            sp.diff(gam[a][b][c], xs[a])
            - sp.diff(gam[a][b][a], xs[c])
            + sum(gam[a][a][d] * gam[d][b][c] - gam[a][c][d] * gam[d][b][a] for d in range(n))
            for a in range(n)
        )
    )


def outgoing_lines() -> dict:
    """Coordinates (u = t - x, x, rho, phi): ds^2 = -du^2 - (1 + k) du dx + drho^2 + rho^2 dphi^2
    with k = k(u, x, rho) arbitrary. The tangent d_x of the lines u, rho, phi = const is null."""
    u, x, rho, phi = sp.symbols("u x rho phi", real=True)
    k = sp.Function("k")(u, x, rho)
    xs = [u, x, rho, phi]
    g = sp.Matrix(
        [[-1, -(1 + k) / 2, 0, 0], [-(1 + k) / 2, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]]
    )
    gam = _christoffel(g, xs)
    acc = [gam[a][1][1] for a in range(4)]  # Gamma^a_xx
    gi = sp.simplify(g.inv())

    def riem_low(a_, b_, c_, d_):  # R_{a b c d} = g_{a m} R^m_{b c d}
        tot = 0
        for m in range(4):
            e = sp.diff(gam[m][b_][d_], xs[c_]) - sp.diff(gam[m][b_][c_], xs[d_])
            e += sum(
                gam[m][c_][n] * gam[n][b_][d_] - gam[m][d_][n] * gam[n][b_][c_] for n in range(4)
            )
            tot += g[a_, m] * e
        return sp.simplify(tot)

    # screen tidal tensor R_{A x B x} for A, B in (rho, phi): its trace-free part sources shear
    screen = {f"{a_}{b_}": str(riem_low(a_, 1, b_, 1)) for a_, b_ in ((2, 2), (2, 3), (3, 3))}
    del gi
    return {
        "screen_tidal_R_AxBx": screen,
        "g_xx": str(g[1, 1]),
        "Gamma_a_xx": [str(a) for a in acc],
        "pregeodesic": all(acc[a] == 0 for a in (0, 2, 3)),  # only the x-component survives
        "R_xx": str(_ricci(gam, xs, 1, 1)),
    }


def flat_where_k_depends_on_x_only() -> dict:
    """Riemann tensor of ER Eq. (13) with k = k(x): all components vanish."""
    t, x, rho, phi = sp.symbols("t x rho phi", real=True)
    k = sp.Function("k")(x)
    xs = [t, x, rho, phi]
    g = sp.Matrix(
        [[-1, (1 - k) / 2, 0, 0], [(1 - k) / 2, k, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]]
    )
    gam = _christoffel(g, xs)
    nonzero = []
    for a, b, c, d in itertools.product(range(4), repeat=4):
        e = sp.diff(gam[a][b][d], xs[c]) - sp.diff(gam[a][b][c], xs[d])
        e += sum(gam[a][c][m] * gam[m][b][d] - gam[a][d][m] * gam[m][b][c] for m in range(4))
        if sp.simplify(e) != 0:
            nonzero.append((a, b, c, d))
    return {"nonzero_riemann_components": nonzero}


def max_outgoing_speed(k: float, n: int = 20000, seed: int = P.seed) -> float:
    """Largest dx/dt over random future causal vectors at a point with the given k (> 0)."""
    rng = np.random.default_rng(seed)
    best = -math.inf
    for _ in range(n):
        dx, dr, dph = rng.normal(size=3)
        # solve -(dt - dx)(dt + k dx) + dr^2 + dph^2 = -s, s >= 0, for the future root dt
        s = rng.exponential()
        q = dr**2 + dph**2 + s
        # dt^2 + (k - 1) dx dt - k dx^2 - q = 0
        bq = (k - 1) * dx
        dt = (-bq + math.sqrt(bq**2 + 4 * (k * dx**2 + q))) / 2
        best = max(best, dx / dt)
    return best


def theta(s, eps):
    """Compact C-infinity step (ER Eq. 7): exactly 0 for s <= 0 and 1 for s >= eps."""
    z = np.clip(np.asarray(s, dtype=float) / eps, 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore"):
        f = np.where(z > 0, np.exp(-1 / np.where(z > 0, z, 1)), 0.0)
        g = np.where(z < 1, np.exp(-1 / np.where(z < 1, 1 - z, 1)), 0.0)
    return f / (f + g)


def theta_tanh(s, eps):
    """ER Eq. (35): (tanh(2 (2 s / eps - 1)) + 1) / 2, with tails on both sides."""
    return (np.tanh(2 * (2 * np.asarray(s, dtype=float) / eps - 1)) + 1) / 2


def er_k(t, x, rho, p: Params = P, step=theta, front: bool = True):
    """ER Eq. (12); front=False gives the eternal idealisation (no front factor)."""
    f = step(t - x - rho, p.eps) if front else 1.0
    return 1 - p.eta * step(p.rho_max - rho, p.eps) * f * (
        step(x, p.eps) - step(x + p.eps - p.d_len, p.eps)
    )


def incoming_rays(p: Params = P) -> dict:
    """Integrate incoming rays dx/dt = -1/k at fixed rho (exact null curves of the (t, x) block),
    entering at x = D + 5 eps after completion, and record along each ray the largest
    |d k / d rho| and |d k / d u| (u = t - x, at fixed x): zero means the ray stays in the region
    where k depends on x alone, which is flat."""
    return _rays(p, launch_offset=2.0)


def _rays(p: Params, launch_offset: float) -> dict:
    """Rays enter at x = D + 5 eps at t = D + rho + launch_offset eps + U(0, 50) eps."""
    rng = np.random.default_rng(p.seed)
    h = 1e-4 * p.eps
    rho = rng.uniform(0, p.rho_max - 1.5 * p.eps, p.n_rays)
    t = p.d_len + rho + launch_offset * p.eps + rng.uniform(0, 50, p.n_rays)
    t0 = t.copy()
    dx = 1e-2 * p.eps
    worst = 0.0
    for x in np.arange(p.d_len + 5 * p.eps, -5 * p.eps, -dx):
        k = er_k(t, x, rho, p)
        dkr = np.abs(er_k(t, x, rho + h, p) - er_k(t, x, rho - h, p)) / (2 * h)
        dku = np.abs(er_k(t + h, x, rho, p) - er_k(t - h, x, rho, p)) / (2 * h)
        worst = max(worst, float(dkr.max()), float(dku.max()))
        t = t + k * dx  # dt = k |dx| along dx/dt = -1/k
    times = t - t0
    return {
        "max_dk_drho_or_du_along_rays": worst,
        "transit_time_min_eps": float(times.min()),
        "transit_time_max_eps": float(times.max()),
        "exterior_light_time_eps": p.d_len + 10 * p.eps,
    }


def early_rays(p: Params = P) -> dict:
    """Positive control: rays launched before completion must cross the front (diagnostic > 0)."""
    return _rays(p, launch_offset=-(p.d_len + 60.0))


@lru_cache(maxsize=1)
def axis_null_ricci() -> dict:
    """Incoming null direction l = (1, -1/k, 0, 0) for k = k(x, rho): R_ab l^a l^b, its axis limit
    with k = k0(x) + k2(x) rho^2 / 2 + ... (valid only if k is a smooth function of rho^2, i.e.
    smooth at the axis; not for the ER tanh profile), and the radial acceleration of the curve."""
    t, x, rho, phi = sp.symbols("t x rho phi", real=True)
    k = sp.Function("k")(x, rho)
    xs = [t, x, rho, phi]
    g = sp.Matrix(
        [[-1, (1 - k) / 2, 0, 0], [(1 - k) / 2, k, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]]
    )
    gam = _christoffel(g, xs)
    lv = [1, -1 / k, 0, 0]
    rll = sp.simplify(
        sum(
            _ricci(gam, xs, a, b) * lv[a] * lv[b]
            for a in range(4)
            for b in range(4)
            if lv[a] and lv[b]
        )
    )
    k0, k2, k4 = (sp.Function(n)(x) for n in ("k0", "k2", "k4"))
    near = sp.simplify(rll.subs(k, k0 + k2 * rho**2 / 2 + k4 * rho**4 / 24).doit())
    acc = sp.simplify(sum(gam[2][a][b] * lv[a] * lv[b] for a in range(4) for b in range(4)))
    return {
        "R_ll": str(sp.factor(rll)),
        "R_ll_on_axis": str(sp.factor(sp.limit(near, rho, 0))),
        "radial_acceleration": str(sp.factor(acc)),
    }


def radial_derivative_at_axis_tanh(eta: float, r_over_eps: float) -> str:
    """k_rho(0) for ER's tanh radial factor, k = 1 - eta theta(R - rho) (mid-tube, B = 1):
    k_rho(0) = eta theta'(R) = eta 2 sech^2(y) / eps with y = 2 (2 R / eps - 1). Nonzero, so
    k(sqrt(y^2 + z^2)) has a kink on the axis: its Cartesian y-derivative jumps by 2 k_rho(0)."""
    import mpmath as mp

    mp.mp.dps = 60
    y = 2 * (2 * mp.mpf(r_over_eps) - 1)
    return mp.nstr(eta * 2 * mp.sech(y) ** 2, 6)


def front_axis_derivative(t_minus_x_over_eps: float) -> float:
    """d/d rho of the front factor theta(t - x - rho) at rho = 0 (compact step, eps = 1): nonzero
    while the front crosses the axis (0 < t - x < eps): a kink there in Cartesian coordinates."""
    h = 1e-6
    return float(
        (theta(t_minus_x_over_eps - h, 1.0) - theta(t_minus_x_over_eps + h, 1.0)) / (2 * h)
    )


@lru_cache(maxsize=1)
def axial_normalisation() -> dict:
    """Eternal tube, axis, k0 = k(x, 0), k2 = d^2k/drho^2(x, 0) (k smooth at the axis). Affine
    tangent of the incoming axial ray with exterior Killing energy E = -g(d_t, l) = 1:
    l = (2 k0 / (1 + k0), -2 / (1 + k0), 0, 0). Returns g(l, l), E, R_ab l^a l^b, d lambda / dx."""
    k0, k2 = sp.symbols("k0 k2", positive=True)
    lt, lx = 2 * k0 / (1 + k0), -2 / (1 + k0)
    g = sp.Matrix([[-1, (1 - k0) / 2], [(1 - k0) / 2, k0]])
    lv = sp.Matrix([lt, lx])
    norm = sp.simplify((lv.T * g * lv)[0])
    energy = sp.simplify(-(sp.Matrix([1, 0]).T * g * lv)[0])
    rll_unit = -(1 + k0) * k2 / k0**2  # R_ll for l^t = 1 (axis_null_ricci)
    rll = sp.simplify(rll_unit * lt**2)
    dlam_dx = sp.simplify(1 / lx)
    return {
        "g_ll": str(norm),
        "killing_energy": str(energy),
        "R_ll_affine": str(rll),
        "dlambda_dx": str(dlam_dx),
        "anec_integrand_per_dx": str(sp.simplify(rll * sp.Abs(dlam_dx))),
    }


def _rll_general(k, k1, k2, rho):
    """General incoming R_ll (l^t = 1) for k = k(x, rho), keeping the k_rho / rho term."""
    return -((1 + k) * k2 - k1**2) / (2 * k**2) - (1 + k) * k1 / (2 * rho * k**2)


def smooth_axis_examples(eta: float = 0.999, r_len: float = 10.0, d_len: float = 400.0) -> dict:
    """Two radial factors that are smooth at the axis (functions of rho^2), eternal tube,
    k = 1 - eta g(rho) B(x) with B the compact box over [0, D]: g2 = exp(-(rho/R)^2) (k2 > 0) and
    g4 = exp(-(rho/R)^4) (k2 = 0 but curved off the axis). Mathematical checks, not source models.
    ANEC along the axial ray with E = 1: -2 int k2 dx, cross-checked against the general R_ll at
    rho = 1e-5 R times the affine weight."""
    xg = np.linspace(-5.0, d_len + 5.0, 200001)
    bx = theta(xg, 1.0) - theta(xg + 1.0 - d_len, 1.0)
    out = {}
    for name, power in (("exp_rho2", 2), ("exp_rho4", 4)):
        rho = 1e-5 * r_len
        gv = math.exp(-((rho / r_len) ** power))
        g1 = -power * rho ** (power - 1) / r_len**power * gv
        g2v = (power * (power - 1) * rho ** (power - 2) / r_len**power) * (-gv) + (
            power * rho ** (power - 1) / r_len**power
        ) ** 2 * gv
        g2_axis = -2 / r_len**2 if power == 2 else 0.0  # g''(0)
        k0 = 1 - eta * bx
        k2_axis = -eta * g2_axis * bx
        anec_formula = float(-2 * np.trapezoid(k2_axis, xg))
        kk = 1 - eta * gv * bx
        rll_small_rho = _rll_general(kk, -eta * g1 * bx, -eta * g2v * bx, rho)
        weight = (2 * k0 / (1 + k0)) ** 2 * (1 + k0) / 2  # (l^t)^2 |d lambda / dx|
        anec_direct = float(np.trapezoid(rll_small_rho * weight, xg))
        # off-axis curvature: R_ll at rho = R / 2, mid-tube
        rr = r_len / 2
        g_r = math.exp(-((rr / r_len) ** power))
        g1r = -power * rr ** (power - 1) / r_len**power * g_r
        g2r = (power * (power - 1) * rr ** (power - 2) / r_len**power) * (-g_r) + (
            power * rr ** (power - 1) / r_len**power
        ) ** 2 * g_r
        k_mid = 1 - eta * g_r
        out[name] = {
            "k_rho_at_axis": 0.0,
            "k2_axis_mid_tube": float(-eta * g2_axis) + 0.0,
            "anec_minus_2_int_k2": anec_formula + 0.0,
            "anec_direct_rho_1e-5R": anec_direct,
            "R_ll_unit_at_rho_R_over_2_mid_tube": float(
                _rll_general(k_mid, -eta * g1r, -eta * g2r, rr)
            ),
            "axis_minimises_k": True,  # g decreasing in rho
        }
    return out


def w_monotone_check(step, n: int = 5000, p: Params = P) -> dict:
    """Random points and future causal vectors in the eternal tube: dw = dt + k(x, 0) dx >= 0."""
    rng = np.random.default_rng(p.seed + 1)
    worst = math.inf
    for _ in range(n):
        x = rng.uniform(-3, p.d_len + 3)
        rho = rng.uniform(0, p.rho_max + 3)
        k = float(er_k(0.0, x, rho, p, step=step, front=False))
        k_ax = float(er_k(0.0, x, 0.0, p, step=step, front=False))
        dx, dr, dph = rng.normal(size=3)
        q = dr**2 + dph**2 + rng.exponential() * rng.integers(0, 2)  # null or timelike
        bq = (k - 1) * dx
        dt = (-bq + math.sqrt(bq**2 + 4 * (k * dx**2 + q))) / 2
        worst = min(worst, (dt + k_ax * dx) / math.hypot(dt, dx))
    return {"min_normalised_dw": worst}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived symbolically from ER Eqs. (12)-(13) (no data); numerical checks use "
        "the stated ASSUMPTION profiles (compact C-infinity step; ER Eq. 35 tanh step)",
        "params": asdict(P),
        "outgoing_lines": outgoing_lines(),
        "flat_where_k_depends_on_x_only": flat_where_k_depends_on_x_only(),
        "max_outgoing_speed": {f"k={k:g}": max_outgoing_speed(k) for k in (1.0, 0.5, 1e-3)},
        "incoming_rays_after_completion": incoming_rays(),
        "incoming_rays_before_completion_control": early_rays(),
        "axis_null_ricci": axis_null_ricci(),
        "axis_regularity": {
            "tanh_radial_k_rho_at_axis_eta0.999_R10eps": radial_derivative_at_axis_tanh(
                0.999, 10.0
            ),
            "compact_radial_k_rho_at_axis": 0.0,
            "front_d_drho_at_axis_t_minus_x_half_eps": front_axis_derivative(0.5),
            "front_d_drho_at_axis_t_minus_x_2eps": front_axis_derivative(2.0),
        },
        "axial_normalisation": axial_normalisation(),
        "smooth_axis_examples": smooth_axis_examples(),
        "w_monotone_check": {
            "compact": w_monotone_check(theta),
            "tanh": w_monotone_check(theta_tanh),
        },
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
