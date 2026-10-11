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
3. Axial incoming ray. For k(x, rho) (eternal tube, or after completion with the compact front),
   w = t + int k(x, 0) dx never decreases along causal curves if the axis minimises k(x, .) for
   every x (true for ER Eq. 12 with any step that is monotone in rho). The axial ray keeps w
   constant, so it is a complete achronal null geodesic. On it R_ll = -(1 + k0) k2 / k0^2 with
   k0 = k(x, 0), k2 = d^2 k / d rho^2 (x, 0) >= 0: its ANEC integral is <= 0, and < 0 unless k2
   vanishes along the whole axis (compact step: 0, a marginal case; tanh step: < 0).
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
    with k = k0(x) + k2(x) rho^2 / 2 + ..., and the radial acceleration of the curve."""
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


def axis_anec_density(eta: float, rho_max_over_eps: float, step: str) -> str:
    """R_ll (l^t = 1) on the axis at mid-tube, -(1 + k0) k2 / k0^2, in units of 1/eps^2, for the
    eternal tube. Compact step with rho_max >= eps: theta is identically 1 near the axis, so k2 = 0.
    Tanh step (ER Eq. 35): theta(z) = (tanh y + 1) / 2 with y = 2 (2 z - 1), so
    theta'' = -16 sech^2(y) tanh(y) and k2 = -eta theta''(rho_max) > 0 (evaluated in mpmath; the
    tail is ~exp(-8 rho_max / eps))."""
    import mpmath as mp

    if step == "compact":
        return "0" if rho_max_over_eps >= 1 else "undefined (rho_max < eps)"
    mp.mp.dps = 60
    y = 2 * (2 * mp.mpf(rho_max_over_eps) - 1)
    k0 = 1 - eta * (mp.tanh(y) + 1) / 2
    k2 = 16 * eta * mp.sech(y) ** 2 * mp.tanh(y)
    return mp.nstr(-(1 + k0) * k2 / k0**2, 6)


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
        "axis_anec_density_mid_tube": {
            f"{step}_rho_max_{r:g}eps_eta_{e:g}": axis_anec_density(e, r, step)
            for step in ("compact", "tanh")
            for r in (3.0, 10.0, 50.0)
            for e in (0.95, 0.999)
        },
        "w_monotone_check": {
            "compact": w_monotone_check(theta),
            "tanh": w_monotone_check(theta_tanh),
        },
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
