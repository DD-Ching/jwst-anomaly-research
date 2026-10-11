"""Gap G2 of goal round 6: do the shortcut null geodesics of the full Krasnikov tube violate ANEC?

Goal round 7 (docs/goal_round7/README.md). Offline, symbolic; no data. Geometry: the complete
Everett & Roman (gr-qc/9702049, "ER") tube, Eqs. (12)-(13), including the front theta(t - x - rho)
laid at light speed and both end caps:

    ds^2 = -dt^2 + (1 - k) dx dt + k dx^2 + drho^2 + rho^2 dphi^2
         = -(dt - dx)(dt + k dx) + drho^2 + rho^2 dphi^2,
    k = 1 - eta theta(rho_max - rho) theta(t - x - rho) [theta(x) - theta(x + eps - D)].

Checked here [derived]:
1. Outgoing lines u = t - x = const, rho = const, phi = const are null pregeodesics for ANY
   k(u, x, rho), and R_ab l^a l^b vanishes identically along them (R_xx = 0 in (u, x, rho, phi)
   coordinates). They are achronal whenever k > 0 (eta < 1): every causal curve has
   dx/dt <= 1, with equality only along these lines.
2. Wherever k depends on x alone the spacetime is flat (Riemann = 0). In the ER tube this holds in
   the core {rho < rho_max - eps, t - x - rho > eps}, end caps included, and outside the tube.
3. An incoming (-x) ray has d(t - x)/dt = 1 + 1/k > 0. So a ray that is in the core when it enters
   at x = D after completion stays in the core: its complete geodesic lies in a flat region
   (ANEC integral 0; the generic condition fails). This is checked numerically on the ER profile.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass

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
    return {
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
    """ER's smooth step: 0 for s < 0, 1 for s > eps (the round-5 tanh form, clipped)."""
    s = np.asarray(s, dtype=float)
    y = (np.tanh(2 * (2 * s / eps - 1)) + 1) / 2
    return np.where(s <= 0, 0.0, np.where(s >= eps, 1.0, y))


def er_k(t, x, rho, p: Params = P):
    """ER Eq. (12)."""
    return 1 - p.eta * theta(p.rho_max - rho, p.eps) * theta(t - x - rho, p.eps) * (
        theta(x, p.eps) - theta(x + p.eps - p.d_len, p.eps)
    )


def incoming_rays(p: Params = P) -> dict:
    """Integrate incoming rays dx/dt = -1/k at fixed rho (exact null curves of the (t, x) block),
    entering at x = D + 5 eps after completion, and record along each ray the largest
    |d k / d rho| and |d k / d u| (u = t - x, at fixed x): zero means the ray stays in the region
    where k depends on x alone, which is flat."""
    rng = np.random.default_rng(p.seed)
    h = 1e-4 * p.eps
    rho = rng.uniform(0, p.rho_max - 1.5 * p.eps, p.n_rays)
    t = p.d_len + rho + 2 * p.eps + rng.uniform(0, 50, p.n_rays)  # after the front passed x = D
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived symbolically from ER Eqs. (12)-(13) (no data); the ray check uses "
        "the stated ASSUMPTION profile",
        "params": asdict(P),
        "outgoing_lines": outgoing_lines(),
        "flat_where_k_depends_on_x_only": flat_where_k_depends_on_x_only(),
        "max_outgoing_speed": {f"k={k:g}": max_outgoing_speed(k) for k in (1.0, 0.5, 1e-3)},
        "incoming_rays": incoming_rays(),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
