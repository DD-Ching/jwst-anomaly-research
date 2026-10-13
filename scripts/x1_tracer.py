"""X1 S1: null-geodesic tracer for the static Simpson-Visser (SV) metric, with the S1 checks.

Metric (units G = c = M = 1): ds^2 = -f dt^2 + dr^2 / f + R^2 dOmega^2, f = 1 - 2/R,
R = sqrt(r^2 + a^2), r in (-inf, inf). For a > 2 the throat r = 0 (R = a) is traversable;
r < 0 is the "far side". a = 0 is exactly Schwarzschild.

Null rays have E = 1 and impact parameter b. With u = 1/R the orbit equation is exact quadrature:
  dphi/du = b / ( sqrt(1 - a^2 u^2) * sqrt(1 - b^2 u^2 (1 - 2u)) )
- near side: a ray from our infinity turns at the smallest root u_t of 1 - b^2 u^2 (1 - 2u)
  if u_t < 1/a;
- otherwise (b < b_c(a)) it crosses the throat at u = 1/a and continues on the far side.
Both endpoint singularities are integrable and removed by u = u_end (1 - s^2) substitutions.

Images of a point emitter at areal radius R_e and unit position n (observer direction o) satisfy
  Phi(b) = gamma + 2 pi k   (image on the side of n's sky projection), or
  Phi(b) = 2 pi - gamma + 2 pi k  (opposite side),  cos gamma = n . o.
Point-source flux per image: g^4 * b |db/dPhi| / (|cos alpha_e| sin(gamma)), with alpha_e the
emission angle to the radial direction (bolometric; g = redshift factor of a Keplerian emitter;
Omega^2 = 1/R^3 holds exactly in the areal radius). Light-travel-time delays are optional
(spot_track(delays=True): each image is shifted to its observed phase).

An independent second route (ODE in the affine parameter, r'' = -W'(r)/2 with W = b^2 f / R^2)
gives the capture/transmission boundary and the far-side swept angle.

  python scripts/x1_tracer.py   -> results/x1/s1_checks.json

All outputs are model computations of a conditional, unestablished geometry
(docs/state/units/X1.md); none is evidence that such a wormhole exists.
"""

from __future__ import annotations

import functools
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy import integrate, optimize, special

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_farside as wf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x1"
QUAD = {"epsabs": 0.0, "epsrel": 1e-10, "limit": 400}
R_START = 2.0e3  # ODE rays start here (M); the regular range beyond is added by quadrature
N_GRID = 1200  # b grid per image curve (converged to 1e-10 against 4,800 points; X1.md)
# ASSUMPTION: S1 pass tolerances, pre-written in docs/state/units/X1.md "Tolerances"
TOL = {"deflection_rel": 1e-6, "b_crit_rel": 1e-4, "swept_angle_rel": 1e-6, "appA_M": 0.5}


def _pot(u: float, b: float) -> float:
    """1 - b^2 u^2 (1 - 2u): positive where the ray can propagate."""
    return 1.0 - b * b * u * u * (1.0 - 2.0 * u)


def turning_u(b: float, a: float) -> float | None:
    """Smallest root u_t of _pot below min(1/3, 1/a) if the ray turns on the near side."""
    u_hi = 1.0 / 3.0  # _pot has its minimum at u = 1/3
    if a > 0:
        u_hi = min(u_hi, 1.0 / a)
    if _pot(u_hi, b) > 0:
        return None
    return optimize.brentq(_pot, 0.0, u_hi, args=(b,), xtol=1e-15, rtol=1e-15)


def _phi_to_turn(b: float, a: float, ut: float, u0: float = 0.0) -> float:
    """Integral of dphi/du from u0 to the turning point ut (u = ut (1 - s^2), regular in s)."""

    def f(s):
        u = ut * (1.0 - s * s)
        # _pot(u) = (ut - u) * q exactly, with ut - u = ut s^2
        q = b * b * ((u + ut) - 2.0 * (u * u + u * ut + ut * ut))
        return 2.0 * math.sqrt(ut) * b / (math.sqrt(1.0 - a * a * u * u) * math.sqrt(q))

    return integrate.quad(f, 0.0, math.sqrt(max(0.0, 1.0 - u0 / ut)), **QUAD)[0]


def _phi_range(b: float, a: float, u0: float, u1: float) -> float:
    """Integral of dphi/du over [u0, u1]; the only singular endpoint allowed is u1 = 1/a."""
    if a > 0 and abs(u1 - 1.0 / a) < 1e-14:
        ua = 1.0 / a
        s0 = math.sqrt(max(0.0, 1.0 - u0 / ua))

        def f(s):
            u = ua * (1.0 - s * s)
            # 1 - a^2 u^2 = s^2 (2 - s^2)
            return 2.0 * ua * b / (math.sqrt(2.0 - s * s) * math.sqrt(_pot(u, b)))

        return integrate.quad(f, 0.0, s0, **QUAD)[0]

    def g(u):
        return b / (math.sqrt(1.0 - a * a * u * u) * math.sqrt(_pot(u, b)))

    return integrate.quad(g, u0, u1, **QUAD)[0]


def deflection_schw_numeric(b: float) -> float:
    """Tracer deflection angle for a = 0 (Schwarzschild): 2 Phi(0 -> u_t) - pi."""
    ut = turning_u(b, 0.0)
    return 2.0 * _phi_to_turn(b, 0.0, ut) - math.pi


def deflection_schw_exact(b: float) -> float:
    """Darwin's elliptic-integral Schwarzschild deflection (periastron P: b^2 = P^3 / (P - 2))."""
    p = optimize.brentq(lambda x: x**3 / (x - 2.0) - b * b, 3.0, 10.0 * b)
    q = math.sqrt((p - 2.0) * (p + 6.0))
    m = (q - p + 6.0) / (2.0 * q)
    phi_inf = math.asin(math.sqrt((q - p + 2.0) / (q - p + 6.0)))
    return -math.pi + 4.0 * math.sqrt(p / q) * (special.ellipk(m) - special.ellipkinc(phi_inf, m))


def phi_to_emitter(b: float, a: float, re: float, side: str) -> list[float]:
    """Swept angles from our infinity to an emitter at areal radius re.

    side 'near': [direct, after the turn] when the ray turns inside re; [direct] when it does
    not turn (captured or transmitted); [] when it turns outside re.
    side 'far': [through the throat] when b < b_c(a) (no turn); [] otherwise.
    """
    if re <= max(a, 3.0):
        raise ValueError(
            f"emitter radius {re} must exceed the throat a = {a} and the photon sphere"
        )
    ue = 1.0 / re
    ut = turning_u(b, a)
    if side == "near":
        if ut is None:  # no near-side turn: the ray crosses R = re once on the way in
            return [_phi_range(b, a, 0.0, ue)]
        if ut < ue:  # turns outside the emitter radius
            return []
        pt = _phi_to_turn(b, a, ut)
        pe = pt - _phi_to_turn(b, a, ut, ue)
        return [pe, 2.0 * pt - pe]
    if ut is not None or a <= 2.0:
        return []
    ua = 1.0 / a
    pin = _phi_range(b, a, 0.0, ua)
    return [pin + _phi_range(b, a, ue, ua)]


def _dtdu_reg(u: float, b: float, a: float) -> float:
    """dt/du - 1/u^2 - 2/u (finite at u = 0); dt/du = 1/(u^2 f sqrt(1 - a^2 u^2) sqrt(_pot))."""
    if u < 1e-6:
        return 4.0 + 0.5 * (a * a + b * b)
    log_h = -math.log1p(-2.0 * u) - 0.5 * math.log1p(-a * a * u * u)
    log_h -= 0.5 * math.log1p(-b * b * u * u * (1.0 - 2.0 * u))
    return (math.expm1(log_h) - 2.0 * u) / (u * u)


def _t_open(b: float, a: float, u1: float) -> float:
    """Coordinate time from our infinity to u1 (no turn), up to a constant common to all rays."""
    return (
        integrate.quad(_dtdu_reg, 0.0, u1, args=(b, a), **QUAD)[0] - 1.0 / u1 + 2.0 * math.log(u1)
    )


def _t_turn(b: float, a: float, ut: float, u0: float) -> float:
    """Coordinate time from u0 (> 0) to the turning point ut (u = ut (1 - s^2), regular in s)."""

    def f(s):
        u = ut * (1.0 - s * s)
        q = b * b * ((u + ut) - 2.0 * (u * u + u * ut + ut * ut))
        fu = 1.0 - 2.0 * u
        return 2.0 * math.sqrt(ut) / (u * u * fu * math.sqrt(1.0 - a * a * u * u) * math.sqrt(q))

    return integrate.quad(f, 0.0, math.sqrt(max(0.0, 1.0 - u0 / ut)), **QUAD)[0]


def _t_throat(b: float, a: float, u0: float) -> float:
    """Coordinate time from u0 (> 0) to the throat u = 1/a (u = (1 - s^2)/a, regular in s)."""
    ua = 1.0 / a

    def f(s):
        u = ua * (1.0 - s * s)
        return 2.0 * ua / (u * u * (1.0 - 2.0 * u) * math.sqrt(2.0 - s * s) * math.sqrt(_pot(u, b)))

    return integrate.quad(f, 0.0, math.sqrt(max(0.0, 1.0 - u0 / ua)), **QUAD)[0]


def time_to_emitter(b: float, a: float, re: float, side: str) -> list[float]:
    """Light travel times from the emitter to the observer, branch by branch as in phi_to_emitter.

    Relative values only: the divergent distance term is dropped identically for every ray.
    """
    ue = 1.0 / re
    ut = turning_u(b, a)
    if side == "near":
        if ut is None:
            return [_t_open(b, a, ue)]
        if ut < ue:
            return []
        um = 0.5 * ue
        t_half = _t_open(b, a, um)
        to_turn = t_half + _t_turn(b, a, ut, um)
        te = to_turn - _t_turn(b, a, ut, ue)
        return [te, 2.0 * to_turn - te]
    if ut is not None or a <= 2.0:
        return []
    ua = 1.0 / a
    um = 0.5 * min(ue, ua)
    to_throat = _t_open(b, a, um) + _t_throat(b, a, um)
    return [to_throat + _t_throat(b, a, ue)]


def b_crit_ode(a: float, r_start: float = R_START, tol: float = 1e-6) -> float:
    """Capture/transmission boundary from ODE rays: largest b whose ray reaches r < 0 (far side)."""

    def w_prime(r, b):
        rr = math.sqrt(r * r + a * a)
        # W = b^2 (1/R^2 - 2/R^3); dW/dr = b^2 (-2/R^3 + 6/R^4) * r / R
        return b * b * (-2.0 / rr**3 + 6.0 / rr**4) * r / rr

    def transmits(b):
        def rhs(_, y):
            return [y[1], -0.5 * w_prime(y[0], b)]

        def crossed(_, y):
            return y[0]

        def turned(_, y):
            return y[1]

        crossed.terminal = turned.terminal = True
        turned.direction = 1.0
        rr = math.sqrt(r_start**2 + a * a)
        v0 = -math.sqrt(max(0.0, 1.0 - b * b * (1.0 - 2.0 / rr) / rr**2))
        sol = integrate.solve_ivp(
            rhs, (0.0, 1e7), [r_start, v0], events=(crossed, turned), rtol=1e-12, atol=1e-12
        )
        return len(sol.t_events[0]) > 0

    lo, hi = 0.5 * wf.b_crit(a), 1.5 * wf.b_crit(a)
    while (hi - lo) / lo > tol:
        mid = 0.5 * (lo + hi)
        if transmits(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def phi_ode(b: float, a: float, re: float, side: str, branch: int = 0) -> float:
    """Second route for the swept angle to an emitter at areal radius re.

    Integrates (r, r', phi) in the affine parameter from R_START inwards (through the throat for
    side 'far') to the branch-th crossing of the emitter radius; the regular range beyond R_START
    is added by quadrature. Near side: branch 0 = inward crossing, 1 = outward after the turn.
    """

    def rhs(_, y):
        r = y[0]
        rr = math.sqrt(r * r + a * a)
        wp = b * b * (-2.0 / rr**3 + 6.0 / rr**4) * r / rr
        return [y[1], -0.5 * wp, b / rr**2]

    r_e = math.sqrt(re * re - a * a) * (1.0 if side == "near" else -1.0)

    def reached(_, y):
        return y[0] - r_e

    reached.terminal = side == "far" or branch == 0
    reached.direction = 1.0 if (side == "near" and branch == 1) else -1.0
    rr0 = math.sqrt(R_START**2 + a * a)
    v0 = -math.sqrt(_pot(1.0 / rr0, b))
    sol = integrate.solve_ivp(
        rhs, (0.0, 1e7), [R_START, v0, 0.0], events=reached, rtol=1e-12, atol=1e-13
    )
    return float(sol.y_events[0][0][2]) + _phi_range(b, a, 0.0, 1.0 / rr0)


def _fold_b(a: float, re: float) -> float:
    """Impact parameter of the ray that turns exactly at R = re: b = re / sqrt(f(re))."""
    return re / math.sqrt(1.0 - 2.0 / re)


@functools.lru_cache(maxsize=32)
def _tables(a: float, re: float, side: str, b_max: float, n: int = N_GRID):
    """Phi along one continuous image curve, as rows (b, branch, Phi).

    Near side: the direct branch for b up to the fold b_f (ray tangent at re), then the
    after-turn branch back down in b; the grid is refined geometrically towards the fold.
    Far side: one monotonic branch.
    """
    if side == "near":
        bf = _fold_b(a, re) * (1.0 - 1e-12)
        g = np.unique(
            np.concatenate([np.linspace(1e-4, bf, n), bf * (1.0 - np.geomspace(1e-11, 1e-2, 300))])
        )
        g = g[g <= bf]
        path = [(b, 0) for b in g] + [(b, 1) for b in g[::-1]]
    else:
        path = [(b, 0) for b in np.linspace(1e-4, b_max, n)]
    rows = []
    for b, br in path:
        ph = phi_to_emitter(b, a, re, side)
        rows.append((b, br, ph[br] if len(ph) > br else np.nan))
    return np.array(rows), (_fold_b(a, re) if side == "near" else b_max)


def _images(a: float, re: float, side: str, gamma: float, tables, kmax: int = 0):
    """Images (b, side_sign, dPhi/db, branch) with Phi(b) = gamma + 2 pi k or 2 pi - gamma + 2 pi k.

    Returns (images, n_dropped); an image is dropped only if its finite difference leaves
    the branch.
    """
    rows, b_end = tables
    targets = []
    for k in range(kmax + 1):
        targets += [(gamma + 2 * math.pi * k, +1.0), (2 * math.pi - gamma + 2 * math.pi * k, -1.0)]
    out, dropped = [], 0
    for tgt, sgn in targets:
        d = rows[:, 2] - tgt
        ok = (
            np.isfinite(d[:-1])
            & np.isfinite(d[1:])
            & (np.sign(d[:-1]) != np.sign(d[1:]))
            & (rows[:-1, 1] == rows[1:, 1])
        )
        for j in np.nonzero(ok)[0]:
            br = int(rows[j, 1])

            def h(b, br=br, tgt=tgt):
                ph = phi_to_emitter(b, a, re, side)
                return ph[br] - tgt if len(ph) > br else math.nan

            lo, hi = sorted((rows[j, 0], rows[j + 1, 0]))
            bs = optimize.brentq(h, lo, hi, xtol=1e-14, rtol=1e-14)
            eps = min(1e-6 * bs, 0.25 * (b_end - bs)) if b_end > bs else 1e-6 * bs
            hp, hm = h(bs + eps), h(bs - eps)
            if not (np.isfinite(hp) and np.isfinite(hm)):
                dropped += 1
                continue
            out.append((bs, sgn, (hp - hm) / (2 * eps), br))
    return out, dropped


def spot_track(
    a: float,
    re: float,
    incl_deg: float,
    side: str,
    n_phase: int = 72,
    kmax: int = 0,
    primary_only: bool = False,
    eta: float = 1.0,
    delays: bool = False,
):
    """Centroid track (M on the sky) and total flux of a point spot on a circular orbit.

    Rows: (phase, x_sky, y_sky, flux, n_dropped_images). The angular velocity is eta times
    Keplerian (eta < 1: sub-Keplerian). Without delays, phase is the emission phase of every image;
    with delays=True each image is shifted to its observed phase (emission phase + Omega * light
    travel time, relative across images) and the images are resampled on a common observed-phase
    grid before they are summed (n_dropped_images is then the per-orbit total in the first row).

    Observer direction o = (sin i, 0, cos i); orbit in the z = 0 plane (of the emitter's own sheet).
    Sky axes: x_sky along (cos i, 0, -sin i), y_sky along (0, 1, 0).
    """
    if delays and kmax > 0:
        raise ValueError("delays=True supports k = 0 images only (image classes keyed by side)")
    inc = math.radians(incl_deg)
    o = np.array([math.sin(inc), 0.0, math.cos(inc)])
    ex = np.array([math.cos(inc), 0.0, -math.sin(inc)])
    ey = np.array([0.0, 1.0, 0.0])
    omega = eta * re**-1.5
    ut_inv = math.sqrt(1.0 - 2.0 / re - eta * eta / re)  # 1/u^t = sqrt(f - R^2 Omega^2)
    b_max = 3.0 * re if side == "near" else wf.b_crit(a) * (1 - 1e-9)
    tables = _tables(a, re, side, b_max)
    rows, per_image = [], {}
    phases = np.linspace(0.0, 2 * math.pi, n_phase, endpoint=False)
    for ph in phases:
        n = np.array([math.cos(ph), math.sin(ph), 0.0])
        cg = float(np.clip(n @ o, -1, 1))
        gamma = math.acos(cg)
        proj = n - cg * o
        pn = np.linalg.norm(proj)
        e = proj / pn if pn > 1e-12 else ex
        cx = cy = ftot = 0.0
        images, dropped = _images(a, re, side, gamma, tables, 0 if primary_only else kmax)
        for bs, sgn, dphidb, br in images:
            if primary_only and sgn < 0:
                continue
            pos = sgn * bs * e  # 3-vector in the sky plane
            # photon axial angular momentum (about z) for E = 1, from the screen position
            lz = float(np.cross(pos, o)[2])
            g = ut_inv / (1.0 - omega * lz)
            # emission solid angle per observer solid angle: sin(al) dal / (sin(g) dg), al = angle
            # of the photon to the radial direction at the emitter (static frame);
            # cos(al) = sqrt(_pot(1/re, b))
            cos_al = math.sqrt(max(_pot(1.0 / re, bs), 1e-12))
            mag = bs / (abs(dphidb) * cos_al * max(math.sin(gamma), 1e-9))
            w = g**4 * mag
            if delays:
                dt = time_to_emitter(bs, a, re, side)[br]
                per_image.setdefault(sgn, []).append(
                    (ph + omega * dt, w * float(pos @ ex), w * float(pos @ ey), w)
                )
            cx += w * float(pos @ ex)
            cy += w * float(pos @ ey)
            ftot += w
        cx, cy = (cx / ftot, cy / ftot) if ftot else (np.nan, np.nan)
        rows.append((ph, cx, cy, ftot, dropped))
    if not delays:
        return np.array(rows)
    # Each image class (same-side / opposite-side) is one continuous curve in emission phase;
    # resample its weighted position and weight periodically on the observed-phase grid.
    sx = np.zeros(n_phase)
    sy = np.zeros(n_phase)
    sw = np.zeros(n_phase)
    for pts in per_image.values():
        # k = 0 images exist at every phase; a few dropped ones (counted) are interpolated over
        if len(pts) < n_phase // 2:
            raise RuntimeError(f"image class present at {len(pts)} of {n_phase} phases")
        arr = np.array(pts)
        t = np.unwrap(arr[:, 0])
        if np.any(np.diff(t) <= 0) or t[-1] - t[0] >= 2 * math.pi:
            raise RuntimeError("observed phase is not monotonic in emission phase over one orbit")
        t0 = t - 2 * math.pi * np.floor(t[0] / (2 * math.pi))
        tt = np.concatenate([t0 - 2 * math.pi, t0, t0 + 2 * math.pi])
        for acc, col in ((sx, 1), (sy, 2), (sw, 3)):
            acc += np.interp(phases, tt, np.tile(arr[:, col], 3))
    n_drop = int(sum(r[4] for r in rows))
    return np.column_stack([phases, sx / sw, sy / sw, sw, np.r_[n_drop, np.zeros(n_phase - 1)]])


def run_checks() -> dict:
    out: dict = {
        "labels": "derived (our computation) of a conditional geometry; tolerances from X1.md"
    }
    c1 = []
    for b in (6.0, 8.0, 20.0):
        num, ex = deflection_schw_numeric(b), deflection_schw_exact(b)
        c1.append({"b_M": b, "tracer": num, "exact": ex, "rel_err": abs(num / ex - 1)})
    out["check1_schwarzschild_deflection"] = {
        "rows": c1,
        "tolerance": TOL["deflection_rel"],
        "pass": all(r["rel_err"] < TOL["deflection_rel"] for r in c1),
    }
    # ODE rays only: the quadrature's turn test reduces algebraically to the analytic b_c
    c2 = []
    for a in (2.1, 3.0, 5.0, 7.0):
        bo, ba = b_crit_ode(a), wf.b_crit(a)
        c2.append({"a_M": a, "analytic": ba, "ode": bo, "rel_err_ode": abs(bo / ba - 1)})
    out["check2_b_crit"] = {
        "rows": c2,
        "tolerance": TOL["b_crit_rel"],
        "pass": all(r["rel_err_ode"] < TOL["b_crit_rel"] for r in c2),
        "note": "this boundary also tests check 4: rays above it never reach r < 0",
    }
    c2b = []
    cases = [(a, "far", 0, fr) for a in (2.1, 3.0, 5.0, 7.0) for fr in (0.3, 0.9, 0.999)]
    cases += [
        (a, "near", br, fr) for a in (0.0, 3.0, 5.0, 7.0) for br in (0, 1) for fr in (1.001, 1.3)
    ]
    for a, side, br, fr in cases:
        b = fr * wf.b_crit(a)
        q = phi_to_emitter(b, a, 9.0, side)
        if len(q) <= br:
            continue
        o = phi_ode(b, a, 9.0, side, br)
        row = {"a_M": a, "side": side, "branch": br, "b_over_bc": fr, "quadrature": q[br], "ode": o}
        c2b.append(row | {"rel_err": abs(q[br] / o - 1)})
    out["check2b_swept_angle_R9"] = {
        "rows": c2b,
        "tolerance": TOL["swept_angle_rel"],
        "pass": all(r["rel_err"] < TOL["swept_angle_rel"] for r in c2b),
        "note": "added in S1 (not pre-written): quadrature vs ODE swept angle, far side below b_c"
        " and near side (both branches) above it",
    }
    # GRAVITY 2023 App. A (after Gralla et al. 2019, Gates et al. 2020): the loop appears at
    # r_o = R + 1 M. That relation describes the direct (primary) image; we report the primary-only
    # loop (the check) and the centroid of all light (primary + secondary), which S2 must use.
    c3, tracks = {}, {}
    for name, inc in (("face_on_i179.5", 179.5), ("i155", 155.0)):
        for prim in (True, False):
            key = f"{name}_{'primary' if prim else 'all_light'}"
            tr = tracks[key] = spot_track(0.0, 9.0, inc, "near", n_phase=72, primary_only=prim)
            rad = np.hypot(tr[:, 1], tr[:, 2])
            c3[key] = {
                "semi_major_M": float(rad.max()),
                "mean_radius_M": float(rad.mean()),
                "min_radius_M": float(rad.min()),
                "flux_max_over_min": float(tr[:, 3].max() / tr[:, 3].min()),
                "dropped_images": int(tr[:, 4].sum()),
            }
    prim_sm = c3["i155_primary"]["semi_major_M"]
    out["check3_schwarzschild_spot_R9_i155"] = {
        "rows": c3,
        "gravity_appA_R_plus_1_M": 10.0,
        "diff_primary_semi_major_M": prim_sm - 10.0,
        "diff_primary_face_on_M": c3["face_on_i179.5_primary"]["mean_radius_M"] - 10.0,
        "tolerance_M": TOL["appA_M"],
        "pass": bool(abs(prim_sm - 10.0) <= TOL["appA_M"]),
        "note": "App. A is an approximation for the direct image, so this is a sanity check, not"
        " an exact benchmark; k = 0 images only; no light-travel-time delays;"
        " bolometric g^4 weights",
    }
    near_flux = float(tracks["i155_all_light"][:, 3].mean())
    c4 = []
    for a in (2.1, 3.0, 5.0, 7.0):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", integrate.IntegrationWarning)
            trf = spot_track(a, 9.0, 155.0, "far", n_phase=72)
        n_quad = sum(issubclass(w.category, integrate.IntegrationWarning) for w in caught)
        radf = np.hypot(trf[:, 1], trf[:, 2])
        c4.append(
            {
                "a_M": a,
                "b_crit_M": wf.b_crit(a),
                "centroid_radius_max_M": float(np.nanmax(radf)),
                "centroid_radius_mean_M": float(np.nanmean(radf)),
                "flux_rel_to_schwarzschild_near_mean": float(np.nanmean(trf[:, 3]) / near_flux),
                "quad_roundoff_warnings": int(n_quad),
                "dropped_images": int(trf[:, 4].sum()),
                "inside_b_crit": bool(np.nanmax(radf) < wf.b_crit(a)),
            }
        )
    out["check4_far_side_inside_b_crit"] = {
        "rows": c4,
        "pass": all(r["inside_b_crit"] for r in c4),
        "note": "far-side tables stop at b_c, so inside_b_crit holds by construction; the"
        " transmission boundary itself is tested by the ODE rays of check 2. The figure"
        " comparison with 2408.07350 Fig. 7 is figure-read only and not done",
    }
    return out


def main() -> int:
    out = run_checks()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s1_checks.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
