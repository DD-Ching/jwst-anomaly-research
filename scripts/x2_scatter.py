"""X2 S1: test-field scattering on Simpson-Visser barriers, with the S1 checks.

Metric (G = c = M = 1): f = 1 - 2/R, R = sqrt(r^2 + a^2), dr* = dr / f, 2 < a < 3 (two barriers,
one on each side of the throat r = 0). Test fields of spin s = 0 (scalar) and s = 1
(electromagnetic) obey psi'' + (w^2 - V(r*)) psi = 0 with (Churilova & Stuchlik, arXiv:1911.11823,
Eq. 6, rewritten with R)
  V = f l(l+1)/R^2 + (1 - s) (1/R) d^2R/dr*^2,   d^2R/dr*^2 = f d/dr (f r / R).
For s = 1 the barrier peaks at R = 3 exactly; for s = 0 the peak is found numerically.

One barrier = the potential on our side of the throat, r* >= 0. A wave coming from the throat
side meets it and is transmitted to our infinity (T) or reflected back towards the throat (R). For
the echo-relevant range a - 2 <~ 1e-2 the potential at the throat V(0) is <~ 1e-2 of the peak; the
throat side is treated as a region of constant V(0).

Two independent routes: (A) Numerov on a uniform r* grid, vectorised over frequency; (B) adaptive
DOP853 on the ODE with V interpolated from r(r*). Route B is used for spot checks only.

  python scripts/x2_scatter.py   -> results/x2/s1_checks.json  (about a minute)

Test-field computations of a conditional geometry. Gravitational waves are NOT shown to obey these
potentials (docs/state/units/X2.md, premise iii).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.special import lambertw

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x2"


# ------------------------------------------------------------------------------------ geometry
def f_of_r(r, a):
    return 1.0 - 2.0 / np.sqrt(r * r + a * a)


def potential_r(r, a, ell, s):
    """V(r): s = 0 or 1 from Eq. 6 of 1911.11823; s = "rw" is the ad hoc proxy
    f [l(l+1)/R^2 - 6/R^3] (Schwarzschild Regge-Wheeler form on the SV areal radius;
    not derived)."""
    R2 = r * r + a * a
    R = np.sqrt(R2)
    f = 1.0 - 2.0 / R
    if s == "rw":
        return f * (ell * (ell + 1) / R2 - 6.0 / R**3)
    fp = 2.0 * r / R**3  # df/dr
    d = fp * r / R + f * (1.0 / R - r * r / R**3)  # d/dr (f r / R)
    return f * ell * (ell + 1) / R2 + (1 - s) * f * d / R


def peak(a, ell, s):
    """(r_peak, R_peak, V_peak) of one barrier (r > 0), found numerically."""
    res = minimize_scalar(
        lambda r: -potential_r(r, a, ell, s),
        bounds=(1e-6, 40.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    r = float(res.x)
    return r, math.hypot(r, a), float(-res.fun)


class Grid:
    """Uniform r* grid on [0, xmax] for one SV barrier, r(r*) from dr/dr* = f."""

    def __init__(self, a, ell, s, xmax=2500.0, h=0.01):
        self.a, self.ell, self.s, self.h = a, ell, s, h
        x = np.arange(0.0, xmax + h / 2, h)
        sol = solve_ivp(
            lambda t, y: [f_of_r(y[0], a)], (0.0, xmax), [0.0], t_eval=x, rtol=1e-12, atol=1e-12
        )
        self.x, self.r = x, sol.y[0]
        self.V = potential_r(self.r, a, ell, s)
        self.v0 = float(self.V[0])

    def V_at(self, xq):
        return np.interp(np.abs(xq), self.x, self.V)


# ----------------------------------------------------------------------- route A: Numerov
def _numerov_inward(V, h, w, psi_end):
    """Integrate psi'' = (V - w^2) psi from the last grid point to the first, vectorised in w.

    psi_end: (psi[N-1], psi[N-2]) arrays over w. Returns psi[0], psi[1].
    """
    w2 = (w * w)[:, None]
    n = V.size
    p1, p0 = psi_end  # psi at N-1 and N-2
    g = lambda i: h * h / 12.0 * (V[i] - w2[:, 0])  # noqa: E731
    gi1, gi0 = g(n - 1), g(n - 2)
    for i in range(n - 3, -1, -1):
        gm = g(i)
        p_new = (2 * p0 * (1 + 5 * gi0) - p1 * (1 - gi1)) / (1 - gm)
        p1, p0 = p0, p_new
        gi1, gi0 = gi0, gm
    return p0, p1  # psi[0], psi[1]


def single_barrier(grid: Grid, w):
    """|R|^2, |T|^2 (route A) for waves from the throat side, w an array."""
    w = np.atleast_1d(np.asarray(w, float))
    if np.any(w * w <= max(grid.v0, grid.V[-1]) * 1.0001):
        raise ValueError("w^2 must exceed V at the throat and at the outer grid edge")
    xe = grid.x[-1]
    k_out = np.sqrt(w * w - grid.V[-1])
    p_end = (np.exp(1j * k_out * xe), np.exp(1j * k_out * (xe - grid.h)))
    psi0, psi1 = _numerov_inward(grid.V, grid.h, w, p_end)
    k0 = np.sqrt(w * w - grid.v0)
    # psi = A e^{i k0 x} + B e^{-i k0 x} at x = 0 and x = h
    e1 = np.exp(1j * k0 * grid.h)
    A = (psi1 - psi0 / e1) / (e1 - 1 / e1)
    B = psi0 - A
    T2 = (k_out / k0) / np.abs(A) ** 2
    R2 = np.abs(B) ** 2 / np.abs(A) ** 2
    return R2, T2


def double_barrier_T2(grid: Grid, w):
    """|T|^2 through the full symmetric double barrier (route A), w an array."""
    w = np.atleast_1d(np.asarray(w, float))
    Vfull = np.concatenate([grid.V[:0:-1], grid.V])  # x from -xmax to +xmax
    xe = grid.x[-1]
    k = np.sqrt(w * w - grid.V[-1])
    p_end = (np.exp(1j * k * xe), np.exp(1j * k * (xe - grid.h)))
    psi0, psi1 = _numerov_inward(Vfull, grid.h, w, p_end)
    x0 = -xe
    e0, e1 = np.exp(1j * k * x0), np.exp(1j * k * (x0 + grid.h))
    # solve psi(x0) = A e0 + B/e0, psi(x0+h) = A e1 + B/e1
    det = e0 / e1 - e1 / e0
    A = (psi0 / e1 - psi1 / e0) / det
    return 1.0 / np.abs(A) ** 2


# ------------------------------------------------------------------------ route B: adaptive ODE
def single_barrier_ode(grid: Grid, w):
    """|R|^2, |T|^2 at one w with DOP853 (independent of the Numerov stepping)."""
    xe = grid.x[-1]
    k_out = math.sqrt(w * w - grid.V[-1])
    y0 = [np.exp(1j * k_out * xe), 1j * k_out * np.exp(1j * k_out * xe)]
    sol = solve_ivp(
        lambda x, y: [y[1], (grid.V_at(x) - w * w) * y[0]],
        (xe, 0.0),
        y0,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    psi, dpsi = sol.y[0, -1], sol.y[1, -1]
    k0 = math.sqrt(w * w - grid.v0)
    A = 0.5 * (psi + dpsi / (1j * k0))
    B = 0.5 * (psi - dpsi / (1j * k0))
    return abs(B) ** 2 / abs(A) ** 2, (k_out / k0) / abs(A) ** 2


# -------------------------------------------------------------------------- check helpers
def schwarzschild_scalar_T2(w, ell=0, xmin=-200.0, xmax=1500.0, h=0.01):
    """|T|^2 for a scalar wave on Schwarzschild (Numerov), for check (2)."""
    w = np.atleast_1d(np.asarray(w, float))
    x = np.arange(xmin, xmax + h / 2, h)
    # r = 2 (1 + W(exp(y))), y = x/2 - 1; for large y, Newton on w + ln w = y
    y = x / 2.0 - 1.0
    big = y >= 300
    small_w = np.real(lambertw(np.exp(np.minimum(y, 300.0))))
    wv = np.where(big, y - np.log(np.maximum(y, 1)), small_w)
    for _ in range(50):
        wv = wv - (wv + np.log(wv) - y) / (1 + 1 / wv) * big
    r = 2.0 * (1.0 + wv)
    V = (1 - 2 / r) * (ell * (ell + 1) / r**2 + 2 / r**3)
    k = w
    p_end = (np.exp(1j * k * x[-1]), np.exp(1j * k * x[-2]))
    psi0, psi1 = _numerov_inward(V, h, w, p_end)
    e0, e1 = np.exp(1j * k * x[0]), np.exp(1j * k * x[1])
    det = e0 / e1 - e1 / e0
    A = (psi0 / e1 - psi1 / e0) / det
    return 1.0 / np.abs(A) ** 2


def wkb_T2(w, vpeak, v2):
    """First-order WKB (Schutz-Will / Iyer-Will): 1 / (1 + exp(2 pi (V0 - w^2) / sqrt(-2 V0'')))."""
    return 1.0 / (1.0 + np.exp(2 * np.pi * (vpeak - w * w) / np.sqrt(-2.0 * v2)))


def vpp_at_peak(grid: Grid, r_peak):
    i = int(np.argmin(np.abs(grid.r - r_peak)))
    hh = grid.h * 20
    xs = grid.x[i]
    return float((grid.V_at(xs + hh) - 2 * grid.V_at(xs) + grid.V_at(xs - hh)) / hh**2)


def turning_round_trip(grid: Grid, w):
    """4 x_t, x_t = r* where V first reaches w^2 going out from the throat (classical cavity)."""
    above = grid.V >= w * w
    if not above.any():
        raise ValueError("w^2 is above the barrier: no turning point")
    return 4.0 * grid.x[int(np.argmax(above))]


# --------------------------------------------------------------------------------------- checks
def main() -> int:
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    import x2_echo_delay as xe

    out = {"label": "derived; test fields s = 0, 1, l = 2; conditional geometry (X2 premises)"}
    peaks = {}
    rows = []
    for s in (0, 1):
        for a in (2.001, 2.01):
            g = Grid(a, 2, s)
            rp, Rp, vp = peak(a, 2, s)
            peaks[f"s{s}_a{a}"] = {"r_peak": rp, "R_peak": Rp, "V_peak": vp, "V_throat": g.v0}
            v2 = vpp_at_peak(g, rp)
            fr = np.array([0.6, 0.9, 1.0, 1.1, 1.4])
            w = np.sqrt(fr * vp)
            R2, T2 = single_barrier(g, w)
            for j, wj in enumerate(w):
                rb, tb = single_barrier_ode(g, float(wj)) if j in (0, 2, 4) else (None, None)
                rows.append(
                    {
                        "s": s,
                        "a": a,
                        "w2_over_Vpeak": float(fr[j]),
                        "R2": float(R2[j]),
                        "T2": float(T2[j]),
                        "flux_err": float(abs(R2[j] + T2[j] - 1)),
                        "T2_routeB": tb,
                        "wkb_T2": float(wkb_T2(wj, vp, v2)),
                    }
                )
    out["peaks_l2"] = peaks
    out["single_barrier_l2"] = rows
    ws = np.array([0.002, 0.004, 0.008])
    ratio = np.pi * schwarzschild_scalar_T2(ws) / ws**2 / (16 * np.pi)
    slope = (ratio[1] - ratio[0]) / (ws[1] - ws[0])
    out["schwarzschild_l0_sigma_over_16pi"] = {
        "w": ws.tolist(),
        "ratio": ratio.tolist(),
        "extrapolated_w0": float(ratio[0] - slope * ws[0]),
    }
    # resonance spacing of the double barrier (s = 1, l = 2) at frequencies below the barrier top
    res = []
    for a in (2.001, 2.0001):
        g = Grid(a, 2, 1, xmax=2500.0, h=0.02)
        _, _, vp = peak(a, 2, 1)
        for lo, hi in ((0.45, 0.50), (0.70, 0.74)):
            wgrid = np.linspace(lo * math.sqrt(vp), hi * math.sqrt(vp), 4000)
            t2 = double_barrier_T2(g, wgrid)
            pk = [i for i in range(1, len(t2) - 1) if t2[i] > t2[i - 1] and t2[i] > t2[i + 1]]
            wp = wgrid[pk]
            dw = np.diff(wp)
            wm = 0.5 * (wp[1:] + wp[:-1])
            res.append(
                {
                    "a": a,
                    "w_over_sqrtVpeak": [lo, hi],
                    "n_resonances": len(wp),
                    "time_2pi_over_dw": float(np.mean(2 * math.pi / dw)) if dw.size else None,
                    "peak_to_peak_4L": xe.echo_delay(a),
                    "turning_point_round_trip": float(
                        np.mean([turning_round_trip(g, x) for x in wm])
                    )
                    if dw.size
                    else None,
                }
            )
    out["resonance_spacing_s1_l2"] = res
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s1_checks.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
