"""Two-way signalling through a collapsing Ellis-Bronnikov wormhole with self-gravitating pulses.

Model (Kain, arXiv:2602.18231v1 = CQG 43, 045014, sections 2-3; same model as Shinkai & Hayward,
gr-qc/0205041, and Calhoun, Fay & Kain, arXiv:2210.04905): 4D, spherical symmetry, the symmetric
massless Ellis-Bronnikov wormhole (throat radius x0 = 1, lengths and times in units of x0,
G = c = 1) supported by a massless ghost scalar phi (eta = -1). ASSUMPTION: the ghost field is a
premise of the model, not a known or manufacturable material. Two ordinary massless scalars carry
the signals: chi1 (outbound, placed on u = 0, travels from the + side to the - side) and chi2
(return, placed on v = 0, travels from the - side to the + side). Both are sources in the Einstein
equations, so the geometry includes the back-reaction of both legs.

Metric ds^2 = -exp(2 sigma) du dv + r^2 dOmega^2 (Kain eq. 1).
Evolution (Kain eqs. 4, 51, with T_uv = 0):
  r_uv     = -(4 r_u r_v + e^{2 sigma}) / (4 r)
  sigma_uv = (4 r_u r_v + e^{2 sigma}) / (4 r^2) - 4 pi (eta phi_u phi_v + sum chi_u chi_v)
  f_uv     = -(r_u f_v + r_v f_u) / r                for f in (phi, chi1, chi2)
Constraints (Kain eq. 3; stress tensor Kain eqs. 53-54), used for initial data and as monitors:
  r_uu = 2 sigma_u r_u - 4 pi r T_uu,  r_vv = 2 sigma_v r_v - 4 pi r T_vv.
Scheme: Kain's second-order predictor-corrector on each (u, v) cell (section 2.3), marched one
anti-diagonal (u + v = const, a spacelike slice) at a time so that each step is one numpy vector op.
Initial data: the static solution on u = 0 and v = 0 with sigma kept static and r re-integrated from
the null constraint where a pulse is present (Kain's "first option", section 3.3); the pulse is
Kain's sin^2 profile in chi_v (eqs. 57-58) or chi_u. Equation numbers follow the arXiv v1 source.

Energy bookkeeping: the Misner-Sharp mass m = (r/2)(1 + 4 e^{-2 sigma} r_u r_v) obeys
  m_u = -8 pi r^2 e^{-2 sigma} r_v T_uu,   m_v = -8 pi r^2 e^{-2 sigma} r_u T_vv   (T_uv = 0),
so along any curve dm splits exactly into per-field, per-direction pieces. The energy a field
carries through a station sphere r = R is the integral of its piece along the station's worldline
(Kodama energy; in the static era the Kodama vector is r'(y) d/dt, so at r = 2 the Kodama energy is
0.866 of the Killing energy, which is what a static observer measures there since sigma = 0).
Stations are constant-r timelike worldlines; their proper time is tau = int e^sigma sqrt(du dv).

  python scripts/wh_twoway.py checks        # static, Kain Fig. 9 and noise checks -> checks.json
  python scripts/wh_twoway.py scan --n 800   # frozen problem: E_ret x u1 grid -> scan_n800.json
  python scripts/wh_twoway.py analyze       # boundaries, Richardson, case table -> boundary.json
  python scripts/wh_twoway.py imprint | figures   (results/wh_twoway/RESULT.md lists every command)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "wh_twoway"
FOUR_PI = 4.0 * math.pi
ETA = -1.0  # ghost sign (Kain eqs. 32, 49)
FIELDS = ("r", "s", "p", "c1", "c2")  # r, sigma, phi (ghost), chi1 (outbound), chi2 (return)


@dataclass(frozen=True)
class Params:
    n: int = 100  # grid points per unit length (du = dv = 1/n)
    u_max: float = 20.0
    v_max: float = 20.0
    a1: float = 0.002  # outbound pulse amplitude (Kain's A), chi1_v = a1 sin^2 on (v1, v2) at u = 0
    v1: float = 1.0
    v2: float = 3.0
    a2: float = 0.0  # return pulse amplitude, chi2_u = a2 sin^2 on (u1, u2) at v = 0
    u1: float = 10.0
    u2: float = 11.0
    r_station: float = 2.0  # areal radius of both stations (A on the + side, B on the - side)
    out_h: float = 0.05  # spacing of the stored output grid (must be a multiple of 1/n)
    r_excise: float = 0.05  # points with r below this (inside the black hole) are dropped


def pulse(x: np.ndarray, amp: float, x1: float, x2: float) -> tuple[np.ndarray, np.ndarray]:
    """Kain eqs. 57-58: derivative amp sin^2(pi (x - x1)/(x2 - x1)) on (x1, x2) and its integral."""
    w = x2 - x1
    s = np.clip((x - x1) / w, 0.0, 1.0)
    deriv = amp * np.sin(math.pi * s) ** 2
    val = amp / (4.0 * math.pi) * (2.0 * math.pi * w * s - w * np.sin(2.0 * math.pi * s))
    return deriv, val


def _surface(x: np.ndarray, amp: float, x1: float, x2: float, sign: float) -> dict:
    """Initial data on one null surface (x = v on u = 0 with sign +1, x = u on v = 0 with sign -1).

    r'' = r [4/(4 + x^2)^2 - 4 pi chi'^2] (null constraint with sigma' = 0, eta = -1), r(0) = 1,
    r'(0) = 0; the transverse derivative q obeys q' = -(4 q r' + 1)/(4 r) (the r_uv equation),
    q(0) = 0. RK4 on the grid itself (fourth-order, so it does not limit the second-order scheme).
    """
    h = x[1] - x[0]

    def rhs(xx, y):
        d, _ = pulse(np.asarray(xx), amp, x1, x2)
        r, rp, q = y
        return np.array(
            [
                rp,
                r * (4.0 / (4.0 + xx * xx) ** 2 - FOUR_PI * float(d) ** 2),
                -(4 * q * rp + 1) / (4 * r),
            ]
        )

    y = np.array([1.0, 0.0, 0.0])
    out = np.empty((len(x), 3))
    out[0] = y
    for k in range(len(x) - 1):
        xx = x[k]
        k1 = rhs(xx, y)
        k2 = rhs(xx + h / 2, y + h / 2 * k1)
        k3 = rhs(xx + h / 2, y + h / 2 * k2)
        k4 = rhs(xx + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[k + 1] = y
    d, val = pulse(x, amp, x1, x2)
    r, rp, q = out[:, 0], out[:, 1], out[:, 2]
    # chi-part of dm along the surface: -8 pi r^2 q chi'^2 (sigma = 0); the pulse's Misner-Sharp
    # energy
    dens = -2.0 * FOUR_PI * r * r * q * d * d
    energy = float(np.sum(0.5 * (dens[1:] + dens[:-1])) * h)
    m = 0.5 * r * (1.0 + 4.0 * q * rp)
    return {
        "r": r,
        "s": np.zeros_like(x),
        "p": sign * np.arctan(x / 2.0) / math.sqrt(FOUR_PI),
        "chi": val,
        "energy": energy,
        "m_end": float(m[-1]),
    }


def rhs_uv(r, s, ru, rv, pu, pv, c1u, c1v, c2u, c2v):
    x = 4.0 * ru * rv + np.exp(2.0 * s)
    r_uv = -x / (4.0 * r)
    s_uv = x / (4.0 * r * r) - FOUR_PI * (ETA * pu * pv + c1u * c1v + c2u * c2v)
    p_uv = -(ru * pv + rv * pu) / r
    c1_uv = -(ru * c1v + rv * c1u) / r
    c2_uv = -(ru * c2v + rv * c2u) / r
    return r_uv, s_uv, p_uv, c1_uv, c2_uv


def _cell(f0: dict, base: dict, du_: dict, dv_: dict, h2: float) -> dict:
    """Kain's cell update f4 = f2 + f3 - f1 + F(f0, f_u, f_v) du dv for every field."""
    src = rhs_uv(
        f0["r"], f0["s"], du_["r"], dv_["r"], du_["p"], dv_["p"],
        du_["c1"], dv_["c1"], du_["c2"], dv_["c2"],
    )  # fmt: skip
    return {f: base[f] + src[n] * h2 for n, f in enumerate(FIELDS)}


class Station:
    """Accumulates the energy fluxes and proper time along a constant-r worldline (r = R)."""

    def __init__(self, side: int, r_station: float):
        self.side = side  # +1: + side (A), -1: - side (B)
        self.r_station = r_station
        self.rows: list[tuple] = []

    def add(self, u, v, s, r, ru, rv, tuu, tvv, m):
        self.rows.append((u, v, s, r, ru, rv, *tuu, *tvv, m))

    def table(self) -> dict:
        if len(self.rows) < 2:
            return {}
        a = np.array(self.rows)
        u, v, s, r, ru, rv = (a[:, i] for i in range(6))
        tuu = a[:, 6:9]  # phi, chi1, chi2
        tvv = a[:, 9:12]
        m = a[:, 12]
        du, dv = np.diff(u), np.diff(v)
        mid = lambda z: 0.5 * (z[1:] + z[:-1])  # noqa: E731
        es = np.exp(s)
        dtau = mid(es) * np.sqrt(np.clip(du * dv, 0.0, None))
        # clock zero: the event where the worldline meets the initial null surface. The first stored
        # row lies
        # slightly later, still in the static era, where tau = t - y_station exactly (sigma = 0, t =
        # (u + v)/2).
        tau0 = 0.5 * (u[0] + v[0]) - math.sqrt(self.r_station**2 - 1.0)
        tau = tau0 + np.concatenate([[0.0], np.cumsum(dtau)])
        pref = 2.0 * FOUR_PI * r * r * np.exp(-2.0 * s)  # 8 pi r^2 e^{-2 sigma}
        # dm pieces: from T_uu (u-moving waves) and T_vv (v-moving waves), per field
        dmu = -mid((pref * rv)[:, None] * tuu) * du[:, None]
        dmv = -mid((pref * ru)[:, None] * tvv) * dv[:, None]
        return {"u": u, "v": v, "tau": tau, "r": r, "m": m, "dm_u": dmu, "dm_v": dmv}


def evolve(p: Params, keep_grid: bool = True) -> dict:
    t0 = time.time()
    h = 1.0 / p.n
    nu = int(round(p.u_max * p.n))
    nv = int(round(p.v_max * p.n))
    step = int(round(p.out_h * p.n))
    assert abs(step * h - p.out_h) < 1e-12 and step >= 1
    uu = np.arange(nu + 1) * h
    vv = np.arange(nv + 1) * h
    su = _surface(vv, p.a1, p.v1, p.v2, +1.0)  # u = 0 surface, parameter v, outbound chi1
    sv = _surface(uu, p.a2, p.u1, p.u2, -1.0)  # v = 0 surface, parameter u, return chi2
    init_u0 = {"r": su["r"], "s": su["s"], "p": su["p"], "c1": su["chi"], "c2": np.zeros(nv + 1)}
    init_v0 = {"r": sv["r"], "s": sv["s"], "p": sv["p"], "c1": np.zeros(nu + 1), "c2": sv["chi"]}

    # three rotating diagonals, each indexed by the u-index i (length nu + 1)
    diag = [{f: np.full(nu + 1, np.nan) for f in FIELDS} for _ in range(3)]
    no, mo = nu // step + 1, nv // step + 1
    grid = {f: np.full((no, mo), np.nan) for f in FIELDS} if keep_grid else None
    sta = {+1: Station(+1, p.r_station), -1: Station(-1, p.r_station)}
    # null-constraint residual |C| and the size of its terms, maximised over out_h x out_h blocks
    # (points with
    # r > 0.5); summarize() reads them over each station's causal past
    cmap = {c: np.zeros((no, mo)) for c in ("u", "v", "su", "sv")}

    def set_boundary(d, k):
        if k <= nv:
            for f in FIELDS:
                d[f][0] = init_u0[f][k]
        if k <= nu:
            for f in FIELDS:
                d[f][k] = init_v0[f][k]

    def store(d, k):
        if grid is None or k % step:
            return
        ilo, ihi = max(0, k - nv), min(k, nu)
        i0 = ilo + (-ilo) % step
        idx = np.arange(i0, ihi + 1, step)
        if len(idx) == 0:
            return
        for f in FIELDS:
            grid[f][idx // step, (k - idx) // step] = d[f][idx]

    d0 = diag[0]
    set_boundary(d0, 0)
    store(d0, 0)
    d1 = diag[1]
    set_boundary(d1, 1)
    store(d1, 1)
    h2 = h * h
    with np.errstate(all="ignore"):
        for k in range(2, nu + nv + 1):
            dm2, dm1, dk = diag[(k - 2) % 3], diag[(k - 1) % 3], diag[k % 3]
            for f in FIELDS:
                dk[f].fill(np.nan)
            set_boundary(dk, k)
            ia, ib = max(1, k - nv), min(k - 1, nu)
            if ia <= ib:
                sl1 = slice(ia - 1, ib)  # points 1 (diag k-2) and 2 (diag k-1): index i-1
                sl3 = slice(ia, ib + 1)  # point 3 (diag k-1): index i
                f1 = {f: dm2[f][sl1] for f in FIELDS}
                f2 = {f: dm1[f][sl1] for f in FIELDS}
                f3 = {f: dm1[f][sl3] for f in FIELDS}
                f0 = {f: 0.5 * (f2[f] + f3[f]) for f in FIELDS}
                base = {f: f2[f] + f3[f] - f1[f] for f in FIELDS}
                du_ = {f: (f3[f] - f1[f]) / h for f in FIELDS}
                dv_ = {f: (f2[f] - f1[f]) / h for f in FIELDS}
                f4 = _cell(f0, base, du_, dv_, h2)  # predictor (one-sided derivatives)
                du_ = {f: (f3[f] - f1[f] + f4[f] - f2[f]) / (2 * h) for f in FIELDS}
                dv_ = {f: (f2[f] - f1[f] + f4[f] - f3[f]) / (2 * h) for f in FIELDS}
                f4 = _cell(f0, base, du_, dv_, h2)  # corrector (centred derivatives)
                bad = ~(f4["r"] > p.r_excise) | ~(np.abs(f4["s"]) < 50.0)
                for f in FIELDS:
                    v_ = f4[f]
                    v_[bad] = np.nan
                    dk[f][sl3] = v_
            store(dk, k)
            # --- diagnostics on cell centres of the cells completed by diagonal k (slice u + v =
            # (k-1) h)
            ca, cb = max(1, k - nv), min(k - 1, nu)
            if ca > cb:
                continue
            i = np.arange(ca, cb + 1)
            g1 = {f: dm2[f][i - 1] for f in FIELDS}
            g2 = {f: dm1[f][i - 1] for f in FIELDS}
            g3 = {f: dm1[f][i] for f in FIELDS}
            g4 = {f: dk[f][i] for f in FIELDS}
            cen = {f: 0.25 * (g1[f] + g2[f] + g3[f] + g4[f]) for f in FIELDS}
            cu = {f: (g3[f] - g1[f] + g4[f] - g2[f]) / (2 * h) for f in FIELDS}
            cv = {f: (g2[f] - g1[f] + g4[f] - g3[f]) / (2 * h) for f in FIELDS}
            rc = cen["r"]
            for side, st in sta.items():
                order = np.arange(len(i)) if side > 0 else np.arange(len(i))[::-1]
                rr = rc[order]
                if not (rr[0] > p.r_station):
                    continue  # outer end of this slice inside the station radius (domain edge)
                below = np.nonzero(~(rr > p.r_station))[0]
                if len(below) == 0:
                    continue
                b = below[0]
                ja, jb = order[b - 1], order[b]
                if not (np.isfinite(rc[ja]) and np.isfinite(rc[jb])):
                    continue
                w = (rc[ja] - p.r_station) / (rc[ja] - rc[jb])

                def lerp(z, ja=ja, jb=jb, w=w):
                    return (1 - w) * z[ja] + w * z[jb]

                u_c = (lerp(i) - 0.5) * h
                v_c = (k - 1) * h - u_c
                tuu = [ETA * lerp(cu["p"]) ** 2, lerp(cu["c1"]) ** 2, lerp(cu["c2"]) ** 2]
                tvv = [ETA * lerp(cv["p"]) ** 2, lerp(cv["c1"]) ** 2, lerp(cv["c2"]) ** 2]
                ru_, rv_, s_ = lerp(cu["r"]), lerp(cv["r"]), lerp(cen["s"])
                m_ = 0.5 * p.r_station * (1 + 4 * np.exp(-2 * s_) * ru_ * rv_)
                st.add(u_c, v_c, s_, p.r_station, ru_, rv_, tuu, tvv, m_)
            # --- null-constraint residuals at the points of diagonal k-1 (needs k-2, k-1, k)
            # C_v at (i, j-1): second difference along v uses index i on diagonals k-2, k-1, k.
            ii = np.arange(max(0, k - nv), min(k - 2, nu) + 1)
            if len(ii):
                rr = dm1["r"][ii]
                r_vv = (dk["r"][ii] - 2 * rr + dm2["r"][ii]) / h2
                r_v = (dk["r"][ii] - dm2["r"][ii]) / (2 * h)
                s_v = (dk["s"][ii] - dm2["s"][ii]) / (2 * h)
                tvv = sum(
                    eta * ((dk[f][ii] - dm2[f][ii]) / (2 * h)) ** 2
                    for f, eta in (("p", ETA), ("c1", 1.0), ("c2", 1.0))
                )
                c = r_vv - 2 * s_v * r_v + FOUR_PI * rr * tvv
                scale = np.abs(r_vv) + np.abs(2 * s_v * r_v) + np.abs(FOUR_PI * rr * tvv)
                ok = np.isfinite(c) & (rr > 0.5)
                if ok.any():
                    bi, bj = ii[ok] // step, np.minimum((k - 1 - ii[ok]) // step, mo - 1)
                    np.maximum.at(cmap["v"], (bi, bj), np.abs(c[ok]))
                    np.maximum.at(cmap["sv"], (bi, bj), scale[ok])
            # C_u at (i-1, j): second difference along u uses i-2, i-1, i on diagonals k-2, k-1, k.
            ii = np.arange(max(2, k - nv), min(k, nu) + 1)
            if len(ii):
                rr = dm1["r"][ii - 1]
                r_uu = (dk["r"][ii] - 2 * rr + dm2["r"][ii - 2]) / h2
                r_u = (dk["r"][ii] - dm2["r"][ii - 2]) / (2 * h)
                s_u = (dk["s"][ii] - dm2["s"][ii - 2]) / (2 * h)
                tuu = sum(
                    eta * ((dk[f][ii] - dm2[f][ii - 2]) / (2 * h)) ** 2
                    for f, eta in (("p", ETA), ("c1", 1.0), ("c2", 1.0))
                )
                c = r_uu - 2 * s_u * r_u + FOUR_PI * rr * tuu
                scale = np.abs(r_uu) + np.abs(2 * s_u * r_u) + np.abs(FOUR_PI * rr * tuu)
                ok = np.isfinite(c) & (rr > 0.5)
                if ok.any():
                    bi, bj = (ii[ok] - 1) // step, np.minimum((k - ii[ok]) // step, mo - 1)
                    np.maximum.at(cmap["u"], (bi, bj), np.abs(c[ok]))
                    np.maximum.at(cmap["su"], (bi, bj), scale[ok])
    return {
        "params": asdict(p),
        "grid": grid,
        "grid_h": p.out_h,
        "stations": {"A": sta[+1].table(), "B": sta[-1].table()},
        "cons": cmap,
        "id_energy": {"chi1": su["energy"], "chi2": sv["energy"]},
        "seconds": time.time() - t0,
    }


# -----------------------------------------------------------------------------------------------
# analysis
FIELD_NAMES = ("ghost", "chi1", "chi2")


def amplitude_for_energy(energy: float, x1: float, x2: float, n: int, x_max: float) -> float:
    """Pulse amplitude whose initial-surface Misner-Sharp energy equals `energy` (secant steps)."""
    if energy <= 0:
        return 0.0
    x = np.arange(int(round(x_max * n)) + 1) / n
    a = 1e-3
    for _ in range(3):
        e = _surface(x, a, x1, x2, 1.0)["energy"]
        a *= math.sqrt(energy / e)
    return a


def crossing_time(tau_mid: np.ndarray, flux: np.ndarray, frac: float) -> float | None:
    """Station proper time at which the cumulative flux reaches `frac` of its total."""
    c = np.cumsum(flux)
    if c[-1] <= 0:
        return None
    c = c / c[-1]
    j = int(np.searchsorted(c, frac))
    if j >= len(c):
        return None
    if j == 0:
        return float(tau_mid[0])
    w = (frac - c[j - 1]) / (c[j] - c[j - 1])
    return float(tau_mid[j - 1] + w * (tau_mid[j] - tau_mid[j - 1]))


def station_summary(tab: dict, side: str, r_station: float, tau_dead: float | None = None) -> dict:
    """Per-field inflow (toward the throat) and outflow energies through one station, with times.

    + side (A): v-moving waves (T_vv) flow in, u-moving waves (T_uu) flow out; - side (B): reverse.
    Kodama energies (dm pieces) and the energies a station observer measures locally,
    dE_loc = dE_K / sqrt(1 - 2 m / R) (the station 4-velocity is the normalised Kodama vector).
    Only flux with station proper time <= tau_dead counts (receiver deadline), if given.
    """
    if not tab:
        return {}
    tau = tab["tau"]
    tau_mid = 0.5 * (tau[1:] + tau[:-1])
    m_mid = 0.5 * (tab["m"][1:] + tab["m"][:-1])
    lapse = np.sqrt(np.clip(1.0 - 2.0 * m_mid / r_station, 1e-12, None))
    keep = np.ones_like(tau_mid, bool) if tau_dead is None else tau_mid <= tau_dead
    if side == "A":
        inflow, outflow = tab["dm_v"], -tab["dm_u"]
    else:
        inflow, outflow = tab["dm_u"], -tab["dm_v"]
    lap_dead = lapse if tau_dead is None else lapse[tau_mid <= tau_dead]
    out = {
        "lapse_min_dead": float(np.min(lap_dead)) if len(lap_dead) else None,
        "tau_end": float(tau[-1]),
        "u_end": float(tab["u"][-1]),
        "v_end": float(tab["v"][-1]),
        "m_min": float(np.min(tab["m"])),
        "m_max": float(np.max(tab["m"])),
        "lapse_min": float(np.min(lapse)),
    }
    for k, name in enumerate(FIELD_NAMES):
        for label, flux in (("in", inflow[:, k]), ("out", outflow[:, k])):
            fk = np.where(keep, flux, 0.0)
            rec = {
                "E_kodama": float(np.sum(fk)),
                "E_local": float(np.sum(fk / lapse)),
                "tau10": crossing_time(tau_mid, np.clip(fk, 0, None), 0.1),
                "tau50": crossing_time(tau_mid, np.clip(fk, 0, None), 0.5),
                "tau90": crossing_time(tau_mid, np.clip(fk, 0, None), 0.9),
            }
            if rec["tau50"] is not None:
                rec["u50"] = float(np.interp(rec["tau50"], tau, tab["u"]))
                rec["v50"] = float(np.interp(rec["tau50"], tau, tab["v"]))
            out[f"{name}_{label}"] = rec
    return out


def passes(e_emit: float, e_recv: float, frac: float = 0.5) -> bool:
    """Frozen success rule: the receiver gets at least `frac` of the emitted Kodama energy."""
    return bool(e_emit > 0 and e_recv >= frac * e_emit)


LAPSE_MIN_VALID = (
    0.5  # stations count as untrapped, timelike observers only above this sqrt(1 - 2m/R)
)


def run_valid(s: dict, tau_dead: float = None) -> bool:
    """A run may be scored only if both stations stay well untrapped (timelike tracks, no clipping),
    A's track reaches the deadline and B's track covers the whole return pulse (review item 3)."""
    td = TAU_DEAD if tau_dead is None else tau_dead
    p = s["params"]
    return bool(
        s["A"].get("tau_end", 0.0) >= td
        and (s["A"].get("lapse_min_dead") or s["A"].get("lapse_min", 0.0)) > LAPSE_MIN_VALID
        and s["B"].get("lapse_min", 0.0) > LAPSE_MIN_VALID
        and s["B"].get("u_end", 0.0) >= p["u2"] + 0.5
    )


def summarize(res: dict, tau_dead: float) -> dict:
    p = res["params"]
    sa = station_summary(res["stations"]["A"], "A", p["r_station"], tau_dead)
    sb = station_summary(res["stations"]["B"], "B", p["r_station"])
    ta, tb = res["stations"]["A"], res["stations"]["B"]
    out = {"params": p, "seconds": res["seconds"], "id_energy": res["id_energy"], "A": sa, "B": sb}
    # radial null rays are the coordinate lines: a ray leaving B inward at u_e reaches A iff A's
    # track
    # reaches u = u_e before the deadline; map u_e to B's and A's proper times
    if ta and tb:
        keep = ta["tau"] <= tau_dead
        out["test_ray"] = {
            "u_A_max": float(np.max(ta["u"][keep])),
            "tauB_at_u_A_max": float(np.interp(np.max(ta["u"][keep]), tb["u"], tb["tau"])),
        }
    e_emit = sb.get("chi2_in", {}).get("E_kodama", 0.0)
    e_recv = sa.get("chi2_out", {}).get("E_kodama", 0.0)
    if p["a2"] > 0:
        out["return"] = {
            "E_emit_B": e_emit,
            "E_recv_A": e_recv,
            "fraction": e_recv / e_emit if e_emit > 0 else None,
            "E_emit_B_local": sb["chi2_in"]["E_local"],
            "E_recv_A_local": sa["chi2_out"]["E_local"],
            "E_reflected_B": sb["chi2_out"]["E_kodama"],
            "tauB_emit50": sb["chi2_in"]["tau50"],
            "uB_emit50": sb["chi2_in"].get("u50"),
            "tauA_recv50": sa["chi2_out"]["tau50"],
            "pass": passes(e_emit, e_recv),
        }
    if p["a1"] > 0:
        out["outbound"] = {
            "E_emit_A": sa["chi1_in"]["E_kodama"],
            "E_recv_B": sb["chi1_out"]["E_kodama"],
            "tauA_emit50": sa["chi1_in"]["tau50"],
            "tauB_recv50": sb["chi1_out"]["tau50"],
            "tauB_recv90": sb["chi1_out"]["tau90"],
        }
    if p["a1"] > 0 and p["a2"] > 0:
        r_ = out["return"]
        o_ = out["outbound"]
        r_["reply_compatible"] = bool(
            r_["tauB_emit50"] is not None
            and o_["tauB_recv90"] is not None
            and sb["chi2_in"]["tau10"] is not None
            and sb["chi2_in"]["tau10"] >= o_["tauB_recv90"]
        )
        if r_["tauA_recv50"] is not None and o_["tauA_emit50"] is not None:
            r_["round_trip_tauA"] = r_["tauA_recv50"] - o_["tauA_emit50"]
        for td in TAU_DEAD_ALT:
            alt = station_summary(ta, "A", p["r_station"], td)["chi2_out"]["E_kodama"]
            out["return"][f"fraction_dead_{td:g}"] = alt / e_emit if e_emit > 0 else None
    # Misner-Sharp bookkeeping: the per-field pieces must add up to the change of m along each track
    out["bookkeeping"] = {}
    for name, tab in (("A", ta), ("B", tb)):
        if tab:
            cum = np.concatenate([[0.0], np.cumsum(tab["dm_u"].sum(1) + tab["dm_v"].sum(1))])
            err = np.abs(cum - (tab["m"] - tab["m"][0]))
            keep = tab["tau"] <= tau_dead
            out["bookkeeping"][name] = float(np.max(err[keep]))
    out["constraint"] = constraint_in_past(res, tau_dead)
    out["constraint_tauA_11"] = constraint_in_past(res, 11.0)
    return out


def constraint_in_past(res: dict, tau_dead: float) -> dict:
    """Largest null-constraint residual (and largest term) in each station's causal past.

    A's past up to the deadline is the rectangle u <= u_A(tau_dead), v <= v_A(tau_dead) (the track
    is future-directed); B's is taken up to the end of its track. Blocks are out_h wide, so the
    rectangle is
    rounded outward by one block.
    """
    cm, h = res["cons"], res["grid_h"]
    out = {}
    for name, tab, td in (("A", res["stations"]["A"], tau_dead), ("B", res["stations"]["B"], None)):
        if not tab:
            continue
        keep = np.ones_like(tab["tau"], bool) if td is None else tab["tau"] <= td
        iu = int(np.max(tab["u"][keep]) // h) + 1
        iv = int(np.max(tab["v"][keep]) // h) + 1
        rec = {}
        for c, sc in (("u", "su"), ("v", "sv")):
            blk, sblk = cm[c][: iu + 1, : iv + 1], cm[sc][: iu + 1, : iv + 1]
            rec[c] = {"max": float(np.max(blk)), "max_term": float(np.max(sblk))}
        out[name] = rec
    return out


# ------------------------------------------------------------------------------------- frozen
# problem
# Fixed 2026-10-28 before any return-pulse run (docs/results/wh_twoway.md, "Fixed problem").
R_STATION = 2.0  # both stations, areal radius in units of the throat radius x0
Y_STATION = math.sqrt(R_STATION**2 - 1.0)  # static-era proper distance of a station from the throat
E_OUT = 1.0e-4  # outbound initial-surface Misner-Sharp energy (x0 units)
V1 = 2.0 * Y_STATION + 0.5  # outbound chi1 on u = 0 for v in (V1, V1 + W_OUT)
W_OUT = 1.0
W_RET = 1.0  # return chi2 on v = 0 for u in (u1, u1 + W_RET); only E_ret and u1 vary
U_MAX, V_MAX = 14.0, 24.0
TAU_DEAD = 12.0  # receiver deadline: station-A proper time (deviation: was 14.0, see the record)
TAU_DEAD_ALT = (11.5, 12.5)  # deadline sensitivity
E_RET = (1e-8, 1e-3, 1e-2, 1e-1)
E_RET_POST_HOC = (
    0.3,
    0.5,
)  # added after the first n = 800 scan showed no window change up to 0.1 (deviation)
U1_GRID = (
    8.0,
    8.5,
    9.0,
    9.5,
    9.75,
    10.0,
    10.25,
    10.5,
    11.0,
)  # 9.75 and 10.25 added post hoc (deviation)
PASS_FRACTION = 0.5


def frozen_params(n: int, e_ret: float = 0.0, u1: float = 10.0, e_out: float = E_OUT) -> Params:
    a1 = amplitude_for_energy(e_out, V1, V1 + W_OUT, n, V_MAX) if e_out > 0 else 0.0
    a2 = amplitude_for_energy(e_ret, u1, u1 + W_RET, n, U_MAX) if e_ret > 0 else 0.0
    return Params(
        n=n, u_max=U_MAX, v_max=V_MAX, a1=a1, v1=V1, v2=V1 + W_OUT, a2=a2, u1=u1, u2=u1 + W_RET,
        r_station=R_STATION,
    )  # fmt: skip


def run_case(args: tuple) -> dict:
    """Pool worker: (n, e_ret, u1, e_out, grid_path) -> summary (grid pickled only if asked)."""
    n, e_ret, u1, e_out, grid_path = args
    res = evolve(frozen_params(n, e_ret, u1, e_out), keep_grid=grid_path is not None)
    if grid_path is not None:
        import pickle

        with open(grid_path, "wb") as fh:
            pickle.dump(res, fh)
    s = summarize(res, TAU_DEAD)
    s["case"] = {"n": n, "e_ret": e_ret, "u1": u1, "e_out": e_out}
    s["A_reaches_deadline"] = bool(s["A"].get("tau_end", 0.0) >= TAU_DEAD)
    s["valid"] = run_valid(s)
    return s


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(f"{float(o):.10g}")
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(obj), indent=1) + "\n")


def pool_map(cases: list, workers: int = 4) -> list:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ProcessPoolExecutor(workers) as ex:
        return list(ex.map(run_case, cases))


# ------------------------------------------------------------------------------------------- checks
def throat_history(res: dict, times) -> list:
    """Minimum of r on the slice u + v = 2t (the throat radius), on the stored output grid."""
    r, h = res["grid"]["r"], res["grid_h"]
    out = []
    for t in times:
        k = int(round(2 * t / h))
        ii = np.arange(max(0, k - (r.shape[1] - 1)), min(k, r.shape[0] - 1) + 1)
        vals = r[ii, k - ii]
        out.append(float(np.nanmin(vals)) if np.isfinite(vals).any() else None)
    return out


def escape_u(res: dict, r_far: float = 3.0) -> float:
    """Largest stored u whose radial ray (u = const) is beyond r_far at v = v_max (escape proxy)."""
    r, h = res["grid"]["r"], res["grid_h"]
    ok = np.isfinite(r[:, -1]) & (r[:, -1] > r_far)
    return float(np.nonzero(ok)[0].max() * h)


def _check_case(args):
    kind, n = args
    if kind == "static":
        res = evolve(Params(n=n, u_max=4.0, v_max=4.0, a1=0.0, a2=0.0))
        r = res["grid"]["r"]
        u = np.arange(r.shape[0]) * res["grid_h"]
        uu, vv = np.meshgrid(u, u, indexing="ij")
        err = np.abs(r - np.sqrt(1.0 + ((vv - uu) / 2.0) ** 2))
        return {"kind": kind, "n": n, "max_abs_dr": float(np.nanmax(err)),
                "max_abs_sigma": float(np.nanmax(np.abs(res["grid"]["s"]))),
                "constraint_max": max(float(np.max(res["cons"][c])) for c in ("u", "v")),
                }  # fmt: skip
    if (
        kind == "kain"
    ):  # Kain 2602.18231v1 Fig. 9 (file fig_3_07): A = 0.002, v in (1, 3), x0 = 1, a1 = 0
        res = evolve(Params(n=n, u_max=20.0, v_max=20.0, a1=0.002, v1=1.0, v2=3.0, a2=0.0))
        t = np.arange(0.0, 9.0, res["grid_h"])
        rth = np.array([np.nan if x is None else x for x in throat_history(res, t)])
        sel = (t >= 5.0) & (t <= 6.5)
        lam = float(np.polyfit(t[sel], np.log(1.0 - rth[sel]), 1)[0])
        return {"kind": kind, "n": n, "u_escape_v20": escape_u(res), "efold_fit_t5_6.5": 1.0 / lam,
                "one_minus_rth_t7": float(1.0 - rth[np.argmin(np.abs(t - 7.0))]),
                "id_energy": res["id_energy"]["chi1"]}  # fmt: skip
    if kind.startswith("noise"):  # outbound only, launched from r = 2, energy E_kain or E_OUT
        e = 5.397872426385404e-05 if kind == "noise_ekain" else E_OUT
        res = evolve(frozen_params(n, 0.0, 10.0, e))
        return {"kind": kind, "n": n, "e_out": e, "throat_t": [6, 8, 10, 11],
                "throat_r": throat_history(res, [6, 8, 10, 11])}  # fmt: skip
    raise ValueError(kind)


def run_checks(workers: int = 4) -> dict:
    cases = [("static", n) for n in (40, 80, 160)] + [("kain", n) for n in (100, 200, 400)]
    cases += [("noise_ekain", n) for n in (100, 200, 400)] + [("noise_eout", n) for n in (200, 400)]
    with ProcessPoolExecutor(workers) as ex:
        rows = list(ex.map(_check_case, cases))
    return {"checks": rows, "lit": {"ggs09_efold_linear": 0.846, "ggs09_efold_code": 0.836,
            "kain26_fig9_escape_u_figure_read": 9.2,
            "kain26_text_event_horizon_u": 9.5}}  # fmt: skip


def load_rows() -> list:
    """All scored runs; a run in validity_recheck.json (same case, re-run with the deadline-limited
    station check) replaces its scan entry."""
    rows = []
    for name in ("scan_n400.json", "scan_n800.json", "converge_n1600.json"):
        if (OUT / name).exists():
            rows += json.loads((OUT / name).read_text())["runs"]
    if (OUT / "validity_recheck.json").exists():
        for new in json.loads((OUT / "validity_recheck.json").read_text())["runs"]:
            key = (new["case"]["n"], new["case"]["e_ret"], new["case"]["u1"])
            rows = [
                r for r in rows if (r["case"]["n"], r["case"]["e_ret"], r["case"]["u1"]) != key
            ] + [new]
    return rows


def richardson(res: dict) -> dict:
    """Boundary per E_ret at n = 400, 800, 1600; observed order and a second-order extrapolation.

    The extrapolation assumes the error ~ n^-2 (the scheme's order); the observed ratio is reported
    so the assumption can be checked. Uncertainty = |b(1600) - b(800)| (the last change), not a
    formal bound.
    """
    out = {}
    for label, key in (
        ("dead_12", "boundary"),
        ("dead_11.5", "boundary_dead_11.5"),
        ("dead_12.5", "boundary_dead_12.5"),
    ):
        for e in E_RET + E_RET_POST_HOC:
            b = {n: res[key].get(f"n{n}_E{e:g}", {}).get("tauB_boundary") for n in (400, 800, 1600)}
            rec = {f"n{n}": b[n] for n in b}
            if None not in b.values():
                d1, d2 = b[800] - b[400], b[1600] - b[800]
                rec["ratio"] = d1 / d2 if d2 != 0 else None
                rec["extrap"] = b[1600] + d2 / 3.0
                rec["last_change"] = abs(d2)
            out[f"{label}_E{e:g}"] = rec
    return out


def case_table(rows: list) -> list:
    """Energies and proper times of selected cases at every resolution (the pass/fail set)."""
    sel = {
        (1e-8, 9.5),
        (1e-8, 10.5),
        (0.1, 9.5),
        (0.1, 10.5),
        (0.5, 8.0),
        (0.5, 9.5),
        (0.5, 10.0),
        (0.5, 10.5),
    }
    keep = ("E_emit_B", "E_recv_A", "E_emit_B_local", "E_recv_A_local", "E_reflected_B", "fraction",
            "tauB_emit50", "tauA_recv50", "round_trip_tauA", "reply_compatible",
            "pass")  # fmt: skip
    out = []
    for r in sorted(rows, key=lambda r: (r["case"]["e_ret"], r["case"]["u1"], r["case"]["n"])):
        c = r["case"]
        if (c["e_ret"], c["u1"]) in sel:
            rec = {"n": c["n"], "e_ret": c["e_ret"], "u1": c["u1"]}
            rec.update({k: r["return"][k] for k in keep})
            rec.update({"outbound_" + k: v for k, v in r["outbound"].items()})
            rec["test_ray_u_A_max"] = r["test_ray"]["u_A_max"]
            rec["test_ray_tauB"] = r["test_ray"]["tauB_at_u_A_max"]
            out.append(rec)
    return out


def _imprint_case(args):
    n, e_ret, u1, with_pulse = args
    p = frozen_params(n, e_ret, u1)
    if not with_pulse:
        p = replace(p, a2=0.0)
    res = evolve(p, keep_grid=False)
    t = res["stations"]["B"]
    return {k: t[k] for k in ("tau", "dm_u", "dm_v")}


def run_imprint(workers: int = 4) -> dict:
    """Review item 2: the pre-loaded return pulse perturbs the ghost field on its way to B, so
    the ghost inflow through B during the chi2 crossing differs from the no-pulse run. A transmitter
    at B would not do this. Imprint = (ghost Kodama inflow with chi2 - without) / chi2 Kodama
    inflow, over the chi2 10-90 % crossing window widened by 0.25 on each side."""
    cases = [(n, e, 9.75, w) for n in (400, 800) for e in (0.1, 0.5) for w in (True, False)]
    with ProcessPoolExecutor(workers) as ex:
        tabs = list(ex.map(_imprint_case, cases))
    out = []
    for i in range(0, len(cases), 2):
        n, e, u1, _ = cases[i]
        on, off = tabs[i], tabs[i + 1]
        tm = 0.5 * (on["tau"][1:] + on["tau"][:-1])
        chi = on["dm_u"][:, 2]
        t10, t90 = crossing_time(tm, chi, 0.1), crossing_time(tm, chi, 0.9)
        win_on = (tm >= t10 - 0.25) & (tm <= t90 + 0.25)
        tm_off = 0.5 * (off["tau"][1:] + off["tau"][:-1])
        win_off = (tm_off >= t10 - 0.25) & (tm_off <= t90 + 0.25)
        g_on, g_off = float(np.sum(on["dm_u"][win_on, 0])), float(np.sum(off["dm_u"][win_off, 0]))
        chi_in = float(np.sum(chi))
        out.append({"n": n, "e_ret": e, "u1": u1, "chi2_in": chi_in, "ghost_in_with": g_on,
                    "ghost_in_without": g_off, "imprint": (g_on - g_off) / chi_in})  # fmt: skip
    return {"imprint": out}


# -------------------------------------------------------------------------------------------
# figures
def _diagram(ax, res, title, lim):
    r, h = res["grid"]["r"], res["grid_h"]
    u = np.arange(r.shape[0]) * h
    v = np.arange(r.shape[1]) * h
    ru, rv = np.gradient(r, h, axis=0), np.gradient(r, h, axis=1)
    cs = ax.contour(
        v, u, r, levels=[0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], colors="0.6", linewidths=0.6
    )
    ax.clabel(cs, fontsize=6)
    ax.contour(v, u, ru, levels=[0.0], colors="k", linewidths=1.4)
    ax.contour(v, u, rv, levels=[0.0], colors="tab:blue", linewidths=1.4)
    ax.contourf(v, u, np.isfinite(r).astype(float), levels=[-0.5, 0.5], colors=["#f4b6b6"])
    for side, col in (("A", "tab:green"), ("B", "tab:purple")):
        t = res["stations"][side]
        ax.plot(t["v"], t["u"], col, lw=1.6, label=f"station {side} (r = 2)")
    p = res["params"]
    ax.axvspan(
        p["v1"], p["v2"], color="tab:orange", alpha=0.25, lw=0, label="outbound chi1 (v band)"
    )
    if p["a2"] > 0:
        ax.axhspan(
            p["u1"], p["u2"], color="tab:red", alpha=0.25, lw=0, label="return chi2 (u band)"
        )
    ax.set_xlim(*lim[:2])
    ax.set_ylim(*lim[2:])
    ax.set_xlabel("v  (ingoing null coordinate, x0)")
    ax.set_ylabel("u  (outgoing null coordinate, x0)")
    ax.set_title(title, fontsize=9)


def make_figures(workers: int = 4) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = [(400, 0.1, 9.5, E_OUT, None), (400, 0.1, 10.5, E_OUT, None)]
    with ProcessPoolExecutor(workers) as ex:
        runs = list(ex.map(grid_case, cases))
        kain = ex.submit(
            evolve, Params(n=400, u_max=20.0, v_max=20.0, a1=0.002, v1=1.0, v2=3.0)
        ).result()
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
    r, h = kain["grid"]["r"], kain["grid_h"]
    u = np.arange(r.shape[0]) * h
    ru, rv = np.gradient(r, h, axis=0), np.gradient(r, h, axis=1)
    cs = axs[0].contour(
        u, u, r, levels=[0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 2.0], colors="0.5", linewidths=0.7
    )
    axs[0].clabel(cs, fontsize=7)
    axs[0].contour(u, u, ru, levels=[0.0], colors="k", linewidths=1.5)
    axs[0].contour(u, u, rv, levels=[0.0], colors="tab:blue", linewidths=1.5)
    axs[0].contour(
        u, u, np.isfinite(r).astype(float), levels=[0.5], colors="tab:red", linewidths=1.5
    )
    axs[0].set_xlim(6, 11)
    axs[0].set_ylim(6, 11)
    axs[0].set_xlabel("v")
    axs[0].set_ylabel("u")
    axs[0].set_title("Kain26 Fig. 9(b) case (A = 0.002), ours at n = 400 (by eye)", fontsize=8)
    t = np.arange(0.0, 9.0, h)
    rth = np.array([np.nan if x is None else x for x in throat_history(kain, t)])
    axs[1].semilogy(t, 1.0 - rth, "k", lw=1.2, label="1 - throat radius (ours, n = 400)")
    t0 = 5.0
    y0 = 1.0 - rth[np.argmin(np.abs(t - t0))]
    axs[1].semilogy(
        t, y0 * np.exp((t - t0) / 0.846), "r--", lw=0.8, label="slope 1/0.846 (GGS09 linear)"
    )
    axs[1].set_ylim(1e-5, 1)
    axs[1].set_xlabel("t on the throat (x0/c)")
    axs[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_benchmark.png", dpi=90)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.0))
    labels = ("E_ret = 0.1, u1 = 9.5 (passes)", "E_ret = 0.1, u1 = 10.5 (fails)")
    for ax, res, lab in zip(axs, runs, labels, strict=True):
        _diagram(ax, res, f"{lab}; n = 400", (0, 24, 0, 14))
    axs[0].legend(fontsize=7, loc="upper left")
    fig.suptitle(
        "r(u, v): black r_u = 0, blue r_v = 0 (apparent horizons); pink: excised (r < 0.05)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig_spacetime.png", dpi=90)
    plt.close(fig)

    rows = load_rows()
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    colors = dict(
        zip(E_RET + E_RET_POST_HOC, plt.cm.viridis(np.linspace(0.0, 0.9, 6)), strict=True)
    )
    for (n, e), rec in ((k, v) for k, v in _series(rows).items()):
        style = {
            400: dict(ls="--", lw=0.8),
            800: dict(ls="-", lw=1.3),
            1600: dict(ls="", marker="o", ms=4),
        }[n]
        ax.plot(rec["tauB_emit50"], rec["fraction"], color=colors.get(e, "k"), **style,
                label=f"E_ret = {e:g}" if n == 800 else None)  # fmt: skip
    ax.axhline(PASS_FRACTION, color="r", lw=0.8)
    b8 = [r for r in rows if r["case"]["n"] == 800 and r["case"]["e_ret"] == 1e-8]
    if b8:
        ax.axvline(b8[0]["test_ray"]["tauB_at_u_A_max"], color="0.3", ls=":", lw=1,
                   label="test null ray: last B time reaching A by the deadline (n = 800)",
                   )  # fmt: skip
        ax.axvline(b8[0]["outbound"]["tauB_recv90"], color="tab:orange", ls=":", lw=1,
                   label="outbound 90 % received at B")  # fmt: skip
    ax.set_xlabel("B proper time of return emission (50 % point), x0/c")
    ax.set_ylabel("received / emitted Kodama energy (A by tau_A = 12)")
    ax.set_title("dashed n = 400, solid n = 800, dots n = 1600", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_fraction.png", dpi=90)
    plt.close(fig)


def _series(rows: list) -> dict:
    out = {}
    for key in sorted({(r["case"]["n"], r["case"]["e_ret"]) for r in rows}):
        ser = sorted(
            (r for r in rows if (r["case"]["n"], r["case"]["e_ret"]) == key),
            key=lambda r: r["case"]["u1"],
        )
        out[key] = {"tauB_emit50": [r["return"]["tauB_emit50"] for r in ser],
                    "fraction": [r["return"]["fraction"] for r in ser]}  # fmt: skip
    return out


def grid_case(args: tuple) -> dict:
    """Pool worker returning the full evolution (with the stored grid) for a frozen-problem case."""
    n, e_ret, u1, e_out, _ = args
    return evolve(frozen_params(n, e_ret, u1, e_out), keep_grid=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mode", choices=["scan", "checks", "figures", "analyze", "imprint"])
    ap.add_argument("--n", type=int, nargs="+", default=[800])
    ap.add_argument("--e-ret", type=float, nargs="+", default=list(E_RET + E_RET_POST_HOC))
    ap.add_argument("--u1", type=float, nargs="+", default=list(U1_GRID))
    ap.add_argument("--e-out", type=float, default=E_OUT)
    ap.add_argument("--out", default="scan.json")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    if a.mode == "imprint":
        write_json(OUT / "imprint.json", run_imprint(a.workers))
        print((OUT / "imprint.json").read_text())
        return 0
    if a.mode == "figures":
        make_figures(a.workers)
        return 0
    if a.mode == "analyze":
        rows = load_rows()
        res = {"tau_dead": TAU_DEAD, "boundary": boundary(rows)}
        for td in TAU_DEAD_ALT:
            res[f"boundary_dead_{td:g}"] = boundary(rows, f"fraction_dead_{td:g}")
        res["richardson"] = richardson(res)
        res["cases"] = case_table(rows)
        write_json(OUT / "boundary.json", res)
        for k, v in res["boundary"].items():
            print(
                k,
                "monotone",
                v["monotone"],
                "tauB_boundary",
                v.get("tauB_boundary"),
                "u1_boundary",
                v.get("u1_boundary"),
            )
        print(json.dumps(_clean(res["richardson"]), indent=1))
        return 0
    if a.mode == "checks":
        write_json(OUT / "checks.json", run_checks(a.workers))
        print((OUT / "checks.json").read_text())
        return 0
    cases = [(n, e, u1, a.e_out, None) for n in a.n for e in a.e_ret for u1 in a.u1]
    t0 = time.time()
    rows = pool_map(cases, a.workers)
    write_json(OUT / a.out, {"tau_dead": TAU_DEAD, "pass_fraction": PASS_FRACTION, "runs": rows})
    print(f"{len(rows)} runs in {time.time() - t0:.0f} s -> {OUT / a.out}")
    for s in rows:
        c, r_ = s["case"], s.get("return", {})
        print(
            f"n={c['n']} E2={c['e_ret']:.0e} u1={c['u1']:.3f} tauB_emit={r_.get('tauB_emit50')} "
            f"frac={r_.get('fraction')} pass={r_.get('pass')} reply={r_.get('reply_compatible')} "
            f"tauA_recv={r_.get('tauA_recv50')} uAmax={s.get('test_ray', {}).get('u_A_max')}"
        )
    return 0


# ------------------------------------------------------------------------------------------
# analysis
def boundary(rows: list, frac_key: str = "fraction") -> dict:
    """Per (n, E_ret): B emission time at which the received fraction crosses PASS_FRACTION.

    Linear interpolation in tauB_emit50 between the last passing and the first failing point (the
    fraction falls monotonically with emission time in every scanned series; checked here).
    """
    out = {}
    for key in sorted({(r["case"]["n"], r["case"]["e_ret"]) for r in rows}):
        ser = sorted((r for r in rows if (r["case"]["n"], r["case"]["e_ret"]) == key
                      and run_valid(r)),
                     key=lambda r: r["case"]["u1"])  # fmt: skip
        tb = [r["return"]["tauB_emit50"] for r in ser]
        fr = [r["return"][frac_key] for r in ser]
        mono = all(b <= a + 1e-12 for a, b in zip(fr, fr[1:], strict=False))
        rec = {
            "monotone": mono,
            "u1": [r["case"]["u1"] for r in ser],
            "tauB_emit50": tb,
            "fraction": fr,
        }
        j = next((i for i in range(1, len(fr)) if fr[i - 1] >= PASS_FRACTION > fr[i]), None)
        if j is not None:
            w = (fr[j - 1] - PASS_FRACTION) / (fr[j - 1] - fr[j])
            rec["tauB_boundary"] = tb[j - 1] + w * (tb[j] - tb[j - 1])
            rec["u1_boundary"] = ser[j - 1]["case"]["u1"] + w * (
                ser[j]["case"]["u1"] - ser[j - 1]["case"]["u1"]
            )
            if (
                mono and len(fr) >= 3
            ):  # shape-preserving cubic through all points: interpolation-error check
                from scipy.interpolate import PchipInterpolator
                from scipy.optimize import brentq

                f = PchipInterpolator(tb, np.array(fr) - PASS_FRACTION)
                rec["tauB_boundary_pchip"] = float(brentq(f, tb[j - 1], tb[j]))
        out[f"n{key[0]}_E{key[1]:g}"] = rec
    return out


if __name__ == "__main__":
    sys.exit(main())
