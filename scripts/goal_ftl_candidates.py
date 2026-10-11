"""Can a receiver read a freely chosen message before exterior light arrives? Two candidates.

Goal round 5 (docs/goal_round5/README.md). Offline, no data.

A. Krasnikov tube (Krasnikov, PRD 57, 4760 (1998), gr-qc/9511068; Everett & Roman, PRD 56, 2100
   (1997), gr-qc/9702049, "ER"). Metric (ER Eq. 13), c = G = hbar = 1:
       ds^2 = -dt^2 + (1 - k) dx dt + k dx^2 + d rho^2 + rho^2 d phi^2,
   with k = 1 - eta theta_eps(rho_max - rho) in the middle of a long tube long after it was built
   (eta = 2 - delta; ER Eq. 12) and theta_eps(s) = (tanh(2 (2 s / eps - 1)) + 1) / 2 (ER Eq. 35).
   Inside, light moving in -x has dx/dt = -1/(1 - eta): a message sent from the far end at time t_s
   reaches the near end at t_s + D (1 - eta), against t_s + D for exterior light, an advance of
   eta D. The tube is laid at <= c by an outbound trip, so it exists only for t > x.
   Computed here [derived]:
   1. Positive controls: the energy density of ER Eq. (14) from an independent Einstein-tensor
      calculation, and ER's T ~ -1/(8 pi eps^2) for delta -> 0 (Eq. 39).
   2. For small eta: the most negative static energy density and the largest orthonormal
      curvature component, both as functions of eta and eps.
   3. The Ford-Roman quantum inequality (ER Eq. 43; massless scalar, flat-space form) applied
      with sampling time tau0 = sigma * (curvature radius): the thickest wall allowed for a given
      advance over D = 20 and 100 light years, and the negative energy that implies (ER Eq. 51).
B. Scharnhorst effect (QED between perfect plates; Liberati, Sonego & Visser, Ann. Phys. 298, 167
   (2002), gr-qc/0107091, "LSV", Eq. 3.3): xi = 11 pi^2 alpha^2 / (4050 (a m_e)^4), light normal
   to the plates moves at c sqrt(1 + xi). Computed here [derived]: the time advance across one gap
   and whether it is resolvable inside the soft-photon regime (omega << m_e) where xi is derived.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

PLANCK_LENGTH_M = 1.616255e-35
C_M_S = 299_792_458.0
LY_M = 9.4607304725808e15
YEAR_S = 365.25 * 86400.0
ALPHA = 1 / 137.035999
HBAR_C_EV_M = 1.973269804e-7  # hbar c in eV m
ME_EV = 0.51099895e6
LAMBDA_C_BAR_M = HBAR_C_EV_M / ME_EV  # reduced Compton wavelength, 3.86e-13 m
M_SUN_KG = 1.98847e30
PLANCK_MASS_KG = 2.176434e-8


@dataclass(frozen=True)
class Params:
    """Numerical choices (ASSUMPTIONs)."""

    sigma: float = 0.01  # sampling time / curvature radius for the flat-space QI (ER: sigma << 1)
    n_wall: float = 100.0  # rho_max / eps
    alpha_band: float = 0.01  # fraction of the wall holding the most negative energy (ER Eq. 49)
    rho_max_m: float = 1.0  # tube radius for the energy estimate (ER's example)
    grid: int = 20001
    distances_ly: tuple = (20.0, 100.0)
    advances_s: tuple = (1.0, 86400.0, YEAR_S)
    scharnhorst_min_gap_compton: float = 1.0  # smallest a / lambda_C considered (EFT edge)


P = Params()

t, x, rho, phi = sp.symbols("t x rho phi", real=True)


@lru_cache(maxsize=1)
def krasnikov_tensors():
    """Einstein tensor and orthonormal Riemann components for ER Eq. (13) with k = k(rho)."""
    k = sp.Function("k")(rho)
    coords = (t, x, rho, phi)
    g = sp.Matrix(
        [
            [-1, (1 - k) / 2, 0, 0],
            [(1 - k) / 2, k, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, rho**2],
        ]
    )
    ginv = sp.simplify(g.inv())
    n = 4
    gam = [
        [
            [
                sp.simplify(
                    sum(
                        ginv[a, d]
                        * (
                            sp.diff(g[d, b], coords[c])
                            + sp.diff(g[d, c], coords[b])
                            - sp.diff(g[b, c], coords[d])
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

    def riem(a, b, c, d):  # R^a_{bcd}
        expr = sp.diff(gam[a][b][d], coords[c]) - sp.diff(gam[a][b][c], coords[d])
        expr += sum(gam[a][c][e] * gam[e][b][d] - gam[a][d][e] * gam[e][b][c] for e in range(n))
        return sp.simplify(expr)

    R = [
        [[[riem(a, b, c, d) for d in range(n)] for c in range(n)] for b in range(n)]
        for a in range(n)
    ]
    ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(R[a][b][a][d] for a in range(n))))
    rs = sp.simplify(sum(ginv[a, b] * ric[a, b] for a in range(n) for b in range(n)))
    ein = sp.simplify(ric - rs * g / 2)
    # lower all indices, then contract with the static orthonormal basis (ER Eqs. 24-27)
    rl = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    rl[(a, b, c, d)] = sum(g[a, e] * R[e][b][c][d] for e in range(n))
    basis = [
        [1, 0, 0, 0],
        [(1 - k) / (1 + k), 2 / (1 + k), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1 / rho],
    ]
    comps = {}
    idx = range(n)
    for A in idx:
        for B in idx:
            if B <= A:
                continue
            for C in idx:
                for D in idx:
                    if D <= C or (C, D) < (A, B):
                        continue
                    val = sum(
                        rl[(a, b, c, d)] * basis[A][a] * basis[B][b] * basis[C][c] * basis[D][d]
                        for a in idx
                        for b in idx
                        for c in idx
                        for d in idx
                        if basis[A][a] != 0
                        and basis[B][b] != 0
                        and basis[C][c] != 0
                        and basis[D][d] != 0
                    )
                    val = sp.simplify(val)
                    if val != 0:
                        comps[(A, B, C, D)] = val
    t_tt = sp.simplify(ein[0, 0] / (8 * sp.pi))  # energy density of the static observer (ER Eq. 31)
    return k, t_tt, comps


def er_eq14(k):
    """ER Eq. (14): T_tt = [32 pi (1 + k)^2]^-1 (-4 (1 + k) k'/rho + 3 k'^2 - 4 (1 + k) k'')."""
    kp, kpp = sp.diff(k, rho), sp.diff(k, rho, 2)
    return (-4 * (1 + k) * kp / rho + 3 * kp**2 - 4 * (1 + k) * kpp) / (32 * sp.pi * (1 + k) ** 2)


def theta(s, eps):
    return (np.tanh(2 * (2 * s / eps - 1)) + 1) / 2


@lru_cache(maxsize=1)
def _numeric():
    k, t_tt, comps = krasnikov_tensors()
    k0, k1, k2 = sp.symbols("k0 k1 k2")
    subs = {sp.diff(k, rho, 2): k2, sp.diff(k, rho): k1, k: k0}

    def lam(e):
        return sp.lambdify((rho, k0, k1, k2), e.subs(subs), "numpy")

    return lam(t_tt), [lam(v) for v in comps.values()], lam(er_eq14(k))


def wall_profile(eta: float, eps: float, n_wall: float = P.n_wall, num: int = P.grid):
    """k, k', k'' across the wall at mid-tube (theta(x) = 1, theta(t - x - rho) = 1)."""
    rho_max = n_wall * eps
    r = np.linspace(rho_max - 3 * eps, rho_max + 2 * eps, num)
    s = rho_max - r
    th = theta(s, eps)
    # d theta / d rho = -d theta / d s
    sech2 = 1 / np.cosh(2 * (2 * s / eps - 1)) ** 2
    dth_ds = 2 * sech2 / eps
    d2th_ds2 = -16 * np.tanh(2 * (2 * s / eps - 1)) * sech2 / eps**2
    k = 1 - eta * th
    k1 = eta * dth_ds
    k2 = -eta * d2th_ds2
    return r, k, k1, k2


def wall_extremes(eta: float, eps: float) -> dict:
    """Most negative static energy density and largest |orthonormal Riemann component|."""
    f_t, f_r, f_er = _numeric()
    r, k, k1, k2 = wall_profile(eta, eps)
    tt = f_t(r, k, k1, k2)
    er = f_er(r, k, k1, k2)
    rmax = np.max([np.max(np.abs(np.asarray(f(r, k, k1, k2)) * np.ones_like(r))) for f in f_r])
    return {
        "T_min": float(tt.min()),
        "T_max": float(tt.max()),
        "max_abs_diff_vs_ER_eq14": float(np.max(np.abs(tt - er))),
        "R_max": float(rmax),
    }


def qi_wall_limit(eta: float, sigma: float = P.sigma) -> dict:
    """Thickest wall (Planck lengths) allowed by the QI |T| <= 3 / (32 pi^2 tau0^4) with
    tau0 = sigma r_c, r_c = R_max^-1/2. Both |T| and R_max scale as eta / eps^2 for small eta
    (checked numerically), so eps_max follows in closed form from one evaluation at eps = 1."""
    w = wall_extremes(eta, 1.0)
    a_t, a_r = -w["T_min"], w["R_max"]  # coefficients of 1 / eps^2
    # a_t / eps^2 <= 3 a_r^2 / (32 pi^2 sigma^4 eps^4)  =>  eps^2 <= 3 a_r^2 / (32 pi^2 sigma^4 a_t)
    eps_max = math.sqrt(3 * a_r**2 / (32 * math.pi**2 * sigma**4 * a_t))
    # the looser prescription tau0 = sigma * eps (ER Eq. 42 with r_c ~ eps)
    eps_max_loose = math.sqrt(3 / (32 * math.pi**2 * sigma**4 * a_t))
    return {
        "eta": eta,
        "T_coeff": a_t,
        "R_coeff": a_r,
        "rc_over_eps": 1 / math.sqrt(a_r),
        "eps_max_planck": eps_max,
        "eps_max_planck_loose_tau0_sigma_eps": eps_max_loose,
    }


def negative_energy_solar(eta: float, eps_planck: float, d_ly: float) -> float:
    """ER Eq. (51): E ~ |T| * 2 pi rho_max (alpha eps) D, with |T| = a_t eta... evaluated with the
    computed coefficient at the given wall thickness. Returned in solar masses (negative)."""
    w = wall_extremes(eta, 1.0)
    t_abs = -w["T_min"] / eps_planck**2  # Planck units
    rho_max = P.rho_max_m / PLANCK_LENGTH_M
    d = d_ly * LY_M / PLANCK_LENGTH_M
    e_planck = t_abs * 2 * math.pi * rho_max * P.alpha_band * eps_planck * d
    return -e_planck * PLANCK_MASS_KG / M_SUN_KG


def krasnikov_task() -> dict:
    """Advance and QI-limited wall for the owner's distances and target advances."""
    rows = []
    for d_ly in P.distances_ly:
        d_s = d_ly * LY_M / C_M_S
        for adv in (*P.advances_s, d_s):
            eta = adv / d_s
            lim = qi_wall_limit(eta)
            rows.append(
                {
                    "D_ly": d_ly,
                    "advance_s": adv,
                    "eta": eta,
                    "eps_max_planck": lim["eps_max_planck"],
                    "eps_max_m": lim["eps_max_planck"] * PLANCK_LENGTH_M,
                    "eps_max_loose_m": lim["eps_max_planck_loose_tau0_sigma_eps"] * PLANCK_LENGTH_M,
                    "negative_energy_solar_masses_at_eps_max": negative_energy_solar(
                        eta, lim["eps_max_planck"], d_ly
                    ),
                }
            )
    return {
        "timing": {
            "deployment": "tube exists only for t > x: built by an outbound trip at <= c, so "
            "nothing can be read early before t = D at the far end",
            "far_to_near_arrival": "t_s + D (1 - eta), exterior light t_s + D: advance eta D",
            "near_to_far": "no advance (the +x branch of the light cone is unchanged)",
            "two_tubes": "a second tube for the other direction allows closed timelike curves "
            "(ER Sec. 4)",
        },
        "rows": rows,
    }


def krasnikov_controls() -> dict:
    """ER Eq. (14) from an independent Einstein tensor; ER Eq. (39) T ~ -1/(8 pi eps^2) as delta
    -> 0; and the small-eta scaling of T and R_max."""
    near = wall_extremes(2 - 1e-3, 1.0)  # mid-tube, delta = 1e-3
    # ER Eq. (38): end cap (theta(x) = 1/2, so eta = 1 - delta/2), observer at rho = rho_max - eps
    delta, n = 1e-3, P.n_wall
    f_t, _, _ = _numeric()
    r, k, k1, k2 = wall_profile(1 - delta / 2, 1.0)
    i = int(np.argmin(np.abs(r - (n - 1.0))))
    t_cap = float(f_t(r[i], k[i], k1[i], k2[i]))
    num = 0.271 * (
        3.504 - 4.034 * n - 0.008 * delta + 0.017 * n * delta - 0.872 * delta**2 + n * delta**2
    )
    t_eq38 = num / (8 * math.pi * (n - 1) * (1.018 + 0.491 * delta) ** 2)
    scal = [qi_wall_limit(e) for e in (1e-2, 1e-4, 1e-6)]
    return {
        "max_abs_diff_vs_ER_eq14_delta_small": near["max_abs_diff_vs_ER_eq14"],
        "endcap_T_at_rho_max_minus_eps": t_cap,
        "ER_eq38_value": t_eq38,
        "small_eta_T_coeff_over_eta": [s["T_coeff"] / s["eta"] for s in scal],
        "small_eta_R_coeff_over_eta": [s["R_coeff"] / s["eta"] for s in scal],
    }


def scharnhorst(gaps_m=(1e-6, 1e-8, 1e-10)) -> dict:
    """LSV Eq. (3.3) and the advance across one gap; resolvability within omega << m_e."""
    omega_c = ME_EV / HBAR_C_EV_M * C_M_S  # m_e c^2 / hbar, s^-1

    def xi(a):
        return 11 * math.pi**2 * ALPHA**2 / (4050 * (a / LAMBDA_C_BAR_M) ** 4)

    rows = []
    a_min = P.scharnhorst_min_gap_compton * LAMBDA_C_BAR_M
    for a in (*gaps_m, a_min):
        z = xi(a)
        # a/c (1 - (1 + xi)^-1/2), written to avoid cancellation
        adv = a / C_M_S * -math.expm1(-0.5 * math.log1p(z))
        rows.append(
            {
                "gap_m": a,
                "xi": z,
                "advance_s": adv,
                "advance_times_mec2_over_hbar": adv * omega_c,
            }
        )
    return {
        "xi_1um": xi(1e-6),
        "LSV_quoted_xi_1um": 4.36e-32,  # LSV below Eq. (3.3); the formula gives 3.17e-32
        "rows": rows,
        "note": "the advance exists only for propagation normal to the plates, across one gap",
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived from published metrics, field-theory results and the flat-space "
        "quantum inequality (no data); estimates under the stated ASSUMPTIONs",
        "params": asdict(P),
        "krasnikov_controls": krasnikov_controls(),
        "krasnikov_task": krasnikov_task(),
        "scharnhorst": scharnhorst(),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
