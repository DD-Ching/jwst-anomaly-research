"""Can a physical field and state source the Krasnikov tube used in goal round 5?

Goal round 6 (docs/goal_round6/README.md). Offline, no data. Reuses the geometry, profile and
quantum-inequality machinery of scripts/goal_ftl_candidates.py (Everett & Roman 1997, "ER").

Computed here [derived]:
1. Weyl versus Ricci on the wall. Krasnikov 2003 (gr-qc/0207057, Sec. III A) notes that the QI
   chain rho ~ l^-4 ~ rho^2 needs max|C| / max|R_ab| <~ 1; where Weyl dominates, the vacuum
   polarisation is not tied to the Ricci part and the Planck-density conclusion lapses. The
   ratio is evaluated in the static orthonormal frame (ER Eqs. 24-27) for the profile actually
   used.
2. The owner's target is a LATENCY (time from sending to reading), not an ADVANCE over light.
   Rows are built for one-way latencies of 1 s, 1 min, 1 day and 1 yr over 20 and 100 ly
   (eta = 1 - latency / D, just below the closed-timelike-curve threshold eta = 1), next to the
   round-5 advance rows.
3. Sensitivity of the wall bound to the sampling-time rule (tau0 = sigma r_c, sigma eps, and
   sigma min(r_c, eps)) and of the negative-energy estimate to the tube radius (a message needs
   only a radius above the carrier wavelength) and to the band fraction alpha (ER Eq. 49).
4. A semiclassical-scale check: the wall thickness and curvature radius in Planck lengths.
5. Null energy on the wall: R_ab l^a l^b for the outgoing (exterior light) and incoming
   (shortcut) null directions, and whether t is a time function (slices t = const spacelike,
   bounded light speeds), which the curved-spacetime QEI theorems assume (global hyperbolicity).
6. The size of the first-order curvature corrections of the Kontou & Olum (arXiv:1410.0665)
   curved-spacetime QI for the static geodesic observer at the sampling time actually used.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp

_here = Path(__file__).resolve().parent
if "goal_ftl_candidates" in sys.modules:
    fc = sys.modules["goal_ftl_candidates"]
else:
    _spec = importlib.util.spec_from_file_location(
        "goal_ftl_candidates", _here / "goal_ftl_candidates.py"
    )
    fc = importlib.util.module_from_spec(_spec)
    sys.modules["goal_ftl_candidates"] = fc
    _spec.loader.exec_module(fc)


@dataclass(frozen=True)
class Params:
    """Numerical choices (ASSUMPTIONs)."""

    sigma: float = 0.01
    latencies_s: tuple = (1.0, 60.0, 86400.0, fc.YEAR_S)
    distances_ly: tuple = (20.0, 100.0)
    radii_m: tuple = (1.0, 1e-3, 1e-6)  # human-size, mm-wave, optical carrier
    alphas: tuple = (0.01, 0.1, 1.0)
    weyl_etas: tuple = (1e-6, 1e-2, 0.5, 1.0, 1.5, 1.9)


P = Params()


@lru_cache(maxsize=1)
def orthonormal_riemann():
    """All orthonormal components R_{ABCD} of ER Eq. (13) with k = k(rho), as numeric functions
    of (rho, k, k', k''), using the static basis of ER Eqs. (24)-(27)."""
    t, x, rho, phi = sp.symbols("t x rho phi", real=True)
    k = sp.Function("k")(rho)
    coords = (t, x, rho, phi)
    g = sp.Matrix(
        [[-1, (1 - k) / 2, 0, 0], [(1 - k) / 2, k, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]]
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

    def riem_up(a, b, c, d):
        e = sp.diff(gam[a][b][d], coords[c]) - sp.diff(gam[a][b][c], coords[d])
        e += sum(gam[a][c][m] * gam[m][b][d] - gam[a][d][m] * gam[m][b][c] for m in range(n))
        return e

    up = {i: riem_up(*i) for i in itertools.product(range(n), repeat=4)}
    low = {
        (a, b, c, d): sp.simplify(sum(g[a, m] * up[(m, b, c, d)] for m in range(n)))
        for a, b, c, d in itertools.product(range(n), repeat=4)
    }
    basis = [
        [1, 0, 0, 0],
        [(1 - k) / (1 + k), 2 / (1 + k), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1 / rho],
    ]
    k0, k1, k2 = sp.symbols("k0 k1 k2")
    subs = {sp.diff(k, rho, 2): k2, sp.diff(k, rho): k1, k: k0}
    funcs = {}
    for A, B, C, D in itertools.product(range(n), repeat=4):
        if A >= B or C >= D:
            continue
        val = sum(
            low[(a, b, c, d)] * basis[A][a] * basis[B][b] * basis[C][c] * basis[D][d]
            for a, b, c, d in itertools.product(range(n), repeat=4)
            if 0 not in (basis[A][a], basis[B][b], basis[C][c], basis[D][d])
        )
        val = sp.simplify(val.subs(subs))
        if val != 0:
            funcs[(A, B, C, D)] = sp.lambdify((rho, k0, k1, k2), val, "numpy")
    return funcs


def riemann_array(r, k, k1, k2) -> np.ndarray:
    """R[A, B, C, D] on the profile points (shape (N, 4, 4, 4, 4)) from the independent part."""
    out = np.zeros((len(r), 4, 4, 4, 4))
    for (A, B, C, D), f in orthonormal_riemann().items():
        v = np.asarray(f(r, k, k1, k2), dtype=float) * np.ones_like(r)
        for (a, b, sab), (c, d, scd) in itertools.product(
            ((A, B, 1), (B, A, -1)), ((C, D, 1), (D, C, -1))
        ):
            out[:, a, b, c, d] = sab * scd * v
    return out


def weyl_ricci(eta: float, eps: float = 1.0) -> dict:
    """max |C_ABCD| / max |R_AB| over the wall (static orthonormal frame, signature -+++)."""
    r, k, k1, k2 = fc.wall_profile(eta, eps)
    R = riemann_array(r, k, k1, k2)
    eta_m = np.diag([-1.0, 1.0, 1.0, 1.0])
    ric = np.einsum("ac,nabcd->nbd", eta_m, R)  # R_bd = eta^{ac} R_abcd
    scal = np.einsum("bd,nbd->n", eta_m, ric)
    g = eta_m
    # Weyl: C = R - (g ^ Ric terms) / 2 + R g g / 6  (4 dimensions)
    term1 = (
        np.einsum("ac,nbd->nabcd", g, ric)
        - np.einsum("ad,nbc->nabcd", g, ric)
        - np.einsum("bc,nad->nabcd", g, ric)
        + np.einsum("bd,nac->nabcd", g, ric)
    )
    term2 = np.einsum("n,ac,bd->nabcd", scal, g, g) - np.einsum("n,ad,bc->nabcd", scal, g, g)
    C = R - term1 / 2 + term2 / 6
    trace_check = float(np.max(np.abs(np.einsum("ac,nabcd->nbd", g, C))))
    c_max = float(np.max(np.abs(C)))
    ric_max = float(np.max(np.abs(ric)))
    # pointwise ratio where the energy density is most negative
    i = int(np.argmin(np.asarray(fc._numeric()[0](r, k, k1, k2))))
    return {
        "eta": eta,
        "max_abs_weyl": c_max,
        "max_abs_ricci": ric_max,
        "ratio_max": c_max / ric_max,
        "ratio_at_most_negative_T": float(np.max(np.abs(C[i])) / np.max(np.abs(ric[i]))),
        "weyl_trace_residual": trace_check,
    }


def kontou_olum_terms(eta: float, sigma: float = P.sigma) -> dict:
    """Size of the first-order curvature corrections in Kontou & Olum (arXiv:1410.0665, Eq. 130)
    for a static (geodesic) observer on the wall, relative to the flat term 3.76, with the
    sampling width t0 = sigma r_c. Order-of-magnitude estimate [derived]: rho-derivatives of the
    static-frame Ricci components stand in for the Fermi-frame derivatives (both are O(1/eps)
    per derivative here), and the unknown curvature-squared couplings (24|a| + 11|b|) are kept
    as a separate factor."""
    r, k, k1, k2 = fc.wall_profile(eta, 1.0)
    R = riemann_array(r, k, k1, k2)
    eta_m = np.diag([-1.0, 1.0, 1.0, 1.0])
    ric = np.einsum("ac,nabcd->nbd", eta_m, R)
    d1 = np.gradient(ric, r, axis=0)
    d2 = np.gradient(d1, r, axis=0)
    d3 = np.gradient(d2, r, axis=0)
    r_c = float(np.max(np.abs(R))) ** -0.5
    t0 = sigma * r_c
    rm, r2, r3 = (float(np.max(np.abs(a))) for a in (ric, d2, d3))
    return {
        "eta": eta,
        "t0_over_eps": t0,
        "ricci_term_rel": 2.63 * rm * t0**2 / 3.76,
        "second_derivative_term_rel": 3.42 * r2 * t0**4 / 3.76,
        "second_derivative_term_rel_per_unit_24a_11b": 197.9 * r2 * t0**4 / 3.76,
        "third_derivative_term_rel": 6.99 * r3 * t0**5 / 3.76,
    }


def wall_energy_integrals(eta: float) -> dict:
    """Integrals of the static energy density across the wall at eps = 1 (they scale as 1/eps):
    I_neg = int min(T, 0) d rho and I_pos = int max(T, 0) d rho. They replace ER's guess
    |T| * alpha * eps for the energy per unit wall area (thin wall, rho ~ rho_max)."""
    r, k, k1, k2 = fc.wall_profile(eta, 1.0)
    tt = np.asarray(fc._numeric()[0](r, k, k1, k2), dtype=float)
    dr = r[1] - r[0]
    return {
        "I_neg": float(np.sum(np.minimum(tt, 0.0)) * dr),
        "I_pos": float(np.sum(np.maximum(tt, 0.0)) * dr),
        "effective_alpha": float(-np.sum(np.minimum(tt, 0.0)) * dr / -tt.min()),
        # proper volume on t = const (sqrt(det h) = sqrt(k) rho): a sensitivity check, since the
        # t-slices are not orthogonal to the static observers
        "I_neg_proper_volume": float(
            np.sum(np.minimum(tt, 0.0) * np.sqrt(np.clip(k, 0, None))) * dr
        ),
    }


@lru_cache(maxsize=1)
def null_ricci():
    """R_ab l^a l^b for the two null directions of the (t, x) block, from the coordinate metric
    (ER Eq. 5 at mid-tube, k = k(rho)), with l^t = 1. Outgoing l = (1, 1): exterior light speed,
    also the fastest +x direction. Incoming l = (1, -1/k): the shortcut direction. Returns the
    two expressions and, for the outgoing curves at fixed rho, the radial geodesic acceleration."""
    t, x, r, ph = sp.symbols("t x rho phi", real=True)
    k = sp.Function("k")(r)
    xs = [t, x, r, ph]
    g = sp.Matrix([[-1, (1 - k) / 2, 0, 0], [(1 - k) / 2, k, 0, 0], [0, 0, 1, 0], [0, 0, 0, r**2]])
    gi = sp.simplify(g.inv())
    gam = [
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
                        for d in range(4)
                    )
                    / 2
                )
                for c in range(4)
            ]
            for b in range(4)
        ]
        for a in range(4)
    ]

    def ric(b, c):
        return sp.simplify(
            sum(
                sp.diff(gam[a][b][c], xs[a])
                - sp.diff(gam[a][b][a], xs[c])
                + sum(gam[a][a][d] * gam[d][b][c] - gam[a][c][d] * gam[d][b][a] for d in range(4))
                for a in range(4)
            )
        )

    ricci = sp.Matrix(4, 4, lambda b, c: ric(b, c))
    out_l = sp.Matrix([1, 1, 0, 0])
    in_l = sp.Matrix([1, -1 / k, 0, 0])
    r_out = sp.simplify((out_l.T * ricci * out_l)[0])
    r_in = sp.simplify((in_l.T * ricci * in_l)[0])
    pairs = itertools.product(range(4), repeat=2)
    acc_out = sp.simplify(sum(gam[2][a][b] * out_l[a] * out_l[b] for a, b in pairs))
    k0, k1, k2 = sp.symbols("k0 k1 k2")
    sub = {sp.Derivative(k, (r, 2)): k2, sp.Derivative(k, r): k1}
    f_in = sp.lambdify((r, k0, k1, k2), r_in.subs(sub).subs(k, k0), "numpy")
    return r_out, r_in, acc_out, f_in, sp.simplify(gi[0, 0])


def null_and_causal(eta: float) -> dict:
    """NEC along both null directions on the wall, and the causal structure of the t-slices."""
    r_out, _, acc_out, f_in, gtt_inv = null_ricci()
    r, k, k1, k2 = fc.wall_profile(eta, 1.0)
    rin = np.asarray(f_in(r, k, k1, k2), dtype=float)
    k_min = float(np.min(k))
    return {
        "outgoing_Rll_symbolic": str(r_out),
        "outgoing_fixed_rho_radial_acceleration": str(acc_out),
        "incoming_Rll_min_eps1": float(rin.min()),
        "incoming_Rll_max_eps1": float(rin.max()),
        "g_inv_tt": str(gtt_inv),
        "k_min": k_min,
        "t_is_time_function": bool(k_min > 0),
        "max_incoming_coordinate_speed": float(1 / k_min) if k_min > 0 else None,
    }


POLARISATION_COEFF = 1 / (2880 * math.pi**2)  # one-loop scale, e.g. the conformal trace anomaly


def self_consistency(eta: float, n_fields: float = 1.0) -> dict:
    """Order-of-magnitude [estimate] of semiclassical self-consistency. Required density
    |T| = a_t / eps^2 (Einstein equations); one-loop vacuum polarisation of N fields in a
    curvature R ~ a_r / eps^2 is ~ N c R^2 with c = 1 / (2880 pi^2) (ASSUMPTION: typical one-loop
    coefficient; the sign and the exact tensor structure are not computed). Equality gives
    eps_sc = sqrt(N c a_r^2 / a_t) Planck lengths. A static negative density cannot come from an
    excited state alone: along a static worldline T is constant, and the flat-space QI with
    unbounded tau0 forces T >= 0, so the curvature itself must produce it."""
    lim = fc.qi_wall_limit(eta, P.sigma)
    a_t, a_r = lim["T_coeff"], lim["R_coeff"]
    eps_sc = math.sqrt(n_fields * POLARISATION_COEFF * a_r**2 / a_t)
    # fields needed for the polarisation to reach the QI-allowed thickness
    n_needed = (lim["eps_max_planck"] / math.sqrt(POLARISATION_COEFF * a_r**2 / a_t)) ** 2
    return {
        "eta": eta,
        "T_coeff": a_t,
        "R_coeff": a_r,
        "eps_self_consistent_planck_N1": eps_sc,
        "eps_qi_max_planck": lim["eps_max_planck"],
        "fields_needed_to_reach_qi_thickness": n_needed,
        "qi_over_polarisation_factor_sigma^-4": P.sigma**-4,
    }


def latency_rows() -> list:
    """Latency L (send to read) over D through one tube: L = D (1 - eta), so eta = 1 - L / D."""
    rows = []
    for d_ly in P.distances_ly:
        d_s = d_ly * fc.LY_M / fc.C_M_S
        for lat in P.latencies_s:
            eta = 1 - lat / d_s
            lim = fc.qi_wall_limit(eta, P.sigma)
            eps = lim["eps_max_planck"]
            rows.append(
                {
                    "D_ly": d_ly,
                    "latency_s": lat,
                    "advance_over_light_s": eta * d_s,
                    "eta": eta,
                    "rc_over_eps": lim["rc_over_eps"],
                    "eps_max_planck_tau_sigma_rc": eps,
                    "eps_max_planck_tau_sigma_min_rc_eps": min(
                        eps,
                        lim["eps_max_planck_shortcut_tau0_sigma_eps"]
                        if lim["eps_max_planck_shortcut_tau0_sigma_eps"] is not None
                        else math.inf,
                    ),
                    "curvature_radius_planck_at_eps_max": eps * lim["rc_over_eps"],
                    "wall_integrals_eps1": (wi := wall_energy_integrals(eta)),
                    "negative_energy_solar_integrated": {
                        f"rho_{rad:g}m": negative_energy(-wi["I_neg"], eps, d_ly, rad, 1.0)
                        for rad in P.radii_m
                    },
                    "negative_energy_solar_ER_band": {
                        f"rho_{rad:g}m_alpha_{al:g}": negative_energy(
                            lim["T_coeff"], eps, d_ly, rad, al
                        )
                        for rad in P.radii_m
                        for al in P.alphas
                    },
                }
            )
    return rows


def negative_energy(t_coeff: float, eps: float, d_ly: float, rho_m: float, alpha: float) -> float:
    """ER Eq. (51) generalised: |T| (2 pi rho_max) (alpha eps) D in solar masses (negative). With
    t_coeff = -I_neg and alpha = 1 it is the integrated negative energy of the wall instead."""
    rho = rho_m / fc.PLANCK_LENGTH_M
    d = d_ly * fc.LY_M / fc.PLANCK_LENGTH_M
    e = t_coeff / eps**2 * 2 * math.pi * rho * alpha * eps * d
    return -e * fc.PLANCK_MASS_KG / fc.M_SUN_KG


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived from the ER metric and the flat-space quantum inequality (no data); "
        "estimates under the stated ASSUMPTIONs",
        "params": asdict(P),
        "weyl_vs_ricci": [weyl_ricci(e) for e in P.weyl_etas],
        "latency_rows": latency_rows(),
        "self_consistency": [self_consistency(e) for e in (1 - 1e-9, 0.95, 1e-2)],
        "kontou_olum_corrections": [
            kontou_olum_terms(e) for e in (1.6e-9, 1e-4, 1e-2, 0.95, 1 - 1e-9)
        ],
        "null_and_causal": [null_and_causal(e) for e in (1e-2, 0.95, 1 - 1e-9)],
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
