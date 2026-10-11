"""Does the beyond-Horndeski wormhole of Mironov, Rubakov & Volkova (2023) meet the complete
stability set of Mironov & Volkova (2024), and what channel would it give?

Goal round 4 (docs/goal_round4/README.md). One fixed model, no data:
- Theory: quadratic beyond Horndeski, action (1) of arXiv:2212.05969 (PRD 107, 104061).
- Solution (their Sec. 4, Eqs. 70, 72, 73, 80, 81, 83): A = 1, B = 1 + sech(r/tau),
  J = tau ln(1 + 2 cosh(r/tau)), pi = tanh(r/tau) - 1, with the reconstructed Lagrangian chosen
  so that H = G = 1, Xi = 0, Gamma = 0, H - 2 F4 B^2 pi'^4 = 1 - w sech(r/tau) (u = 0),
  c_s2^2 = cosh(r/tau) / (1 + cosh(r/tau)); tau = 100, w = 2. F follows from their Eq. (71)
  with A = 1: F = J^2 [B' J' / (2 J) + B (J'^2 / J^2 + J'' / J)].
- Stability set: arXiv:2404.06297v2 (Mironov & Volkova, IJMPA), Eqs. (5)-(8), (21)-(45) and
  Appendix B, Eqs. (129)-(133). It is a SUFFICIENT set; failing it is not a proof of instability.

Because Xi = Gamma = 0 and H = G = 1 everywhere, every quantity in that set is fixed by the metric
and the functions above, so the reconstructed Lagrangian functions are not needed.

Computed here [derived]:
1. Positive controls: F > 1 (their claim below Eq. 72); the Fig. 2 positivity of K11, det K, G11,
   det G, M11^(l^2), det M^(l^2); the leading order in l of Eq. (33) equals Appendix B.
2. Angular propagation speeds of the two even-parity modes for l >> 1, from the full matrices
   K (Eqs. 22, 129, 130) and M^(l^2) (Eqs. 41, 132, 133).
3. The low-l conditions (39b) M22 > 0 and (39c) det G det M > Q12^2 M22 G11 / 4 (Eqs. 33, 37),
   and the odd-parity tachyon condition (6d).
4. Event timing: the light-crossing time between the two asymptotic regions.

Units: reduced Planck units (G4 -> 1/2 asymptotically, so M_Pl = 1; an ASSUMPTION read off the
action normalisation), radial coordinate x = r / tau.
"""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import quad

TAU = 100
W_PARAM = 2
REDUCED_PLANCK_LENGTH_M = 8.1e-35  # sqrt(8 pi G hbar / c^3), rounded
C_M_S = 299_792_458.0

x, L = sp.symbols("x L", real=True)


@lru_cache(maxsize=4)
def expressions(tau_value: int = TAU) -> dict:
    """Every quantity of the stability set as a sympy expression in x = r / tau and L = l(l+1)."""
    tau, w = sp.Integer(tau_value), sp.Integer(W_PARAM)
    sech = 1 / sp.cosh(x)

    def d(f):
        return sp.diff(f, x) / tau

    B = 1 + sech
    J = tau * sp.log(1 + 2 * sp.cosh(x))
    # H = G = 1 throughout (their choice below Eq. 71), so they are dropped from the formulas
    jp, jpp, bp = d(J), d(d(J)), d(B)
    sb = sp.sqrt(B)
    F = J**2 * (bp * jp / (2 * J) + B * (jp**2 / J**2 + jpp / J))  # Eq. (71), A = 1
    h = 1 - w * sech  # H - 2 F4 B^2 pi'^4, Eq. (80) with u = 0
    xi = J * h / (2 * jp)  # Eq. (24c) with Xi = 0, H = 1
    P1 = sb * d(xi / sb)
    cs2 = sp.cosh(x) / (1 + sp.cosh(x))  # Eq. (83)
    two_jj = 2 * J * jp
    N = cs2 * two_jj**2 * (2 * P1 - F)
    Theta = 2 * L * J * h - 2 * B * jp * two_jj  # Eq. (24a), A = 1, Xi = 0
    T = 2 * jp + J / (L - 2) * d(Theta / (F * J))  # Eq. (35)
    D = d(B * J**2) / (J**2 * F)  # Eq. (36)
    dp = d(D)
    W = 6 * B * J * jp**2 - 2 * J * (2 * L - 1) * h  # Eq. (34), Gamma1 = 0
    brk = 2 * L * J * jp - two_jj * (2 + J**2 * dp)
    M11 = (
        sb * N / (4 * J**2)
        + ((L - 2) * sb * T**2 - d(sb * Theta * T / F)) / (4 * L)
        - Theta * (Theta + W) / (4 * sb * J**2 * (L - 2) * F)
    )  # Eq. (33a)
    M12 = (
        Theta / (F * J) * d(sb * (2 + J**2 * dp))
        + (L * J * (D * W + 2 * B * Theta * jp / (F * J)) + 2 * (L - 2) * B * brk) / (2 * sb * J**2)
    ) / (2 * (L - 2))  # Eq. (33b)
    M22 = -(L / (L - 2)) * (
        sb / J**2 * (L * (L - 2) - (2 + J**2 * dp - L) ** 2)
        + d(sb / (J * Theta) * (2 * (L - 2) * B * brk + J * D * (L * W + Theta * (2 + J**2 * dp))))
    )  # Eq. (33c)
    G11 = sb * Theta**2 / (4 * L * (L - 2) * F**2)  # Eq. (27)
    detG = L * B**3 * N / ((L - 2) * F**2)
    b32 = B ** sp.Rational(3, 2)
    Q12 = (
        2 * d(b32 * J / F) * T
        + b32 * J / F * d(T)
        - sb / (2 * J * (L - 2)) / F * (L * W + Theta * (2 + J**2 * dp))
        - 4 * L * B * J / Theta * M11
    )  # Eq. (37)
    detM = M11 * M22 - M12**2
    # high-l matrices (Appendix B) for the angular speeds
    P4 = 2 * jp / h  # Eq. (43), Gamma = 0, A' = 0
    K11 = Theta**2 / (4 * L * (L - 2) * sb * F)
    K12 = sb / (2 * (L - 2)) * (2 * (L - 2) * two_jj - d(B * J**2) / (B * J) * Theta / F)
    K22 = 4 * L * sb / Theta**2 * (L * B * (2 * P1 - F) * two_jj**2 + (L - 2) * F * K12**2)
    m11 = h**2 / (sb * F)  # Eq. (41) / l^2
    m12 = -h / (2 * F * J * sb) * (2 * d(B * J**2) / J - 2 * B * jp - B * F * P4)  # (132) / l^2
    m22 = (
        -(1 / F)
        * sb
        / J**2
        * (2 * F + sb * J**4 * F * d(sb / J**3 * P4) + (d(B * J**2) / J) ** 2 / B)
        - 2 / h * d(B * J**2) / J**2 * m12
    )  # (133) / l^2
    # odd parity, Eq. (5) with A = 1, H = 1, and condition (6d)
    fp = d(F)
    V = (
        B
        / (2 * F)
        * (
            fp / F * (bp / B + 4 * jp / J)
            - (d(bp) / B + 4 * bp * jp / (B * J) - 4 * jp**2 / J**2 + 4 * jpp / J)
            - 4 * F / (J**2 * B)
        )
    )
    return {
        "J": J, "B": B, "F": F, "h": h, "P1m": 2 * P1 - F, "N": N, "Theta": Theta,
        "M11": M11, "M12": M12, "M22": M22, "detM": detM, "G11": G11, "detG": detG, "Q12": Q12,
        "cond39c": detG * detM - Q12**2 / 4 * M22 * G11,
        "K11": K11, "K12": K12, "K22": K22, "m11": m11, "m12": m12, "m22": m22,
        "odd_tachyon": V + 6 / J**2,
    }  # fmt: skip


@lru_cache(maxsize=4)
def numeric(tau_value: int = TAU) -> dict:
    return {k: sp.lambdify((x, L), v, "numpy") for k, v in expressions(tau_value).items()}


def grid(n: int = 4001, xmax: float = 10.0) -> np.ndarray:
    xs = np.linspace(-xmax, xmax, n)
    return xs[np.abs(xs) > 1e-6]  # J' = 0 at the throat: coefficients are evaluated off it


def ev(name: str, xs, ell: int, tau_value: int = TAU):
    vals = numeric(tau_value)[name](xs, float(ell * (ell + 1)))
    return np.asarray(vals, dtype=float) * np.ones_like(xs)


def positive_controls(xs) -> dict:
    """Claims of arXiv:2212.05969 Sec. 4 that this module must reproduce."""
    out = {"F_min": float(ev("F", xs, 10).min())}
    for k in ("K11", "G11", "detG", "m11"):
        out[f"{k}_min_l10"] = float(ev(k, xs, 10).min())
    detk = ev("K11", xs, 10) * ev("K22", xs, 10) - ev("K12", xs, 10) ** 2
    detm = ev("m11", xs, 10) * ev("m22", xs, 10) - ev("m12", xs, 10) ** 2
    out["detK_min_l10"] = float(detk.min())
    out["detM_l2_min"] = float(detm.min())
    out["no_ghost_2P1_minus_F_min"] = float(ev("P1m", xs, 10).min())
    out["N_min"] = float(ev("N", xs, 10).min())
    # leading order of Eq. (33) equals Appendix B (Eqs. 41, 132, 133)
    big = 10**4
    lead = [
        float(np.max(np.abs(ev(a, xs, big) / (big * (big + 1)) / ev(b, xs, big) - 1)))
        for a, b in (("M11", "m11"), ("M12", "m12"), ("M22", "m22"))
    ]
    out["eq33_vs_appendixB_max_rel_diff_l1e4"] = lead
    return out


def angular_speeds(xs, ell: int = 10**4, tau_value: int = TAU) -> dict:
    """c_a^2 = eigenvalues of (J^2 / A) K^-1 M^(l^2) / l^2 for l >> 1 (2404.06297 below Eq. 44)."""
    k11, k12, k22 = (ev(n, xs, ell, tau_value) for n in ("K11", "K12", "K22"))
    m11, m12, m22 = (ev(n, xs, ell, tau_value) for n in ("m11", "m12", "m22"))
    j2 = ev("J", xs, ell, tau_value) ** 2
    c2 = []
    for i in range(len(xs)):
        k = np.array([[k11[i], k12[i]], [k12[i], k22[i]]])
        m = np.array([[m11[i], m12[i]], [m12[i], m22[i]]])
        c2.append(np.sort(np.linalg.eigvals(j2[i] * np.linalg.solve(k, m)).real))
    c2 = np.array(c2)
    i_max = int(np.argmax(c2[:, 1]))
    return {
        "max_c2": float(c2[i_max, 1]),
        "r_over_tau_at_max": float(xs[i_max]),
        "fraction_of_grid_superluminal": float(np.mean(c2[:, 1] > 1 + 1e-9)),
        "min_c2": float(c2[:, 0].min()),
        "edge_values_c2": [float(c2[0, 0]), float(c2[0, 1])],
    }


def low_l_conditions(xs, ells=(2, 3, 4, 5, 6, 8, 10, 20, 50)) -> list:
    """2404.06297 Eqs. (39b), (39c) and (6d) for each l on the grid: where they fail."""
    rows = []
    for ell in ells:
        m22 = ev("M22", xs, ell)
        c39 = ev("cond39c", xs, ell)
        theta = ev("Theta", xs, ell)
        odd = ev("odd_tachyon", xs, ell)

        def span(mask):
            if not mask.any():
                return None
            return [float(xs[mask].min()), float(xs[mask].max()), int(mask.sum())]

        inner = np.abs(xs) < 3.0
        bad = c39 <= 0
        rows.append(
            {
                "l": ell,
                "M22_le_0": span(m22 <= 0),
                "cond39c_fail_fraction_inner_|r|<3tau": float(bad[inner].mean()),
                "cond39c_fail_fraction_outer": float(bad[~inner].mean()),
                "theta_zero_near": [float(v) for v in xs[1:][np.diff(np.sign(theta)) != 0]],
                "odd_tachyon_le_0": span(odd <= 0),
                "grid_points": int(len(xs)),
            }
        )
    return rows


def crossing_time() -> dict:
    """Radial light from r = -R to r = +R through the throat takes 2R - tau * I (A = 1, so Killing
    time is proper time of static clocks at both ends), I = int (1 - 1/sqrt(1 + sech x)) dx."""
    val, err = quad(lambda s: 1 - 1 / math.sqrt(1 + 1 / math.cosh(min(abs(s), 700))), -50, 50)
    throat = TAU * math.log(3.0)
    return {
        "I": val,
        "quad_error": err,
        "throat_radius_planck": throat,
        "throat_radius_m": throat * REDUCED_PLANCK_LENGTH_M,
        "crossing_deficit_planck": TAU * val,
        "crossing_deficit_s": TAU * val * REDUCED_PLANCK_LENGTH_M / C_M_S,
        "clock_rate_ratio_between_ends": 1.0,  # A = 1 everywhere
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=4001)
    a = ap.parse_args(argv)
    xs = grid(a.n)
    out = {
        "provenance": "derived from the published solution and stability set (no data); "
        "sufficient-condition checks, not proofs of instability",
        "model": {"tau": TAU, "w": W_PARAM, "u": 0, "source": "arXiv:2212.05969v1 Sec. 4"},
        "positive_controls": positive_controls(xs),
        "angular_speeds_high_l": angular_speeds(xs),
        "low_l_conditions": low_l_conditions(xs),
        "timing": crossing_time(),
        "scale_check_tau_1e4": angular_speeds(grid(401), tau_value=10**4),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
