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
from dataclasses import asdict, dataclass
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import quad

TAU = 100
W_PARAM = 2
REDUCED_PLANCK_LENGTH_M = 8.1e-35  # sqrt(8 pi G hbar / c^3), rounded
C_M_S = 299_792_458.0

x, L = sp.symbols("x L", real=True)


@dataclass(frozen=True)
class Params:
    """Numerical choices (ASSUMPTIONs); every headline number depends on them."""

    x_max: float = 10.0  # grid half-width in units of tau
    n_grid: int = 4001
    throat_gap: float = 1e-6  # |r/tau| excluded at the throat, where J' = 0
    inner_zone: float = 3.0  # |r/tau| below which a point counts as "near the throat"
    superluminal_margin: float = 1e-9  # c^2 > 1 + margin counts as superluminal
    high_l: int = 10**4  # l used for the l >> 1 angular speeds and the Eq. (33) limit
    control_l: int = 10  # l of the Fig. 2 positivity control (their choice)
    scale_tau: int = 10**4  # second throat size for the scale check
    low_ls: tuple = (2, 3, 4, 5, 6, 8, 10, 20, 50)


P = Params()


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
        "cr1_sq": 1 / F, "cr2_sq": cs2,  # Eq. (31): radial speeds, G / F and c_s2 (Eq. 83)
        "detK_eq22": L * B / (4 * (L - 2)) * two_jj**2 * (2 * P1 - F) / F,
        "detM_eq42": -B * h**2 / (4 * J**2 * F**2)
        * (J**4 / sb * d(sb / J**3 * P4) + 2 / B + (F * P4 + 2 * jp) ** 2 / (4 * F)),  # / l^4
    }  # fmt: skip


@lru_cache(maxsize=4)
def numeric(tau_value: int = TAU) -> dict:
    return {k: sp.lambdify((x, L), v, "numpy") for k, v in expressions(tau_value).items()}


def grid(n: int = P.n_grid, xmax: float = P.x_max) -> np.ndarray:
    xs = np.linspace(-xmax, xmax, n)
    return xs[np.abs(xs) > P.throat_gap]  # J' = 0 at the throat: coefficients evaluated off it


def ev(name: str, xs, ell: int, tau_value: int = TAU):
    vals = numeric(tau_value)[name](xs, float(ell * (ell + 1)))
    return np.asarray(vals, dtype=float) * np.ones_like(xs)


def segments(xs, mask) -> list:
    """Contiguous runs of True as [first x, last x, points]."""
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if start is not None and (not m or i == len(mask) - 1):
            end = i if m else i - 1
            out.append([float(xs[start]), float(xs[end]), end - start + 1])
            start = None
    return out


def positive_controls(xs, tau_value: int = TAU) -> dict:
    """Claims of arXiv:2212.05969 Sec. 4 that this module must reproduce, plus internal checks of
    2404.06297 (its Eq. 22 against Appendix B, its Eq. 42 against Eqs. 41, 132, 133)."""
    lc = P.control_l

    def e(name, ell=lc):
        return ev(name, xs, ell, tau_value)

    out = {"F_min": float(e("F").min()), "F_max": float(e("F").max())}
    for k in ("K11", "G11", "detG", "m11"):
        out[f"{k}_min_l{lc}"] = float(e(k).min())
    detk = e("K11") * e("K22") - e("K12") ** 2
    detm = e("m11") * e("m22") - e("m12") ** 2
    out[f"detK_min_l{lc}"] = float(detk.min())
    out["detM_l2_min"] = float(detm.min())
    out["no_ghost_2P1_minus_F_min"] = float(e("P1m").min())
    out["N_min"] = float(e("N").min())
    out["radial_c2_tensor_range"] = [float(e("cr1_sq").min()), float(e("cr1_sq").max())]
    out["radial_c2_scalar_range"] = [float(e("cr2_sq").min()), float(e("cr2_sq").max())]
    ratio_k = detk / e("detK_eq22")
    ratio_m = detm / e("detM_eq42")
    out["detK_appendixB_over_eq22_range"] = [float(ratio_k.min()), float(ratio_k.max())]
    out["detM_appendixB_over_eq42_range"] = [float(ratio_m.min()), float(ratio_m.max())]
    big = P.high_l
    lead = [
        float(np.max(np.abs(e(a, big) / (big * (big + 1)) / e(b, big) - 1)))
        for a, b in (("M11", "m11"), ("M12", "m12"), ("M22", "m22"))
    ]
    out[f"eq33_vs_appendixB_max_rel_diff_l{big}"] = lead
    return out


def angular_speeds(xs, ell: int = P.high_l, tau_value: int = TAU) -> dict:
    """c_a^2 = eigenvalues of (J^2 / A) K^-1 M^(l^2) / l^2 for l >> 1 (2404.06297 below Eq. 44)."""
    k = np.stack(
        [
            np.stack([ev("K11", xs, ell, tau_value), ev("K12", xs, ell, tau_value)], -1),
            np.stack([ev("K12", xs, ell, tau_value), ev("K22", xs, ell, tau_value)], -1),
        ],
        -2,
    )
    m = np.stack(
        [
            np.stack([ev("m11", xs, ell, tau_value), ev("m12", xs, ell, tau_value)], -1),
            np.stack([ev("m12", xs, ell, tau_value), ev("m22", xs, ell, tau_value)], -1),
        ],
        -2,
    )
    j2 = ev("J", xs, ell, tau_value) ** 2
    lam = np.linalg.eigvals(j2[:, None, None] * np.linalg.solve(k, m))
    c2 = np.sort(lam.real, axis=1)
    i_max = int(np.argmax(c2[:, 1]))
    return {
        "max_c2": float(c2[i_max, 1]),
        "r_over_tau_at_max": float(xs[i_max]),
        "fraction_of_grid_superluminal": float(np.mean(c2[:, 1] > 1 + P.superluminal_margin)),
        "slower_mode_c2_range": [float(c2[:, 0].min()), float(c2[:, 0].max())],
        "edge_values_c2": [float(c2[0, 0]), float(c2[0, 1])],
        "max_abs_imag": float(np.abs(lam.imag).max()),  # complex c^2 would signal an instability
        "non_finite_points": int((~np.isfinite(lam)).any(axis=1).sum()),
    }


def low_l_conditions(xs, ells=P.low_ls, tau_value: int = TAU) -> list:
    """2404.06297 Eqs. (39b), (39c) and (6d) for each l on the grid: where they fail. Non-finite
    values (poles at Theta = 0) are counted separately, never as passes."""
    rows = []
    inner = np.abs(xs) < P.inner_zone
    for ell in ells:
        m22 = ev("M22", xs, ell, tau_value)
        c39 = ev("cond39c", xs, ell, tau_value)
        theta = ev("Theta", xs, ell, tau_value)
        odd = ev("odd_tachyon", xs, ell, tau_value)
        finite = np.isfinite(m22) & np.isfinite(c39) & np.isfinite(odd)
        bad = ~(c39 > 0)
        sign_change = np.diff(np.sign(theta)) != 0
        rows.append(
            {
                "l": ell,
                "M22_le_0_segments": segments(xs, ~(m22 > 0)),
                "cond39c_fail_segments": segments(xs, bad),
                "cond39c_fail_fraction_inner": float(bad[inner].mean()),
                "cond39c_fail_fraction_outer": float(bad[~inner].mean()),
                "theta_sign_changes_at": [float(v) for v in xs[1:][sign_change]],
                "odd_tachyon_fail_segments": segments(xs, ~(odd > 0)),
                "non_finite_points": int((~finite).sum()),
                "grid_points": int(len(xs)),
            }
        )
    return rows


def mp_crosscheck(points=((2, 0.2), (2, 3.0), (2, 8.0), (10, 4.0)), dps: int = 40) -> list:
    """Evaluate the low-l conditions at a few points with mpmath (dps digits) against float64."""
    import mpmath as mp

    mp.mp.dps = dps
    exprs = expressions()
    out = []
    for ell, xv in points:
        lv = ell * (ell + 1)
        row = {"l": ell, "x": xv}
        for name in ("M22", "cond39c"):
            hi = sp.lambdify((x, L), exprs[name], "mpmath")(mp.mpf(xv), mp.mpf(lv))
            lo = float(ev(name, np.array([xv]), ell)[0])
            row[name] = [float(hi), lo, abs(float(hi) - lo) / max(abs(float(hi)), 1e-300)]
        out.append(row)
    return out


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
    ap.add_argument("--n", type=int, default=P.n_grid)
    ap.add_argument("--skip-mp", action="store_true", help="skip the mpmath cross-check (~1 min)")
    a = ap.parse_args(argv)
    xs = grid(a.n)
    out = {
        "provenance": "derived from the published solution and stability set (no data); "
        "sufficient-condition checks, not proofs of instability",
        "model": {"tau": TAU, "w": W_PARAM, "u": 0, "source": "arXiv:2212.05969v1 Sec. 4"},
        "params": asdict(P),
        "positive_controls": positive_controls(xs),
        "angular_speeds_high_l": angular_speeds(xs),
        "low_l_conditions": low_l_conditions(xs),
        "timing": crossing_time(),
        "scale_check": {
            "tau": P.scale_tau,
            "angular_speeds_high_l": angular_speeds(xs, tau_value=P.scale_tau),
            "low_l_conditions": low_l_conditions(xs, tau_value=P.scale_tau),
        },
    }
    if not a.skip_mp:
        out["mpmath_crosscheck"] = mp_crosscheck()
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
