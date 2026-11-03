"""FGM wormhole opening window versus finite quantum wave packets (analytic round, 2026-10-29).

Fu, Grado-White & Marolf, arXiv:1908.03273v2 (FGM). Units: hbar = c = 1, lengths in units of r+,
G in units of r+^2, i.e. G = (l_P / r+)^2. Kruskal normalisation of FGM (3.3) with x*0 = d/2 (2.21).
Frozen plan: docs/results/fgm_wavepacket_window.md section 2 (committed before this script).

Window (FGM (3.8)-(3.9), (3.13)-(3.14), (2.20)): an ingoing ray of advanced time v at mouth A
reaches mouth B iff V = exp(kappa v)/kappa < |DeltaV|, with
  DeltaV = C_RN * 8 pi G * (H(theta)/r+^2) * Int T_UU dU,
  Int T_UU dU = -exp(-kappa x*0) c kappa / 16.
Closure by positive energy falling into mouth B at retarded label u_B (der from (2.20), (3.3)):
  dE_crit(u_B) = (c kappa / 16) exp(kappa (u_B - x*0)).

  python scripts/fgm_window.py   # -> results/fgm_window/window.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "fgm_window"

L_P = 1.616255e-35  # m
HBAR_C_EV_M = 1.973269804e-7  # eV m
C_LIGHT = 299_792_458.0  # m/s

# frozen parameter point (plan section 2)
R_PLUS_M = 1.0e3
R_MINUS = 0.5  # in units of r+
D = 1.0e4  # mouth separation in units of r+
N_STRINGS = 10
G_MU = 1e-7
MARGIN = 100.0  # packet centre placed this far (in r+) before the window edge
F_GEO = 1.0 / D**2  # fraction of reflected energy reaching mouth B, (r+/d)^2 [assumption]


def kappa(r_minus: float = R_MINUS) -> float:
    """Surface gravity of the outer horizon of RN with r+ = 1."""
    return (1.0 - r_minus) / 2.0


def c_rn(r_minus: float = R_MINUS) -> float:
    """FGM (3.9) prefactor: (r+/r-) ((r+ - r-)/r-)^(-1-(r-/r+)^2) exp(2 kappa r+)."""
    a = r_minus**2
    return (
        (1.0 / r_minus) * ((1.0 - r_minus) / r_minus) ** (-1.0 - a) * math.exp(2 * kappa(r_minus))
    )


def green_h(theta: float, r_minus: float = R_MINUS) -> float:
    """FGM (3.14): H(theta)/r+^2 = -P_lambda(-cos theta) / (2 sin(pi lambda))."""
    lam = -0.5 * (1 + mp.sqrt(1 - 8 * kappa(r_minus)))
    val = -mp.legenp(lam, 0, -mp.cos(theta), type=2) / (2 * mp.sin(mp.pi * lam))
    return float(mp.re(val))


def green_h_series(theta: float, r_minus: float = R_MINUS, jmax: int = 4000) -> float:
    """FGM (3.15) summed: sum_j (2j+1)/(4 pi) P_j(cos theta) * 2 / (2 kappa + j(j+1))."""
    x = math.cos(theta)
    p0, p1, tot = 1.0, x, 0.0
    k2 = 2 * kappa(r_minus)
    for j in range(jmax + 1):
        pj = p0 if j == 0 else p1
        tot += (2 * j + 1) / (4 * math.pi) * pj * 2.0 / (k2 + j * (j + 1))
        if j >= 1:
            p0, p1 = p1, ((2 * j + 1) * x * p1 - j * p0) / (j + 1)
    return tot


def window(theta: float, g_tilde: float, c: float, x0: float = D / 2, r_minus: float = R_MINUS):
    """|DeltaV|, v_max and the minimal transit time (FGM (3.10)) in units of r+."""
    k = kappa(r_minus)
    q = k**2 * c_rn(r_minus) * math.pi * g_tilde * green_h(theta, r_minus) * c / 2  # k|DV|e^{k x0}
    v_max = -x0 + math.log(q) / k
    t_min = 2 * x0 - (2 / k) * math.log(q / 2)
    return {"Q": q, "v_max": v_max, "t_min": t_min}


def t_min_fgm_317(theta: float, g_tilde: float, c: float, x0: float = D / 2, r_minus=R_MINUS):
    """FGM (3.17) with H read as H/r+^2 (dimensionally consistent form)."""
    k, a = kappa(r_minus), r_minus**2
    arg = (math.pi * g_tilde * c / 16) * r_minus * ((1 - r_minus) / r_minus) ** (1 - a)
    return 2 * x0 - 4 - (2 / k) * math.log(arg * green_h(theta, r_minus))


def packet_gaussian(w0: float, sw: float) -> dict:
    """One quantum, |f|^2 Gaussian (mean w0, sd sw) truncated to w > 0; analytic moments."""
    z = w0 / sw
    norm = 0.5 * math.erfc(-z / math.sqrt(2))
    sub_w = 0.5 * math.erfc((w0 - 1.0) / (sw * math.sqrt(2))) - 0.5 * math.erfc(z / math.sqrt(2))
    # energy below w = 1 (sub-barrier): int_0^1 w |f|^2 dw, |f|^2 ~ exp(-(w-w0)^2/2sw^2)/sqrt(2pi)sw
    w = np.linspace(0.0, 1.0, 2001)
    dens = np.exp(-((w - w0) ** 2) / (2 * sw**2)) / (math.sqrt(2 * math.pi) * sw) / norm
    e_sub = float(np.trapezoid(w * dens, w))
    above5 = 0.5 * math.erfc((5.0 - w0) / (sw * math.sqrt(2))) / norm
    return {"mean_E": w0, "sd_v": 1 / (2 * sw), "weight_below_1": sub_w / norm,
            "E_below_1": e_sub, "weight_above_5": above5}  # fmt: skip


def packet_bump(w0: float, half: float, n: int = 40001) -> dict:
    """One quantum, |f|^2 proportional to exp(-1/(1-x^2)), x = (w - w0)/half, compact support."""
    x = np.linspace(-1, 1, n)[1:-1]
    w = w0 + half * x
    dens = np.exp(-1.0 / (1 - x**2))
    dens /= np.trapezoid(dens, w)
    mean_e = float(np.trapezoid(w * dens, w))
    sd_w = float(np.sqrt(np.trapezoid((w - mean_e) ** 2 * dens, w)))
    # arrival density near the centre: P(t) = |(2pi)^-1/2 int sqrt(dens) e^{-iwt} dw|^2
    amp = np.sqrt(dens)
    t = np.linspace(-0.5, 0.5, 2001)
    ft = np.array([np.trapezoid(amp * np.exp(-1j * (w - w0) * ti), w) for ti in t]) / math.sqrt(
        2 * math.pi
    )
    p_t = np.abs(ft) ** 2
    within = float(np.trapezoid(p_t, t))
    sd_t = float(np.sqrt(np.trapezoid(t**2 * p_t, t) / within))
    return {"mean_E": mean_e, "sd_w": sd_w, "sd_v": sd_t, "prob_within_0.5": within,
            "weight_below_1": 0.0, "E_below_1": 0.0, "weight_above_5": 1.0}  # fmt: skip


def de_crit(u_b_minus_x0: float, c: float, r_minus: float = R_MINUS) -> float:
    """Killing energy falling into mouth B at u_B that cancels the FGM integrated null energy."""
    k = kappa(r_minus)
    return (c * k / 16) * math.exp(k * u_b_minus_x0)


def run() -> dict:
    g_tilde = (L_P / R_PLUS_M) ** 2
    c = 2 * N_STRINGS
    k = kappa()
    x0 = D / 2
    res = {"params": {"r_plus_m": R_PLUS_M, "r_minus_over_r_plus": R_MINUS, "d_over_r_plus": D,
                      "N": N_STRINGS, "c": c, "G_mu": G_MU, "G_over_r_plus2": g_tilde,
                      "kappa_r_plus": k, "margin_r_plus": MARGIN, "f_geo": F_GEO,
                      "t_merge_r_plus": D**1.5},  # fmt: skip
           "green": {}, "window": {}, "packets": {}}  # fmt: skip
    res["validity"] = {"d_over_r_plus": D, "mu_r_plus2": G_MU / g_tilde, "G_mu": G_MU,
                       "N_G_mu": N_STRINGS * G_MU, "kappa_d": k * D}  # fmt: skip
    for name, th in (("pi", math.pi), ("pi_over_2", math.pi / 2)):
        res["green"][name] = {"closed_form": green_h(th), "series_4000": green_h_series(th)}
        w = window(th, g_tilde, c)
        w["t_min_fgm_3_17"] = t_min_fgm_317(th, g_tilde, c)
        w["t_min_over_d"] = w["t_min"] / D
        w["t_min_ms"] = w["t_min"] * R_PLUS_M / C_LIGHT * 1e3
        w["exterior_light_ms"] = D * R_PLUS_M / C_LIGHT * 1e3
        w["v_max_ms"] = w["v_max"] * R_PLUS_M / C_LIGHT * 1e3
        res["window"][name] = w
    wpi = res["window"]["pi"]
    v_c = wpi["v_max"] - MARGIN
    u_b = v_c + D  # arrival of A's reflection at mouth B (flat propagation over d)
    crit_centre = de_crit(u_b - x0, c)
    res["closure_into_B"] = {"u_B_minus_x0": u_b - x0, "dE_crit_centre": crit_centre,
                             "dE_crit_centre_eV": crit_centre * HBAR_C_EV_M / R_PLUS_M}  # fmt: skip

    packets = {"P1_gaussian": packet_gaussian(500.0, 100.0), "P2_bump": packet_bump(500.0, 200.0)}
    for name, p in packets.items():
        e = p["mean_E"]
        p["mean_E_eV"] = e * HBAR_C_EV_M / R_PLUS_M
        p["R1_prob_in_window"] = (
            0.5 * math.erfc(-MARGIN / (p["sd_v"] * math.sqrt(2)))
            if name.startswith("P1")
            else p["prob_within_0.5"]  # lower bound: mass within +-0.5 r+ of the centre
        )
        p["R2_weight_above_5"] = p["weight_above_5"]
        p["R3_G_E_over_r_plus"] = g_tilde * e
        alpha_crit = math.exp(k * MARGIN) - 1
        p["R5_alpha_crit"] = alpha_crit
        p["R5_shape_factor_needed"] = alpha_crit / (g_tilde * e)
        p["R6_E_into_B_subbarrier"] = p["E_below_1"] * F_GEO
        p["R6_ratio_to_crit"] = p["R6_E_into_B_subbarrier"] / crit_centre
        p["R6_reflection_bound_above_barrier"] = crit_centre / (e * F_GEO)
        packets[name] = p
    res["packets"] = packets
    # standard above-barrier reflection estimate at the P2 band edge (assumption, comparison only)
    w_min = 300.0
    res["reflection_estimate"] = {"omega_min_r_plus": w_min,
                                  "ln_R_estimate": -math.pi * w_min / k}  # fmt: skip
    # post hoc: P2 past tail. Saddle point of the bump's Fourier transform gives
    # |F(dt)|^2 ~ exp(-sqrt(2 W dt)); the closure threshold falls as exp(-k dt). Crossover where
    # k dt - sqrt(2 W dt) = ln(dE_crit / (R E f_geo)), with R from the estimate above.
    big_w = 200.0
    rhs = math.log(crit_centre) + math.pi * w_min / k - math.log(500.0 * F_GEO)  # in logs
    s = (math.sqrt(2 * big_w) + math.sqrt(2 * big_w + 4 * k * rhs)) / (2 * k)
    res["post_hoc_P2_tail"] = {"ln_threshold_over_reflected": rhs, "dt_crossover_r_plus": s**2,
                               "inside_background_lifetime": s**2 < D**1.5}  # fmt: skip
    # post hoc: a narrower Gaussian (not frozen) and one signal-energy quantum falling into B
    nb = packet_gaussian(500.0, 25.0)
    res["post_hoc_narrowband_gaussian"] = {
        "omega0_over_sigma": 20.0,
        "E_into_B_subbarrier": nb["E_below_1"] * F_GEO,
        "ratio_to_crit": nb["E_below_1"] * F_GEO / crit_centre,
        "sd_v": nb["sd_v"],
    }
    res["closure_into_B"]["one_signal_quantum_over_crit"] = 500.0 / crit_centre
    # two-way consequences (separate from the verdict)
    v_out_min = -wpi["v_max"]  # earliest emergence label at B
    res["two_way"] = {
        "reverse_window_edge": wpi["v_max"],
        "earliest_receipt": v_out_min,
        "reply_shortfall_r_plus": v_out_min - wpi["v_max"],
        "log10_block_ratio_reverse": (k * (x0 - v_c)) / math.log(10)
        + math.log10(500.0 / (c * k / 16)),
    }
    return res


def main() -> int:
    res = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "window.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
