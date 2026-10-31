"""Can P1's own reflection into mouth B be bounded from below? (analytic follow-up, 2026-10-31).

Follow-up to scripts/fgm_window.py (#619). Same units (hbar = c = 1, lengths in r+), same frozen
FGM parameter point and the same P1 spectrum (Gaussian, omega0 r+ = 500, sigma r+ = 100,
truncated to omega > 0). Record: docs/results/fgm_p1_reflection_bound.md.

Closure functional with FGM-consistent weights [der from FGM (2.20), (3.3), (3.14)]: positive
Killing energy falling into B with flux F(u, Omega') closes the window on the generator theta_s iff
  J = int du e^{-kappa (u - x0)} int dOmega' F(u, Omega') H(gamma(theta_s, Omega')) / H(theta_s)
    > c kappa / 16,
where the strings sit at theta' = 0 and gamma is the angular separation. Jensen's inequality
(e^{-kappa u} is convex) gives J >= E_B * a * exp(-kappa (<u> - x0)) with <u> the energy-weighted
mean arrival time and a the deposit-averaged angular factor.

  python scripts/fgm_p1_reflection.py   # -> results/fgm_p1_reflection/bound.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import fgm_window as fw

ROOT = Path(__file__).resolve().parents[1]
WINDOW_JSON = ROOT / "results" / "fgm_window" / "window.json"
OUT = ROOT / "results" / "fgm_p1_reflection"

W0, SW = 500.0, 100.0  # P1 spectrum (frozen in #619)
BAND = (0.02, 0.1)  # band used for the setting-S bound, in 1/r+ [choice; see the record]
K_T, K_SIGMA = 2.0, 2.0  # allowances on the leading low-frequency s-wave results [assumption]
DELTA_MEAN = 200.0  # allowed energy-weighted mean delay of the reflected flux, in r+ [assumption]
B_C = 3.149  # critical impact parameter of RN with r- = 0.5 r+, in r+ [der, photon sphere]


def gauss_norm(w0: float = W0, sw: float = SW) -> float:
    return 0.5 * math.erfc(-w0 / (sw * math.sqrt(2)))


def density(w: float, w0: float = W0, sw: float = SW) -> float:
    """|f(omega)|^2 per unit omega for the truncated Gaussian, normalised on omega > 0."""
    norm = math.sqrt(2 * math.pi) * sw * gauss_norm(w0, sw)
    return math.exp(-((w - w0) ** 2) / (2 * sw**2)) / norm


def band_energy(a: float, b: float, weight=None, n: int = 4001, w0=W0, sw=SW) -> float:
    """int_a^b omega |f|^2 weight(omega) d omega (Simpson)."""
    h = (b - a) / (n - 1)
    tot = 0.0
    for i in range(n):
        w = a + i * h
        c = 1 if i in (0, n - 1) else (4 if i % 2 else 2)
        tot += c * w * density(w, w0, sw) * (weight(w) if weight else 1.0)
    return tot * h / 3


def h_mean_over_r2(r_minus: float = fw.R_MINUS) -> float:
    """Sphere average of H/r+^2: the j = 0 term of FGM (3.15), 1 / (4 pi kappa)."""
    return 1.0 / (4 * math.pi * fw.kappa(r_minus))


def angular_factors(theta_s: float) -> dict:
    """Weight of a deposit relative to the string's, for the generator theta_s."""
    h_s = fw.green_h(theta_s)
    grid = [math.pi * i / 400 for i in range(1, 401)]
    h_vals = [fw.green_h(g) for g in grid]
    monotone = all(h_vals[i] >= h_vals[i + 1] for i in range(len(h_vals) - 1))
    return {"H_over_r2": h_s, "uniform_deposit": h_mean_over_r2() / h_s,
            "any_deposit_min": fw.green_h(math.pi) / h_s,
            "equatorial_deposit": fw.green_h(math.pi / 2) / h_s,
            "H_decreasing_in_separation": monotone}  # fmt: skip


def double_factorial_log10(n: int) -> float:
    return sum(math.log10(k) for k in range(n, 0, -2))


def vortex_transfer_log10(m: int, d: float = fw.D) -> float:
    """Upper estimate of the A->B transfer of the sub-barrier part for a vortex completion.

    Sectoral reflected pattern Y_{|m| m} (node on the A-B axis), B absorbing only l_B >= |m|
    with Gamma_l <= (omega r+)^(2l+2) [assumption]: T_m <= [(2|m|+1)!!]^2 (r+/d)^(2|m|+2) / 4.
    """
    return 2 * double_factorial_log10(2 * m + 1) - (2 * m + 2) * math.log10(d) - math.log10(4)


def tail_flux_coeff(w0: float = W0, sw: float = SW) -> float:
    """Energy flux of the one-quantum packet far from its centre: F(dv) ~ coeff / |dv|^3.

    The truncation at omega = 0 leaves |f(0)|^2 = p0 > 0; the endpoint term of
    G(v) = (2 pi)^-1/2 int sqrt(omega) f e^{-i omega v} d omega gives |G|^2 = p0 Gamma(3/2)^2 /
    (2 pi |dv|^3) = p0 / (8 |dv|^3).
    """
    return density(0.0, w0, sw) / 8.0


def log10_tail_ratio(lead: float, transfer: float, crit: float, k: float, coeff: float) -> float:
    """log10 of [int_{lead/2}^{lead} transfer coeff s^-3 e^{k s} ds] / crit (dominant upper end)."""
    return (math.log(transfer * coeff) + k * lead - 3 * math.log(lead) - math.log(k)) / math.log(
        10
    ) - math.log10(crit)


def break_even_lead(transfer: float, crit: float, k: float, coeff: float) -> float:
    """Smallest lead (in r+) at which the reflected early tail alone reaches the threshold."""
    lo, hi = 1.0, 1.0e7
    if log10_tail_ratio(lo, transfer, crit, k, coeff) >= 0:
        return lo
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if log10_tail_ratio(mid, transfer, crit, k, coeff) >= 0:
            hi = mid
        else:
            lo = mid
    return hi


def future_tail_fraction(e_band: float, lag: float) -> float:
    """Reflected-tail energy arriving later than `lag` after u_c, per unit band energy.

    Arrivals after the signal's emergence label (v_out >= -v_max) cross the signal far from B's
    horizon and do not close it; the power-law tail beyond that lag carries coeff / (2 lag^2)
    (before the common transfer factor).
    """
    return tail_flux_coeff() / (2 * max(lag, 1.0) ** 2) / e_band


def m_cross_d(crit: float, k: float, coeff: float, d: float) -> int:
    """Smallest |m| whose early-tail crossing lead exceeds d (beyond it the estimate stops)."""
    for m in range(13, 2000):
        if log10_tail_ratio(d, 10.0, crit, k, coeff) + vortex_transfer_log10(m, d) - 1.0 < 0:
            return m
    return -1


def run() -> dict:
    stored = json.loads(WINDOW_JSON.read_text())
    k = stored["params"]["kappa_r_plus"]
    d = stored["params"]["d_over_r_plus"]
    c = stored["params"]["c"]
    crit = stored["closure_into_B"]["dE_crit_centre"]
    u_rel = stored["closure_into_B"]["u_B_minus_x0"]
    u_c = u_rel + d / 2  # B's label of the reflection of the packet centre (x0 = d/2)
    lag_emerge = -stored["window"]["pi"]["v_max"] - u_c  # earliest emergence label minus u_c
    e_low = stored["packets"]["P1_gaussian"]["E_below_1"]
    f_geo = stored["params"]["f_geo"]
    res: dict = {"inputs_from_window_json": {"dE_crit_centre": crit, "u_B_minus_x0": u_rel,
                                             "E_below_1": e_low, "f_geo": f_geo, "kappa": k,
                                             "d": d, "c": c}}  # fmt: skip
    # 1. reproduction of the frozen comparison
    res["reproduction"] = {
        "critical_fraction": crit / e_low,
        "frozen_R6_ratio": e_low * f_geo / crit,
        "recomputed_E_below_1": band_energy(0.0, 1.0),
    }
    # 2. angular weights (FGM (3.14)-(3.15))
    res["angular"] = {"pi": angular_factors(math.pi), "pi_over_2": angular_factors(math.pi / 2),
                      "H_mean_over_r2": h_mean_over_r2()}  # fmt: skip
    # 3. setting S: scalar field, sub-barrier part in l = 0 about A (the completion that the
    # frozen f_geo = (r+/d)^2 corresponds to). Lower bound on the energy captured by B.
    a, b = BAND
    e_band = band_energy(a, b)
    e_band_refl = band_energy(a, b, weight=lambda w: 1.0 - 4.0 * K_T * w**2)
    e_b_low = e_band_refl / (K_SIGMA * d**2)
    a_pi = res["angular"]["pi"]["uniform_deposit"]
    ratio_centre = e_b_low * a_pi / crit
    res["setting_S"] = {
        "band_r_plus": list(BAND), "E_band": e_band, "E_band_reflected_min": e_band_refl,
        "reflection_min_at_band_top": 1.0 - 4.0 * K_T * b**2,
        "capture_fraction_min": 1.0 / (K_SIGMA * d**2),
        "E_into_B_min": e_b_low, "transfer_fraction_of_E_low_min": e_b_low / e_low,
        "ratio_at_centre": ratio_centre,
        "ratio_with_mean_delay": ratio_centre * math.exp(-k * DELTA_MEAN),
        "mean_delay_allowed": DELTA_MEAN, "break_even_mean_delay": math.log(ratio_centre) / k,
        "log10_ratio_with_mean_delay": math.log10(ratio_centre) - k * DELTA_MEAN / math.log(10),
        "lag_to_emergence": lag_emerge,
        "future_tail_fraction_after_emergence": future_tail_fraction(e_band, lag_emerge),
    }  # fmt: skip
    # single-run quantities in S (order of magnitude)
    p_refl = stored["packets"]["P1_gaussian"]["weight_below_1"]
    p_b = p_refl * 2.0 / d**2
    res["single_run_S"] = {"p_not_entering_A_max": p_refl, "p_into_B_max": p_b,
                           "relative_fluctuation_of_E_B": math.sqrt((1 - p_b) / p_b),
                           "one_quantum_at_band_bottom_over_crit": a / crit}  # fmt: skip
    # 4. completions below the threshold: vortex about the A-B axis (prompt part, centre time)
    a_eq = res["angular"]["pi"]["equatorial_deposit"]
    rows, m_min = [], None
    for m in (1, 5, 8, 10, 12, 13, 15, 20, 50):
        lt = vortex_transfer_log10(m, d)
        lr = lt + math.log10(e_low * a_eq / crit)
        rows.append({"m": m, "log10_transfer_max": lt, "log10_ratio_to_threshold_max": lr})
        if m_min is None and lr < 0:
            m_min = m
    res["vortex_completion"] = {"rows": rows, "smallest_listed_m_below_threshold": m_min,
                                "m_max_for_main_part_capture": W0 * B_C, "b_c_r_plus": B_C,
                                "spin_factor_s1": (b) ** 2, "spin_factor_s2": (b) ** 4}  # fmt: skip
    # 5. early tail (one-quantum packet, power law from the truncation at omega = 0)
    coeff = tail_flux_coeff()
    t_life = stored["params"]["t_merge_r_plus"]
    lead_max = u_c + t_life
    t13 = 10 ** vortex_transfer_log10(13, d)
    nb_coeff = tail_flux_coeff(500.0, 25.0)
    leads = {}
    for m in (13, 50, 200, 1000, 1574):
        lt = 10 ** vortex_transfer_log10(m, d) if m < 300 else 0.0
        if lt > 0:
            leads[str(m)] = break_even_lead(lt, crit, k, coeff)
        else:  # transfer below double range: solve in logs (k X - 3 ln X = ln(crit/(T c)) + ln k)
            rhs = (math.log(crit) - vortex_transfer_log10(m, d) * math.log(10)
                   - math.log(coeff) + math.log(k))  # fmt: skip
            x = rhs / k
            for _ in range(50):
                x = (rhs + 3 * math.log(x)) / k
            leads[str(m)] = x
    res["early_tail"] = {
        "p0": density(0.0), "flux_coeff": coeff,
        "lead_max_to_background_start": lead_max,
        "S_break_even_lead": break_even_lead(f_geo, crit, k, coeff),
        "S_log10_ratio_at_lead_max": log10_tail_ratio(lead_max, f_geo, crit, k, coeff),
        "S_log10_ratio_at_lead_d": log10_tail_ratio(d, f_geo, crit, k, coeff),
        "vortex13_break_even_lead": break_even_lead(t13, crit, k, coeff),
        "vortex_break_even_leads": leads,
        "vortex_m_where_crossing_lead_reaches_d": m_cross_d(crit, k, coeff, d),
        "narrowband_S_break_even_lead": break_even_lead(f_geo, crit, k, nb_coeff),
        "narrowband_p0": density(0.0, 500.0, 25.0),
    }  # fmt: skip
    return res


def main() -> int:
    res = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "bound.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
