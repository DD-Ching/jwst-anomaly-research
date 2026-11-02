"""One-way bit through the FGM window: necessary-condition comparison (analytic round, 2026-11-02).

New model; record docs/results/fgm_oneway_bit.md. Same FGM parameter point as #619
(results/fgm_window/window.json) and the closure functional with FGM's time and angle weights of
docs/results/fgm_p1_reflection_bound.md section 3. Units hbar = c = 1, lengths in r+.

Source S* (fixed): a thin disc aperture of radius A_SRC on the A-B axis, L_SRC from A's centre on
the far side from B, amplitude apodisation bump(rho/a) with a focusing phase on A's centre; current
  J = s * bump(rho / a) * bump((t - t_c) / T) * cos(omega0 (t - t_c)),  bump(x) = exp(-1/(1-x^2)),
with s = +1 or -1 chosen at t_c - T (both choices are zero before). A classical source acting on
the background state displaces it: the two hypotheses are coherent states |+alpha>, |-alpha>.

Closure (mean field, FGM (2.20), (3.3), (3.14)): N_A omega0 * eta_ext * a * W <= c kappa / 16, with
eta_ext the share of the radiated energy that crosses B's horizon via external space, a the FGM
angular weight of its deposit, W the FGM time weight averaged over its arrival profile.
Read-out (Helstrom): N_A eta_WH >= N_th.

  python scripts/fgm_oneway_bit.py   # -> results/fgm_oneway_bit/bound.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import fgm_p1_reflection as pr
import numpy as np
from scipy.special import exp1, j0

ROOT = Path(__file__).resolve().parents[1]
WINDOW_JSON = ROOT / "results" / "fgm_window" / "window.json"
OUT = ROOT / "results" / "fgm_oneway_bit"

# fixed source S* [design choices]
OMEGA0 = 500.0  # carrier, 1/r+ (P1's carrier)
T_HALF = 20.0  # half-duration of the bump envelope, r+
MARGIN = 100.0  # envelope centre enters A this far before the window edge, r+ (as #619)
L_SRC = 10.0  # aperture distance from A's centre, r+
A_SRC = 5.0  # aperture radius, r+
OMEGA_FOCUS = 10.0  # below this frequency the aperture cannot focus inside b_c [choice]
DELAY_MAX = 20.0  # extra arrival delay allowed for leaked rays (winding around A), r+ [assumption]
P_ERR = 0.1  # required error probability
B_C = pr.B_C  # critical impact parameter of RN with r- = 0.5 r+, r+


def helstrom_error(n: float) -> float:
    """Minimum error for |+alpha>, |-alpha>, equal priors: (1 - sqrt(1 - e^{-4n})) / 2."""
    return 0.5 * (1 - math.sqrt(1 - math.exp(-4 * n)))


def helstrom_threshold(p: float = P_ERR) -> float:
    return -math.log(1 - (1 - 2 * p) ** 2) / 4


def pilot_threshold(p: float = P_ERR) -> float:
    """Per-pulse photons for a reference pulse plus a signed pulse with unknown common phase.

    Interfering the two pulses leaves one output port empty; with counting the only error is the
    joint vacuum, P_e = e^{-2n} / 2 (optimal for this pair of phase-averaged states).
    """
    return math.log(1 / (2 * p)) / 2


def thermal_requirement_min(nbar: float, p: float = P_ERR) -> float:
    """Necessary photon number for displaced thermal states (fidelity bound, equal covariance)."""
    return (2 * nbar + 1) * helstrom_threshold(p)


def bump(x):
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    m = np.abs(x) < 1
    out[m] = np.exp(-1.0 / (1.0 - x[m] ** 2))
    return out


def hankel_bump(q: np.ndarray, panels: int = 80, order: int = 120) -> np.ndarray:
    """F(Q) = int_0^1 bump(x) J0(Q x) x dx by composite Gauss-Legendre (double precision)."""
    t, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(0.0, 1.0, panels + 1)
    xs = (edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * t
    ws = (edges[1:, None] - edges[:-1, None]) / 2 * w
    xs, ws = xs.ravel(), ws.ravel()
    f = bump(xs) * xs * ws
    return np.array([np.dot(f, j0(qi * xs)) for qi in np.atleast_1d(q)])


def hankel_envelope(q):
    """Saddle-point envelope of |F(Q)|: 2^{1/4} e^{-1/4} Q^{-5/4} e^{-sqrt Q} [der]."""
    q = np.asarray(q, dtype=float)
    return 2**0.25 * math.exp(-0.25) * q ** (-1.25) * np.exp(-np.sqrt(q))


def aperture_norm() -> float:
    """int_0^1 bump(x)^2 x dx = (e^{-2} - 2 E1(2)) / 2."""
    return 0.5 * (math.exp(-2.0) - 2.0 * exp1(2.0))


def halo_fraction_log(q_c: float) -> float:
    """ln of the beam energy fraction at impact parameter > b_c (paraxial focal field).

    Oscillation-averaged envelope, int_{Q_c}^inf |F|^2 Q dQ ~ 2^{1/2} e^{-1/2} e^{-2 sqrt Q_c} /
    (2 Q_c) (leading order), divided by the aperture norm (Hankel-Parseval).
    """
    c2 = 2**0.5 * math.exp(-0.5) / 2  # envelope^2 * (1/2 for the oscillation average)
    return math.log(c2) - 2 * math.sqrt(q_c) - math.log(q_c) - math.log(aperture_norm())


def halo_check(q_c: float) -> float:
    """Direct quadrature of the halo integral over the asymptotic formula (validity check)."""
    q = np.arange(q_c, q_c + 700.0, 0.02)
    direct = np.trapezoid(hankel_bump(q) ** 2 * q, q) / aperture_norm()
    return float(direct / math.exp(halo_fraction_log(q_c)))


def time_weight(u_rel_centre: float, k: float, half: float = T_HALF, n: int = 20001) -> float:
    """<exp(-kappa (u - x0))> over an arrival profile bump^2((u - u_c)/T), u_c - x0 given."""
    x = np.linspace(-1, 1, n)
    p = bump(x) ** 2
    wts = np.exp(-k * (u_rel_centre + half * x))
    return float(np.trapezoid(p * wts, x) / np.trapezoid(p, x))


def emergence_offset(v_minus_vmax: float, k: float) -> float:
    """v_out + v_max for a ray entering at v (FGM (3.8)-(3.9) with #619's normalisation)."""
    return -math.log(1 - math.exp(k * v_minus_vmax)) / k


def run() -> dict:
    s = json.loads(WINDOW_JSON.read_text())
    k = s["params"]["kappa_r_plus"]
    d = s["params"]["d_over_r_plus"]
    c = s["params"]["c"]
    x0 = d / 2
    v_max = s["window"]["pi"]["v_max"]
    q_close = c * k / 16
    res: dict = {"inputs": {"kappa": k, "d": d, "c": c, "v_max_pi": v_max, "x0": x0,
                            "closure_constant_c_kappa_over_16": q_close},
                 "source": {"omega0": OMEGA0, "T_half": T_HALF, "margin": MARGIN,
                            "L_src": L_SRC, "a_src": A_SRC, "b_c": B_C,
                            "numerical_aperture_paraxial": A_SRC / L_SRC}}  # fmt: skip
    # 1. read-out thresholds
    n_th = helstrom_threshold()
    res["readout"] = {
        "helstrom_threshold": n_th, "helstrom_error_at_threshold": helstrom_error(n_th),
        "pilot_threshold_per_pulse": pilot_threshold(),
        "pilot_threshold_total": 2 * pilot_threshold(),
        "error_without_phase_reference": 0.5,
        "thermal_nbar_at_omega0_log10": -2 * math.pi * OMEGA0 / k / math.log(10),
        "thermal_nbar_at_0p04": 1 / math.expm1(2 * math.pi * 0.04 / k),
        "requirement_at_nbar_0p58": thermal_requirement_min(1 / math.expm1(2 * math.pi * 0.04 / k)),
        "bogoliubov_mixing_log10": -2 * math.pi * OMEGA0 / k / math.log(10),
    }  # fmt: skip
    # 2. timing, paths, emergence map, phase reference
    v_c = v_max - MARGIN
    u_ext = v_c + d  # external arrival at B of the leak from the envelope centre
    early, late = -(MARGIN + T_HALF), -(MARGIN - T_HALF)
    res["paths"] = {
        "v_c": v_c, "external_arrival_minus_x0": u_ext - x0,
        "wormhole_emergence_min": -v_max,
        "wormhole_minus_external_at_B": -v_max - u_ext,
        "emerging_frequency_range": [OMEGA0 * math.expm1(-k * late),
                                     OMEGA0 * math.expm1(-k * early)],
        "emerging_duration": emergence_offset(late, k) - emergence_offset(early, k),
        "phase_sensitivity_rad_per_dln_dV": OMEGA0 * math.exp(k * MARGIN) / k,
        "dln_dV_needed_for_0p1_rad": 0.1 * k / (OMEGA0 * math.exp(k * MARGIN)),
        "dln_dV_from_O_ln_d_term": k * math.log(d),
    }  # fmt: skip
    # 3. leakage terms
    q_c = OMEGA0 * A_SRC * B_C / L_SRC
    ln_h = halo_fraction_log(q_c)
    ln_spec = -2 * math.sqrt((OMEGA0 - OMEGA_FOCUS) * T_HALF)
    tau_iso = (B_C / (2 * d)) ** 2
    res["leakage"] = {
        "Q_c": q_c, "aperture_norm": aperture_norm(),
        "halo_fraction": math.exp(ln_h), "log10_halo_fraction": ln_h / math.log(10),
        "log10_spectral_leak_below_focus": ln_spec / math.log(10),
        "tau_isotropic_estimate": tau_iso, "tau_upper_bound": 1.0,
        "log10_reflection_at_A": -math.pi * OMEGA0 / k / math.log(10),
        "halo_formula_check_ratio_at_Qc200": halo_check(200.0),
    }  # fmt: skip
    # 4. FGM weights and the two intervals
    w_c = time_weight(u_ext - x0, k)
    w_low = w_c * math.exp(-k * DELAY_MAX)
    a_min = 1.0
    a_unif = pr.angular_factors(math.pi)["uniform_deposit"]
    eta_wh = 1.0  # upper bound; achievable ~ 1 - h [assumption]
    log10 = math.log10

    def n_max_log10(tau: float, a: float, w: float) -> float:
        return log10(q_close) - log10(OMEGA0) - ln_h / math.log(10) - log10(tau * a * w)

    res["weights"] = {"W_centre_profile": w_c, "log10_W": log10(w_c), "log10_W_low": log10(w_low),
                      "a_min": a_min, "a_uniform": a_unif}  # fmt: skip
    tau_crit_log10 = (log10(q_close * eta_wh) - log10(n_th * OMEGA0) - ln_h / math.log(10)
                      - log10(a_min * w_low))  # fmt: skip
    res["intervals"] = {
        "N_required_min": n_th / eta_wh, "N_required_min_pilot": 2 * pilot_threshold() / eta_wh,
        "log10_N_allowed_max_tau1_strict": n_max_log10(1.0, a_unif, w_c),
        "log10_N_allowed_max_tau_iso": n_max_log10(tau_iso, a_min, w_low),
        "log10_tau_critical": tau_crit_log10,
        "log10_leak_share_critical": tau_crit_log10 + ln_h / math.log(10),
        "radiated_energy_at_threshold": n_th * OMEGA0,
        "log10_gap_tau_iso_over_critical": log10(tau_iso) - tau_crit_log10,
        "log10_spectral_term_over_halo_term": ln_spec / math.log(10) + log10(tau_iso)
        - (ln_h / math.log(10) + log10(tau_iso)),
        "self_coupling_shape_factor_needed": s["packets"]["P1_gaussian"]["R5_alpha_crit"]
        / (n_th * s["packets"]["P1_gaussian"]["R3_G_E_over_r_plus"]),
    }  # fmt: skip

    # 5. what would change the verdict: carrier frequency at which h * tau_iso reaches critical
    def excess(om: float) -> float:
        qc = om * A_SRC * B_C / L_SRC
        return (halo_fraction_log(qc) + math.log(tau_iso * a_min * w_low * n_th * om)
                - math.log(q_close))  # fmt: skip

    lo, hi = OMEGA0, 1e7
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if excess(mid) > 0 else (lo, mid)
    res["crossover"] = {"omega0_where_tau_iso_reaches_critical": hi,
                        "ratio_to_fixed_carrier": hi / OMEGA0}  # fmt: skip
    return res


def main() -> int:
    res = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "bound.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
