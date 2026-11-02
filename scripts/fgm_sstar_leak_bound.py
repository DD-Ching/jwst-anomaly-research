"""Direct exclusion test for S*: error-controlled lower bound on its weighted leak K (2026-11-02).

Record: docs/results/fgm_sstar_leak_bound.md. Same units (hbar = c = r+ = 1), FGM point, source S*
and read window as scripts/fgm_oneway_bit.py. Closure: N_A K <= c kappa / 16 = 0.3125; read-out
(ideal reference, eta = 1): N_A >= 0.2554; S* is excluded if K > 0.3125 / 0.2554 = 1.2235.

Method (lower bound on K from one specified receiving mode and one early time window):
1. Receiving mode: B's s-wave channel. Absorption by B is a sum over B's partial waves of
   non-negative terms, so the s-wave share is a lower bound; it deposits uniformly on B's horizon
   (FGM angular weight <H>/H(pi), computed).
2. Reciprocity of the static frequency-domain Green's function: the field at B's centre is
   e^{i w d}/(4 pi d) * int J Psi, with Psi the field of a point source at B seen at the aperture
   (B absent). All paths add coherently inside this one amplitude.
3. Psi at the aperture edge by geometric optics on the RN throat metric of FGM section 3 (exact
   null-geodesic quadrature): opposite-side primary path (family 2), same-side path (family 1),
   relativistic orbits (n >= 1).
4. The overlap with the bump-apodised aperture has no stationary point; it is the edge saddle of
   the bump, |int_0^1 x bump(x) e^{iqx} dx| ~ e^{-1/4} sqrt(pi) 2^{-1/4} q^{-3/4} e^{-sqrt q},
   q = sigma a with sigma the overlap phase gradient at the edge (checked against quadrature).
5. Time gate: the first relativistic orbit reaches B 18.8 r+ after family 2; inside that window only
   families 1 and 2 contribute, and family 1 is bounded by its edge factor.

  python scripts/fgm_sstar_leak_bound.py   # -> results/fgm_sstar_leak_bound/bound.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import fgm_oneway_bit as ob
import fgm_p1_reflection as pr
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
WINDOW_JSON = ROOT / "results" / "fgm_window" / "window.json"
ONEWAY_JSON = ROOT / "results" / "fgm_oneway_bit" / "bound.json"
OUT = ROOT / "results" / "fgm_sstar_leak_bound"

RM = 0.5  # r- / r+
M, Q2 = (1 + RM) / 2, RM  # RN mass and charge^2 with r+ = 1
# allowances [assumption], each in the direction that lowers the bound
GAMMA0_MIN = 0.9  # B's s-wave transmission at omega r+ = 500 (above the s-wave barrier)
NORM_FACTOR_MIN = 0.5  # curved-space source normalisation and coupling vs the flat formula
SADDLE_ERR = 0.3  # relative amplitude error of the edge-saddle evaluation
DELTA_CONV = 100.0  # r+: B-side tortoise and x*0 O(ln d/r+) conventions, transmission delay
GATE_SAFETY = 1.0  # r+: kept between the window end and the first relativistic arrival
C_THETA = 1.0 / math.cos(math.atan(0.5))  # <1/cos theta> bound for the sheet's radiation


def f(r: float) -> float:
    return (r - 1) * (r - RM) / r**2


def rstar(r: float) -> float:
    """Tortoise coordinate with r*(0) = 0 (FGM (3.2))."""
    c = 1 / (1 - RM)
    c2 = RM**2 / (1 - RM)
    r0 = c * math.log(1.0) - c2 * math.log(RM)
    return r + c * math.log(abs(r - 1)) - c2 * math.log(abs(r - RM)) - r0


R_PH = (3 * M + math.sqrt(9 * M * M - 8 * Q2)) / 2
B_C = R_PH / math.sqrt(f(R_PH))


def g_fun(u: float, b: float) -> float:
    return 1 / b**2 - u**2 * (1 - 2 * M * u + Q2 * u**2)


def u_turn(b: float) -> float:
    return brentq(lambda u: g_fun(u, b), 1e-9, 1 / R_PH, xtol=1e-16)


def _p(h: float, ut: float) -> float:
    """G(ut - h) = h * P(h), expanded exactly so the turning point carries no cancellation."""
    return (2 * ut - h - 6 * M * ut**2 + 6 * M * ut * h - 2 * M * h * h + 4 * Q2 * ut**3
            - 6 * Q2 * ut**2 * h + 4 * Q2 * ut * h * h - Q2 * h**3)  # fmt: skip


def leg_sweep(r1: float, b: float) -> float:
    ut = u_turn(b)
    smax = math.sqrt(ut - 1 / r1)
    return quad(lambda s: 2 / math.sqrt(_p(s * s, ut)), 0, smax, limit=400,
                epsabs=1e-14, epsrel=1e-13)[0]  # fmt: skip


def leg_time(r1: float, b: float) -> float:
    ut = u_turn(b)
    rt = 1 / ut
    smax = math.sqrt(r1 - rt)

    def g(s: float) -> float:
        r = rt + s * s
        u = 1 / r
        return 2 / (f(r) * b * math.sqrt(ut * u * _p(s * s * ut * u, ut)))

    pts = [x for x in (0.5, 1, 2, 4, 8, 16, 32, 64) if x < smax]
    return quad(g, 0, smax, points=pts, limit=800, epsabs=1e-10, epsrel=1e-13)[0]


def solve_b(r1: float, r2: float, target: float) -> float:
    hi = r1 / math.sqrt(f(r1)) * (1 - 1e-12)
    return brentq(lambda b: leg_sweep(r1, b) + leg_sweep(r2, b) - target, B_C * (1 + 1e-13), hi,
                  xtol=1e-15)  # fmt: skip


def plane_rho(b: float, d: float, length: float, phi_side: int, r_guess: float) -> float:
    """Transverse radius where the ray from B (impact parameter b) crosses z = -L.

    phi_side = +1: the ray reaches the aperture after passing the opposite side (total sweep
    pi + psi); -1: same side (sweep pi - psi).
    """
    s0 = leg_sweep(d, b)

    def h(r: float) -> float:
        phi = s0 + leg_sweep(r, b)
        return r * math.cos(phi) + length

    r = brentq(h, r_guess * 0.97, r_guess * 1.03, xtol=1e-14)
    phi = s0 + leg_sweep(r, b)
    return abs(r * math.sin(phi)) if phi_side else float("nan")


def edge_factor(q: float) -> float:
    """Asymptotic |int_0^1 x bump(x) e^{iqx} dx| (edge saddle of the bump)."""
    return math.exp(-0.25) * math.sqrt(math.pi) * 2**-0.25 * q**-0.75 * math.exp(-math.sqrt(q))


def two_edge_numeric(q: float, panels: int = 400, order: int = 64) -> float:
    """|int_{-1}^{1} bump(x) e^{iqx} dx| by composite Gauss-Legendre (no interior endpoint)."""
    t, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(-1.0, 1.0, panels + 1)
    xs = (edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * t
    ws = (edges[1:, None] - edges[:-1, None]) / 2 * w
    xs, ws = xs.ravel(), ws.ravel()
    return float(abs(np.sum(ws * ob.bump(xs) * np.exp(1j * q * xs))))


def two_edge_mp(q: float, dps: int = 45, panels: int = 600) -> float:
    """Same integral with mpmath at high precision (for q where double precision is too short)."""
    import mpmath as mp

    mp.mp.dps = dps
    one = mp.mpf(1)
    nodes = [-one + 2 * one * i / panels for i in range(panels + 1)]
    val = mp.quad(lambda x: mp.exp(-1 / (1 - x * x)) * mp.cos(q * x), nodes)
    return float(abs(val))


def envelope_ratio(q: float, use_mp: bool = False, span: int = 40) -> float:
    """max over q' in [q, q + 2 pi] of the two-edge integral, over 2 * edge_factor(q')."""
    qs = [q + 2 * math.pi * i / span for i in range(span + 1)]
    fun = two_edge_mp if use_mp else two_edge_numeric
    return max(fun(qq) / (2 * edge_factor(qq)) for qq in qs)


def family_at_edge(b: float, r_e: float, length: float, rho: float, side: int, d: float):
    """Overlap phase gradient sigma/omega, |n_z| and ray-tube amplitude of a family at the edge."""
    fr = f(r_e)
    radial = math.sqrt(1 - b * b * fr / r_e**2)
    sigma = (rho / r_e) * (radial / fr - 1) + side * b * length / r_e**2
    phi = math.pi + side * math.atan(rho / length)
    n_r, n_phi = radial, b * math.sqrt(fr) / r_e
    n_z = abs(n_r * math.cos(phi) - n_phi * math.sin(phi))
    db = 1e-5 * b
    rp = plane_rho(b + db, d, length, side, r_e)
    rm = plane_rho(b - db, d, length, side, r_e)
    drho_db = (rp - rm) / (2 * db)
    amp = math.sqrt(b / (abs(drho_db) * rho * n_z))
    return {"b": b, "sigma_over_omega": abs(sigma), "n_z": n_z, "drho_db": drho_db,
            "psi_amplitude": amp}  # fmt: skip


def run() -> dict:
    s_win = json.loads(WINDOW_JSON.read_text())
    s_one = json.loads(ONEWAY_JSON.read_text())
    k = s_win["params"]["kappa_r_plus"]
    d = s_win["params"]["d_over_r_plus"]
    v_max = s_win["window"]["pi"]["v_max"]
    q_close = s_one["inputs"]["closure_constant_c_kappa_over_16"]
    n_th = s_one["readout"]["helstrom_threshold"]
    omega, a, length, t_half = ob.OMEGA0, ob.A_SRC, ob.L_SRC, ob.T_HALF
    res: dict = {"thresholds": {"N_th": n_th, "closure": q_close, "K_critical": q_close / n_th}}
    # geometry of the paths from the aperture edge to B
    r_e = math.hypot(a, length)
    psi = math.atan(a / length)
    targets = {"family1_same_side": math.pi - psi, "family2_opposite_side": math.pi + psi,
               "orbit1_same_side": 3 * math.pi - psi, "orbit1_opposite_side": 3 * math.pi + psi,
               "orbit2_same_side": 5 * math.pi - psi}  # fmt: skip
    paths = {}
    for name, tgt in targets.items():
        b = solve_b(r_e, d, tgt)
        paths[name] = {"b": b, "b_minus_bc": b - B_C, "travel_time": leg_time(r_e, b)
                       + leg_time(d, b)}  # fmt: skip
    t2 = paths["family2_opposite_side"]["travel_time"]
    for p in paths.values():
        p["minus_family2"] = p["travel_time"] - t2
    gate = min(p["minus_family2"] for n, p in paths.items() if n.startswith("orbit"))
    res["geometry"] = {"r_ph": R_PH, "b_c": B_C, "r_edge": r_e, "paths": paths,
                       "first_orbit_after_family2": gate}  # fmt: skip
    # reciprocal field at the edge, families 2 and 1
    fam2 = family_at_edge(paths["family2_opposite_side"]["b"], r_e, length, a, +1, d)
    fam1 = family_at_edge(paths["family1_same_side"]["b"], r_e, length, a, -1, d)
    q2 = fam2["sigma_over_omega"] * omega * a
    q1 = fam1["sigma_over_omega"] * omega * a
    e2, e1 = edge_factor(q2), edge_factor(q1)
    ratio_1_2 = fam1["psi_amplitude"] * e1 / (fam2["psi_amplitude"] * e2)
    res["edge"] = {"family2": fam2 | {"q": q2, "edge_factor": e2},
                   "family1": fam1 | {"q": q1, "edge_factor": e1},
                   "family1_over_family2_amplitude": ratio_1_2,
                   "envelope_ratio_mp_at_q2": envelope_ratio(q2, use_mp=True, span=16),
                   "envelope_ratio_numeric": {str(qq): envelope_ratio(qq)
                                              for qq in (100.0, 300.0, 600.0)}}  # fmt: skip
    # energy ratio into B's s-wave channel per radiated energy (flat-source normalisation)
    n_ap = ob.aperture_norm()
    i_rho = fam2["psi_amplitude"] * a * a * e2 * (1 - SADDLE_ERR)
    ratio = GAMMA0_MIN * NORM_FACTOR_MIN * i_rho**2 / (4 * C_THETA * d * d * a * a * n_ap)
    # time window and weights
    x = np.linspace(-1.0, 1.0, 400001)
    e_sq = ob.bump(x) ** 2
    win_len = gate - GATE_SAFETY
    m = x <= -1 + win_len / t_half
    f_win = float(np.trapezoid(e_sq[m], x[m]) / np.trapezoid(e_sq, x))
    v_c = v_max - ob.MARGIN
    t_edge = -(r_e - length)  # the edge leads the centre (flat focusing phase)
    delta2 = t_edge + t2 - (rstar(length) + d)  # family-2 arrival minus #619's u_c = v_c + d
    u_end_rel = (v_c + d - d / 2) + delta2 - t_half + win_len + DELTA_CONV
    a_unif = pr.angular_factors(math.pi)["uniform_deposit"]
    k_lower = omega * ratio * a_unif * max(f_win - 2 * ratio_1_2, 0.0) * math.exp(-k * u_end_rel)
    log10_k = math.log10(k_lower)
    k_crit = q_close / n_th
    res["bound"] = {
        "s_wave_energy_ratio_min": ratio, "window_length": win_len, "window_energy_fraction": f_win,
        "delta2_vs_619_reference": delta2, "u_end_minus_x0_max": u_end_rel,
        "angular_weight_uniform": a_unif, "K_lower": k_lower, "log10_K_lower": log10_k,
        "log10_margin_over_K_critical": log10_k - math.log10(k_crit),
        "extra_offset_break_even_r_plus": (log10_k - math.log10(k_crit)) * math.log(10) / k,
        "allowances": {"Gamma0_min": GAMMA0_MIN, "normalisation_min": NORM_FACTOR_MIN,
                       "saddle_rel_err": SADDLE_ERR, "delta_conv": DELTA_CONV,
                       "gate_safety": GATE_SAFETY, "c_theta": C_THETA},
    }  # fmt: skip
    iso = s_one["intervals"]
    n_iso = iso["log10_N_allowed_max_tau_iso"]
    res["comparison_isotropic"] = {"log10_N_allowed_max_tau_iso": n_iso,
                                   "log10_K_iso": math.log10(q_close) - n_iso}  # fmt: skip
    scan = []
    for rho in (0.5, 1.0, 2.0, 3.0, 4.0, a):
        r = math.hypot(rho, length)
        ps = math.atan(rho / length)
        te = -(r - length)
        b2 = solve_b(r, d, math.pi + ps)
        rad = math.sqrt(1 - b2 * b2 * f(r) / r**2)
        sig = (rho / r) * (rad / f(r) - 1) + b2 * length / r**2
        bo = solve_b(r, d, 3 * math.pi - ps)
        scan.append({"rho": rho, "family2_sigma_over_omega": sig,
                     "family2_arrival": te + leg_time(r, b2) + leg_time(d, b2),
                     "orbit1_same_arrival": te + leg_time(r, bo) + leg_time(d, bo)})  # fmt: skip
    t2_edge = scan[-1]["family2_arrival"]
    sig_min = min(x["family2_sigma_over_omega"] for x in scan)
    orbit_gap = min(x["orbit1_same_arrival"] for x in scan) - t2_edge
    res["aperture_scan"] = {"rows": scan, "min_family2_sigma_over_omega": sig_min,
                            "earliest_orbit_minus_family2_edge": orbit_gap}  # fmt: skip
    res["verdict"] = "excluded" if k_lower > k_crit else "C"
    return res


def main() -> int:
    res = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "bound.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
