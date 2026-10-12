"""Far-side hot spots through a Simpson-Visser wormhole throat at Sgr A*: period-astrometry check.

Conditional model (docs/wh_farside/README.md): Sgr A* is a Simpson-Visser wormhole
  ds^2 = -f dt^2 + dr^2 / f + R(r)^2 dOmega^2,  f = 1 - 2M/R,  R = sqrt(r^2 + a^2),  a > 2M,
and a NIR flare is a hot spot on a circular orbit on the *other* side of the throat.

Shared parameters: M (fixed by stellar orbits, not by flares), the throat parameter a, and the
distance R0. Inputs used to fix the orbit: the polarisation period P only, assuming it equals
the orbital period. Prediction for data not used in that step: the astrometric loop radius.
- Far side: every photon that crosses the throat has impact parameter b < b_c(a), so the true images
  and their own centroid lie inside theta_c = b_c * theta_g, whatever the orbit radius (2408.07350
  shows this for ray-traced examples; here it is the general potential-barrier bound). The observed
  loop is a pipeline estimate of the total-light centroid; off-centre static light and astrometric
  errors can make it exceed theta_c. Gaps against it are summary-level standardised differences.
- Same side (ordinary control): loop radius ~ (R_e + 1) theta_g (GRAVITY 2023, Appendix A).
Kepler's law holds exactly in the areal radius: Omega^2 = M / R^3 for this metric.

  python scripts/wh_farside.py   -> results/wh_farside/check.json

Every number is a model_prediction of a conditional, unestablished geometry or a published
measurement; nothing here is evidence that such a wormhole exists.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "wh_farside"

G, C, MSUN, PC = 6.67430e-11, 2.99792458e8, 1.98847e30, 3.0856775814913673e16
UAS = math.pi / 180 / 3600 * 1e-6

# Published inputs (GRAVITY Collaboration 2023, A&A 677, L10, arXiv:2307.11821, Sect. 3-4)
M_S2 = (4.297e6, 0.012e6)  # Msun, stellar orbits (GRAVITY 2022b), independent of flares
R0 = 8277.0  # pc, the distance GRAVITY 2023 assumes
P_POL = (60.0, 3.0)  # min, polarisation period (not astrometry)
R_LOOP = (61.0, 9.0)  # uas, mean astrometric loop radius of the four flares
# EHT 2022 Sgr A* paper VI (ApJL 930, L17; arXiv:2311.09484 Sect. 5): shadow-size deviation delta
EHT_DELTA = {"Keck": (-0.04, 0.09), "VLTI": (-0.08, 0.09)}  # (value, upper 1-sigma)


def theta_g(m_sun: float, r0_pc: float = R0) -> float:
    """GM/c^2 / D in micro-arcseconds."""
    return G * m_sun * MSUN / C**2 / (r0_pc * PC) / UAS


def kepler_radius(p_min: float, m_sun: float) -> float:
    """Areal orbit radius in units of M for period p (Omega^2 = M / R^3, exact here)."""
    tg = G * m_sun * MSUN / C**3
    return (p_min * 60.0 / (2 * math.pi * tg)) ** (2.0 / 3.0)


def b_crit(a: float) -> float:
    """Largest impact parameter (units of M) that crosses the throat: 1/sqrt(max V), V = f/R^2."""
    if a <= 3.0:
        return math.sqrt(27.0)
    return a / math.sqrt(1.0 - 2.0 / a)


def b_crit_numeric(a: float, emitter_r: float = 1e6) -> float:
    """Numerical check of b_crit: maximise V(R) = (1 - 2/R) / R^2 on [a, emitter_r]."""
    rr = np.geomspace(a, emitter_r, 400001)
    v = (1 - 2 / rr) / rr**2
    return float(1 / math.sqrt(v.max()))


def shadow_delta(a: float) -> float:
    """Fractional shadow-size deviation from Schwarzschild for the symmetric SV metric."""
    return b_crit(a) / math.sqrt(27.0) - 1.0


def a_max_eht(delta_up: float) -> float:
    """Largest a whose shadow deviation stays below delta_up (bisection on a > 3)."""
    lo, hi = 3.0, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if shadow_delta(mid) < delta_up:
            lo = mid
        else:
            hi = mid
    return lo


def far_side_bound_uas(a: float, m_sun: float) -> float:
    return b_crit(a) * theta_g(m_sun)


def a_star(r_obs: float, sig: float, m_sun: float, z: float = 2.0) -> float:
    """Largest a for which the far-side bound lies z sigma below the observed loop radius."""
    target = r_obs - z * sig
    lo, hi = 2.0001, 60.0
    if far_side_bound_uas(lo, m_sun) >= target:
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if far_side_bound_uas(mid, m_sun) < target:
            lo = mid
        else:
            hi = mid
    return lo


def gravity_mass_estimator(r_o_uas: float, p_min: float, r0_pc: float = R0) -> float:
    """GRAVITY 2023 Appendix A: R_g = (2 pi / P c)^2 (r_o - R_g)^3, solved for M (Msun)."""
    ro = r_o_uas * UAS * r0_pc * PC  # metres
    k = (2 * math.pi / (p_min * 60.0 * C)) ** 2
    lo, hi = 1e-6 * ro, ro * 0.999
    for _ in range(200):
        rg = 0.5 * (lo + hi)
        if rg - k * (ro - rg) ** 3 < 0:
            lo = rg
        else:
            hi = rg
    return rg * C**2 / G / MSUN


def main() -> int:
    rng = np.random.default_rng(20261012)
    m, p = M_S2[0], P_POL[0]
    tg = theta_g(m)
    re = kepler_radius(p, m)
    ms = rng.normal(M_S2[0], M_S2[1], 20000)
    ps = rng.normal(P_POL[0], P_POL[1], 20000)
    re_mc = np.array([kepler_radius(pp, mm) for pp, mm in zip(ps, ms, strict=True)])
    same_side = (re_mc + 1.0) * np.array([theta_g(mm) for mm in ms])
    ctrl = float(np.mean(same_side)), float(np.std(same_side))
    out = {
        "inputs": {
            "M_S2_Msun": M_S2,
            "R0_pc": R0,
            "P_pol_min": P_POL,
            "loop_radius_uas": R_LOOP,
            "source": "GRAVITY Collaboration 2023, A&A 677, L10 (arXiv:2307.11821)",
        },
        "theta_g_uas": round(tg, 3),
        "kepler_radius_M": {"value": round(re, 2), "mc_sd": round(float(re_mc.std()), 2)},
        "control_same_side": {
            "predicted_loop_uas": [round(ctrl[0], 1), round(ctrl[1], 1)],
            "gap_summary_obs_minus_pred": round(
                (R_LOOP[0] - ctrl[0]) / math.hypot(R_LOOP[1], ctrl[1]), 2
            ),
        },
        "far_side": [],
    }
    for a in (2.1, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0):
        bnd = far_side_bound_uas(a, m)
        orbit_ok = re >= a  # the far-side orbit must lie outside the throat (R >= a)
        out["far_side"].append(
            {
                "a_M": a,
                "b_crit_M": round(b_crit(a), 3),
                "bound_uas": round(bnd, 1),
                "gap_summary_obs_above_bound": round((R_LOOP[0] - bnd) / R_LOOP[1], 2),
                "orbit_outside_throat_if_Ppol_is_Porb": bool(orbit_ok),
                "gravity_estimator_mass_at_bound_Msun": float(
                    f"{gravity_mass_estimator(bnd, p):.3g}"
                ),
            }
        )
    out["a_gap_2_summary_M"] = round(a_star(*R_LOOP, m), 2)
    out["a_gap_3_summary_M"] = round(a_star(*R_LOOP, m, z=3.0), 2)
    out["labels"] = {
        "gap_summary": (
            "standardised gap (obs - value) / published sigma, using GRAVITY's single averaged,"
            " phase-wrapped loop radius 61 +- 9 uas; not a full exclusion significance (no"
            " per-flare data, flare correlations, phase alignment or model systematics propagated)"
        ),
        "bound": (
            "geometric bound on the true far-side image positions (b < b_c); the observed loop is a"
            " pipeline estimate of the total-light centroid, which can exceed it through off-centre"
            " static emission or astrometric errors"
        ),
        "orbit": "a > R_e excludes a far-side Kepler orbit only if P_pol equals the orbital period",
    }
    out["eht_shadow_a_max_M_2sigma_if_calibration_transfers"] = {
        k: round(a_max_eht(v + 2 * u), 2) for k, (v, u) in EHT_DELTA.items()
    }
    out["note_eht"] = (
        "our application of the EHT delta bound to the SV critical curve; EHT calibrated the"
        " ring-to-shadow relation on Kerr and non-Kerr GRMHD libraries, not on this wormhole"
    )
    out["gravity_estimator_mass_at_observed_loop_Msun"] = float(
        f"{gravity_mass_estimator(R_LOOP[0], p):.3g}"
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
