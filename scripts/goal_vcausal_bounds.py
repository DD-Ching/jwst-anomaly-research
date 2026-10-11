"""Preferred-frame ("v-causal") influence-speed bounds: calibration, data adequacy, requirements.

Goal-round-1 unit (docs/goal_round1/README.md section 3). Three questions, no data scanning:
1. Calibration: does the published bipartite bound formula, with each experiment's stated
   parameters, reproduce the published bounds (Salart et al. 2008, Yin et al. 2013, Cocciaro et
   al. 2018)? This is the positive control.
2. Data adequacy: what is the best bound the public NIST 2015 loophole-free Bell data could give,
   from the timing geometry in Shalm et al. 2015 (Table I, Fig. 3)?
3. Requirements: lab simultaneity, integration window and coincidences needed to reach the
   speeds a 20 / 100 light-year conversation would need.

Formula (Salart et al. 2008 Eqs. 2-3; Yin et al. 2013 Eq. 2; Cocciaro et al. 2018 Eq. 5):
    (V/c)^2 >= 1 + (1 - beta^2) (1 - rho^2) / (rho + beta_par)^2,
rho = c t_AB / r_AB (lab time offset of the two measurements over the light time between them) and
beta_par the preferred-frame velocity component along the baseline. Over an integration window T
the Earth's rotation sweeps beta_par; for a baseline perpendicular to the frame's projected
velocity the worst-case residual is beta sin(theta) omega T / 2 (Yin Eq. 5; Cocciaro Eq. 6 uses the
doubled window). These bounds exclude influences SLOWER than V in the frames covered. They do not
exclude faster ones and say nothing about controllable signalling (section 1.1 of the README).
"""

from __future__ import annotations

import argparse
import json
import math

OMEGA_SIDEREAL = 7.2921159e-5  # rad / s
C_M_S = 299_792_458.0
BETA_CMB = 369.82e3 / C_M_S  # Planck 2018 dipole velocity / c (ASSUMPTION: candidate frame)
LY_LIGHT_SECONDS = 365.25 * 86400.0


def speed_bound(rho: float, beta: float, beta_par: float) -> float:
    """Lower bound V/c from a run with lab offset rho and residual frame component beta_par."""
    if not (0 <= rho < 1 and 0 <= beta < 1):
        raise ValueError("need 0 <= rho < 1 and 0 <= beta < 1")
    den = rho + abs(beta_par)
    if den == 0:
        return math.inf
    return math.sqrt(1.0 + (1.0 - beta**2) * (1.0 - rho**2) / den**2)


def rotation_residual(beta: float, window_s: float, sin_theta: float = 1.0) -> float:
    """Worst-case |beta_par| inside one integration window (Yin et al. 2013 Eq. 5)."""
    return abs(beta * sin_theta * OMEGA_SIDEREAL * window_s / 2.0)


# Published parameters (literature audit, docs/goal_round1/README.md section 1.2).
PUBLISHED = {
    "Salart2008": {
        "rho": 5.4e-6,
        "window_s": 360.0,
        "beta": 1e-3,
        "published_bound": 5.4e4,
        "where": "Nature 454, 861: rho-bar Eq. 3 text, fringe period 360 s, Fig. 5 (chi = 90 deg)",
    },
    "Yin2013": {
        "rho": 6.84e-6,
        "window_s": 1800.0,
        "beta": 1e-3,
        "published_bound": 1.38e4,
        "where": "PRL 110, 260407: Eq. 2 and Eq. 5, T = 1800 s, theta = pi/2",
    },
    "Cocciaro2018": {
        "rho": 1.83e-7,
        "window_s": 0.494,
        "beta": 1e-3,
        "published_bound": 5e6,
        "where": "PRA 97, 052124: Eq. 5-7, delta t = 2 Delta t = 0.494 s, Fig. 7 (chi = pi/2)",
    },
}

# NIST 2015 (Shalm et al., PRL 115, 250402; arXiv:1511.03189). Alice-Bob baseline read off Fig. 3b
# (~185 m; ASSUMPTION within the paper's 1 m position uncertainty); a measurement is spread over
# the aggregated 12.6 ns pump slots inside a ~200 ns Pockels window (text below Fig. 3; Table I),
# so even a perfectly centred pair has rho >= (slot spread) / (light time).
NIST2015 = {"baseline_m": 185.0, "slot_ns": 12.6, "aggregate_pulses": (1, 3, 5, 7)}


def calibration() -> list[dict]:
    rows = []
    for name, p in PUBLISHED.items():
        res = rotation_residual(p["beta"], p["window_s"])
        v = speed_bound(p["rho"], p["beta"], res)
        rows.append(
            {
                "experiment": name,
                "reproduced_bound": v,
                "published_bound": p["published_bound"],
                "ratio": v / p["published_bound"],
                "source": p["where"],
            }
        )
    return rows


def nist_best_case() -> list[dict]:
    """Optimistic upper envelope of what NIST 2015 data could bound: lab offset 0 (the real
    geometry is asymmetric, so the true rho is larger) and a frame aligned at every instant."""
    light_ns = NIST2015["baseline_m"] / C_M_S * 1e9
    out = []
    for n in NIST2015["aggregate_pulses"]:
        rho = n * NIST2015["slot_ns"] / light_ns
        out.append(
            {"aggregate_pulses": n, "rho_min": rho, "best_bound": speed_bound(rho, 0.0, 0.0)}
        )
    return out


def required_speed(distance_ly: float, one_way_latency_s: float) -> float:
    return distance_ly * LY_LIGHT_SECONDS / one_way_latency_s


def requirements(target: float, beta: float = BETA_CMB, baselines_m=(1.2e3, 18e3, 1e5)) -> dict:
    """Budget split equally between the lab offset and the rotation residual (ASSUMPTION):
    rho <= 1/(2 V), T <= 1/(V beta omega). Timing precision sigma_t = rho r / c per baseline."""
    rho_max = 1.0 / (2.0 * target)
    t_max = 1.0 / (target * beta * OMEGA_SIDEREAL)
    return {
        "target_V_over_c": target,
        "rho_max": rho_max,
        "window_max_s": t_max,
        "sigma_t_fs": {f"{r / 1e3:g} km": rho_max * r / C_M_S * 1e15 for r in baselines_m},
        "check_bound_at_budget": speed_bound(rho_max, beta, rotation_residual(beta, t_max)),
    }


def chsh_events_per_window(s_quantum: float = 2.6, s_local: float = 2.0, z: float = 5.0) -> float:
    """Coincidences per window to tell the quantum CHSH value from the local bound at z sigma, with
    sigma_S ~ 2 / sqrt(N) (ASSUMPTION: order of magnitude for four equally sampled settings)."""
    return (z * 2.0 / (s_quantum - s_local)) ** 2


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    targets = {
        "20 ly, 1 s one-way": required_speed(20, 1.0),
        "20 ly, 0.3 s": required_speed(20, 0.3),
        "100 ly, 1 s": required_speed(100, 1.0),
        "100 ly, 0.3 s": required_speed(100, 0.3),
    }
    out = {
        "provenance": "model_prediction (formula) + literature parameters; no new measurement",
        "calibration": calibration(),
        "nist2015_best_case": nist_best_case(),
        "targets": {k: {"V_over_c": v, **requirements(v)} for k, v in targets.items()},
        "chsh_coincidences_per_window_5sigma": chsh_events_per_window(),
        "bancal_gap_ratio_vs_chsh": (7.2014 - 7.0) / 7.0 / ((2 * math.sqrt(2) - 2) / 2),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
