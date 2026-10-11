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
BETA_PUBLISHED = 1e-3  # the frame speed the three calibration papers quote their bounds at
BANCAL_S_LOCAL, BANCAL_S_QUANTUM = (
    7.0,
    7.2014,
)  # Bancal et al. 2012 Lemma 1 (Eq. 4) and 4-qubit state
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

# NIST 2015 (Shalm et al., PRL 115, 250402; arXiv:1511.03189): Alice-Bob 184.9 m (Fig. 4a), timing
# uncertainty 3.7 ns (Table I margins), time tags 78.125 ps (File_Folder_Descriptions.pdf), about
# two days of runs in total, and a Bell violation of only p = 2.5e-3 (1 pulse) to 5.9e-9 (5 pulses)
# over the whole dataset (Table I). The Alice-Bob detection offset t_AB is a fixed fibre/geometry
# delay the paper does not give (both photons come from one pump pulse, so aggregating slots does
# not change it).
NIST2015 = {
    "baseline_m": 184.9,
    "timing_sigma_ns": 3.7,
    "tag_ns": 0.078125,
    # Table I: 1.77e8 trials at 99,100 trials/s, i.e. one ~30 min run gives p = 5.9e-9 (5 pulses);
    # ASSUMPTION: no shorter window carries a usable violation
    "window_s_for_significance": 1.77e8 / 99_100.0,
}


def calibration() -> list[dict]:
    rows = []
    for name, p in PUBLISHED.items():
        assert p["beta"] == BETA_PUBLISHED
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
    """Envelopes of what the NIST 2015 data could bound, each assuming a zero lab offset t_AB
    (unknown; any real offset only weakens them): (i) the stated 3.7 ns timing uncertainty,
    (ii) the 78 ps time-tag floor, (iii) the rotation residual of the shortest window in which the
    violation is significant (one ~30 min run; CMB frame), with the tag floor as rho."""
    light_ns = NIST2015["baseline_m"] / C_M_S * 1e9
    rho_t = NIST2015["timing_sigma_ns"] / light_ns
    rho_tag = NIST2015["tag_ns"] / light_ns
    res = rotation_residual(BETA_CMB, NIST2015["window_s_for_significance"])
    return [
        {"case": "timing uncertainty 3.7 ns", "rho": rho_t, "bound": speed_bound(rho_t, 0, 0)},
        {"case": "time-tag floor 78 ps", "rho": rho_tag, "bound": speed_bound(rho_tag, 0, 0)},
        {
            "case": "tag floor + 30 min window rotation (CMB frame)",
            "rho": rho_tag,
            "bound": speed_bound(rho_tag, BETA_CMB, res),
        },
    ]


def required_speed(distance_ly: float, one_way_latency_s: float) -> float:
    return distance_ly * LY_LIGHT_SECONDS / one_way_latency_s


def requirements(target: float, beta: float = BETA_CMB, baselines_m=(1.2e3, 18e3, 1e5)) -> dict:
    """Budget split equally between the lab offset and the rotation residual (ASSUMPTION): solve
    target^2 = 1 + (1 - beta^2)(1 - rho^2) / (2 eps)^2 for eps, then rho = eps and the window T
    whose residual beta omega T / 2 = eps. Timing precision sigma_t = rho r / c per baseline."""
    eps = 0.5 * math.sqrt((1.0 - beta**2) / (target**2 - 1.0))
    eps = 0.5 * math.sqrt((1.0 - beta**2) * (1.0 - eps**2) / (target**2 - 1.0))  # rho^2 term
    t_max = 2.0 * eps / (beta * OMEGA_SIDEREAL)
    return {
        "target_V_over_c": target,
        "rho_max": eps,
        "window_max_s": t_max,
        "sigma_t_fs": {f"{r / 1e3:g} km": eps * r / C_M_S * 1e15 for r in baselines_m},
        "check_bound_at_budget": speed_bound(eps, beta, rotation_residual(beta, t_max)),
    }


def chsh_events_per_window(s_quantum: float = 2.6, s_local: float = 2.0, z: float = 5.0) -> float:
    """Coincidences per window to tell the quantum CHSH value from the local bound at z sigma.
    Four settings, N/4 each, |E| = S/4: var(S) = 4 * (1 - E^2) / (N/4) = 16 (1 - E^2) / N."""
    e = s_quantum / 4.0
    sigma_1 = math.sqrt(16.0 * (1.0 - e * e))  # sigma_S * sqrt(N)
    return (z * sigma_1 / (s_quantum - s_local)) ** 2


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
        "bancal_gap_ratio_vs_chsh": (BANCAL_S_QUANTUM - BANCAL_S_LOCAL)
        / BANCAL_S_LOCAL
        / ((2 * math.sqrt(2) - 2) / 2),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
