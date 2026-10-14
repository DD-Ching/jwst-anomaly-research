"""S0 screen 2026-10-14, direction A: order-of-magnitude reach of pulsar timing for exotic lenses.

Lens family of D-047 (Kitamura, Nakajima & Asada 2013): deflection alpha(b) = eps / b^n, eps of
either sign. Its lens potential gives the excess light-travel time of a pulse whose line of sight
passes the lens at impact parameter b (time-variable part; geometric term second order):
  n = 1:  dt(b) = -(eps / c) ln b          (point mass: eps = 4GM/c^2, eps < 0 for negative mass)
  n > 1:  dt(b) = (eps / c) b^(1-n) / (n - 1)   (Ellis wormhole: n = 2, eps = pi a^2 / 4, Abe 2010)
A lens moving at transverse speed v gives b(t) = sqrt(b0^2 + v^2 (t - t0)^2).

Event count for one pulsar (straight passages through the line of sight within the data span T):
  N = n_lens * 2 b_max * L * v * T
with b_max the largest impact parameter whose signal exceeds the timing threshold delta AND whose
passage time b/v is shorter than T (slower bumps are absorbed by the spin-down fit).

All numbers here are [derived] order-of-magnitude estimates under the [assumption]s in ASSUMPTIONS;
none is a measurement. Writes results/s0/pta_shapiro_reach.json.

  python scripts/s0_pta_shapiro.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "s0"

G, C = 6.674e-11, 2.998e8
MSUN, PC, AU, YR = 1.989e30, 3.086e16, 1.496e11, 3.156e7

ASSUMPTIONS = {
    "timing_threshold_s": 3e-7,  # a bump of 0.3 us after the timing-model fit (assumption)
    "transverse_speed_m_s": 2.0e5,  # 200 km/s (assumption)
    "span_yr": 15.0,  # NANOGrav 15-year-like span (assumption)
    "path_length_pc": 1000.0,  # typical pulsar distance (assumption)
    "n_pulsars": 68,  # NANOGrav 15-year count (to be checked against the release)
    "rho_dm_msun_pc3": 0.01,  # local dark-matter density (assumption)
    "fit_absorption": 0.5,  # fraction of a bump surviving the quadratic spin fit (assumption)
}


def delay(b, eps, n):
    """Excess delay (s) at impact parameter b (m) for alpha = eps / b^n."""
    if n == 1:
        return -(eps / C) * np.log(b)
    return (eps / C) * b ** (1 - n) / (n - 1)


def bump_amplitude(b0, eps, n, v, span):
    """Delay at closest approach relative to the delay at b = sqrt(b0^2 + (v span/2)^2)."""
    b_edge = np.hypot(b0, v * span / 2)
    return delay(b0, eps, n) - delay(b_edge, eps, n)


def b_max(eps, n, a=ASSUMPTIONS):
    """Largest b0 with |bump| > threshold (after fit absorption) and passage time b0 / v < span."""
    v, span = a["transverse_speed_m_s"], a["span_yr"] * YR
    thr = a["timing_threshold_s"] / a["fit_absorption"]
    grid = np.logspace(5, np.log10(v * span), 4000)  # 100 km to v*T
    ok = np.abs(bump_amplitude(grid, eps, n, v, span)) > thr
    return float(grid[ok].max()) if ok.any() else 0.0


def events_per_unit_density(eps, n, a=ASSUMPTIONS):
    """Expected events summed over pulsars for one lens per pc^3."""
    bm = b_max(eps, n, a)
    tube = 2 * bm * a["path_length_pc"] * PC * a["transverse_speed_m_s"] * a["span_yr"] * YR
    return a["n_pulsars"] * tube / PC**3


def main() -> int:
    a = ASSUMPTIONS
    rows_mass, rows_ellis = [], []
    for m in np.logspace(-6, 2, 17):
        eps = 4 * G * m * MSUN / C**2
        per_unit = events_per_unit_density(eps, 1)
        n_dm = a["rho_dm_msun_pc3"] / m  # number density if the objects were all the dark matter
        rows_mass.append(
            {
                "abs_mass_msun": float(m),
                "b_max_au": b_max(eps, 1) / AU,
                "events_if_all_dm": per_unit * n_dm,
                "f_dm_for_3_events": 3 / (per_unit * n_dm) if per_unit > 0 else None,
            }
        )
    for a_km in np.logspace(2, 7, 11):
        eps = np.pi * (a_km * 1e3) ** 2 / 4
        per_unit = events_per_unit_density(eps, 2)
        rows_ellis.append(
            {
                "throat_radius_km": float(a_km),
                "b_max_au": b_max(eps, 2) / AU,
                "n95_per_pc3_if_null": 3 / per_unit if per_unit > 0 else None,
            }
        )
    rec = {
        "label": "derived (order of magnitude) under stated assumptions; not a measurement",
        "assumptions": a,
        "point_mass_either_sign": rows_mass,
        "ellis_n2": rows_ellis,
        "check_ellis_example": {
            "a_km": 1e5,
            "b_au": 7.0,
            "delay_us": float(delay(7 * AU, np.pi * 1e16 / 4, 2) * 1e6),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pta_shapiro_reach.json").write_text(json.dumps(rec, indent=1) + "\n")

    def fmt(x):
        return f"{x:.3g}" if x else "-"

    for r in rows_mass:
        m, bm, ev = r["abs_mass_msun"], r["b_max_au"], r["events_if_all_dm"]
        print(f"M={m:.0e} b_max={bm:.2g} AU N(all DM)={ev:.3g} f(3)={fmt(r['f_dm_for_3_events'])}")
    for r in rows_ellis:
        a_km, bm = r["throat_radius_km"], r["b_max_au"]
        print(f"a={a_km:.0e} km  b_max={bm:.2g} AU  n95={fmt(r['n95_per_pc3_if_null'])} pc^-3")
    print(rec["check_ellis_example"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
