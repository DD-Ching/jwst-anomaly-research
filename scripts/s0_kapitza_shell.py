"""S0-2026-10-15b (Kapitza screen): does a confined high-frequency field (averaged: a 2D radiation
component) stabilise a thin-shell wormhole throat?

Unit file: docs/state/units/S0-2026-10-15b.md. `derived` (exact junction conditions), toy;
no data.

Model: Schwarzschild thin-shell wormhole (Visser 1989; Poisson & Visser 1995, gr-qc/9506083),
G = c = M = 1.
  sigma = -(1/2 pi a) sqrt(1 - 2M/a + adot^2); static p0 = (1/4 pi a0)(1 - M/a0)/sqrt(1 - 2M/a0);
  equation of motion adot^2 + V(a) = 0, V(a) = 1 - 2M/a - (2 pi a sigma(a))^2;
  each component conserved on the shell: d sigma_i/da = -2 (sigma_i + p_i)/a.
Components:
  exotic, linear EoS p_x = p_x0 + beta2_x (sigma_x - sigma_x0);
  radiation (the Brill-Hartle / Isaacson average of fields confined to the shell), p_r = sigma_r/2,
  so sigma_r ~ a^-3.
The throat a0 is fixed. The radiation share is f = sigma_r0/|sigma0|; the exotic part absorbs
the rest, sigma_x0 = sigma0 - sigma_r0 (more negative). Stable iff V''(a0) > 0.

ASSUMPTIONS (new hypotheses, not derived here): the radiation stays confined to the shell
(no leakage); its average obeys a 2D radiation EoS; the exotic EoS is linear near a0.

  python scripts/s0_kapitza_shell.py      # -> results/s0/kapitza_shell.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parents[1] / "results" / "s0"
M = 1.0


def static(a0):
    s0 = -np.sqrt(1 - 2 * M / a0) / (2 * np.pi * a0)
    p0 = (1 - M / a0) / np.sqrt(1 - 2 * M / a0) / (4 * np.pi * a0)
    return s0, p0


def potential(a0, b2, f, agrid):
    s0, p0 = static(a0)
    sr0 = f * abs(s0)
    sx0, px0 = s0 - sr0, p0 - sr0 / 2

    def rhs(a, y):
        return [-2 * (y[0] + px0 + b2 * (y[0] - sx0)) / a]

    out = []
    for a in agrid:
        sx = sx0 if a == a0 else solve_ivp(rhs, [a0, a], [sx0], rtol=1e-12, atol=1e-15).y[0, -1]
        out.append(1 - 2 * M / a - (2 * np.pi * a * (sx + sr0 * (a0 / a) ** 3)) ** 2)
    return np.array(out)


def derivatives(a0, b2, f, h=1e-3):
    g = a0 * np.array([1 - 2 * h, 1 - h, 1.0, 1 + h, 1 + 2 * h])
    v = potential(a0, b2, f, g)
    d = a0 * h
    return (
        v[2],
        (v[3] - v[1]) / (2 * d),
        (-v[4] + 16 * v[3] - 30 * v[2] + 16 * v[1] - v[0]) / (12 * d * d),
    )


def needed_share(a0, b2):
    fs = np.concatenate([[0.0], np.logspace(-3, 2.5, 60)])
    vals = [derivatives(a0, b2, f)[2] for f in fs]
    for i in range(1, len(fs)):
        if vals[i - 1] <= 0 < vals[i]:
            return brentq(lambda f: derivatives(a0, b2, f)[2], fs[i - 1], fs[i], xtol=1e-10)
    return None


def main() -> int:
    rows = []
    for a0 in (2.5, 4.0, 6.0, 10.0, 20.0):
        s0, p0 = static(a0)
        for b2 in (0.0, 0.25, 0.75, 1.0):
            v0, v1, v2 = derivatives(a0, b2, 0.0)
            f = needed_share(a0, b2)
            rows.append(
                {
                    "a0_over_M": a0,
                    "beta2_exotic": b2,
                    "shell_sigma_plus_p_sign": int(np.sign(s0 + p0)),
                    "static_check_V_Vprime": [float(v0), float(v1)],
                    "Vpp_without_radiation": float(v2),
                    "dVpp_df_at_0": float((derivatives(a0, b2, 1e-4)[2] - v2) / 1e-4),
                    "f_needed": None if f is None else float(f),
                    "exotic_negative_energy_factor": None if f is None else float(1 + f),
                }
            )
    OUT.mkdir(parents=True, exist_ok=True)
    rec = {"label": "derived (thin-shell junction conditions), toy model; no data", "rows": rows}
    (OUT / "kapitza_shell.json").write_text(json.dumps(rec, indent=1) + "\n")
    for r in rows:
        print(r["a0_over_M"], r["beta2_exotic"], r["f_needed"], r["exotic_negative_energy_factor"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
