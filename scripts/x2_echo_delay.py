"""X2 S0 minimal check: echo delay of a symmetric Simpson-Visser throat with two potential peaks.

Metric (G = c = M = 1): ds^2 = -f dt^2 + dr^2 / f + R^2 dOmega^2, f = 1 - 2/R, R = sqrt(r^2 + a^2).
For 2 < a < 3 the wave potential of every massless field peaks at R = 3 on each side of the throat
(r = +-sqrt(9 - a^2)). A ringdown pulse that leaks inward bounces between the two peaks, so the
echo spacing is the round trip between them in tortoise distance (dr* = dr / f):
  dt_echo = 4 * L(a),  L(a) = int_0^{sqrt(9 - a^2)} dr / f(R(r)).
Near a -> 2 (eps = 1 - 2/a -> 0): f ~ eps + r^2 / a^3, so L ~ (pi / 2) sqrt(a^3 / eps), a power
law in eps (not the logarithm of horizon-scale ECO models).

The check maps the post-merger windows searched by Miani et al. 2023 (arXiv:2302.12158: echo time
50-350 ms or 0.2-1.2 s after coalescence) onto a - 2 for a given remnant mass.

  python scripts/x2_echo_delay.py   -> results/x2/s0_echo_delay.json

Model computation of a conditional, unestablished geometry; not evidence that any remnant is a
wormhole.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x2"
T_SUN = 4.925490947e-6  # G Msun / c^3 in s


def half_cavity(a: float) -> float:
    """Tortoise distance from the throat to one potential peak (R = 3), in M."""
    if not 2.0 < a < 3.0:
        raise ValueError("two peaks need 2 < a < 3")
    eps = 1.0 - 2.0 / a
    rp = math.sqrt(9.0 - a * a)
    # split off the near-throat Lorentzian part analytically for accuracy at small eps:
    # 1/f = 1/(eps + q r^2) + [1/f - 1/(eps + q r^2)], q = 1/a^3 (the bracket is regular)
    q = 1.0 / a**3

    def f(r):
        return 1.0 - 2.0 / math.sqrt(r * r + a * a)

    def rest(r):
        return 1.0 / f(r) - 1.0 / (eps + q * r * r)

    lor = math.atan(rp * math.sqrt(q / eps)) / math.sqrt(q * eps)
    # the bracket is O(r^2)/eps^2-free but cancels two large terms near r = 0 for tiny eps;
    # an absolute tolerance of 1e-9 M is ample next to lor ~ (pi/2) sqrt(a^3/eps)
    return lor + quad(rest, 0.0, rp, limit=400, epsabs=1e-9, epsrel=1e-10)[0]


def echo_delay(a: float) -> float:
    """Round-trip time between the two peaks, in M."""
    return 4.0 * half_cavity(a)


def a_for_delay(dt_m: float) -> float:
    """a - 2 (in M) giving an echo spacing of dt_m (in M)."""
    if dt_m <= echo_delay(2.999):
        raise ValueError("delay below the two-peak range")
    return brentq(lambda x: echo_delay(2.0 + x) - dt_m, 1e-14, 0.999, xtol=1e-16, rtol=1e-12)


def main() -> int:
    rows = []
    for d in (1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99):
        a = 2.0 + d
        eps = 1.0 - 2.0 / a
        rows.append(
            {
                "a_minus_2M": d,
                "dt_echo_M": round(echo_delay(a), 2),
                "small_eps_formula_M": round(2 * math.pi * math.sqrt(a**3 / eps), 2),
            }
        )
    windows = {}
    for m_final in (30.0, 60.0, 90.0):
        tm = m_final * T_SUN
        windows[str(m_final)] = {
            f"{lo}-{hi}s": [
                float(f"{a_for_delay(hi / tm):.2e}"),
                float(f"{a_for_delay(lo / tm):.2e}"),
            ]
            for lo, hi in ((0.05, 0.35), (0.2, 1.2))
        }
    out = {
        "model": "symmetric Simpson-Visser, 2M < a < 3M, peaks at R = 3M on both sides",
        "table": rows,
        "a_minus_2M_range_for_window_by_final_mass_Msun": windows,
        "window_source": "Miani et al. 2023, arXiv:2302.12158 Table I (post-merger windows)",
        "label": "derived (model computation); premises in docs/state/units/X2.md",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s0_echo_delay.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
