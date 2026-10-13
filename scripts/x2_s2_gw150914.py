"""X2 S2: one-event demonstration (GW150914) under phenomenological assumptions P1-P4.

Design frozen in docs/state/units/X2_s2_design.md before this script was run. Observational input:
Miani et al. 2023 (arXiv:2302.12158v2) Table II for GW150914:
- h_MR = 16.0e-23 Hz^-1/2 (merger-ringdown h_rss);
- h_UL = 3.4e-23 Hz^-1/2 (95 % coverage, PMW [0.2, 1.2] s);
- ratio 0.21.

Echo train in the PMW: pulse k (k = 1, 2, ...) arrives at t_k = k dt(a) and carries T R^(2k-1).
The source is eta times the merger-ringdown spectrum (eta free). Then
  F(a)^2 = sum_{k in PMW} int W |T|^2 |R|^(2(2k-1)) df / int W df,   eta_max(a) = ratio_UL / F(a).

  python scripts/x2_s2_gw150914.py   -> results/x2/s2_gw150914.json

A limit on a phenomenological echo model (test-field barrier as a proxy for gravitational waves,
non-rotating SV, free source amplitude). It is not an observational exclusion of the SV geometry.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import x2_echo_delay as xe  # noqa: E402
import x2_scatter as xs  # noqa: E402

OUT = ROOT / "results" / "x2"
T_SUN = 4.925490947e-6
M_DET = 68.8  # Msun, GWTC-1 GW150914: M_f = 63.1 (source), z = 0.09
CHI = 0.69
RATIO_UL = 0.21
PMW = (0.2, 1.2)


def ringdown_mode(m_det=M_DET, chi=CHI):
    """(f, tau) of the l = m = 2, n = 0 Kerr mode (Berti, Cardoso & Will 2006 fits)."""
    tm = m_det * T_SUN
    mw = 1.5251 - 1.1568 * (1 - chi) ** 0.1292
    q = 0.7 + 1.4187 * (1 - chi) ** (-0.4990)
    f = mw / (2 * math.pi * tm)
    return f, q / (math.pi * f)


def weights(fgrid):
    f0, tau = ringdown_mode()
    g = 1.0 / (2 * math.pi * tau)
    lor = 1.0 / ((fgrid - f0) ** 2 + g**2) + 1.0 / ((fgrid + f0) ** 2 + g**2)
    flat = ((fgrid >= 35.0) & (fgrid <= 350.0)).astype(float)
    return {"ringdown_lorentzian": lor, "flat_35_350Hz": flat}


def analyse(a, s, fgrid, W):
    tm = M_DET * T_SUN
    dt_s = xe.echo_delay(a) * tm
    w = 2 * math.pi * fgrid * tm  # dimensionless M omega
    g = xs.Grid(a, 2, s, xmax=2500.0, h=0.01)
    ok = w * w > max(g.v0, g.V[-1]) * 1.0001
    R2 = np.ones_like(w)
    T2 = np.zeros_like(w)
    R2[ok], T2[ok] = xs.single_barrier(g, w[ok])
    ks = [k for k in range(1, 50) if PMW[0] <= k * dt_s <= PMW[1]]
    res = {}
    for name, wt in W.items():
        norm = np.trapezoid(wt, fgrid)
        e = [np.trapezoid(wt * T2 * R2 ** (2 * k - 1), fgrid) / norm for k in ks]
        F = math.sqrt(sum(e)) if e else 0.0
        spec1 = wt * T2 * R2  # first pulse spectrum (power)
        fmean = float(np.trapezoid(fgrid * spec1, fgrid) / np.trapezoid(spec1, fgrid))
        # rms spectral width -> pulse duration scale ~ 1 / (2 pi sigma_f)
        fvar = float(np.trapezoid((fgrid - fmean) ** 2 * spec1, fgrid) / np.trapezoid(spec1, fgrid))
        dur = 1.0 / (2 * math.pi * math.sqrt(fvar)) if fvar > 0 else None
        frac2 = (sum(e[:2]) / sum(e)) if e else None
        res[name] = {
            "F": F,
            "eta_max": RATIO_UL / F if F > 0 else None,
            "pulses_in_PMW": ks,
            "first_pulse_mean_f_Hz": fmean,
            "pulse_duration_scale_s": dur,
            "first_two_pulses_fraction_of_hrss2": frac2,
            "applicable": bool(
                ks
                and ks[0] == 1
                and len(ks) >= 2
                and 80 <= fmean <= 400
                and (frac2 or 0) >= 0.9
                and (dur or 1) <= 0.1
            ),
        }
    return {"a_minus_2M": a - 2, "dt_s": dt_s, "proxy": str(s), **res}


def main() -> int:
    tm = M_DET * T_SUN
    fgrid = np.linspace(10.0, 1200.0, 2400)
    W = weights(fgrid)
    # a such that dt spans [0.2, 1.2] s (two or more pulses only for dt <= 0.6 s)
    dts = [0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00]
    avals = [2 + xe.a_for_delay(d / tm) for d in dts]
    rows = []
    for s in (1, 0, "rw"):
        for a in avals:
            rows.append(analyse(a, s, fgrid, W))
    f0, tau = ringdown_mode()
    out = {
        "design": "docs/state/units/X2_s2_design.md (frozen before this run)",
        "observed": {
            "source": "Miani et al. 2023 arXiv:2302.12158v2 Table II (GW150914)",
            "h_MR": 16.0e-23,
            "h_UL_95": 3.4e-23,
            "ratio": RATIO_UL,
            "PMW_s": PMW,
        },
        "remnant": {"M_det_Msun": M_DET, "chi": CHI, "f_220_Hz": f0, "tau_220_s": tau},
        "rows": rows,
        "label": "derived; limit on phenomenological model P1-P4, not on the SV geometry",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s2_gw150914.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for r in rows:
        print(
            r["proxy"],
            f"a-2={r['a_minus_2M']:.2e}",
            f"dt={r['dt_s']:.2f}s",
            {
                k: (
                    round(v["eta_max"], 3) if v["eta_max"] else None,
                    v["applicable"],
                    round(v["first_pulse_mean_f_Hz"]),
                    v["pulses_in_PMW"],
                )
                for k, v in r.items()
                if isinstance(v, dict)
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
