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
# Applicability ASSUMPTIONs (design conditions a-d; docs/state/units/X2_s2_design.md)
PARAMS = {
    # P4 timing offset (<~10 ms) plus half a pulse: a pulse must sit this far inside the PMW
    "edge_margin_s": 0.015,
    "f_validated_Hz": (80.0, 400.0),  # Miani et al. sec. IV.1 morphology range
    "min_frac_first_two": 0.9,  # validated two-pulse morphology
    # transform-limited duration (a lower bound; S3 measures the envelope)
    "max_duration_lower_bound_s": 0.1,
}


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
    m = PARAMS["edge_margin_s"]
    ks = [k for k in range(1, 50) if PMW[0] + m <= k * dt_s <= PMW[1] - m]
    edge = [k for k in range(1, 50) if abs(k * dt_s - PMW[0]) < m or abs(k * dt_s - PMW[1]) < m]
    res = {}
    for name, wt in W.items():
        norm = np.trapezoid(wt, fgrid)
        e = [np.trapezoid(wt * T2 * R2 ** (2 * k - 1), fgrid) / norm for k in ks]
        F = math.sqrt(sum(e)) if e else 0.0
        # echo-weighted spectrum over all in-window pulses (design condition b)
        spec = sum(wt * T2 * R2 ** (2 * k - 1) for k in ks) if ks else wt * T2 * R2
        fmean = float(np.trapezoid(fgrid * spec, fgrid) / np.trapezoid(spec, fgrid))
        # rms spectral width -> transform-limited duration 1 / (2 pi sigma_f): a LOWER bound on the
        # pulse duration (the reflection phase can only stretch it; S3 measures the real envelope)
        fvar = float(np.trapezoid((fgrid - fmean) ** 2 * spec, fgrid) / np.trapezoid(spec, fgrid))
        dur = 1.0 / (2 * math.pi * math.sqrt(fvar)) if fvar > 0 else None
        frac2 = (sum(e[:2]) / sum(e)) if e else None
        res[name] = {
            "F": F,
            "eta_max": RATIO_UL / F if F > 0 else None,
            "pulses_in_PMW": ks,
            "pulses_on_PMW_edge_excluded": edge,
            "echo_weighted_mean_f_Hz": fmean,
            "duration_lower_bound_s": dur,
            "first_two_pulses_fraction_of_hrss2": frac2,
            "applicable": bool(
                ks
                and ks[0] == 1
                and len(ks) >= 2
                and not edge
                and PARAMS["f_validated_Hz"][0] <= fmean <= PARAMS["f_validated_Hz"][1]
                and (frac2 or 0) >= PARAMS["min_frac_first_two"]
                and (dur or 1) <= PARAMS["max_duration_lower_bound_s"]
            ),
        }
    return {"a_minus_2M": a - 2, "dt_s": dt_s, "proxy": str(s), **res}


def main() -> int:
    tm = M_DET * T_SUN
    fgrid = np.linspace(10.0, 1200.0, 2400)
    W = weights(fgrid)
    # a such that dt spans [0.2, 1.2] s (two or more pulses only for dt <= 0.6 s)
    # grid kept off the PMW edges (0.2, 1.2 s and their sub-multiples 0.3, 0.4, 0.6 s)
    dts = [0.22, 0.27, 0.32, 0.37, 0.45, 0.52, 0.58, 0.80, 1.00]
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
        "params_assumptions": PARAMS,
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
                    round(v["echo_weighted_mean_f_Hz"]),
                    v["pulses_in_PMW"],
                )
                for k, v in r.items()
                if isinstance(v, dict)
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
