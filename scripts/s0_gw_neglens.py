"""S0-2026-10-15 direction A, minimal step: wave-optics GW lensing by a negative-mass point lens.

Unit file: docs/state/units/S0-2026-10-15.md. Everything here is `derived` or `simulated`; no strain
data is used.

Amplification factor of a point lens of either sign (s = +1 ordinary mass, s = -1 negative mass),
Fermat potential T(x) = |x - y|^2 / 2 - s ln|x| in Einstein units, dimensionless frequency
w = 8 pi G M_Lz f / c^3:

    F_s(w, y) = exp[s pi w / 4 + i (s w / 2) ln(w / 2)]
                * Gamma(1 - i s w / 2) 1F1(i s w / 2; 1; i w y^2 / 2)

The s = +1 case is the standard point-mass result (e.g. Takahashi & Nakamura 2003; the overall time
constant is dropped). The s = -1 case is checked against direct numerical integration of the
diffraction integral (contour rotated by pi/4) in ``tests/test_s0_gw_neglens.py``, and its
geometric-optics limit against the W3 magnification.

Distinguishability (noise-weighted overlap; amplitude, phase and arrival time free; intrinsic
parameters FIXED, so the mismatches are upper bounds on what a full parameter estimation would see):
  1 - O(h F_neg, h)            negative lens vs unlensed;
  min 1 - O(h F_neg, h F_pos)  negative lens vs the best point-mass lens on a (M_Lz, y) grid.
Required SNR for Delta chi^2 = 2 rho^2 (1 - O) >= 9:  rho_req = sqrt(4.5 / (1 - O)).

ASSUMPTIONS: inspiral-only amplitude |h| ~ f^-7/6 from 20 Hz to f_ISCO = c^3 / (6^1.5 pi G M)
(merger and ringdown omitted, so high-frequency power is underweighted); an aLIGO design PSD fit,
S(x) = 1e-49 [x^-4.14 - 5 x^-2 + 111 (1 - x^2 + x^4/2) / (1 + x^2/2)], x = f / 215 Hz,
attributed to Ajith 2011 (PRD 84, 084037, arXiv:1107.1267) but not checked against the paper
body; a flat PSD is run as a check.

  python scripts/s0_gw_neglens.py      # writes results/s0/gw_neglens.json
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "s0"
W_PER_MSUN_HZ = 8.0 * np.pi * 4.925490947e-6  # 8 pi G M_sun / c^3 in s
F_LOW = 20.0
DF = 1.0
SOURCES = {"M20": 20.0, "M60": 60.0}  # detector-frame total mass of the binary (Msun)
NEG_CASES = [(m, y) for m in (30.0, 100.0, 300.0, 1000.0) for y in (0.3, 1.0, 1.5, 2.5)]
POS_M = np.logspace(0.5, 4.0, 36)
POS_Y = np.linspace(0.05, 3.0, 30)


def amp_factor(w: float, y: float, s: int) -> complex:
    if w <= 0:
        return 1.0 + 0.0j
    w_, y_ = mp.mpf(w), mp.mpf(y)
    val = (
        mp.exp(s * mp.pi * w_ / 4 + 1j * (s * w_ / 2) * mp.log(w_ / 2))
        * mp.gamma(1 - 1j * s * w_ / 2)
        * mp.hyp1f1(1j * s * w_ / 2, 1, 1j * w_ * y_ * y_ / 2)
    )
    return complex(val)


def lens_curve(freqs: np.ndarray, m_lz: float, y: float, s: int) -> np.ndarray:
    return np.array([amp_factor(W_PER_MSUN_HZ * m_lz * f, y, s) for f in freqs])


def psd(freqs: np.ndarray, kind: str) -> np.ndarray:
    if kind == "flat":
        return np.ones_like(freqs)
    x = freqs / 215.0
    return 1e-49 * (x**-4.14 - 5 * x**-2 + 111 * (1 - x**2 + x**4 / 2) / (1 + x**2 / 2))


def overlap(a: np.ndarray, b: np.ndarray, weight: np.ndarray, pad: int = 16) -> float:
    """max over time shift and phase of |<a, b>| / sqrt(<a,a><b,b>), <a,b> = sum a b* weight."""
    prod = a * np.conj(b) * weight
    n = pad * prod.size
    corr = np.fft.ifft(prod, n) * n
    return float(
        np.max(np.abs(corr))
        / np.sqrt(np.sum(np.abs(a) ** 2 * weight) * np.sum(np.abs(b) ** 2 * weight))
    )


def main() -> int:
    out = {
        "label": "derived (closed-form wave optics) / simulated (fiducial inspiral); no data",
        "assumptions": __doc__.split("ASSUMPTIONS:")[1].split("python scripts")[0].strip(),
        "criterion": (
            "rho_req = sqrt(4.5 / (1 - O)) for Delta chi^2 = 9; intrinsic parameters fixed"
        ),
        "sources": {},
    }
    for name, mtot in SOURCES.items():
        f_isco = 1.0 / (6**1.5 * np.pi * mtot * 4.925490947e-6)
        freqs = np.arange(F_LOW, f_isco, DF)
        h = freqs ** (-7.0 / 6.0)
        pos = {(m, y): lens_curve(freqs, m, y, +1) for m in POS_M for y in POS_Y}
        rows = []
        for kind in ("aLIGO", "flat"):
            wgt = 1.0 / psd(freqs, kind)
            for m, y in NEG_CASES:
                fneg = lens_curve(freqs, m, y, -1)
                a = h * fneg
                mm_unl = 1.0 - overlap(a, h, wgt)
                best, arg = 1.0, None
                for key, fpos in pos.items():
                    mm = 1.0 - overlap(a, h * fpos, wgt)
                    if mm < best:
                        best, arg = mm, key
                gain = np.abs(fneg)
                rows.append(
                    {
                        "psd": kind,
                        "M_Lz": m,
                        "y": y,
                        "w_range": [
                            float(W_PER_MSUN_HZ * m * freqs[0]),
                            float(W_PER_MSUN_HZ * m * freqs[-1]),
                        ],
                        "abs_F_range": [float(gain.min()), float(gain.max())],
                        "mismatch_vs_unlensed": mm_unl,
                        "rho_req_vs_unlensed": float(np.sqrt(4.5 / max(mm_unl, 1e-12))),
                        "mismatch_vs_best_positive": best,
                        "rho_req_vs_positive": float(np.sqrt(4.5 / max(best, 1e-12))),
                        "best_positive": {"M_Lz": float(arg[0]), "y": float(arg[1])},
                    }
                )
                r = rows[-1]
                lo, hi = r["abs_F_range"]
                print(
                    f"{name} {kind:5s} M_L={m:6.0f} y={y:3.1f} |F| {lo:.3f}-{hi:.3f}"
                    f"  vs unlensed 1-O={mm_unl:.2e} (rho {r['rho_req_vs_unlensed']:.1f})"
                    f"  vs +mass 1-O={best:.2e} (rho {r['rho_req_vs_positive']:.1f})"
                    f" at M={arg[0]:.0f} y={arg[1]:.2f}"
                )
        out["sources"][name] = {
            "M_total_det": mtot,
            "f_range_hz": [F_LOW, float(f_isco)],
            "cases": rows,
        }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gw_neglens.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
