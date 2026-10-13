"""X1 S2: far-side vs ordinary hot-spot discriminant at Sgr A* (docs/state/units/X1.md).

Hypothesis tested (conditional geometry, premises in X1.md): a NIR flare is a Keplerian hot spot
on the far side of a Simpson-Visser throat (a in {2.1, 3, 5, 7, 9} M). Ordinary model: a
same-side Schwarzschild hot spot (a = 0) with eta = v / v_K in [0.6, 1]. Both orbits are fixed by
the polarisation period P = 60 +- 3 min (premise P_pol = P_orb), i = 155 deg (GRAVITY 2023), and
both models get the same nuisance static component (fraction 0-0.5 of the mean total flux, offset
0-20 uas in any direction). Light-travel-time delays are included; images k = 0 (primary +
secondary); bolometric g^4 weights.

Observables, each with an ASSUMPTION for the achievable precision:
- loop radius L: mean distance of the total-light centroid from its own mean (free centre, as the
  centre is degenerate in GRAVITY fits); precision 9 uas (the 4-flare averaged loop, GRAVITY 2023);
- flux contrast C = ln(F_max / F_min) of the total light; precision 0.05 (few-% light curves).
- flux peak count per orbit (maxima with prominence > 2 % of the mean flux; ASSUMPTION); static
  light does not change it.
Gate (pre-written in X1.md): some observable's far-side vs ordinary gap, minimised over both
models' nuisance ranges, exceeds its precision.

  python scripts/x1_s2.py   -> results/x1/s2_discriminant.json

Every number is a model computation (derived, ours) of an unestablished geometry; nothing here is
evidence for or against a wormhole at Sgr A*.
"""

from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_farside as wf  # noqa: E402
import x1_tracer as xt  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results" / "x1"
M_SUN = wf.M_S2[0]
THETA_G = wf.theta_g(M_SUN)  # uas per M
INCL = 155.0  # deg, GRAVITY 2023
PERIODS = (57.0, 60.0, 63.0)  # min, P_pol = 60 +- 3 (GRAVITY 2023)
A_FAR = (2.1, 3.0, 5.0, 7.0, 9.0)
ETAS = (0.6, 0.7, 0.8, 0.9, 1.0)
# ASSUMPTION: nuisance ranges and precisions, pre-written in X1.md (S2 row, "Tolerances")
STATIC_FRAC = np.linspace(0.0, 0.5, 6)
STATIC_OFF_UAS = (0.0, 5.0, 10.0, 15.0, 20.0)
STATIC_DIR = np.linspace(0.0, 2 * math.pi, 8, endpoint=False)
PREC = {"loop_uas": 9.0, "ln_flux_contrast": 0.05}
N_PHASE = 72


def observables(track: np.ndarray) -> np.ndarray:
    """(L [uas], C) for every static-component nuisance setting; rows over the nuisance grid."""
    f = track[:, 3] / track[:, 3].mean()
    c = track[:, 1:3] * THETA_G
    out = []
    for fs in STATIC_FRAC:
        s = fs / (1.0 - fs)
        for d in STATIC_OFF_UAS:
            for psi in STATIC_DIR if d > 0 else (0.0,):
                off = d * np.array([math.cos(psi), math.sin(psi)])
                tot = f + s
                cen = (f[:, None] * c + s * off) / tot[:, None]
                loop = np.hypot(*(cen - cen.mean(0)).T).mean()
                out.append((loop, math.log(tot.max() / tot.min())))
    return np.array(out)


def flux_peaks(flux: np.ndarray, rel: float = 0.02) -> int:
    """Local maxima of the periodic flux curve with prominence above rel x the mean flux."""
    n, m, count = len(flux), flux.mean(), 0
    for i in range(n):
        if flux[i] >= flux[i - 1] and flux[i] > flux[(i + 1) % n]:
            left = min(flux[(i - k) % n] for k in range(1, n // 2))
            right = min(flux[(i + k) % n] for k in range(1, n // 2))
            count += flux[i] - max(left, right) > rel * m
    return int(count)


def track(a: float, re: float, side: str, eta: float) -> tuple[np.ndarray, int]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tr = xt.spot_track(a, re, INCL, side, n_phase=N_PHASE, eta=eta, delays=True)
    return tr, len(caught)


def ln_contrast(tr: np.ndarray, fs: float) -> float:
    f = tr[:, 3] / tr[:, 3].mean()
    st = fs / (1.0 - fs)
    return math.log((f.max() + st) / (f.min() + st))


def fragility(a: float) -> dict:
    """Sensitivity of a flux-contrast pass: static-fraction cap that closes the gap; inclination."""
    trs = {}
    for p in PERIODS:
        rk = wf.kepler_radius(p, M_SUN)
        trs[p] = (
            [track(0.0, e ** (2.0 / 3.0) * rk, "near", e)[0] for e in ETAS],
            track(a, rk, "far", 1.0)[0],
        )
    far_max = max(ln_contrast(f, 0.0) for _, f in trs.values())
    closes = None
    for fs in np.arange(0.50, 0.95, 0.01):
        o_min = min(ln_contrast(t, fs) for o, _ in trs.values() for t in o)
        if o_min - far_max <= PREC["ln_flux_contrast"]:
            closes = round(float(fs), 2)
            break
    incl = []
    rk = wf.kepler_radius(60.0, M_SUN)
    for i in (135.0, 145.0, 155.0, 165.0):
        far_c = ln_contrast(xt.spot_track(a, rk, i, "far", n_phase=N_PHASE, delays=True), 0.0)
        ord_c = min(
            ln_contrast(
                xt.spot_track(
                    0.0, e ** (2.0 / 3.0) * rk, i, "near", n_phase=N_PHASE, eta=e, delays=True
                ),
                float(STATIC_FRAC.max()),
            )
            for e in ETAS
        )
        incl.append({"incl_deg": i, "far_C_fs0": far_c, "ordinary_min_C_fs05": ord_c})
    return {
        "far_max_C": far_max,
        "gap_closes_at_static_fraction": closes,
        "incl_P60": incl,
        "note": "a flux-contrast pass needs the orbital modulation separated from intrinsic flare"
        " evolution to the assumed 0.05 precision",
    }


def main() -> int:
    ordinary, far = [], []
    for p in PERIODS:
        rk = wf.kepler_radius(p, M_SUN)
        for eta in ETAS:
            re = eta ** (2.0 / 3.0) * rk  # same observed period, slower orbit -> smaller R
            tr, nw = track(0.0, re, "near", eta)
            ordinary.append(
                {"P_min": p, "eta": eta, "R_M": re, "obs": observables(tr), "warn": nw}
                | {"flux_peaks": flux_peaks(tr[:, 3])}
            )
        for a in A_FAR:
            if rk <= a:
                far.append({"P_min": p, "a_M": a, "R_M": rk, "skipped": "orbit inside throat"})
                continue
            tr, nw = track(a, rk, "far", 1.0)
            far.append(
                {"P_min": p, "a_M": a, "R_M": rk, "obs": observables(tr), "warn": nw}
                | {"flux_peaks": flux_peaks(tr[:, 3])}
            )
    o_all = np.vstack([o["obs"] for o in ordinary])
    rows = []
    for a in A_FAR:
        fa = [f for f in far if f["a_M"] == a and "obs" in f]
        if not fa:
            rows.append({"a_M": a, "status": "no far-side orbit outside the throat"})
            continue
        f_all = np.vstack([f["obs"] for f in fa])
        row = {"a_M": a, "periods_used": [f["P_min"] for f in fa]}
        gate = False
        for j, key in enumerate(("loop_uas", "ln_flux_contrast")):
            lo_o, hi_o = float(o_all[:, j].min()), float(o_all[:, j].max())
            lo_f, hi_f = float(f_all[:, j].min()), float(f_all[:, j].max())
            gap = max(0.0, lo_o - hi_f, lo_f - hi_o)  # interval separation over all nuisances
            row[key] = {
                "ordinary_range": [lo_o, hi_o],
                "far_range": [lo_f, hi_f],
                "gap": gap,
                "gap_over_precision": gap / PREC[key],
            }
            gate |= gap > PREC[key]
        o_pk = sorted({o["flux_peaks"] for o in ordinary})
        f_pk = sorted({f["flux_peaks"] for f in fa})
        row["flux_peaks"] = {"ordinary": o_pk, "far": f_pk, "disjoint": not set(o_pk) & set(f_pk)}
        gate |= row["flux_peaks"]["disjoint"]
        row["gate_pass"] = bool(gate)
        if gate:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                row["fragility"] = fragility(a)
        rows.append(row)

    def summ(items):
        return [
            {k: v for k, v in it.items() if k != "obs"}
            | (
                {
                    "loop_uas_range": [float(it["obs"][:, 0].min()), float(it["obs"][:, 0].max())],
                    "ln_contrast_range": [
                        float(it["obs"][:, 1].min()),
                        float(it["obs"][:, 1].max()),
                    ],
                }
                if "obs" in it
                else {}
            )
            for it in items
        ]

    out = {
        "labels": "derived (our computation) of a conditional, unestablished geometry; precisions"
        " and nuisance ranges are ASSUMPTIONs pre-written in docs/state/units/X1.md",
        "inputs": {
            "M_sun": M_SUN,
            "theta_g_uas": THETA_G,
            "incl_deg": INCL,
            "periods_min": PERIODS,
            "precision": PREC,
            "static_fraction": STATIC_FRAC.tolist(),
            "static_offset_uas": STATIC_OFF_UAS,
            "n_phase": N_PHASE,
            "model": "point spot, k = 0 images, bolometric g^4, light-travel delays included",
        },
        "discriminant": rows,
        "ordinary_tracks": summ(ordinary),
        "far_tracks": summ(far),
        "gate_pass_any_a": any(r.get("gate_pass", False) for r in rows),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s2_discriminant.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(rows, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
