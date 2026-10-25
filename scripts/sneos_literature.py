"""SN Eos 101.2 / 101.1 flux ratios from the published photometry and lens-model magnifications.

Literature check named in ``results/transients/sneos_ratio/PLAN.md`` (run after that measurement
was frozen). Inputs are transcribed from Coulter et al. 2026 (arXiv:2601.04156v2) Tables S1 / S2
and Allingham et al. 2026 (arXiv:2602.14074) Table 3; they are published observed values and
model_predictions, not our measurements.
Output: ``results/transients/sneos_ratio/literature.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

# Coulter et al. 2026 Table S1: (program, MJD, band, m_101.1, err, m_101.2, err), AB mag.
PHOT = [
    (17504, 60380.1, "F814W", 26.40, 0.16, 26.08, 0.12),
    (17504, 60380.1, "F110W", 25.51, 0.05, 24.81, 0.04),
    (17504, 60382.1, "F814W", 26.61, 0.34, 26.30, 0.27),
    (17504, 60382.1, "F110W", 25.37, 0.06, 24.83, 0.04),
    (17504, 60384.1, "F814W", 26.23, 0.25, 25.91, 0.20),
    (17504, 60384.1, "F110W", 25.59, 0.12, 25.04, 0.07),
    (17504, 60386.2, "F814W", 25.63, 0.27, 26.04, 0.35),
    (17504, 60386.2, "F110W", 25.45, 0.15, 24.74, 0.09),
    (17504, 60387.0, "F110W", 25.51, 0.18, 24.57, 0.08),
    (6882, 60919.1, "F277W", 26.80, 0.01, 26.63, 0.03),
    (6882, 60919.1, "F300M", 26.27, 0.01, 26.25, 0.03),
    (6882, 60919.1, "F356W", 25.94, 0.01, 25.81, 0.02),
    (6882, 60919.0, "F410M", 25.27, 0.01, 25.26, 0.01),
    (6882, 60919.1, "F444W", 25.44, 0.01, 25.25, 0.01),
    (9493, 60956.6, "F277W", 27.21, 0.01, 27.11, 0.05),
    (9493, 60956.6, "F356W", 26.20, 0.01, 26.07, 0.03),
    (9493, 60956.6, "F444W", 25.63, 0.01, 25.51, 0.02),
]
# |μ| of 101.1, 101.2 (model_prediction; each quoted with a ~20 % systematic by its authors).
MU_PUBLISHED = {
    "Allingham+2026 Table 3": (25.79, 27.49),
    "Coulter+2026 Table S2": (25.46, 29.88),
}
R_PINNED_MAP = 1.162  # our pinned AstroLensPy v2.1 maps (sneos/summary.json)
R_OURS = (1.382, 0.083)  # R_obs, σ_ln from sneos_ratio/summary.json (difference images)
FLOOR_LN = 0.30  # ASSUMPTION, PLAN.md rule 6
LN10_04 = 0.4 * np.log(10.0)


def ratio(m1: float, e1: float, m2: float, e2: float) -> tuple[float, float]:
    """F(101.2) / F(101.1) from AB magnitudes and its σ in ln."""
    return float(10 ** (-0.4 * (m2 - m1))), float(LN10_04 * np.hypot(e1, e2))


def weighted(ln_r: np.ndarray, s: np.ndarray) -> dict:
    """Inverse-variance mean of ln R, its σ and χ² about the mean."""
    w = 1.0 / s**2
    mu = float((ln_r * w).sum() / w.sum())
    return {
        "r": float(np.exp(mu)),
        "sigma_ln": float(w.sum() ** -0.5),
        "chi2": float(((ln_r - mu) ** 2 * w).sum()),
        "dof": int(len(ln_r) - 1),
    }


def difference_ratio(early: tuple, late: tuple) -> tuple[float, float]:
    """Ratio of the 101.2 and 101.1 flux drops between two epochs (as in our difference images)."""

    def flux(m):
        return 10 ** (-0.4 * m)

    d1 = flux(early[3]) - flux(late[3])
    d2 = flux(early[5]) - flux(late[5])
    s1 = LN10_04 * np.hypot(flux(early[3]) * early[4], flux(late[3]) * late[4])
    s2 = LN10_04 * np.hypot(flux(early[5]) * early[6], flux(late[5]) * late[6])
    return float(d2 / d1), float(np.hypot(s1 / d1, s2 / d2))


def run(out: Path) -> dict:
    groups: dict[str, list] = {}
    entries = []
    for row in PHOT:
        r, s = ratio(*row[3:])
        entries.append({"program": row[0], "mjd": row[1], "band": row[2], "r": r, "sigma_ln": s})
        key = f"HST {row[2]}" if row[0] == 17504 else "JWST 2025"
        groups.setdefault(key, []).append((np.log(r), s))
    combined = {k: weighted(*map(np.array, zip(*v, strict=True))) for k, v in groups.items()}
    by_band = {(r[0], r[2]): r for r in PHOT}
    diff = {}
    for band in ("F277W", "F356W", "F444W"):
        r, s = difference_ratio(by_band[(6882, band)], by_band[(9493, band)])
        diff[band] = {"r": r, "sigma_rel": s}
    r_model = {k: v[1] / v[0] for k, v in MU_PUBLISHED.items()}
    r_model["pinned map"] = R_PINNED_MAP
    lo, hi = min(r_model.values()), max(r_model.values())
    offsets = {
        k: {
            "ln_offset_from_nearest_model": float(np.log(v["r"]) - np.log(np.clip(v["r"], lo, hi))),
            "beyond_floor": bool(abs(np.log(v["r"]) - np.log(np.clip(v["r"], lo, hi))) > FLOOR_LN),
        }
        for k, v in combined.items()
    }
    summary = {
        "sources": ["arXiv:2601.04156v2 Tables S1, S2", "arXiv:2602.14074 Table 3"],
        "provenance": {
            "photometry": "observed (published)",
            "magnifications": "model_prediction (published)",
            "ratios": "derived",
        },
        "entries": entries,
        "combined": combined,
        "difference_ratio_6882_minus_9493": diff,
        "r_model": r_model,
        "r_model_range": [lo, hi],
        "offset_vs_model_range": offsets,
        "ours_difference_images": {"r": R_OURS[0], "sigma_ln": R_OURS[1]},
        "floor_ln_assumption": FLOOR_LN,
    }
    out.mkdir(parents=True, exist_ok=True)
    Path(out, "literature.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


if __name__ == "__main__":
    s = run(Path(__file__).resolve().parents[1] / "results" / "transients" / "sneos_ratio")
    print(
        json.dumps(
            {"combined": s["combined"], "diff": s["difference_ratio_6882_minus_9493"]}, indent=1
        )
    )
