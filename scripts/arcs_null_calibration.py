"""`arcs` shape-noise null calibrated on the weak-shear background rows of the stored runs.

Plan: results/arcs/null-calibration/PLAN.md (frozen before the weak-row orientations were read).
Fits σ_ε, or an isotropic fraction f at σ_ε 0.26, to the weak-shear (``strong_shear`` False)
background rows of the six #406 fields, then re-computes the strong-shear expectations of the
pooled fields and Abell S1063 under those calibrations. No network, no new data.

    python scripts/arcs_null_calibration.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_pooled as ap  # noqa: E402

OUT = ap.ROOT / "null-calibration"
SIGMA0 = 0.26
SIGMA_GRID = np.round(np.arange(0.15, 0.6001, 0.01), 2)  # plan: [0.15, 0.60], step 0.01
F_GRID = np.round(np.arange(0.0, 1.0001, 0.01), 2)  # plan: [0, 1], step 0.01
CHI2_95 = 3.84  # Δ(−2 ln L) for a 95 % profile interval, one parameter
P_ISO = 1.0 / 3.0  # P(axis ≥ 60° from the shear) for a uniform axis angle


def rows(path: Path, strong: bool) -> Table:
    """Background rows of one run, strong-shear or weak-shear (same cuts otherwise)."""
    t = Table.read(path)
    t = t[t["population"] == "background"]
    return t[t["strong_shear"] if strong else ~t["strong_shear"]]


def p_anti(t: Table, sigma: float, e_min: float, anti_deg: float) -> np.ndarray:
    """Per-row shape-noise P(anti) at round(g_min, 2), as `lens_consistency.anti_expected`."""
    return np.array([ap._p_anti(round(float(g), 2), sigma, e_min, anti_deg) for g in t["g_min"]])


def with_iso(p: np.ndarray, f: float) -> np.ndarray:
    """Mixture: a fraction f of rows has an orientation unrelated to the model shear."""
    return (1.0 - f) * p + f * P_ISO


def neg2lnl(p: np.ndarray, a: np.ndarray) -> float:
    """−2 ln Bernoulli likelihood of the anti flags ``a`` given per-row probabilities ``p``."""
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return float(-2.0 * np.sum(np.where(a, np.log(p), np.log1p(-p))))


def profile(grid: np.ndarray, curve: np.ndarray) -> dict:
    """Grid MLE and the 95 % interval (Δ(−2 ln L) ≤ 3.84); a bound at the grid edge is open."""
    i = int(np.argmin(curve))
    ok = grid[curve - curve[i] <= CHI2_95]
    return {
        "best": float(grid[i]),
        "lo95": float(ok.min()),
        "hi95": float(ok.max()),
        "lo_at_edge": bool(ok.min() == grid[0]),
        "hi_at_edge": bool(ok.max() == grid[-1]),
    }


def two_sided_p(k: int, lam: float) -> float:
    """Two-sided Poisson P: twice the smaller tail, capped at 1."""
    return float(min(1.0, 2 * min(poisson.cdf(k, lam), poisson.sf(k - 1, lam))))


def tally(k: int, lam: float) -> dict:
    return {"anti": int(k), "expected": float(lam), "p_ge_observed": float(poisson.sf(k - 1, lam))}


def outcome(p_b_best: float, p_b_upper: float, alpha: float = ap.ALPHA) -> str:
    """Frozen plan rule. ``p_b_best``: pooled B's larger P at σ̂_ε or f̂; ``p_b_upper``: its larger
    P at the 95 % upper bounds of σ_ε and f."""
    if p_b_best >= alpha:
        return "calibration explains the pooled excess"
    if p_b_upper < alpha:
        return "null calibrated; excess stands"
    return "weak rows lack power"


def main() -> None:
    meta = json.loads((ap.ROOT / ap.S1063).parent.joinpath("arcs.json").read_text())["assumptions"]
    e_min, anti_deg = float(meta["min_ellipticity"]), float(meta["anti_deg"])
    paths = {**{f: ap.ROOT / p for f, p in ap.FIELDS.items()}, "S1063": ap.ROOT / ap.S1063}
    for f, p in paths.items():
        if json.loads(p.parent.joinpath("arcs.json").read_text())["assumptions"] != meta:
            raise SystemExit(f"error: {f} uses different `arcs` assumptions")
    weak = {f: rows(p, strong=False) for f, p in paths.items()}
    strong = {f: rows(p, strong=True) for f, p in paths.items()}

    def pooled_rows(tabs: dict[str, Table], fields: tuple[str, ...]) -> Table:
        cols = [tabs[f]["g_min", "orientation_class"] for f in fields]
        return vstack(cols, metadata_conflicts="silent")

    def p_at(t: Table, sigma: float, f: float = 0.0) -> np.ndarray:
        return with_iso(p_anti(t, sigma, e_min, anti_deg), f)

    calib = {}
    for v, fields in ap.VARIANTS.items():
        w = pooled_rows(weak, fields)
        a = ap.anti_mask(w)
        p0 = p_at(w, SIGMA0)
        sig = profile(SIGMA_GRID, np.array([neg2lnl(p_at(w, s), a) for s in SIGMA_GRID]))
        fis = profile(F_GRID, np.array([neg2lnl(with_iso(p0, f), a) for f in F_GRID]))
        calib[v] = {
            "weak_n": len(w),
            "weak_anti": int(a.sum()),
            "weak_expected_sigma0": float(p0.sum()),
            "weak_p_two_sided_sigma0": two_sided_p(int(a.sum()), float(p0.sum())),
            "sigma_eps": sig,
            "f_iso_at_sigma0": fis,
        }

    def strong_under(t: Table, sigma: float, f: float) -> dict:
        return tally(int(ap.anti_mask(t).sum()), float(p_at(t, sigma, f).sum()))

    propagated = {}
    for v, fields in ap.VARIANTS.items():
        c = calib[v]
        settings = {
            "sigma0": (SIGMA0, 0.0),
            "sigma_best": (c["sigma_eps"]["best"], 0.0),
            "sigma_hi95": (c["sigma_eps"]["hi95"], 0.0),
            "f_best": (SIGMA0, c["f_iso_at_sigma0"]["best"]),
            "f_hi95": (SIGMA0, c["f_iso_at_sigma0"]["hi95"]),
        }
        s = pooled_rows(strong, fields)
        propagated[v] = {
            name: {
                "sigma_eps": sg,
                "f_iso": f,
                "pooled": strong_under(s, sg, f),
                "S1063": strong_under(strong["S1063"], sg, f),
            }
            for name, (sg, f) in settings.items()
        }
    pb = propagated["B"]
    best = max(pb["sigma_best"]["pooled"]["p_ge_observed"], pb["f_best"]["pooled"]["p_ge_observed"])
    upper = max(
        pb["sigma_hi95"]["pooled"]["p_ge_observed"], pb["f_hi95"]["pooled"]["p_ge_observed"]
    )
    w63 = weak["S1063"]
    summary = {
        "plan": "results/arcs/null-calibration/PLAN.md",
        "label": "derived; pre-registered for the weak rows (the strong-shear excess was known)",
        "fields": ap.FIELDS,
        "s1063": ap.S1063,
        "assumptions": meta,
        "calibration": calib,
        "propagated": propagated,
        "outcome": outcome(best, upper),
        "s1063_weak": {
            "n": len(w63),
            "anti": int(ap.anti_mask(w63).sum()),
            "expected_sigma0": float(p_at(w63, SIGMA0).sum()),
        },
        "weak_by_field": {
            f: {
                "n": len(t),
                "anti": int(ap.anti_mask(t).sum()),
                "expected_sigma0": float(p_at(t, SIGMA0).sum()),
            }
            for f, t in weak.items()
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("calibration", "outcome", "s1063_weak")}, indent=1))
    for v, d in propagated.items():
        for name, r in d.items():
            pool, s63 = r["pooled"], r["S1063"]
            print(
                f"{v} {name:10s} σ {r['sigma_eps']:.2f} f {r['f_iso']:.2f}: pooled "
                f"{pool['anti']} vs {pool['expected']:.2f} P {pool['p_ge_observed']:.4f} | "
                f"S1063 {s63['anti']} vs {s63['expected']:.2f} P {s63['p_ge_observed']:.4f}"
            )


if __name__ == "__main__":
    main()
