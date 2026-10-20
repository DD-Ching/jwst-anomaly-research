"""Pooled `arcs` anti count over the Abell S1063 replication fields (exploratory, post hoc).

Plan: results/arcs/pooled-tally/PLAN.md. Reads each field's stored deciding-run ``arcs.ecsv``,
recomputes the shape-noise expectation with ``lens_consistency.anti_expected`` (checked against the
stored summary), pools the counts and scans σ_ε. No network, no new data.

    python scripts/arcs_pooled.py
"""

from __future__ import annotations

import json
import sys
from functools import cache
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import brentq
from scipy.stats import binom, chi2, poisson

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1] / "results" / "arcs"
OUT = ROOT / "pooled-tally"
# field -> deciding run's arcs.ecsv (each field's RESULT.md / summary.json)
FIELDS = {
    "MACS0417": "macs0417-relics/run2-canucs-v1/arcs.ecsv",
    "MACS1423": "macs1423-canucs/macs1423-canucs-maps/arcs.ecsv",
    "WHL0137": "whl0137-relics/whl0137-relics-lenstool/arcs.ecsv",
    "RXC J0600": "rxc0600-relics/rxc0600-relics-glafic/arcs.ecsv",
    "MACS0308": "macs0308-relics/macs0308-relics-glafic/arcs.ecsv",
    "RXC J0018": "rxc0018-relics/rxc0018-relics-lenstool/arcs.ecsv",
}
# variant B drops MACS1423: the clump-pair census explains its screen excess (CHANGELOG 2026-10-16)
VARIANTS = {"A": tuple(FIELDS), "B": tuple(f for f in FIELDS if f != "MACS1423")}
S1063 = "abells1063-cats/arcs.ecsv"
ALPHA = 0.05  # ASSUMPTION (plan)
SIGMAS = (0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40)


def background(path: Path) -> Table:
    """Strong-shear background rows of one run, as `lens_consistency` counts them."""
    t = Table.read(path)
    t = t[t["strong_shear"]]
    return t[t["population"] == "background"]


@cache
def _p_anti(g: float, sigma_eps: float, min_ellipticity: float, anti_deg: float) -> float:
    return lc.anti_probability(g, sigma_eps, min_ellipticity, anti_deg)


def field_stats(rows: Table, sigma_eps: float, min_ellipticity: float, anti_deg: float) -> dict:
    """Like `lens_consistency.anti_expected` (same rounding of g_min), but λ is kept unrounded."""
    lam = sum(
        _p_anti(round(float(g), 2), sigma_eps, min_ellipticity, anti_deg)
        for g in np.asarray(rows["g_min"], float)
    )
    if not len(rows) or not np.isfinite(lam) or lam <= 0:
        raise ValueError(f"no usable shape-noise expectation (n {len(rows)}, λ {lam})")
    k = int(np.sum(anti_mask(rows)))
    return {
        "n": len(rows),
        "anti": k,
        "expected": lam,
        "p_ge_observed": float(poisson.sf(k - 1, lam)),
    }


def pooled(stats: dict[str, dict], fields: tuple[str, ...]) -> dict:
    n = sum(stats[f]["anti"] for f in fields)
    lam = sum(stats[f]["expected"] for f in fields)
    return {"anti": n, "expected": lam, "p_ge_observed": float(poisson.sf(n - 1, lam))}


def outcome(p_a: float, p_b: float, alpha: float = ALPHA) -> str:
    """Frozen plan rule: both variants below alpha, both at or above it, or mixed."""
    if p_a < alpha and p_b < alpha:
        return "pooled excess"
    if p_a >= alpha and p_b >= alpha:
        return "no pooled excess"
    return "mixed"


def ratio_interval(n: int, lam: float, cl: float = 0.95) -> tuple[float, float]:
    """Central Poisson (Garwood) interval on n over λ: range of the observed / expected ratio."""
    a = 1 - cl
    lo = chi2.ppf(a / 2, 2 * n) / 2 if n > 0 else 0.0
    hi = chi2.ppf(1 - a / 2, 2 * n + 2) / 2
    return float(lo / lam), float(hi / lam)


def s1063_vs_pooled(s1063: dict, pool: dict) -> dict:
    """Post hoc: is S1063's anti count in line with the pooled fields' anti rate per expected count?

    ``scaled``: S1063's λ times the pooled observed / expected ratio (a plug-in point value; its
    uncertainty is not propagated), P(>= observed), also at the ratio's 95 % bounds.
    ``conditional``: given the combined count, P(S1063 holds >= its share) with binomial
    p = λ_S1063 / (λ_S1063 + λ_pool). This cancels a common factor in λ, but a change of σ_ε does
    not scale every field's λ by the same factor (see the σ_ε scan).
    """
    ratio = pool["anti"] / pool["expected"]
    lo, hi = ratio_interval(pool["anti"], pool["expected"])
    k = s1063["anti"]
    total = k + pool["anti"]
    share = s1063["expected"] / (s1063["expected"] + pool["expected"])
    return {
        "pool_ratio": ratio,
        "pool_ratio_95": [lo, hi],
        "scaled_expected": ratio * s1063["expected"],
        "p_scaled": float(poisson.sf(k - 1, ratio * s1063["expected"])),
        "p_scaled_at_ratio_95": [float(poisson.sf(k - 1, r * s1063["expected"])) for r in (lo, hi)],
        "p_conditional": float(binom.sf(k - 1, total, share)),
    }


def anti_mask(rows: Table) -> np.ndarray:
    return np.asarray(rows["orientation_class"] == "anti")


def inside(rows: Table) -> np.ndarray:
    """Inside the critical curve as `anti_probability` treats it: |round(g_min, 2)| > 1."""
    return np.abs(np.round(np.asarray(rows["g_min"], float), 2)) > 1


def sigma_crossing(p_of_sigma, lo: float, hi: float, alpha: float = ALPHA) -> float | None:
    """σ_ε in [lo, hi] at which a pooled P reaches alpha (brentq), or None without a sign change."""

    def f(s: float) -> float:
        return p_of_sigma(s) - alpha

    if f(lo) * f(hi) > 0:
        return None
    return float(brentq(f, lo, hi, xtol=1e-3))


def load_run(rel: str) -> tuple[Table, dict]:
    """Background rows plus the run's stored `arcs.json` (assumptions and shape-noise summary)."""
    path = ROOT / rel
    return background(path), json.loads((path.parent / "arcs.json").read_text())


def stored_check(stats: dict, stored: dict) -> dict:
    """Recomputed vs stored shape-noise summary (stored `expected` is rounded to 3 decimals)."""
    same = (
        stats["n"] == stored["n"]
        and stats["anti"] == stored["anti"]
        and abs(stats["expected"] - stored["expected"]) <= 5e-4
    )
    return {"stored": stored, "matches": bool(same)}


def main() -> None:
    runs = {f: load_run(p) for f, p in FIELDS.items()}
    s1063_rows, s1063_meta = load_run(S1063)
    assumptions = {name: meta["assumptions"] for name, (_, meta) in runs.items()}
    assumptions["S1063"] = s1063_meta["assumptions"]
    ref = next(iter(assumptions.values()))
    if any(a != ref for a in assumptions.values()):
        raise SystemExit(f"error: runs use different `arcs` assumptions: {assumptions}")
    e_min, anti_deg = float(ref["min_ellipticity"]), float(ref["anti_deg"])
    if float(ref["sigma_eps"]) != SIGMAS[0]:
        raise SystemExit(f"error: runs used σ_ε {ref['sigma_eps']}, the scan starts at {SIGMAS[0]}")

    def stats_at(sigma: float, fields: tuple[str, ...] = tuple(FIELDS)) -> dict[str, dict]:
        return {f: field_stats(runs[f][0], sigma, e_min, anti_deg) for f in fields}

    per_field = stats_at(SIGMAS[0])
    checks = {
        f: stored_check(per_field[f], meta["anti_expected_shape_noise_background"])
        for f, (_, meta) in runs.items()
    }
    s1063_ref = field_stats(s1063_rows, SIGMAS[0], e_min, anti_deg)
    checks["S1063"] = stored_check(s1063_ref, s1063_meta["anti_expected_shape_noise_background"])
    if not all(c["matches"] for c in checks.values()):
        raise SystemExit(f"error: recomputed summaries differ from the stored ones: {checks}")
    gsplit = {
        f: {
            "|g_min|<=1": int(np.sum(anti_mask(t) & ~inside(t))),
            "|g_min|>1": int(np.sum(anti_mask(t) & inside(t))),
        }
        for f, (t, _) in runs.items()
    }
    res = {v: pooled(per_field, fs) for v, fs in VARIANTS.items()}
    scan = []
    for s in SIGMAS:
        pf = per_field if s == SIGMAS[0] else stats_at(s)
        pv = {v: pooled(pf, fs) for v, fs in VARIANTS.items()}
        s1063 = field_stats(s1063_rows, s, e_min, anti_deg)
        scan.append(
            {
                "sigma_eps": s,
                **{f"expected_{v}": r["expected"] for v, r in pv.items()},
                **{f"p_{v}": r["p_ge_observed"] for v, r in pv.items()},
                "s1063": s1063,
                **{
                    f"p_conditional_{v}": s1063_vs_pooled(s1063, r)["p_conditional"]
                    for v, r in pv.items()
                },
            }
        )
    crossing = {
        v: sigma_crossing(
            lambda x, fs=fs: pooled(stats_at(x, fs), fs)["p_ge_observed"], SIGMAS[0], SIGMAS[-1]
        )
        for v, fs in VARIANTS.items()
    }
    summary = {
        "plan": "results/arcs/pooled-tally/PLAN.md",
        "label": "exploratory (post hoc)",
        "fields": FIELDS,
        "assumptions": ref,
        "stored_check": checks,
        f"per_field_sigma_{SIGMAS[0]}": per_field,
        "anti_by_gmin": gsplit,
        "anti_by_gmin_pooled": {
            v: {k: sum(gsplit[f][k] for f in fs) for k in ("|g_min|<=1", "|g_min|>1")}
            for v, fs in VARIANTS.items()
        },
        "pooled": res,
        "outcome": outcome(res["A"]["p_ge_observed"], res["B"]["p_ge_observed"]),
        "sigma_scan": scan,
        "sigma_at_alpha": crossing,
        "s1063_at_sigma_at_alpha": {
            v: None if x is None else field_stats(s1063_rows, x, e_min, anti_deg)
            for v, x in crossing.items()
        },
        "s1063_vs_pooled_posthoc": {v: s1063_vs_pooled(s1063_ref, res[v]) for v in VARIANTS},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    keys = ("pooled", "outcome", "sigma_at_alpha", "s1063_at_sigma_at_alpha", "anti_by_gmin_pooled")
    print(json.dumps({k: summary[k] for k in keys}, indent=1))
    for r in scan:
        print(
            f"sigma {r['sigma_eps']:.2f}: A {r['expected_A']:.2f} P {r['p_A']:.4f} | "
            f"B {r['expected_B']:.2f} P {r['p_B']:.4f} | S1063 {r['s1063']['anti']} vs "
            f"{r['s1063']['expected']:.2f} P {r['s1063']['p_ge_observed']:.4f} | "
            f"cond A {r['p_conditional_A']:.3f} B {r['p_conditional_B']:.3f}"
        )


if __name__ == "__main__":
    main()
