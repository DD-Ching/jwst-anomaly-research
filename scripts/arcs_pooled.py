"""Pooled `arcs` anti count over the Abell S1063 replication fields (exploratory, post hoc).

Plan: results/arcs/pooled-tally/PLAN.md. Reads each field's stored deciding-run ``arcs.ecsv``,
recomputes the shape-noise expectation with ``lens_consistency.anti_expected`` (checked against the
stored summary), pools the counts and scans σ_ε. No network, no new data.

    python scripts/arcs_pooled.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.stats import binom, poisson

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
    t = Table.read(path)
    t = t[t["strong_shear"]]
    return t[t["population"] == "background"]


def pooled(stats: dict[str, dict], fields: tuple[str, ...]) -> dict:
    n = sum(stats[f]["anti"] for f in fields)
    lam = sum(stats[f]["expected"] for f in fields)
    return {"anti": n, "expected": round(lam, 3), "p_ge_observed": float(poisson.sf(n - 1, lam))}


def outcome(p_a: float, p_b: float, alpha: float = ALPHA) -> str:
    """Frozen plan rule: both variants below alpha, both at or above it, or mixed."""
    if p_a < alpha and p_b < alpha:
        return "pooled excess"
    if p_a >= alpha and p_b >= alpha:
        return "no pooled excess"
    return "mixed"


def s1063_vs_pooled(s1063: dict, pool: dict) -> dict:
    """Post hoc: is S1063's anti count in line with the pooled fields' anti rate per expected count?

    ``scaled``: S1063's λ times the pooled observed/expected ratio, P(>= observed).
    ``conditional``: given the combined count, P(S1063 holds >= its share) with binomial
    p = λ_S1063 / (λ_S1063 + λ_pool); this cancels a common scale error in σ_ε (first order).
    """
    ratio = pool["anti"] / pool["expected"]
    lam = ratio * s1063["expected"]
    total = s1063["anti"] + pool["anti"]
    share = s1063["expected"] / (s1063["expected"] + pool["expected"])
    return {
        "pool_ratio": round(ratio, 3),
        "scaled_expected": round(lam, 3),
        "p_scaled": float(poisson.sf(s1063["anti"] - 1, lam)),
        "p_conditional": float(binom.sf(s1063["anti"] - 1, total, share)),
    }


def main() -> None:
    tables = {f: background(ROOT / p) for f, p in FIELDS.items()}
    args = (0.5, 60.0)  # D-024 defaults: min ellipticity, anti_deg
    per_field = {f: lc.anti_expected(t, 0.26, *args) for f, t in tables.items()}
    gsplit = {
        f: {
            "g_min<1": int(np.sum((t["orientation_class"] == "anti") & (np.abs(t["g_min"]) < 1))),
            "g_min>=1": int(np.sum((t["orientation_class"] == "anti") & (np.abs(t["g_min"]) >= 1))),
        }
        for f, t in tables.items()
    }
    res = {v: pooled(per_field, fs) for v, fs in VARIANTS.items()}
    s1063 = background(ROOT / S1063)
    scan = []
    for s in SIGMAS:
        pf = {f: lc.anti_expected(t, s, *args) for f, t in tables.items()}
        scan.append(
            {
                "sigma_eps": s,
                **{f"p_{v}": pooled(pf, fs)["p_ge_observed"] for v, fs in VARIANTS.items()},
                **{f"expected_{v}": pooled(pf, fs)["expected"] for v, fs in VARIANTS.items()},
                "s1063": lc.anti_expected(s1063, s, *args),
            }
        )
    summary = {
        "plan": "results/arcs/pooled-tally/PLAN.md",
        "label": "exploratory (post hoc)",
        "fields": FIELDS,
        "per_field_sigma_0.26": per_field,
        "anti_by_gmin": gsplit,
        "pooled": res,
        "outcome": outcome(res["A"]["p_ge_observed"], res["B"]["p_ge_observed"]),
        "sigma_scan": scan,
        "s1063_vs_pooled_posthoc": {v: s1063_vs_pooled(scan[0]["s1063"], res[v]) for v in VARIANTS},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {k: summary[k] for k in ("pooled", "outcome", "s1063_vs_pooled_posthoc")}, indent=1
        )
    )
    for r in scan:
        print(
            f"sigma {r['sigma_eps']:.2f}: A {r['expected_A']:.2f} P {r['p_A']:.4f} | "
            f"B {r['expected_B']:.2f} P {r['p_B']:.4f} | S1063 {r['s1063']['anti']} vs "
            f"{r['s1063']['expected']:.2f} P {r['s1063']['p_ge_observed']:.4f}"
        )


if __name__ == "__main__":
    main()
