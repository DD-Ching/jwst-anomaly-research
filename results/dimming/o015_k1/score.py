"""PLAN.md scoring: G0 (finite k = 0 fluxes unchanged under k = 1), `on_sky` / n of the o015 pairs by
exposure class, k = 0 beside k = 1, and the frozen rules R1-R3.
Usage (repo root): score.py  -> results/dimming/o015_k1/score.json"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "results/dimming/nan_k1")
from score import EPS_GT4, EPS_K0, EPS_K1, g0, wilson  # noqa: E402

D = Path("results/dimming/o015_k1")
N_MIN, PERSIST_MAX = 10, 0.02  # PLAN.md


def classes(v: Table) -> dict:
    n_exp = np.asarray(v["n_exposures"], int)
    on = np.asarray(v["verdict"]) == "on_sky"
    out = {}
    for name, m in (("le4", n_exp <= 4), ("gt4", n_exp > 4)):
        k, n = int((on & m).sum()), int(m.sum())
        out[name] = {"on_sky": k, "n": n, "wilson": wilson(k, n)}
    k = int((np.asarray(v["verdict"]) == "persistence").sum())
    out["persistence"] = {"k": k, "n": len(v), "wilson": wilson(k, len(v))}
    out["verdicts"] = {
        str(a): int(b) for a, b in zip(*np.unique(v["verdict"], return_counts=True), strict=True)
    }
    return out


def decide(g0_pass: bool, k1: dict) -> str:
    """Frozen R1-R3 (PLAN.md); every condition of R1 must hold."""
    le4 = k1["le4"]
    if le4["n"] < N_MIN:
        return "R3"
    lo, hi = le4["wilson"]
    contains_or_above = hi >= EPS_K1
    r1 = (
        g0_pass
        and contains_or_above
        and lo > EPS_K0
        and k1["persistence"]["wilson"][0] <= PERSIST_MAX
    )
    return "R1" if r1 else "R2"


def main() -> None:
    ve0, ve1 = Table.read(D / "k0/verdicts.ecsv"), Table.read(D / "k1/verdicts.ecsv")
    ex0, ex1 = Table.read(D / "k0/exposures.ecsv"), Table.read(D / "k1/exposures.ecsv")
    gate = g0(ex0, ve0, ex1, ve1)
    c0, c1 = classes(ve0), classes(ve1)
    gt4 = c0["gt4"]
    res = {
        "G0": gate,
        "k0": c0,
        "k1": c1,
        "decision": decide(gate["pass"], c1),
        "gt4_holds": None if gt4["n"] == 0 else gt4["wilson"][1] >= EPS_GT4,
        "assumptions": {"eps_k1": EPS_K1, "eps_k0": EPS_K0, "eps_gt4": EPS_GT4, "n_min": N_MIN},
        "provenance": "derived",
    }
    (D / "score.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("decision", "gt4_holds")}), c0, c1, gate["pass"])


if __name__ == "__main__":
    main()
