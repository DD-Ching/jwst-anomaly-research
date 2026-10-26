"""SN Requiem (MACS J0138.0-2155): measured time delays vs the seven blind lens models.

Frozen plan: ``results/lens_consistency/requiem_dt/PLAN.md``. Published table values only:

- Suyu et al. (arXiv:2509.12319, e-print of 2026-03-03) Table "tab:snrequiem": blind model delays
  of 2b and 2c relative to 2a (median and 68 % interval, flat LCDM, H0 70). The seven blind models
  are the ones of the SN Encore check (``encore_mu_check.MODELS``).
- Rodney et al. 2021 (arXiv:2106.08935v2) Table 1 "tab:time_delays": Delta t_obs of SN2 and SN3
  relative to SN1 (colour curve combined with their own lens-model magnifications). By position
  (Table 1 coordinates vs Suyu et al. offsets from the BCG) SN1 = 2b, SN2 = 2a, SN3 = 2c.

    python scripts/requiem_dt_check.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

OUT = Path("results/lens_consistency/requiem_dt")
Z_FLAG = 3.0  # ASSUMPTION: plan rule
N_DRAW = 200_000
SEED = 20261026

# Suyu et al. Table tab:snrequiem (model_prediction): (median, +err, -err), days
MODELS = {
    "glafic": {"ba": (-49.8, 4.1, 4.4), "ca": (-80.4, 5.4, 6.0)},
    "GLEE": {"ba": (-47.9, 3.5, 3.6), "ca": (-78.4, 7.2, 8.7)},
    "Lenstool I": {"ba": (-74.0, 52.0, 57.0), "ca": (-83.0, 51.0, 84.0)},
    "Lenstool II": {"ba": (-52.2, 3.8, 5.5), "ca": (-87.7, 7.9, 8.2)},
    "MrMARTIAN": {"ba": (-29.8, 2.4, 2.2), "ca": (-70.6, 36.4, 6.2)},
    "WSLAP+": {"ba": (-61.0, 14.0, 14.0), "ca": (-95.0, 43.0, 43.0)},
    "Zitrin-analytic": {"ba": (-35.0, 11.2, 14.2), "ca": (-58.3, 25.2, 33.9)},
}
# Rodney et al. Table 1 (observed with their lens-model magnification prior), relative to SN1 (= 2b)
RODNEY = {"SN2": (114.0, 28.0, 31.0), "SN3": (-17.0, 19.0, 16.0)}


def ln_sigma(value: tuple[float, float, float]) -> float:
    """Half the 16-84 % width in ln space of |value| (+ and - errors swap for a negative median)."""
    med, hi, lo = value
    if med < 0:
        med, hi, lo = -med, lo, hi
    return 0.5 * (math.log(med + hi) - math.log(med - lo))


def split_normal(rng: np.random.Generator, value: tuple[float, float, float], n: int) -> np.ndarray:
    """Draws with the published + / - half-widths on either side of the median."""
    med, hi, lo = value
    z = rng.standard_normal(n)
    return med + np.where(z > 0, hi, lo) * z


def observed_draws(seed: int = SEED, n: int = N_DRAW) -> dict[str, np.ndarray]:
    """Delays relative to 2a: dt_2b,2a = -dt(SN2 - SN1); dt_2c,2a = dt(SN3 - SN1) - dt(SN2 - SN1).
    The two published delays are drawn independently (no covariance is published)."""
    rng = np.random.default_rng(seed)
    d2 = split_normal(rng, RODNEY["SN2"], n)
    d3 = split_normal(rng, RODNEY["SN3"], n)
    return {"ba": -d2, "ca": d3 - d2}


def summarize(x: np.ndarray) -> tuple[float, float, float]:
    lo, med, hi = np.percentile(x, [16, 50, 84])
    return float(med), float(hi - med), float(med - lo)


def ratio_obs(draws: dict[str, np.ndarray]) -> dict:
    """R = dt_2c,2a / dt_2b,2a; undefined in ln when either delay's 68 % interval crosses zero."""
    crosses = any(s[0] - s[2] < 0 < s[0] + s[1] for s in map(summarize, draws.values()))
    ok = draws["ba"] != 0
    r = draws["ca"][ok] / draws["ba"][ok]
    pos = r[r > 0]
    lo, med, hi = np.percentile(pos, [16, 50, 84])
    return {
        "defined": not crosses,
        "median": float(med),
        "sigma_ln": float(0.5 * (math.log(hi) - math.log(lo))),
        "fraction_positive": float(pos.size / r.size),
    }


def ratio_z(model: dict, robs: dict) -> float:
    pred = model["ca"][0] / model["ba"][0]
    s_m = math.hypot(ln_sigma(model["ca"]), ln_sigma(model["ba"]))
    return (math.log(robs["median"]) - math.log(pred)) / math.hypot(robs["sigma_ln"], s_m)


def delay_z(model: dict, key: str, obs: tuple[float, float, float]) -> float:
    """Linear z on one delay at H0 70, symmetrised half-widths in quadrature."""
    s = math.hypot(0.5 * (obs[1] + obs[2]), 0.5 * (model[key][1] + model[key][2]))
    return (obs[0] - model[key][0]) / s


def outcome(rows: dict, ratio_defined: bool) -> tuple[list[str], str]:
    """Flag a quantity only when |z| > 3 in every model; U: ratio undefined and no delay flag."""
    keys = ["ratio", "ba", "ca"] if ratio_defined else ["ba", "ca"]
    flags = [k for k in keys if all(abs(r[k]) > Z_FLAG for r in rows.values())]
    if flags:
        return flags, "F"
    return flags, "N" if ratio_defined else "U"


def run() -> dict:
    draws = observed_draws()
    obs = {k: summarize(v) for k, v in draws.items()}
    robs = ratio_obs(draws)
    rows = {}
    for name, m in MODELS.items():
        rows[name] = {
            "model_ratio": m["ca"][0] / m["ba"][0],
            "ratio": ratio_z(m, robs),
            "ba": delay_z(m, "ba", obs["ba"]),
            "ca": delay_z(m, "ca", obs["ca"]),
        }
    flags, out = outcome(rows, robs["defined"])
    return {
        "rows": rows,
        "observed_delays": obs,
        "observed_ratio": robs,
        "flag_reach_factor_ratio": math.exp(Z_FLAG * robs["sigma_ln"]),
        "flags": flags,
        "outcome": out,
        "inputs": {"suyu_models": MODELS, "rodney_table1": RODNEY, "seed": SEED, "n_draw": N_DRAW},
        "assumptions": {"z_flag": Z_FLAG},
    }


def main() -> None:
    result = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print("obs", result["observed_delays"], result["observed_ratio"])
    for name, r in result["rows"].items():
        print(
            f"{name:16s} R_m {r['model_ratio']:.2f} z_R {r['ratio']:+.2f} "
            f"z_ba {r['ba']:+.2f} z_ca {r['ca']:+.2f}"
        )
    print("outcome", result["outcome"], result["flags"], result["flag_reach_factor_ratio"])


if __name__ == "__main__":
    main()
