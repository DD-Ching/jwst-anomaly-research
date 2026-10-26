"""SN Refsdal (MACS J1149.5+2223): measured SX delay and SX / S1 μ ratio vs blind models.

Frozen plan: ``results/lens_consistency/refsdal_sx/PLAN.md``. Published table values only:

- Treu et al. 2016 (ApJ 817, 60; arXiv:1510.05750 source) Table "tab:predictions": blind model
  predictions of Delta t_X1 (days, each model's fiducial cosmology) and mu(X) / mu(1).
- Kelly et al. 2023 (Science 380, eabh1322; arXiv:2305.06367 source) Table "tab:measurements":
  SX - S1 delay and mu(X) / mu(1), maximum-likelihood value and 16 / 50 / 84 percentiles.
- Post hoc only (not in the outcome): Kelly et al. Tables "tab:mod_time_delays" /
  "tab:mod_mag_ratios", models revised or first computed after the reappearance (asterisked there),
  and the blind Jauzac15.1 prediction, which is not in the Treu et al. table.

    python scripts/refsdal_sx_check.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

OUT = Path("results/lens_consistency/refsdal_sx")
Z_FLAG = 3.0  # ASSUMPTION: plan rule

# Treu et al. 2016 tab:predictions (model_prediction): (median, +err, -err)
BLIND = {
    "Die-a": {"dt": (262.0, 55.0, 55.0), "mu": (0.31, 0.10, 0.10)},
    "Gri-g": {"dt": (361.0, 19.0, 27.0), "mu": (0.36, 0.11, 0.09)},
    "Ogu-g": {"dt": (311.0, 24.0, 24.0), "mu": (0.27, 0.05, 0.05)},
    "Ogu-a": {"dt": (336.0, 21.0, 21.0), "mu": (0.27, 0.03, 0.03)},
    "Sha-g": {"dt": (277.0, 11.0, 21.0), "mu": (0.25, 0.05, 0.02)},
    "Sha-a": {"dt": (233.0, 46.0, 13.0), "mu": (0.19, 0.01, 0.04)},
    "Zit-g": {"dt": (224.0, 262.0, 262.0), "mu": (0.31, 0.05, 0.05)},
}
# Kelly et al. 2023 tab:mod_time_delays / tab:mod_mag_ratios rows not identical to BLIND (post hoc).
# Their Zitrin-c* mu(X) lower error is printed without a sign (0.023) and is read as -0.023.
POST_HOC = {
    "Jauzac15.1 (blind, not in Treu table)": {"dt": (449.0, 45.0, 45.0), "mu": (0.19, 0.02, 0.02)},
    "Jauzac15.2*": {"dt": (361.0, 42.0, 42.0), "mu": (0.19, 0.02, 0.02)},
    "Oguri-g*": {"dt": (330.7, 21.9, 17.4), "mu": (0.29, 0.03, 0.02)},
    "Oguri-a*": {"dt": (353.8, 20.7, 17.3), "mu": (0.28, 0.04, 0.02)},
    "Sharon-a*": {"dt": (298.0, 58.0, 17.0), "mu": (0.19, 0.01, 0.04)},
    "Sharon-g*": {"dt": (331.0, 30.0, 35.0), "mu": (0.25, 0.05, 0.02)},
    "Zitrin-c*": {"dt": (36.74, 122.0, 11.0), "mu": (0.2921, 0.0557, 0.023)},
    "Keeton*": {"dt": (340.01, 20.5, 22.2), "mu": (0.28, 0.07, 0.05)},
    "Chen2020*": {"dt": (332.35, 9.32, 9.32), "mu": None},
    "Zitrin2020-ltm*": {"dt": (267.4, 37.9, 4.3), "mu": None},
    "Zitrin2020-p*": {"dt": (325.8, 19.8, 14.2), "mu": None},
}
# Kelly et al. 2023 tab:measurements (observed): (ML value, p16, p50, p84)
OBS_DT = (376.02, 370.50, 376.03, 381.65)
OBS_MU = (0.299, 0.272, 0.305, 0.354)


def half_width(value: tuple[float, float, float]) -> float:
    return 0.5 * (value[1] + value[2])


def ln_half_width(value: tuple[float, float, float]) -> float:
    med, hi, lo = value
    return 0.5 * (math.log(med + hi) - math.log(med - lo))


def obs_sigma() -> tuple[float, float]:
    """Delay sigma (days) and ratio sigma_ln from the published 16-84 percentiles."""
    return 0.5 * (OBS_DT[3] - OBS_DT[1]), 0.5 * (math.log(OBS_MU[3]) - math.log(OBS_MU[1]))


def delay_z(model: tuple[float, float, float]) -> float:
    s_obs, _ = obs_sigma()
    return (OBS_DT[0] - model[0]) / math.hypot(s_obs, half_width(model))


def ratio_z(model: tuple[float, float, float]) -> float:
    _, s_obs = obs_sigma()
    return (math.log(OBS_MU[0]) - math.log(model[0])) / math.hypot(s_obs, ln_half_width(model))


def row(m: dict) -> dict:
    zd = delay_z(m["dt"])
    zr = ratio_z(m["mu"]) if m["mu"] else None
    flagged = abs(zd) > Z_FLAG or (zr is not None and abs(zr) > Z_FLAG)
    return {"z_dt": zd, "z_mu": zr, "flagged": flagged, "both_given": zr is not None}


def outcome(rows: dict) -> str:
    """F only if every model with both quantities is flagged; U if none gives both."""
    both = [r for r in rows.values() if r["both_given"]]
    if not both:
        return "U"
    return "F" if all(r["flagged"] for r in both) else "N"


def run() -> dict:
    blind = {k: row(m) for k, m in BLIND.items()}
    post = {k: row(m) for k, m in POST_HOC.items()}
    s_dt, s_ln = obs_sigma()
    return {
        "outcome": outcome(blind),
        "blind": blind,
        "n_blind_flagged": sum(r["flagged"] for r in blind.values()),
        "post_hoc": post,
        "observed": {"dt": OBS_DT, "mu": OBS_MU, "sigma_dt": s_dt, "sigma_ln_mu": s_ln},
        "reach": {"dt_days_obs_only": Z_FLAG * s_dt, "mu_factor_obs_only": math.exp(Z_FLAG * s_ln)},
        "inputs": {"treu2016": BLIND, "kelly2023_post_hoc": POST_HOC},
        "assumptions": {"z_flag": Z_FLAG},
    }


def main() -> None:
    result = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    for group in ("blind", "post_hoc"):
        print(group)
        for name, r in result[group].items():
            zm = "  —  " if r["z_mu"] is None else f"{r['z_mu']:+.2f}"
            print(f"  {name:40s} z_dt {r['z_dt']:+.2f} z_mu {zm} {'FLAG' if r['flagged'] else ''}")
    print("outcome", result["outcome"], result["n_blind_flagged"], result["reach"])


if __name__ == "__main__":
    main()
