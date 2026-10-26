"""SN Encore (MACS J0138.0-2155): standard-candle magnifications vs the seven blind lens models.

Frozen plan: ``results/lens_consistency/encore_mu/PLAN.md``. Inputs are published table values only
(no images):

- Suyu et al. (arXiv:2509.12319, e-print of 2026-03-03) Table "tab:snencore": blind model
  predictions for images 1a and 1b (median and 68 % interval, flat LCDM, H0 70, Om 0.3). The seven
  blind models are the ones Pierel et al. Table "tab:h0" combine; the GLEE-baseline and post-blind
  columns are not used.
- Pierel et al. (arXiv:2509.12301) Table "tab:final_delay_mag" (final mu and delay, systematics
  included), Sect. "Magnifications" (statistical-only mu), Table "tab:systematics" (simulated micro-
  and millilensing) and Table "tab:fitting" (directly fitted mu_1b / mu_1a, statistical only).

    python scripts/encore_mu_check.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

OUT = Path("results/lens_consistency/encore_mu")
IMAGES = ("1a", "1b")
Z_FLAG = 3.0  # ASSUMPTION: plan rule
H0_RANGE = (64.0, 80.0)  # ASSUMPTION: plan rule
H0_FID = 70.0  # fiducial H0 of every model (Suyu et al. table footnote)
SAME_SENSE_MIN = 5  # ASSUMPTION: plan rule, models of 7 with a positive offset

# Suyu et al. Table tab:snencore (model_prediction): (median, +err, -err); mu signed as published
MODELS = {
    "glafic": {"1a": (-27.5, 2.1, 2.3), "1b": (48.4, 7.6, 5.5), "dt": (-32.4, 2.2, 2.5)},
    "GLEE": {"1a": (-28.0, 2.5, 3.7), "1b": (40.0, 7.6, 6.3), "dt": (-37.1, 2.7, 2.6)},
    "Lenstool I": {"1a": (-19.9, 1.6, 1.7), "1b": (21.7, 5.6, 3.6), "dt": (-75.0, 55.0, 54.0)},
    "Lenstool II": {"1a": (-26.1, 3.3, 4.4), "1b": (31.8, 6.1, 5.9), "dt": (-35.6, 3.1, 4.3)},
    "MrMARTIAN": {"1a": (-31.9, 1.4, 7.0), "1b": (27.5, 14.7, 2.3), "dt": (-40.6, 9.8, 2.9)},
    "WSLAP+": {"1a": (-53.0, 14.0, 14.0), "1b": (31.8, 8.8, 8.8), "dt": (-112.0, 32.0, 32.0)},
    "Zitrin-analytic": {"1a": (-18.3, 4.4, 4.4), "1b": (18.3, 3.9, 4.1), "dt": (-40.2, 7.5, 11.1)},
}
# Pierel et al. (observed / derived from the SN Ia standard candle; |mu|)
OBS = {"1a": (21.8, 9.3, 8.6), "1b": (32.4, 11.1, 11.0), "dt": (-39.8, 3.9, 3.3)}  # final table
OBS_STAT = {"1a": (21.8, 1.2, 1.3), "1b": (32.4, 2.4, 2.2)}  # Sect. Magnifications
SYS_MICRO = {"1a": (9.3, 4.2), "1b": (7.1, 5.9)}  # tab:systematics (+, -), simulated
SYS_MILLI = {"1a": (6.2, 6.2), "1b": (9.2, 9.2)}  # tab:systematics (+, -), simulated, f_sub 5 %
FITTED_RATIO = (1.5, 0.5, 0.5)  # tab:fitting, SNTD phot_starred, statistical only


def ln_sigma(value: tuple[float, float, float]) -> float:
    """Half the 16-84 % width in ln space of |value| (+ and - errors swap for a negative median)."""
    med, hi, lo = value
    if med < 0:
        med, hi, lo = -med, lo, hi
    return 0.5 * (math.log(med + hi) - math.log(med - lo))


def mag(x: float) -> float:
    return 2.5 * math.log10(x)


def sigma_obs_ln(image: str, final: bool = True) -> float:
    """Observed ln-mu error. ``final``: the final-table interval, which already contains the
    published micro- and millilensing terms (so none is added again). Otherwise (post hoc check):
    statistical interval with both simulated terms added in quadrature."""
    if final:
        return ln_sigma(OBS[image])
    med = OBS_STAT[image][0]
    extra = [0.5 * sum(SYS_MICRO[image]) / med, 0.5 * sum(SYS_MILLI[image]) / med]
    return math.sqrt(ln_sigma(OBS_STAT[image]) ** 2 + sum(e * e for e in extra))


def image_z(model: dict, image: str, final: bool = True) -> float:
    """(ln |mu_obs| - ln |mu_model|) / combined ln sigma."""
    s = math.hypot(sigma_obs_ln(image, final), ln_sigma(model[image]))
    return (math.log(OBS[image][0]) - math.log(abs(model[image][0]))) / s


def ratio_z(model: dict, final: bool = True) -> float:
    """Same for |mu_1b / mu_1a|. The paper does not separate a term shared by both images, so all
    errors add in quadrature (the plan's conservative branch)."""
    s2 = sum(sigma_obs_ln(k, final) ** 2 + ln_sigma(model[k]) ** 2 for k in IMAGES)
    obs = math.log(OBS["1b"][0] / OBS["1a"][0])
    pred = math.log(abs(model["1b"][0] / model["1a"][0]))
    return (obs - pred) / math.sqrt(s2)


def fitted_ratio_z(model: dict) -> float:
    """Post hoc check: the directly fitted ratio (statistical only) with the model ratio error."""
    s2 = ln_sigma(FITTED_RATIO) ** 2 + sum(ln_sigma(model[k]) ** 2 for k in IMAGES)
    pred = math.log(abs(model["1b"][0] / model["1a"][0]))
    return (math.log(FITTED_RATIO[0]) - pred) / math.sqrt(s2)


def offset_mag(model: dict) -> float:
    """Mean over images of 2.5 log10(|mu_model| / mu_obs) (positive: the model over-predicts)."""
    return sum(mag(abs(model[i][0]) / OBS[i][0]) for i in IMAGES) / len(IMAGES)


def h0_mst(model: dict) -> dict:
    """H0 at which one mass sheet (mu / lambda^2, lambda dt) fits both mu and the delay."""
    mean_ln = sum(math.log(abs(model[i][0]) / OBS[i][0]) for i in IMAGES) / len(IMAGES)
    lam_mu = math.exp(0.5 * mean_ln)
    lam_t = OBS["dt"][0] / model["dt"][0]
    h0 = H0_FID * lam_mu / lam_t
    # post hoc 1 sigma: observed mu errors fully correlated across the two images (conservative),
    # model terms independent, delay terms independent
    n = len(IMAGES)
    s_obs = sum(sigma_obs_ln(i) for i in IMAGES) / n
    s_mod2 = sum(ln_sigma(model[i]) ** 2 for i in IMAGES) / n**2
    s_mu = math.sqrt(s_obs**2 + s_mod2)
    s_t = math.hypot(ln_sigma(OBS["dt"]), ln_sigma(model["dt"]))
    s_ln = math.hypot(0.5 * s_mu, s_t)
    return {
        "lambda_mu": lam_mu,
        "lambda_t": lam_t,
        "h0_mst": h0,
        "in_range": H0_RANGE[0] <= h0 <= H0_RANGE[1],
        "h0_mst_1sigma": [h0 * math.exp(-s_ln), h0 * math.exp(s_ln)],
    }


def outcome(rows: dict) -> tuple[list[str], list[str], str]:
    """Plan rule: an image or the ratio is flagged only when |z| > 3 in every model."""
    flags_image = [i for i in IMAGES if all(abs(r["image_z"][i]) > Z_FLAG for r in rows.values())]
    flags_ratio = ["1b/1a"] if all(abs(r["ratio_z"]) > Z_FLAG for r in rows.values()) else []
    return flags_image, flags_ratio, "F" if flags_image or flags_ratio else "N"


def same_sense(offsets: list[float]) -> bool:
    """Plan rule: pooled offset > 0 and at least SAME_SENSE_MIN models with a positive offset."""
    return sum(offsets) / len(offsets) > 0 and sum(o > 0 for o in offsets) >= SAME_SENSE_MIN


def run() -> tuple[dict, dict]:
    rows, posthoc_rows = {}, {}
    for name, m in MODELS.items():
        mst = h0_mst(m)
        rows[name] = {
            "image_z": {i: image_z(m, i) for i in IMAGES},
            "ratio_z": ratio_z(m),
            "model_ratio": abs(m["1b"][0] / m["1a"][0]),
            "offset_mag": offset_mag(m),
            **{k: mst[k] for k in ("lambda_mu", "lambda_t", "h0_mst", "in_range")},
        }
        posthoc_rows[name] = {
            "h0_mst_1sigma": mst["h0_mst_1sigma"],
            "image_z_stat_plus_sys_quadrature": {i: image_z(m, i, final=False) for i in IMAGES},
            "ratio_z_stat_plus_sys_quadrature": ratio_z(m, final=False),
            "ratio_z_fitted_stat_only": fitted_ratio_z(m),
        }
    flags_image, flags_ratio, out = outcome(rows)
    offsets = [r["offset_mag"] for r in rows.values()]
    result = {
        "rows": rows,
        "flags_image": flags_image,
        "flags_ratio": flags_ratio,
        "outcome": out,
        "observed_ratio": OBS["1b"][0] / OBS["1a"][0],
        "pooled_offset_mag": sum(offsets) / len(offsets),
        "offset_range_mag": [min(offsets), max(offsets)],
        "n_models_offset_positive": sum(o > 0 for o in offsets),
        "same_sense_as_h0pe": same_sense(offsets),
        "n_models_h0_mst_in_range": sum(r["in_range"] for r in rows.values()),
        # derived: smallest |mu_obs / mu_model| factor a flag needs with a zero model error
        "flag_reach_factor": {i: math.exp(Z_FLAG * sigma_obs_ln(i)) for i in IMAGES},
        "inputs": {
            "suyu_models": MODELS,
            "pierel_final": OBS,
            "pierel_stat": OBS_STAT,
            "pierel_sys_micro": SYS_MICRO,
            "pierel_sys_milli": SYS_MILLI,
            "pierel_fitted_ratio": FITTED_RATIO,
            "sigma_obs_ln": {i: sigma_obs_ln(i) for i in IMAGES},
        },
        "assumptions": {
            "z_flag": Z_FLAG,
            "h0_range": list(H0_RANGE),
            "same_sense_min": SAME_SENSE_MIN,
        },
    }
    ranges = [r["h0_mst_1sigma"] for r in posthoc_rows.values()]
    posthoc = {
        "rows": posthoc_rows,
        "n_models_1sigma_reaching_h0_range": sum(
            lo <= H0_RANGE[1] and hi >= H0_RANGE[0] for lo, hi in ranges
        ),
    }
    return result, posthoc


def main() -> None:
    result, posthoc = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    (OUT / "posthoc.json").write_text(json.dumps(posthoc, indent=2) + "\n")
    for name, r in result["rows"].items():
        zi = " ".join(f"{k} {v:+.2f}" for k, v in r["image_z"].items())
        p = posthoc["rows"][name]
        print(
            f"{name:16s} z[{zi}] ratio {r['model_ratio']:.2f} z {r['ratio_z']:+.2f} "
            f"(fit {p['ratio_z_fitted_stat_only']:+.2f}) offset {r['offset_mag']:+.3f} mag "
            f"H0_MST {r['h0_mst']:.1f} ({p['h0_mst_1sigma'][0]:.1f}-{p['h0_mst_1sigma'][1]:.1f})"
        )
    print(
        f"outcome {result['outcome']} pooled {result['pooled_offset_mag']:+.3f} "
        f"range {result['offset_range_mag']} same_sense {result['same_sense_as_h0pe']} "
        f"positive {result['n_models_offset_positive']}/7"
    )


if __name__ == "__main__":
    main()
