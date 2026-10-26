"""SN H0pe (PLCK G165.7+67.0): standard-candle magnifications vs the seven blind lens models.

Frozen plan: ``results/lens_consistency/h0pe_mu/PLAN.md`` (rules 1-2 and addendum A). Inputs are
published table values only (no images):

- Pascale et al. 2025 (arXiv:2403.18902v2) Table 1: unblinded model predictions and the
  photometric measurements of Pierel et al. 2024 (arXiv:2403.18954v2), 68 % intervals.
- Agrawal et al. (arXiv:2510.07637v2) Sect. 4: per-model Delta-mu range and pooled value
  (comparison only).

    python scripts/h0pe_mu_check.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

OUT = Path("results/lens_consistency/h0pe_mu")
IMAGES = ("a", "b", "c")
SIGMA_ML_MAG = (
    0.15  # ASSUMPTION: plan placeholder, replaced by Pierel et al. Table 4 (post hoc check)
)
Z_FLAG = 3.0  # ASSUMPTION: plan rule
REPRODUCE_MAG = 0.3  # ASSUMPTION: addendum A rule 1
H0_RANGE = (64.0, 80.0)  # ASSUMPTION: addendum A rule 2
H0_FID = 70.0  # fiducial H0 of every model (Pascale et al. Sect. 2)

# (median, +err, -err); observed for "phot", model_prediction otherwise. Time delays in days.
TABLE1 = {
    "1 GLAFIC": {
        "dt_ab": (-105.18, 5.16, 7.89),
        "dt_cb": (-50.71, 3.37, 5.09),
        "a": (8.02, 0.64, 0.57),
        "b": (12.23, 1.30, 1.51),
        "c": (9.32, 0.80, 0.87),
    },
    "2 Zitrin-analytic": {
        "dt_ab": (-105.52, 5.09, 6.32),
        "dt_cb": (-41.06, 14.27, 13.19),
        "a": (11.25, 1.05, 0.90),
        "b": (16.03, 1.77, 1.81),
        "c": (14.48, 1.91, 1.65),
    },
    "3 LENSTOOL": {
        "dt_ab": (-102.68, 6.47, 6.80),
        "dt_cb": (-54.09, 3.18, 3.76),
        "a": (8.75, 1.01, 0.60),
        "b": (11.32, 1.29, 1.14),
        "c": (10.47, 1.14, 0.74),
    },
    "4 MARS": {
        "dt_ab": (-136.27, 20.23, 20.32),
        "dt_cb": (-63.72, 29.24, 28.46),
        "a": (6.82, 0.50, 0.44),
        "b": (9.17, 0.80, 0.72),
        "c": (7.55, 0.86, 0.70),
    },
    "5 Chen2020": {
        "dt_ab": (-112.25, 6.43, 6.57),
        "dt_cb": (-53.35, 2.71, 2.99),
        "a": (6.52, 0.24, 0.22),
        "b": (10.35, 0.46, 0.41),
        "c": (6.68, 0.24, 0.22),
    },
    "6 WSLAP+": {
        "dt_ab": (-273.35, 95.30, 95.30),
        "dt_cb": (342.75, 92.54, 92.54),
        "a": (16.38, 0.71, 0.71),
        "b": (47.54, 3.23, 3.23),
        "c": (33.35, 0.74, 0.74),
    },
    "7 Zitrin-LTM": {
        "dt_ab": (-96.45, 5.41, 5.85),
        "dt_cb": (-27.64, 6.52, 4.80),
        "a": (5.66, 0.15, 0.14),
        "b": (9.77, 0.47, 0.51),
        "c": (8.80, 0.67, 0.46),
    },
}
PHOT = {
    "dt_ab": (-116.60, 10.83, 9.29),
    "dt_cb": (-48.30, 3.58, 3.97),
    "a": (4.43, 1.52, 1.60),
    "b": (8.00, 3.42, 2.34),
    "c": (6.43, 1.25, 1.13),
}
# Agrawal et al. Sect. 4 (SALT3-NIR, Pantheon+ flat LCDM): pooled mean and per-model range, mag
AGRAWAL_POOLED_MAG = 1.073
AGRAWAL_RANGE_MAG = (0.622, 2.391)
# Pierel et al. 2024 (arXiv:2403.18954v2): host-dust A_V of the joint BayeSN fit (Table 3, observed;
# context only) and Table 4 (simulated) delta-mu rows, (median, +err, -err) in mu units
PIEREL_AV = (1.21, 0.11, 0.11)
PIEREL_MICROLENSING = {"a": (0.22, 0.19, 0.20), "b": (-0.63, 2.51, 1.49), "c": (0.30, 0.20, 0.47)}
# shared by all images (standard-candle calibration): intrinsic SN Ia scatter and m_B interpolation
PIEREL_COMMON = {
    "a": (0.40, 0.56),
    "b": (0.76, 0.91),
    "c": (0.59, 0.83),
}  # b m_B: mean of +0.99 / -0.83


def ln_sigma(value: tuple[float, float, float]) -> float:
    """Half the 16-84 % width in ln space for a positive quantity."""
    med, hi, lo = value
    return 0.5 * (math.log(med + hi) - math.log(med - lo))


def mag(x: float) -> float:
    return 2.5 * math.log10(x)


def sigma_ml_ln(image: str, published: bool = True) -> float:
    """Microlensing term in ln mu: Pierel et al. Table 4 half-width over mu_phot (the published
    value the plan asks for), or the 0.15 mag placeholder."""
    if not published:
        return SIGMA_ML_MAG * math.log(10) / 2.5
    _, hi, lo = PIEREL_MICROLENSING[image]
    return 0.5 * (hi + lo) / PHOT[image][0]


def sigma_common_ln(image: str) -> float:
    """Shared standard-candle calibration part of the photometric ln-mu error (Table 4 rows)."""
    return math.hypot(*PIEREL_COMMON[image]) / PHOT[image][0]


def image_z(model: dict, image: str, published_ml: bool = True) -> float:
    """(ln mu_obs - ln mu_model) / combined ln sigma, microlensing term included."""
    s_ml = sigma_ml_ln(image, published_ml)
    s = math.hypot(math.hypot(ln_sigma(PHOT[image]), ln_sigma(model[image])), s_ml)
    return (math.log(PHOT[image][0]) - math.log(model[image][0])) / s


def ratio_z(
    model: dict, image: str, published_ml: bool = True, remove_common: bool = True
) -> float:
    """Same for mu_image / mu_a. With remove_common the shared calibration term cancels (plan: the
    ratio "removes the absolute-calibration term"); without it, all errors add in quadrature."""
    s2 = sum(
        ln_sigma(PHOT[k]) ** 2 + ln_sigma(model[k]) ** 2 + sigma_ml_ln(k, published_ml) ** 2
        for k in (image, "a")
    )
    if remove_common:
        s2 -= 2 * sigma_common_ln(image) * sigma_common_ln("a")
    obs = math.log(PHOT[image][0] / PHOT["a"][0])
    pred = math.log(model[image][0] / model["a"][0])
    return (obs - pred) / math.sqrt(s2)


def offset_mag(model: dict) -> float:
    """Rule 1: mean over images of 2.5 log10(mu_model / mu_phot)."""
    return sum(mag(model[i][0] / PHOT[i][0]) for i in IMAGES) / len(IMAGES)


def h0_mst(model: dict) -> dict:
    """Rule 2: H0 at which one mass sheet (mu/lambda^2, lambda dt) fits both mu and delays."""
    mean_ln = sum(math.log(model[i][0] / PHOT[i][0]) for i in IMAGES) / len(IMAGES)
    lam_mu = math.exp(0.5 * mean_ln)
    lam_t = 0.5 * sum(PHOT[k][0] / model[k][0] for k in ("dt_ab", "dt_cb"))
    h0 = H0_FID * lam_mu / lam_t
    # post hoc: 1-sigma in ln H0_MST of the image-mean ln(mu ratio): photometric errors fully
    # correlated across images (conservative), model and microlensing terms independent; the two
    # delays share image b but no covariance is published, so they are taken as independent
    n = len(IMAGES)
    s_phot = sum(ln_sigma(PHOT[i]) for i in IMAGES) / n
    s_ind2 = sum(ln_sigma(model[i]) ** 2 + sigma_ml_ln(i) ** 2 for i in IMAGES) / n**2
    s_mu = math.sqrt(s_phot**2 + s_ind2)
    s_t = 0.5 * math.sqrt(
        sum(ln_sigma_abs(PHOT[k]) ** 2 + ln_sigma_abs(model[k]) ** 2 for k in ("dt_ab", "dt_cb"))
    )
    s_ln = math.hypot(0.5 * s_mu, s_t)
    return {
        "lambda_mu": lam_mu,
        "lambda_t": lam_t,
        "h0_mst": h0,
        "in_range": H0_RANGE[0] <= h0 <= H0_RANGE[1],
        "h0_mst_1sigma": [h0 * math.exp(-s_ln), h0 * math.exp(s_ln)],
    }


def ln_sigma_abs(value: tuple[float, float, float]) -> float:
    """ln_sigma of |value| for a negative delay (the + and - errors swap roles)."""
    med, hi, lo = value
    return ln_sigma((abs(med), lo, hi) if med < 0 else value)


def run() -> tuple[dict, dict]:
    """Frozen-rule result and the post hoc block (written to separate files)."""
    rows = {}
    posthoc_rows = {}
    for name, m in TABLE1.items():
        mst = h0_mst(m)
        rows[name] = {
            "image_z": {i: image_z(m, i) for i in IMAGES},
            "ratio_z": {f"{i}/a": ratio_z(m, i) for i in ("b", "c")},
            "offset_mag": offset_mag(m),
            **{k: mst[k] for k in ("lambda_mu", "lambda_t", "h0_mst", "in_range")},
        }
        posthoc_rows[name] = {
            "h0_mst_1sigma": mst["h0_mst_1sigma"],
            "image_z_placeholder_ml": {i: image_z(m, i, published_ml=False) for i in IMAGES},
            "ratio_z_quadrature": {
                f"{i}/a": ratio_z(m, i, remove_common=False) for i in ("b", "c")
            },
        }
    flags_image = [i for i in IMAGES if all(abs(r["image_z"][i]) > Z_FLAG for r in rows.values())]
    flags_ratio = [
        k for k in ("b/a", "c/a") if all(abs(r["ratio_z"][k]) > Z_FLAG for r in rows.values())
    ]
    offsets = [r["offset_mag"] for r in rows.values()]
    pooled = sum(offsets) / len(offsets)
    no6 = [r["offset_mag"] for n, r in rows.items() if not n.startswith("6 ")]
    posthoc = {
        "rows": posthoc_rows,
        "pooled_offset_mag_without_model6": sum(no6) / len(no6),
        "n_usable_models_1sigma_reaching_h0_range": sum(
            r["h0_mst_1sigma"][0] <= H0_RANGE[1] and r["h0_mst_1sigma"][1] >= H0_RANGE[0]
            for n, r in posthoc_rows.items()
            if not n.startswith("6 ")
        ),
    }
    result = {
        "rows": rows,
        "flags_image": flags_image,
        "flags_ratio": flags_ratio,
        "outcome": "F" if flags_image or flags_ratio else "N",
        "pooled_offset_mag": pooled,
        "offset_range_mag": [min(offsets), max(offsets)],
        "n_models_offset_positive": sum(o > 0 for o in offsets),
        "n_models_h0_mst_in_range": sum(r["in_range"] for r in rows.values()),
        "agrawal_pooled_mag": AGRAWAL_POOLED_MAG,
        "agrawal_range_mag": list(AGRAWAL_RANGE_MAG),
        "rule1": "reproduces"
        if abs(pooled - AGRAWAL_POOLED_MAG) <= REPRODUCE_MAG
        else "method-dependent",
        "inputs": {
            "pascale2025_table1_models": TABLE1,
            "pierel2024_photometry": PHOT,
            "pierel2024_av": PIEREL_AV,
            "pierel2024_microlensing": PIEREL_MICROLENSING,
            "pierel2024_common_calibration": PIEREL_COMMON,
            "sigma_ml_ln": {i: sigma_ml_ln(i) for i in IMAGES},
            "agrawal_reference_mag": "SALT3-NIR, Pantheon+ flat LCDM",
        },
        "assumptions": {
            "sigma_ml_mag": SIGMA_ML_MAG,
            "z_flag": Z_FLAG,
            "reproduce_mag": REPRODUCE_MAG,
            "h0_range": list(H0_RANGE),
        },
    }
    return result, posthoc


def main() -> None:
    result, posthoc = run()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    (OUT / "posthoc.json").write_text(json.dumps(posthoc, indent=2) + "\n")
    for name, r in result["rows"].items():
        zi = " ".join(f"{k} {v:+.2f}" for k, v in r["image_z"].items())
        zr = " ".join(f"{k} {v:+.2f}" for k, v in r["ratio_z"].items())
        print(
            f"{name:18s} z[{zi}] ratio[{zr}] offset {r['offset_mag']:+.3f} mag "
            f"lam_mu {r['lambda_mu']:.3f} lam_t {r['lambda_t']:.3f} H0_MST {r['h0_mst']:.1f} "
            f"(post hoc {posthoc['rows'][name]['h0_mst_1sigma'][0]:.1f}-"
            f"{posthoc['rows'][name]['h0_mst_1sigma'][1]:.1f})"
        )
    print(
        f"outcome {result['outcome']} pooled {result['pooled_offset_mag']:.3f} "
        f"range {result['offset_range_mag']} rule1 {result['rule1']}"
    )


if __name__ == "__main__":
    main()
