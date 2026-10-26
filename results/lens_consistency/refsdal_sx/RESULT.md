# SN Refsdal (MACS J1149.5+2223): measured SX delay and magnification ratio vs the blind models: result

Plan `PLAN.md` (frozen before either table was read). Published tables only, no images (no contact sheet).
Numbers: `summary.json`, written by `scripts/refsdal_sx_check.py`. Labels: observed (Kelly et al. 2023 Table
"tab:measurements": SX − S1 delay 376.02 d, 16–84 % 370.50–381.65; μ(X)/μ(1) 0.299, 0.272–0.354; maximum-likelihood
values used), model_prediction (Treu et al. 2016 Table "tab:predictions": seven blind models, random errors only,
each model's fiducial cosmology), derived (z), assumption (|z| > 3).

## Frozen outcome: **N** (5 of 7 blind models within |z| ≤ 3 on both quantities)
| Model | Δt_X1 (d) | z_Δt | μX/μ1 | z_μ |
|---|---|---|---|---|
| Die-a | 262 ± 55 | +2.06 | 0.31 ± 0.10 | −0.10 |
| Gri-g | 361 +19 / −27 | +0.63 | 0.36 +0.11 / −0.09 | −0.61 |
| Ogu-g | 311 ± 24 | +2.64 | 0.27 ± 0.05 | +0.45 |
| Ogu-a | 336 ± 21 | +1.84 | 0.27 ± 0.03 | +0.59 |
| Sha-g | 277 +11 / −21 | **+5.84** | 0.25 +0.05 / −0.02 | +0.96 |
| Sha-a | 233 +46 / −13 | **+4.76** | 0.19 +0.01 / −0.04 | +2.32 |
| Zit-g | 224 ± 262 | +0.58 | 0.31 ± 0.05 | −0.17 |

Reach (observation alone): a delay offset ≥ 16.7 d or a ratio factor ≥ 1.48 is needed for a flag.

## Reading (no exotic interpretation)
- The two flags (Sharon) have a documented ordinary cause: Kelly et al. 2023 report a clerical error in the redshifts
  entering the Sharon time-delay distance (delays revised by ≈ 27 %). Post hoc with their recalculated delays,
  Sharon-a* 298 d and Sharon-g* 331 d give z +2.06 and +1.37. Model / calculation error, not a lensing anomaly.
- The magnification ratio agrees with every blind model (|z| ≤ 2.32).
- Post hoc pattern: the measured delay is longer than all seven blind predictions (z +0.58 to +5.84). At fixed lens
  model Δt ∝ 1/H0, so this is the H0 < 70 preference Kelly et al. infer from the same delay, not an anomaly in the
  lens geometry (the ratio, H0-free, agrees).
- Post hoc, models outside the frozen outcome (Kelly et al. tables): Jauzac15.1 (blind, not in Treu et al.'s table)
  z −1.61 / +2.68; Zitrin-c* (+5.08), Chen2020* (+4.02) and Zitrin2020-ltm* (+4.98) exceed 3 in delay; the others
  (Oguri-g*, Oguri-a*, Jauzac15.2*, Keeton*, Zitrin2020-p*) are within 3. Their Zitrin-c* μ lower error is printed
  without a sign and read as −0.023.

Not a candidate. SN Refsdal SX is the one cluster-scale case where a blind delay prediction is tested at the
per-cent level; the best blind models pass.

## Deviations
None from the plan. The ML value (not the median) is used as the measured point, as Kelly et al. tabulate it first.

## What would decide more
Nothing pending: no further reappearance is predicted in the future (SY is predicted at −4,521 to −7,665 d, in the past; model_prediction).
Rerun only if revised blind-table errata appear.

## Cost
Two arXiv e-print sources (≈ 15 MB, scratch, deleted after use), < 1 s compute.
