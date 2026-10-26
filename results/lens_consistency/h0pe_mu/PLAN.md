# SN H0pe (PLCK G165.7+67.0): measured image magnifications vs the published blind lens models: plan

Frozen before any published magnification table is read in this session (2026-10-26).

## Question
Do the standard-candle magnifications of SN H0pe images A, B, C (SN Ia; Pierel et al. 2024) agree with the
blind lens-model predictions compiled before unblinding (Pascale et al. 2025), within the stated errors?

## Hypotheses and ordinary explanations
- H0 (expected): every image agrees with the model ensemble.
- Ordinary explanations of a disagreement, tested in this order: SN Ia standardisation scatter and dust
  (the paper's own error budget), stellar microlensing (typically 0.1–0.3 mag at these μ; ASSUMPTION until a
  published estimate is found), model-to-model systematics (spread of the blind models), substructure.

## Frozen rule (ASSUMPTION)
For each image i and each blind model m: z_im = (ln μ_obs,i − ln μ_m,i) / sqrt(σ_obs² + σ_m²) in ln space,
using the published 68 % intervals (symmetrised in ln as half the 16–84 width).
- Image flag: |z_im| > 3 for every model with a prediction for that image, after adding a 0.15 mag
  (σ_ln 0.138) microlensing term in quadrature (ASSUMPTION; replaced by a published value if one exists).
- Ratio flag: the same on μ_B/μ_A and μ_C/μ_A (removes the absolute-calibration term).
- Outcome N (no flag), F (flag: vet before any interpretation), U (undetermined: inputs missing).

## Inputs
Published tables only (arXiv e-prints, sha256 recorded). No images are read in this step.

## What would change the plan
F: vet (microlensing amplitude, dust, model spread); N: record and close; U: record the missing input.
