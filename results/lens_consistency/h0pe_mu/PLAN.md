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

## Addendum A (2026-10-26, written after reading the abstracts and Pascale et al. 2025 Table 1, before any computation)
Inputs fixed: Pascale et al. 2025 (arXiv:2403.18902v2) Table 1, unblinded rows (models 1–7, photometry). The
photometric μ of Pierel et al. 2024 are standard-candle magnifications under a fiducial cosmology, with host dust
fitted (A_V 1.21 ± 0.11; Pierel et al. Table 3), so they are the reference.
1. Offset O_m = mean over images of 2.5 log10(μ_m,i / μ_phot,i) (mag; derived). Compare with Agrawal et al.
   (arXiv:2510.07637v2) per-model Δμ 0.62–2.39 mag (SALT, Pantheon+) and pooled ≈ 1.07 mag. Rule (ASSUMPTION): the
   table-level offset "reproduces" the paper if |O_m − Δμ_m| ≤ 0.3 mag for the pooled value; otherwise
   "method-dependent" (the > 1 mag figure depends on that paper's light-curve re-fit, not on the tabulated μ).
2. Mass-sheet test (derived): λ_μ,m = sqrt(exp(mean_i ln(μ_m,i / μ_phot,i))), λ_t,m = mean of Δt_obs / Δt_m over
   ab and cb (fiducial H0 70). One mass sheet explains both iff H0_MST,m = 70 λ_μ,m / λ_t,m lies in 64–80 km/s/Mpc
   (ASSUMPTION: Planck to SH0ES ± ~3σ). Report per model; model 6 (WSLAP+, Δt_cb sign opposite) is reported but
   not interpreted.

## Deviations
- Post hoc (after the run): the pooled offset without WSLAP+ and 1σ intervals on H0_MST (`posthoc.json`).
- Deviation (after the review, before reading any changed outcome): the first run kept the 0.15 mag placeholder
  although Pierel et al. Table 4 publishes microlensing; the rule now uses Table 4 as written above. The first
  ratio errors added the shared calibration twice; it now cancels as the ratio rule intends. Both first-run
  variants stay in `posthoc.json`; the outcome (N) is the same in each.
