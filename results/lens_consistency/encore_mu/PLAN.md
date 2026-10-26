# SN Encore (MACS J0138.0−2155): measured image magnifications vs the seven blind lens models: plan

Frozen before any published magnification or delay table is read in this session (2026-10-26, 03:25Z). Only the
search-engine abstract snippets have been seen (they quote two μ values; no per-model table, no errors).

## Question
Do the standard-candle magnifications of SN Encore images 1a and 1b (SN Ia; Pierel et al., arXiv:2509.12301) and
their ratio agree with the blind lens-model predictions (Suyu et al., arXiv:2509.12319) within the stated errors?
Same test as SN H0pe (D-034 addendum bn, `results/lens_consistency/h0pe_mu/`), on an independent cluster.

## Hypotheses and ordinary explanations
- H0 (expected): both images and the ratio agree with the model ensemble.
- Ordinary explanations of a disagreement, in order: SN Ia standardisation scatter and host dust (the paper's error
  budget), stellar micro- and millilensing (the paper's own simulations), model-to-model systematics (spread of the
  seven blind models), substructure near image 1b.
- What would change the plan: a ratio flag (image-specific, not removed by calibration) leads to `/vet-candidate`;
  a common offset again in the same sense as H0pe makes a two-cluster "shared over-prediction" question
  (calibration vs model slope/sheet) worth a record.

## Frozen rule (ASSUMPTION; identical to addendum bn where it applies)
For image i and blind model m: z_im = (ln μ_obs,i − ln μ_m,i) / sqrt(σ_obs² + σ_m² + σ_ml,i²), ln space, published
68 % intervals symmetrised as half the 16–84 width in ln.
- σ_ml,i: the paper's published micro/millilensing term per image if it exists; else 0.15 mag (σ_ln 0.138).
- Image flag: |z_im| > 3 for every model with a prediction for that image.
- Ratio flag: the same on μ_1b/μ_1a, with terms shared by both images (SN Ia luminosity, cosmology) cancelled
  where the paper separates them; otherwise all terms in quadrature (conservative) and say so.
- Offset O_m = mean over images of 2.5 log10(μ_m,i / μ_obs,i) (mag, derived); report per model and the unweighted
  pool. Reading "same sense as H0pe" iff the pooled O > 0 and ≥ 5 of 7 models have O_m > 0 (ASSUMPTION).
- Mass sheet (derived): λ_μ,m = sqrt(exp(mean_i ln(μ_m,i / μ_obs,i))), λ_t,m = Δt_obs / Δt_m(1b,1a) at the models'
  fiducial H0; H0_MST,m = H0_fid λ_μ,m / λ_t,m; "consistent with one sheet and a normal H0" iff in 64–80 km/s/Mpc.
  Point estimates are the frozen statistic; 1σ intervals are post hoc.
- Outcome N (no flag), F (flag: vet first), U (inputs missing).

## Inputs
arXiv e-prints only (source tarballs or HTML, sha256 recorded in the result). No images, so no contact sheet.
Script: `scripts/encore_mu_check.py`, values transcribed into it with the table and row they come from.
