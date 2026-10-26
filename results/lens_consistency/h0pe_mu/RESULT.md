# SN H0pe (PLCK G165.7+67.0): standard-candle magnifications vs the seven blind lens models: result

Plan `PLAN.md` (rules frozen before any table was read; addendum A before any computation). Published tables only,
no images (so there is no contact sheet; the tables were read in the e-prints). Numbers: `summary.json`, written by
`scripts/h0pe_mu_check.py`. Labels: observed (Pierel et al. 2024 photometric μ and delays, SN Ia standard candle with
host dust fitted, A_V 1.21 ± 0.11), model_prediction (Pascale et al. 2025 Table 1, unblinded), derived (z, offsets,
λ, H0_MST), assumption (thresholds below).

## Frozen outcome: **N** (no image flag, no ratio flag)
- Microlensing term: Pierel et al. Table 4 (simulated), σ_ln 0.044 / 0.250 / 0.052 for a / b / c, as the plan
  requires when a published value exists (the 0.15 mag placeholder is kept only as a check in `posthoc.json`).
- Per-image |z| > 3 in every model: no image. Largest |z| outside WSLAP+: Zitrin-analytic image c, −3.56; the rule
  needs all seven.
- Ratios μ_b/μ_a and μ_c/μ_a, with the shared standard-candle terms (Table 4 intrinsic scatter and m_B interpolation)
  cancelled as the plan intends: |z| ≤ 0.96 in all 14 model-ratio pairs. With every error in quadrature instead
  (check, `posthoc.json`) |z| ≤ 0.82. The image-to-image pattern of every model agrees with the SN.

## Rule 1: magnification offset (addendum A)
- All 7 models predict more magnification than the SN: offset +0.247 (Chen 2020) to +0.883 mag (Zitrin-analytic),
  +1.714 mag for WSLAP+.
- Pooled over the 7 models (unweighted mean): 0.633 mag; post hoc without WSLAP+: 0.453 mag. Agrawal et al.
  (arXiv:2510.07637v2) report 1.073 mag pooled (SALT3-NIR, Pantheon+ flat ΛCDM), per-model 0.622–2.391 mag.
- Frozen reading: **method-dependent** (|0.633 − 1.073| > 0.3 mag; also without WSLAP+). The tabulated μ against the
  photometric μ give a smaller offset than that paper's Δμ. Whether the gap comes from its light-curve re-fit, its
  pooling (per light-curve point) or its cosmology baseline is not tested here.

## Rule 2: one mass sheet for μ and delays (addendum A)
- H0_MST (one λ scaling μ by 1/λ² and delays by λ, fiducial H0 70): GLAFIC 85.6, Zitrin-analytic 92.1, LENSTOOL
  88.8, MARS 98.0, Chen 2020 80.7, Zitrin-LTM 53.7 km/s/Mpc; WSLAP+ 1079.2 (its Δt_cb has the opposite sign; not
  interpreted).
- Frozen rule (point estimate in 64–80): **0 of 7** models in range.
- **Post hoc** 1σ intervals (`posthoc.json`; photometric μ errors fully correlated across images, model and
  microlensing terms independent, the two delays independent although they share image b): 67.9–95.9 (Chen 2020)
  to 72.0–133.3 (MARS); Zitrin-LTM 43.9–65.8 (post hoc). All 6 usable intervals reach 64–80, so these errors cannot
  reject a mass sheet plus a normal H0. The point estimates lean to a high H0 for the five models that over-predict
  μ but match the delays.

## Reading (no exotic interpretation)
- The disagreement is one factor common to all three images (ratios agree). That is what a shared calibration term
  (SN Ia luminosity and the assumed cosmology, host-light subtraction without a template), a mass-sheet or radial-slope
  term in the models, or a uniform demagnification would all give. On these data they cannot be separated; the
  docs/exotic_lensing.md "demagnified images" signature is not distinguishable from those ordinary terms here.
- Not a candidate. Microlensing (simulated in Pierel et al. Table 4, e.g. δμ_a 0.22 +0.19 / −0.20) and host dust
  (fitted in the photometric μ) are not tested here as the source of the shared factor.

## What would decide it
The template-subtracted photometry (GO 4744; Agrawal et al. 2026 and Pascale et al. 2026, both in prep.) or the
revised lens models. Rerun `scripts/h0pe_mu_check.py` with their tables.

## Cost
Four arXiv e-prints (< 30 MB, scratch, deleted), < 1 s compute.
