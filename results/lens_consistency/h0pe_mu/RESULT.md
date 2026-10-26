# SN H0pe (PLCK G165.7+67.0): standard-candle magnifications vs the seven blind lens models: result

Plan `PLAN.md` (rules frozen before any table was read; addendum A before any computation). Published tables only,
no images (so there is no contact sheet; the tables were read in the e-prints). Numbers: `summary.json`, written by
`scripts/h0pe_mu_check.py`. Labels: observed (Pierel et al. 2024 photometric μ and delays, SN Ia standard candle with
host dust fitted, A_V 1.21 ± 0.11), model_prediction (Pascale et al. 2025 Table 1, unblinded), derived (z, offsets,
λ, H0_MST), assumption (thresholds below).

## Frozen outcome: **N** (no image flag, no ratio flag)
- Per-image |z| > 3 in every model (0.15 mag microlensing term added): no image. Largest |z| outside WSLAP+: Zitrin-
  analytic image c, −3.10; the rule needs all seven.
- Ratios μ_b/μ_a and μ_c/μ_a: |z| ≤ 0.86 in all 14 model-ratio pairs. The image-to-image pattern of every model
  agrees with the SN.

## Rule 1: magnification offset (addendum A)
- All 7 models predict more magnification than the SN: offset +0.247 (Chen 2020) to +0.883 mag (Zitrin-analytic),
  +1.714 mag for WSLAP+.
- Pooled over the 7 models: 0.633 mag; post hoc without WSLAP+: 0.453 mag. Agrawal et al. (arXiv:2510.07637v2)
  report 1.073 mag pooled (SALT3-NIR, Pantheon+), per-model 0.622–2.391 mag.
- Frozen reading: **method-dependent** (|0.633 − 1.073| > 0.3 mag). The tabulated μ give about 0.4 mag less offset than
  that paper's light-curve re-fit; its "> 1 mag" figure is not a property of the published μ alone. Which choices of
  that re-fit add the difference is not tested here.

## Rule 2: one mass sheet for μ and delays (addendum A)
- H0_MST (one λ scaling μ by 1/λ² and delays by λ, fiducial H0 70): GLAFIC 85.6, Zitrin-analytic 92.1, LENSTOOL
  88.8, MARS 98.0, Chen 2020 80.7, Zitrin-LTM 53.7 km/s/Mpc; WSLAP+ 1079.2 (its Δt_cb has the opposite sign; not
  interpreted).
- Frozen rule (point estimate in 64–80): **0 of 7** models in range.
- **Post hoc** 1σ intervals (photometric μ errors taken as fully correlated): 68.2–95.4 (Chen 2020) up to 72.0–133.2
  (MARS); Zitrin-LTM 44.1–65.5 (post hoc). Every interval except WSLAP+'s reaches 64–80, so these errors cannot reject a mass sheet
  plus a normal H0. The point estimates lean to a high H0 for the five models that over-predict μ but match the
  delays.

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
