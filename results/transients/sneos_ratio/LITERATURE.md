# SN Eos 101.2 / 101.1 flux ratio: literature check (2026-10-25)

The check named in `PLAN.md` ("whether Allingham et al. 2026 or the SN Eos discovery paper already report this
ratio"), done after `RESULT.md` was frozen and merged (#541). Tool `scripts/sneos_literature.py`, tests
`tests/test_sneos_literature.py`; every number below is in `literature.json`. Published photometry is observed
(published); published |μ| are model_predictions; ratios are derived.

## What the papers report
- Neither paper states a 101.2 / 101.1 flux ratio or discusses microlensing of SN Eos (HTML full texts searched for
  "flux ratio", "microlens"). Coulter et al. 2026 (arXiv:2601.04156v2) Table S1 gives per-image AB magnitudes at 17
  band × epoch entries: HST PID 17504 (WFC3 F814W / F110W, 2024-03-11 to 03-18, difference against 2012 templates),
  JWST 6882 (2025-09-01) and 9493 (2025-10-08). Ratios here are computed from that table.
- Published |μ| (both quoted with a ~20 % systematic per image, `mu_systematic_fraction_per_image`): Allingham et al. 2026
  (arXiv:2602.14074) Table 3, 25.79 / 27.49 → R_model 1.066; Coulter et al. Table S2, 25.46 / 29.88 → 1.174. Our
  pinned v2.1 maps give 1.162. Published model range used below: 1.066–1.174.
- Delay (`tau_published_days`): our pinned map has 101.2 arriving 0.32 d **before** 101.1 (`sneos/summary.json`,
  Δt −0.324 d), matching Allingham et al. (−0.37 d, −5.00 / +0.37; "101.2 must arrive before 101.1"). Coulter et al.
  Table S2 lists 101.2 0.99 ± 0.33 d after 101.1. The 2026-10-25 CHANGELOG line "101.2 trails by 0.32 d" was a wording error (corrected).

## Ratios from the published photometry (derived)
| Group | R = F(101.2)/F(101.1) | σ_ln (stat) | σ_ln scaled by √(χ²/dof) | χ² / dof |
|---|---|---|---|---|
| JWST 6882 (2025-09-01), 5 bands | 1.100 | 0.008 | 0.038 | 93.9 / 4 |
| JWST 9493 (2025-10-08), 3 bands | 1.118 | 0.016 | 0.016 | 0.25 / 2 |
| HST F814W 2024, 4 epochs | 1.243 | 0.137 | 0.137 | 2.4 / 3 |
| HST F110W 2024, 5 epochs | 1.805 | 0.039 | 0.047 | 5.7 / 4 |

- Both JWST epochs sit inside the model range. The 6882 entries span 1.009 (F410M) to 1.191 (F444W) (`entries`);
  their χ² 93.9 / 4 means the scatter exceeds the quoted errors (cause not tested).
- **HST F110W 2024 is 0.43 in ln above the top of the model range** (`offset_vs_model_range`), beyond the plan's
  0.30 floor (an ASSUMPTION). F814W 2024 is within the floor (σ_ln 0.137).
- Flux-drop ratio 6882 → 9493 from the same table (comparable to our difference images): F277W 1.33, F356W 1.13,
  F444W 1.58 (σ_rel 0.12, 0.14, 0.11). Our 9493 pairs: F356W 1.37 ± 0.21, F444W 1.93 ± 0.41 (`summary.json`), each
  within about 1σ of these (not stored). Our combined 1.382 is above the published direct-image ratios 1.100 / 1.118.

## Reading (post hoc; not a candidate)
- The published model range is wider than the single R_model in `PLAN.md`; with it, our "agrees" stands.
- The 2024 F110W ratio is a flux-ratio offset in published photometry at an epoch our plan never covered.
  Intrinsic variability through the delay is too small to make it: the week-long F110W light curves have slopes
  −0.004 ± 0.018 (101.1) and −0.014 ± 0.011 mag/d (101.2) (`hst_f110w_slopes`); at the 2σ bound 0.040 mag/d times the
  largest published delay (1.32 d Coulter; 5.37 d Allingham with its systematic) the ratio changes by ≤ 0.05 or
  ≤ 0.20 in ln (`delay_times_slope`), against 0.43. This holds only if the light curve did not change faster between
  the HST visits than across the week. A ratio
  that differs between 2024 (rest-frame FUV, days after explosion) and 2025 (end of plateau) points to an effect
  that varies in time or wavelength. Ordinary candidates, untested: microlensing by intracluster stars at |μ| ≈ 25–30
  (stronger for the compact early photosphere), HST difference-photometry systematics (template, PSF, host
  residuals), and macro-model error (would not change with time). A difference between direct and flux-drop ratios
  in 2025 also hints at a time-varying ratio or host light in the direct photometry (not tested).
- Not run through `/vet-candidate`: it is a published measurement, not a screen hit of ours, and microlensing is
  the expected first explanation. An own measurement on the HST 17504 / 2012 images (frozen plan, injections) would
  test it.

## Limits
Table values transcribed by hand from the arXiv HTML (v2 of 2601.04156; 2602.14074 as served 2026-10-25); quoted
errors taken at face value; the two μ tables are model averages from different runs of the same group.
