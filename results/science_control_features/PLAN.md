# Plan: science vs control feature distributions (TASKS "Then" 8, D-022)

Written 2026-11-06 before any fetch. Samples: `smacs0723_nircam` (science) and `ceers_t021` (control) from
`configs/reference_sample.yaml`, both on jwst 3.0.0 catalogs and DJA v7.4 `dja05` colours since 2026-11-06.

## Question
D-022 treats CEERS t021 as the control whose top-k statistics are compared with SMACS. That comparison is fair only
if the ranked (gated, galaxy-stratum) feature distributions are comparable, or if their differences are understood.
D-022 addendum: "revisit if the two fields' size distributions differ strongly".

## Hypotheses and ordinary explanations (all expected, none exotic)
- H1 depth: SMACS ERO and CEERS t021 have different exposure times, so `ref_mag` and the S/N-gated morphology
  shift (fainter sources in the deeper field; small sizes near the limit).
- H2 cluster: SMACS holds a z = 0.39 cluster (bright, large, red members; ICL), so `ref_log_isophotal_area`,
  `ref_log_nn_dist` (crowding) and the F150W−F200W / F200W−F277W colours differ at bright magnitudes.
- H3 processing: the same pipeline version (3.0.0) and DJA recipe in both fields, so no processing difference is
  expected; any difference in `ref_log_ci_*` at fixed magnitude would point to PSF / background differences.

## Measurement
- Run the pipeline on both samples (one run), read the saved `features` and gate tables of the ranked galaxy strata.
- Per shared feature: median, 16/84 % quantiles, two-sample KS D and the Anderson–Darling p, on all ranked galaxies and
  in `ref_mag` bins (≤ 24, 24–26, 26–28) to separate depth (H1) from population (H2).
- Shared features are those both fields compute: `ref_mag`, colours for band pairs present in both
  (F150W−F200W, F200W−F277W, F277W−F356W), and the reference-band morphology (F200W in both).

## Decision rule (ASSUMPTION thresholds)
- KS D ≤ 0.1 within the magnitude bins for every shared feature: comparable; no change.
- D > 0.2 in a magnitude bin for any feature: "differ strongly" → the control comparison needs magnitude (and, for
  size, area) matching before any science vs control top-k statement; record as D-022 addendum and a TASKS item.
- In between: record, no rule change.
- Feature differences explained by H2 at bright magnitudes only (cluster members) do not count against the control.

## What changes the plan
A strong difference at fixed magnitude in a morphology feature (H3) would reopen the jwst 3.0.0 / DJA comparability
question for both fields' rankings.

## Data
Pipeline catalogs (small) and both DJA catalogs (SMACS small; CEERS 250.5 MB, D-022 stated reason) to the scratch data
root; deleted after the run (Cloud disk decision).
