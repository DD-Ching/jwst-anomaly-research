# Plan: SMACS point-source screen before top-k vetting (TASKS Next M2)

Written 2026-11-07 ~13:15Z, before the pipeline run and before any `flux_radius` of a star, image or target was
read. Follows `results/smacs_jwst3/zphot_nofit/RESULT.md` (D-022 addendum 5: eazy no-fit rows are compact).

## Question and hypotheses
- Can DJA v7.4 `flux_radius` (detection image, pixels) separate point sources from lensed galaxies in SMACS well
  enough to label top-k rows as point-like before vetting?
- H-pt: `f200w_1816` (rank 14) and `f200w_384` (rank 19) are point sources (faint stars or compact non-lensed
  objects), so they are not lens-model residual candidates.
- Ordinary alternatives: compact high-z galaxies or compact lensed knots also have small `flux_radius`; at
  24–28 mag the radius is noisy. The validation below bounds the second case only through catalogued images.

## Screen (ASSUMPTION thresholds)
- Locus stars: the run's star-stratum sources (Gaia-only split) matched to DJA phot at 0.2″, with ref_mag in
  [19, 24] (ASSUMPTION: brighter stars saturate in F200W and inflate the radius); if fewer than 10, use [18, 24].
- `r_psf` = median `flux_radius` of locus stars; `s` = 1.4826 × MAD; point-like ⇔ `flux_radius ≤ r_psf + 3 s`.

## Validation, frozen before reading the numbers
- V1 star recovery (held out): split locus stars by DJA `id` parity; derive the cut on one half, test on the
  other, both ways. Pass if recovery ≥ 0.90 in both directions.
- V2 lensed images kept: Mahler+2022 ICLv2 catalogued multiple images (`arcs.dat`, `lens_consistency.image_list`)
  matched to DJA phot at 0.2″; pass if ≥ 10 match and ≤ 10 % are point-like.
- Fewer than 10 locus stars in [18, 24] → screen undetermined.

## Decision rule
- V1 and V2 pass → screen validated on this field. Apply to `f200w_1816` and `f200w_384`: point-like → labelled
  point-like (not a lens-residual candidate); otherwise kept for vetting.
- Either fails → screen not validated; record numbers; the two rows stay for normal vetting.
- Descriptive only (no decision): point-like fraction of galaxy-stratum rows per ref_mag bin, and of the top 20.

## What changes the plan
A validated screen becomes a proposed pre-vetting label next to D-019–D-021 (adoption into the screen stack is a
separate PR); a failed V2 means compact lensed knots look stellar and size cannot be used for lensed fields.
