# Third epoch for SMACS `f200w_384`: does program 12774 reproduce the 2736 → 6882 shift?

Frozen in the claim commit, before any 12774 centroid is measured. Follows `../epochs_1816_384_img/RESULT.md`:
`384` was undetermined (F150W residual 25.03 mas above all 16 peers; F444W 22.61 mas below the peer p95 27.60;
post hoc, the two band residuals point 12.3° apart).

## New data (observed, MAST query 2026-11-07)
- Program 12774 o002 t032 (target `SMACS0723-GLSN-1`), NIRCam F115W/F150W/F300M/F444W, observed 2026-07-06,
  public 2026-07-07; level-3 `_i2d.fits` and `_cat.ecsv` on `s3://stpubdata/jwst/public/jw12774/L3/t/o002/`.
  It covers `384` (2″ cone). `1816` is not covered by 12774 (its record is closed anyway).
- 12774 is 31 days after 6882 o057 (2026-06-05), so it adds no proper-motion leverage; it is an independent
  re-measurement of the 2026 position with a different program, dither and drizzle grid.

## Hypotheses
- H-move: the 2736 → 6882 residual is real source motion (≈ 6 mas/yr; model_prediction from that residual).
  Predicts 2736 → 12774 ≈ 2736 → 6882 (expected change over 31 d ≲ 0.5 mas).
- H-6882: the residual is a 6882-specific systematic (drizzle grid, dither, a detector artefact at that
  position). Predicts 2736 → 12774 consistent with the peer null.
- Not separable here: a 2736-specific systematic predicts the same as H-move (both use epoch 1).

## Method (same as `../epochs_1816_384_img/PLAN.md`, epoch 2 replaced by 12774)
- `scripts/smacs_epochs_img.py` functions unchanged: 1.5″ S3 byte-range cutouts, `centroid_2dg` on a 7 × 7 px
  box, peers = 2736 F150W CAT sources within 60″, ±1.0 mag, ±50 % `semimajor_sigma`, isolated (1″), with a
  mutual 0.3″ F150W match in the 12774 CAT; frame tie = median peer offset; leave-one-out peer residuals.
- Bands F150W and F444W. Synthetic +20 mas injection on the 12774 cutouts first (pass: median within 5 mas).
- Also reported (derived, not part of the rule): 6882 → 12774 residual of `384` with the same peers that are
  in all three catalogs.

## Frozen outcome rules (ASSUMPTIONs)
- ≥ 10 peers with valid centroids per band, else "undetermined".
- "shift reproduced": in both bands the 2736 → 12774 residual exceeds the peer 95th percentile and points
  within 45° of the same band's 2736 → 6882 residual vector (`../epochs_1816_384_img/summary.json`).
- "shift not reproduced": in both bands the 2736 → 12774 residual is ≤ the peer median.
- Anything else: "undetermined".

## What would change the plan
"reproduced" removes H-6882; the record becomes "compact, moves ≈ 6 mas/yr between 2022 and 2026 in three
images; likely a Galactic star (needs spectroscopy or a later epoch to exclude a 2736 systematic)".
"not reproduced" retires the 2736 → 6882 shift as a 6882 systematic; record "compact, unresolved; no motion
established; needs spectroscopy". "undetermined" retires the astrometric route for `384`.
