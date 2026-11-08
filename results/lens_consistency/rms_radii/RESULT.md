# Match / search radii scaled to the image-plane rms (result)

Plan: `PLAN.md` (frozen before the run). Numbers: `summary.json`, `census.ecsv`, `posthoc.json`.
Script: `scripts/rms_radii_census.py` (offline; functions `rms_scaled_radii`, `rms_outlier_driven`,
`image_plane_rms` in `scripts/lens_consistency.py`, unit-tested before the run).

## Frozen outcome: "rms scale outlier-driven"
- 55 stored `validate` image-plane tables read, 55 used (n ≥ 3 finite offsets): 30 with rms ≤ 1″, 10 "scalable"
  (1″ < rms ≤ 2″), 15 "gate-K failure, not scalable" (rms > 2″).
- In 6 of the 10 scalable tables the rms drops to ≤ 1″ once the single largest offset is removed
  (`rms_outlier_driven`), so the frozen reading is "outlier-driven" (≥ half). One image, not the model as a
  whole, sets the rms there; scaling every match radius by that rms would widen it for all images of the field.
- Decision per the plan: rms scaling is **not** wired into `images`. No candidate, no flag (a census of stored
  derived tables; no new photometry).

## Numbers (from `summary.json` / `census.ecsv`)
- Scalable tables: 82 catalogued images beyond the fixed 1.5″ match radius, 52 beyond the scaled radius.
- None of the 15 gate-K failures is single-outlier-driven (their rms stays > 1″ after one removal).

## Post hoc (not in PLAN.md; `posthoc.json`)
- 12 of the 25 tables with rms > 1″ have a robust scale (1.4826 × median |dθ|) ≤ 1″: 5 / 10 scalable and 7 / 15
  gate-K failures. So in about half of them a few images, not only one, set the rms; a leave-one-out test misses
  those. This motivates a robust (median-based) scale, but this run does not gate one.

## Limitations
- The 55 tables are not independent: 25 rms > 1″ tables come from 18 result directories, and several are the
  same field and model with a different input list (held-out fits, fixed vs fitted z, e.g. MACS0553 twice).
  The tally is over tables, not fields.
- Offsets are image-plane distances to the nearest prediction from `validate`. The census does not re-run
  `forced_check`, so it says nothing about how forced classes would change.
