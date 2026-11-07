# Plan: segment-overlap MIRI band association (2026-11-07)

Trigger: `PLAN_sizemerge.md` outcome (D-014 addendum 5). A size-scaled radius that passes the chance rule merges only
1 / 4 of the vetted "single band" F770W rows (`_221`, `_239`, `_118`, `_240`), whose other-band centroids lie
0.18–0.66″ away.

**Hypothesis (ordinary).** These rows are the same objects as their other-band detections; the centroids differ
because isophotal centroids of extended or blended MIRI sources shift between bands, but the footprints overlap.
Footprint overlap should then associate them at a lower chance rate than a centroid radius large enough to reach
them.

**Method (frozen before measuring).**
- Data: MAST L3 `jw02736-o002_t001_miri_{f770w,f1000w,f1500w,f1800w}_segm.fits` and `_cat.ecsv` (jwst 3.0.0),
  streamed to scratch and deleted after use.
- Each non-reference band's segmentation map is resampled onto the F770W pixel grid by nearest neighbour through
  the two WCSs. For an F770W segment i and a band-b segment j: `n_ij` = shared F770W-grid pixels, `A_i`, `A_j` =
  their pixel counts on that grid, overlap `o_ij = n_ij / min(A_i, A_j)`.
- Pairs are F770W row ↔ band-b row only (reference-anchored; no chaining through other bands). Accept pairs with
  `o_ij ≥ f`, one-to-one, greedy by decreasing `o_ij`. Grid f ∈ {0.1, 0.25, 0.5, 0.75}.
- "Extra" pairs: accepted pairs whose catalog centroids are more than 0.1″ apart (pairs the current radius cannot
  make). Chance: repeat with each band-b map shifted on the sky by the same four 11.3″ vectors as
  `scripts/miri_sizemerge.py` (base (±8″, ±8″), rotated 120° per band); chance extra = mean extra pairs over shifts.
- **Rule (function `adopt_f`, unit-tested before the run):** chance fraction = chance extra / real extra; adopt the
  smallest f whose chance fraction ≤ 0.10 (ASSUMPTION, as in the size-aware plan). None passes → keep 0.1″ and
  record the limit.
- Check, not tuning: at the adopted f, the four vetted rows should each gain ≥ 1 band (expected 4 / 4); `_176`
  stays single (no other-band detection within 5″).

**What would change the plan.** ≥ 3 / 4 vetted merges at a passing f: segment overlap explains the "single band"
rows; next is wiring it as a `merge_bands` option and re-ranking. < 3 / 4 at every passing f: footprint
association does not explain them either; record it, close the MIRI association line and keep 0.1″ with the
stated limitation. If the chance fraction is above 0.10 already at f = 0.75: overlap is no better than radius.
