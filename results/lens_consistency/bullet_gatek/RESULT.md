# Bullet Cluster: gate K fails as written; gate K-WLS (post-failure deviation) passes (2026-10-24)

Plan `PLAN.md` (frozen in 3613c8f, before any map value was read). Script `scripts/bullet_gatek.py`; tests
`tests/test_bullet_gatek.py`. Numbers from `summary.json` (plan), `posthoc_compare.json`, `posthoc_k1_sources.json` and
`posthoc_grid010/summary.json` (post hoc), and `summary_wls.json` (`PLAN_1b.md`, deviation). Every image-plane number is a model_prediction; the image list is observed.

## Result (frozen rules)
- **Gate K fails:** 135 gold images at their system z, 0.25″ grid: rms **1.196″** (rule ≤ 1.0″), 97.0 % within 3″
  (rule ≥ 80 %). The rms condition fails, so the combined gate fails.
- The failure comes from two systems, K1b and K1e (knots of galaxy K1, z 3.24): the solver predicts **one** image for
  each, where three are catalogued; their four unmatched images lie 5.3–7.1″ from the one predicted image.
- Lenstool's own `chires.dat` for this fit: rms 0.41″ over the same 135 images (χ² 92.1).
- Per the plan, no flag is read and the extra-image test does not run on these maps.

## Post hoc checks (not in the plan; they do not change the gate outcome)
- **Comparison with Lenstool** (`posthoc_compare.json`): over all 135 images our and Lenstool's per-image residuals
  have the same median (0.335″ vs 0.334″); outside K1b / K1e they correlate at 0.88. Without those two systems: 129
  images, rms 0.435″ (Lenstool 0.419″), 100 % within 3″. This agreement makes a loading error (orientation, units,
  D_LS / D_S scaling) unlikely; it is not a frame test against the JWST catalogue (frame check F is still to run).
- **Grid** (`posthoc_grid010/`): a 0.1″ grid gives the same rms (1.196″) and the same single images for K1b / K1e, so
  grid resolution is not the cause.
- **Caustic proximity** (`posthoc_k1_sources.json`): the back-traced sources of K1b and K1e scatter by only 0.061″ and
  0.026″ about their means, yet the predicted multiplicity switches between 1 and 3 depending on which image's own
  source is used (K1b: 1, 3, 3; K1e: 2, 3, 1; the reproduced neighbour K1d also gives 3, 3, 1). Since the
  own-image sources of one family lie within twice that scatter of each other (≤ 0.12″, derived, not stored) and give
  different multiplicities, a caustic passes within that distance of these sources, and the number of predicted
  images depends on the source-position estimator. Lenstool's estimator was not reproduced; this is a reading of the
  numbers above, not a separate test.

## Gate K-WLS (`PLAN_1b.md`, post-failure deviation; frozen in 3774ed9 before its run)
- The Abell 2744 precedent (D-034 addendum k) applied unchanged: only the family source changes, to D-041's
  magnification-weighted least-squares source.
- **Passes:** rms **0.391″**, 100 % of the 135 images within 3″. K1e gets two predicted images for three catalogued
  (its rms 0.369″: two catalogued images share one predicted image). 18 families have more predicted than
  catalogued images (`summary_wls.json`).
- The pass is own-constraint validation only and is labelled post-failure wherever quoted.

## Reading
- As written (`PLAN.md`) gate K fails. With the WLS source (`PLAN_1b.md`, deviation) the maps reproduce the gold
  list, so the next stage runs **with the WLS source**, as on Abell 2744 (D-034 addendum y).
- Next: a frozen gate C2 / extra-image plan with the CANUCS Bullet catalogue, including frame check F against
  that catalogue. The 18 families with extra predicted images are its targets; none is read here.
- Data: the two maps (219 / 205 MiB gzip, the > 200 MB reason in `PLAN.md`) were held in the session cache only and
  deleted at the end of the run.
