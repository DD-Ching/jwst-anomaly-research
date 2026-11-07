# Plan: size-aware MIRI band association (2026-11-07)

Trigger: `docs/candidates/smacs_miri_f770w_top20_singles.md`. Of the five new top-20 F770W "singles", 4 have
counterparts in other MIRI bands 0.18–0.66″ away, outside the fixed 0.1″ merge radius (TASKS Then 7).

**Hypothesis (ordinary).** The fixed 0.1″ radius is smaller than the band-to-band centroid scatter of resolved MIRI
sources (FWHM 0.27–0.59″; isophotal centroids of extended sources shift between bands), so real multi-band
sources are split into single-band rows. Those rows then fail `require_multiband` or rank as unusual single-band
detections.

**Method (frozen before measuring).**
- Per pair, accept if `sep ≤ max(r0, α · r_eq)`, where r0 = 0.1″ (the current radius), `r_eq = sqrt(A_min / π) ·
  pixel scale` and `A_min` is the smaller isophotal area of the two detections (the row's anchor and the candidate).
  Matching stays one-to-one and greedy by separation.
- α grid {0, 0.25, 0.5, 0.75, 1.0} on the SMACS 0723 MIRI L3 catalogs (program 2736, jwst 3.0.0).
- Chance rate: re-merge with each non-reference band catalog shifted by (±8″, ±8″), 4 shifts. For each α,
  chance fraction = (mean shifted matches beyond those at α = 0) / (real matches beyond α = 0).
- **Rule:** adopt the largest α whose chance fraction is ≤ 0.10 (ASSUMPTION). If none of α > 0 passes, keep 0.1″.
- The five vetted sources are a check, not a tuning input: expected merged at the adopted α are `_221`, `_239`,
  `_118`, `_240`; `_176` stays single (nearest counterpart ≥ 5″).
- If an α is adopted: re-run `configs/smacs_jwst3_miri_confirm.yaml` with it and compare the gate count and the
  top 20 with run `20261107T022631Z-036f249e`; view the contact sheet.

**What would change the plan.** If fewer than 3 / 4 expected merges happen at the adopted α, the size-aware radius
is not the explanation; record it and look at segmentation-overlap association instead. If the chance fraction
exceeds 0.10 already at α = 0.25, keep the fixed radius and record the limit.
