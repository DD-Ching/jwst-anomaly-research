# SMACS eazy `z_phot = -1` rows are compact: the frozen rule supports H-star (point-like sources), not a coverage flag

Plan: `PLAN.md` (pushed before the run; definitions added before any ranked-source number was seen). Run
`20261107T121058Z-86f3f464` of `configs/smacs_jwst3_gaia_only.yaml` at main `84e9096` (2938 ranked, as #764).
Measurement: `scripts/smacs_zphot_nofit.py <run> --zout <DJA v7.4 zout> --phot <DJA v7.4 phot> --out .` →
`summary.json` (both files sha256-checked). Contact sheet of the run looked at.

## Numbers (derived; ASSUMPTION thresholds from the plan)
- 2340 / 2938 ranked sources match the zout at 0.2″: 156 no-fit (`z_phot = -1`), 2184 fitted.
- **Coverage (H-cov): not supported.** `nusefilt` of no-fit rows: median 7, 5th percentile 6; fitted rows: median 8,
  5th percentile 6. Share of no-fit rows below the fitted 5th percentile: 0.0 (rule needs ≥ 0.8). In the ranked
  stratum every no-fit row has ≥ 6 usable filters.
- **Size (H-star): supported by the frozen rule.** DJA `flux_radius` (detection image, pixels): no-fit median 1.972
  vs fitted 10th percentile 3.105 (fitted median 5.441). The 133 matched star-stratum sources have median 2.222 and
  90th percentile 5.389 (the "compact" cut). Compact: 0.9359 of no-fit vs 0.4895 of fitted rows.
- Top-20 no-fit rows: ranks 3 (`f200w_1395`, 4.483 px), 14 (`f200w_1816`, 2.276 px), 19 (`f200w_384`, 1.906 px) are
  compact; rank 13 (`f200w_3884`, 11.734 px) is not. 3 ≥ 2 → `decide` returns "H-star supported".
- **Post hoc (deviation: added after the rule ran), size at fixed magnitude**, since faint galaxies are also small.
  Median `flux_radius` per `ref_mag` bin, no-fit vs fitted: 24–26 2.014 vs 5.917 (58 / 155 rows); 26–27 1.977 vs
  5.987 (44 / 212); 27–28 1.885 vs 5.361 (36 / 402). No-fit rows are about a third of the fitted size in every
  populated bin, close to the star median, so magnitude alone does not explain the result.

## Reading
- The eazy no-fit flag, among SMACS ranked galaxy-stratum rows, marks point-like sources at 24–28 mag that the
  Gaia-only star split does not remove. Which population they belong to is not tested here: faint Milky Way stars
  below the Gaia limit, compact sources around the cluster, and compact high-z galaxies are ordinary candidates.
- Why eazy did not fit them is not established (the flag is not a filter-count cut here); a stellar-template or
  point-source mask upstream in DJA is a hypothesis, not checked.
- Of the 4 top-20 no-fit rows, 2 (`f200w_1395`, `f200w_3884`) were already screened as spiked stellar-colour sources
  (D-019); `f200w_1816` and `f200w_384` look point-like in the contact sheet (#14, #19). Neither is a lens-model
  residual candidate on this evidence.

## What this changes
- Proposed, not adopted: a point-source screen for SMACS top-k vetting (DJA `flux_radius` against the stellar locus,
  validated on Gaia stars and catalogued multiple images first). TASKS Next (M2).
- D-022 addendum 5 records the result; the addendum 2 control rule is unchanged.

## Limitations
- One field, one eazy run; 2938 − 2340 ranked rows lack a zout match within 0.2″.
- `flux_radius` is measured on DJA's detection image, not per band; saturated stars inflate the star 90th percentile,
  so the "compact" cut is loose (0.4895 of fitted rows pass it).
- Data: zout member and phot catalog streamed to scratch and deleted after the run (Cloud disk decision).
