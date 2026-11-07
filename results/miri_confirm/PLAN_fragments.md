# Plan: deblend-fragment check for newly admitted MIRI F770W singles (2026-11-07)

Written before measuring (TASKS Then 7; follows `PLAN_rerank.md` and RESULT.md § Re-rank).

**Question.** The five single-band F770W sources that entered the SMACS MIRI top 20 after D-014 addendum 3
(`f770w_221` #2, `_239` #4, `_176` #13, `_118` #16, `_240` #18) lead with high CI_50_30 / CI_70_50. Is that score a
deblending effect (a child segment cut out of a bright neighbour's footprint, whose light gradient sets the
concentration), or a property of the source?

**Hypotheses.** H0 (ordinary, processing): the source's segment is a deblended child adjacent to a much brighter
segment, and adjacency alone produces the high CI. H1: the CI excess is intrinsic (to be vetted).

**Data (observed).** MAST L3 `jw02736-o002_t001_miri_f770w_cat.ecsv` and `_segm.fits` (S3 `stpubdata`, the same
products the pipeline run used; CAL_VER checked), the run's `quality` table (gate pass / reason) reconstructed from
the D-014 addendum 3 lists in `summary.json` / `summary_bands.json`.

**Measures (derived).**
1. Adjacency: a source is *bright-adjacent* when its segment shares a 4-connected pixel boundary with a segment
   whose F770W `isophotal_flux` is at least 5 × its own. ASSUMPTION: the 5× factor; 4-connectivity.
2. Population test of H0's mechanism, over all F770W catalog sources with finite CI and S/N ≥ 5: CI_50_30 of
   bright-adjacent vs other sources, one-sided Mann–Whitney; repeated in the source's flux tertiles.
3. Base rate of bright-adjacency among (a) the 113 baseline gate-passing sources, (b) the 55 newly admitted,
   (c) the five new top-20 entries.

**Frozen rule.**
- Mechanism supported when measure 2 gives p < 0.01 and the bright-adjacent median CI_50_30 is higher, in all flux
  tertiles that hold ≥ 10 sources per group.
- A new top-20 single is *explained by deblending* (ordinary processing; not an anomaly; no vetting) when it is
  bright-adjacent and the mechanism is supported. Otherwise it goes to `/vet-candidate`.
- Plan change trigger: if more than half of the 55 newly admitted are bright-adjacent while fewer than a quarter
  of the 113 baseline sources are, propose (not apply) a bright-adjacent flag for single-band confirmation
  (needs a new plan; post hoc here).

**What would change the plan.** Mechanism not supported → the CI excess is not explained by adjacency, and all
five go to vetting.
