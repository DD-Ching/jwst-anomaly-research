# SMACS galaxy ranking scored within DJA z_phot strata: members are under-represented in the top 20; no rule change

Plan: `PLAN.md` (pushed before the run). Run `20261107T111137Z-3978d5de` of `configs/smacs_jwst3_gaia_only.yaml`
at main `da1d1be`. Re-score: `scripts/smacs_zphot_stratify.py <run> --zout <SMACS DJA v7.4 zout> --out .`
(the script checks the zout sha256 against SOURCES.md and stores it in `summary.json`) → `summary.json`. Contact sheet of the run looked at.

## Tool check (observed)
- The baseline re-score reproduces all 2938 pipeline galaxy-stratum ranks exactly (the script stops otherwise).
- CDS XMatch still refuses `vizier:I/355/gaiadr3`; the Gaia `field` fallback carried the star split (as #763).

## Strata (derived; ASSUMPTION ranges as D-022 addendum 2, 0.2″ zout match)
| Stratum | Ranked | Overall top 20 |
|---|---|---|
| member (z_phot 0.25–0.55) | 614 (0.209) | 2 (0.10) |
| background (z_phot ≥ 0.8) | 1274 | 12 |
| other | 1050 | 6 |

## Frozen rule (`decide`, unit-tested)
- Members in the top 20: 2 (< 4) and share 0.10 (< 2 × 0.209). **Verdict: no rule change.** H0 (bright members
  ranked anomalous as a minority population) is not supported; the overall top 20 is mostly background (H1).
- Background scored alone: its top 9 (its share of 20) holds 7 of the overall top 20 (overall ranks 1, 2, 6, 7, 9,
  11, 18, 40, 80). Members scored alone: their top 4 have overall ranks 38, 58, 97, 213 (none in the top 20).
- **"What changes the plan" trigger, as written (`image_trigger`): met.** The 14 ICLv2 multiple images that match a
  ranked background source (0.3″; 25 / 60 match any ranked source) have median rank 1420.5 overall and 611.5 within the
  background stratum, a raw ratio of 2.323 (> 2). **Deviation (post hoc reading, not frozen):** the background list is
  shorter by 2938 / 1274 = 2.306, so the raw ratio is the list-length ratio; as percentiles the images sit at 0.4835
  overall and 0.48 within the background. The plan's consequence (wire the option before new top-k vetting) is
  therefore not taken; this deviation is stated in the PR for the owner. The rule should have been written in
  percentiles.

## Post hoc (deviation: added after the first pass, not a rule)
- The "other" stratum splits into 156 eazy `z_phot = -1` rows (no fit), 598 without a zout match and 296 between the
  ranges. `z_phot = -1` rows are 4 / 20 of the top 20 (0.20) against 156 / 2938 of the ranked stratum: 2 of them
  (`f200w_1395`, `f200w_3884`) are already screened as spiked stellar-colour sources (D-019), and `f200w_1816`,
  `f200w_384` look compact and point-like in the contact sheet (visual, not stored). Hypothesis (untested): an eazy no-fit flag marks
  uncatalogued stars or compact non-galaxy sources missed by the star split.

## What this changes
- TASKS Next (M2) "DJA photo-z as a ranking feature or stratifier": stratifier not adopted for SMACS rankings
  (D-022 addendum 4). The D-022 addendum 2 control rule (science-vs-control comparisons use z_phot ≥ 0.8 in `ref_mag`
  24–28) is unchanged; it concerns feature distributions, not the ranking.

## Limitations
- One field, one score (the pipeline's robust-z / isolation forest / LOF ensemble), z_phot from one eazy run.
- The top 20 is taken before the D-019–D-021 screen (5 of the 20 are screened; `summary.json` lists them).
- Data: pipeline catalogs and the zout member streamed to scratch, deleted after the run (Cloud disk decision).
