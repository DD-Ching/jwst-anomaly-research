# W3 ZTF v11 — choose the area by observing conditions first, then validate and (only if passed) search

Plan fixed 2026-10-19 before any metadata of this round is read. Owner authorization: the round after #359.
The W3 model, the code (`scripts/w3_ztf.py`, DR24 light curves; psfcat fill not used unless section 2 shows a DR24
product gap) and every frozen threshold are unchanged. No abundance limit is derived.

## 1. Candidate areas (fixed now; ≤ 3)
- The three DR24 light-curve partitions at Norder 5 under the order-4 objects pixel 959, the area with the densest
  catalogue sampling seen so far (v7: every month to 2024-06 has detections in Npix 3838):
  - **Npix 3837:** never read;
  - **Npix 3838:** v4 calibration sample of 1,000 groups (used);
  - **Npix 3839:** v1–v3 sample of 1,000 groups (used; its holdout was searched in v5).
- **Eligible sources:** groups from the v4 selection rule (`select()`: r band, refmag 14.5–17.5, |sharp| ≤ 0.3,
  ngoodobs ≥ 300, non-transient, isolated, not ambiguous) that are not within 2 × 1.5″ of any v1–v3 or v4 sample
  group. Previously used or searched sources are never reused.
- Npix 7134 is not a candidate: all its groups were used in v6/v10.

## 2. Observing-condition check (official metadata only; no light curve read)
- **Positions:** per area, 200 eligible groups drawn with seed 11 (all if fewer). Each group's quadrants are those
  of its r-band objectids (all fields).
- **Exposures:** IRSA TAP `ztf.ztf_current_meta_sci`, r band, those quadrants. An exposure counts for a group if
  infobits = 0 and the group position is ≥ 10″ inside its footprint (`in_quad`). Exposures shared by two
  objectids of the group count once (expid); times are the heliocentric mid-times (`exposure_times`).
- **Period:** each group's full span of such exposures (first to last), seasonal gaps included. Nothing is removed
  afterwards.
- **Cells:** the 10 calibrated cells (ρ, blend, t_E): (0.01, 1, 10/30/100), (0.1, 1, 10/30), (0.01, 0.7,
  10/30/100), (0.1, 0.7, 10/30). W3 model `w3_dmag`, n = 1, ε < 0.
- **Event draws:** per cell 2,000 draws, seed 2100 + cell index: group uniform over the 200, t0 uniform over its
  span, u0 uniform 0–2 (as `stage_inject`).
- **Necessary-sampling indicator:** the noise-free model is evaluated at the counted exposures. A draw counts if
  `flagged_runs` (the frozen pre-screen: ≥ 3 consecutive epochs with Δm ≥ 1 on ≥ 3 nights, ≥ 20 non-dim epochs
  before and after) finds a run overlapping the umbra t0 ± t_E √(4 − u0²). This is an upper bound on recovery
  from sampling alone, not a recovery rate. It ignores noise, the good-exposure rule, DR24 product inclusion and
  vetting.
- **Exposures without a usable product, listed separately per area:**
  - r-band exposures with infobits ≠ 0 in those quadrants (excluded by the frozen rule);
  - psfcat existence for 20 counted exposures per area (fixed by seed 12; HTTP HEAD). DR24 light-curve inclusion
    cannot be known before reading light curves; it is checked after the read (section 4).

## 3. Choice and stop rule (computed by code)
- Score of an area: its lowest necessary-sampling fraction over the 10 cells. The area with the highest score is
  chosen; ties go to the area with more eligible groups.
- **If no area scores ≥ 0.5, the round stops:** no light curve is downloaded and no threshold is changed. The report
  names the missing sampling (umbra epochs, nights, or bracket) per cell.

## 4. If an area is chosen
- **Sample:** `stage_sample` with VERSION 11 in that partition: ≤ 1,000 eligible groups (seed 7, as before),
  excluding v1–v3 and v4 sample positions. Split: search if crc32("v11:" + key) is even, else calibration. Control
  group: comparison stars of the area without any sample member. All objectids of a source (fields, quadrants,
  filters) within 1.5″ form one group.
- **Data:** DR24 light-curve partition (7 columns) and TAP metadata, as in v6. `context()` logs every load; the
  search group gives 0 rows until the validation has passed.
- **Same-area caveat:** stars of one area share fields, exposures and systematics with the v1–v5 samples. This is
  not an independent-area validation.
- **Realised sampling:** after the read, the same indicator on the calibration carriers' actual series (frozen
  `series()`), reported next to section 2.

## 5. Validation (unchanged from v6 section 5; computed by code)
- `stage_inject --cells calibrated --n 200`: each of the 10 cells recovers ≥ 0.5; each control class passes
  vetting in ≤ 10 %; the real calibration light curves leave ≤ 5 survivors.
- Any failure: record and stop. No second area, no threshold change, and no partial pass from the cells that
  reached 0.5.

## 6. Search (only if section 5 passes, the review leaves nothing that changes the verdict, and ≥ 35 min remain)
- One search of the search group (`search --half holdout`), nothing added.
- Image check: `../v6/image_script.py.txt` configured for v11 (order, ≤ 5 flags, ≤ 12 exposures each, ≤ 64 in
  total, ≤ 50 MB; 2 before, up to 8 inside, 2 after). Classes: data problem, real sky change, ordinary
  explanation, undetermined. A missing record is never zero flux; a real dimming is not called a W3 or wormhole
  candidate.

## 7. Budget
- Downloads ≤ 1.5 GB (light-curve partition ≈ 0.85 GB, objects ≈ 0.09 GB, TAP); cutouts ≤ 50 MB. Each product is
  downloaded once and cached; the cache is deleted after the round. 2 h in total, with ≥ 20 min for review and
  records.
