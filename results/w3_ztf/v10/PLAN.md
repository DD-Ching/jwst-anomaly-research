# W3 ZTF v10 — official psfcat in the W3 chain, calibration re-run, then (only if it passes) one fixed search

Plan fixed 2026-10-19 before any psfcat flux of this round is read. Owner authorization: the round after #355
(psfcat path, calibration re-run, and the 539-group search only if every condition below passes; 2 h in total).
No new photometry code, no threshold change, no new source. Labels: `observed` for psfcat, TAP and DR24;
`derived` for counts; `simulated` for injections.

## 1. Data and sets (unchanged from v6)
- Area Npix 7134, DR24 light-curve partition re-read (the session cache had been deleted; same 12,209 rows,
  `../v6/fetch_info.json`), TAP exposure metadata re-queried (same 21,632 rows).
- Sets from `../v6/sample.json`: calibration group 533, search group 539 (`holdout`), control group = the v6
  comparison stars (no sample member). `context()` logs every load to `source_usage.jsonl` in this folder;
  search-member rows must be 0 until section 6 starts.

## 2. Fill set (fixed by `fillset_script.py.txt` → `fillset.json`, before any psfcat is read)
- **Rule (product gap, from metadata and DR24 coverage only):** for every field / quadrant / filter (zg, zr) that
  holds control stars, an official exposure is filled if infobits = 0, at least one control star is in play inside
  its footprint (10″ margin), and DR24 has **no** detection (any catflags) of any control star of that quadrant.
  It does not use injections, event windows or which case failed.
- **Result:** 724 exposures in 33 field / quadrant / filter combinations, MJD 58243–60839 (mostly 2022–2023 and
  2025; per year in `fillset.json`). Exposures with partial DR24 coverage are not filled.
- **Overlap exposures:** per combination with a fill, the 2 infobits-0 exposures with the most clean DR24 control
  detections (ties: earliest), 66 in total.
- **Cost:** 790 distinct psfcat files, about 173 MB at the v9 mean of 214 kB; each downloaded once, cached in the
  scratchpad, used for every star of its quadrant, deleted after the round. Total downloads ≤ 1.5 GB with the
  387 MB partition.

## 3. Merging rules
- **Magnitude:** psfcat `mag` + header MAGZP (v9-verified definition), error `sigmag`; no colour term (as v9; DR24
  assumes zero colour). Matching by sky position to the DR24 object position: nearest source ≤ 1.5″; a second
  source ≤ 3″ = ambiguous.
- **Per combination validation (v9 criteria, fixed):** pairs = control star × overlap exposure with a DR24
  catflags-0 detection matched to that exposure and a unique psfcat match with flags 0. Pass if n ≥ 8,
  |median Δ| ≤ 0.05 mag and robust scatter ≤ max(0.03, 2 × median combined error). A combination that fails gets no
  fill and is reported; no product enters without passing.
- **One entry per star and exposure; DR24 is never overwritten.** psfcat is used only for fill exposures, where DR24
  has no detection of the star. Exposures with infobits ≠ 0 stay excluded; a DR24 flagged detection keeps the
  exposure unknown for the group (unchanged rule).
- **psfcat states, per star in a fill exposure:** unique match flags 0 → detection (clean or dim by the frozen
  DIM_MAG rule against the star's median of clean DR24 + psfcat magnitudes); flags ≠ 0 or ambiguous → flagged
  (unknown, as DR24 flagged); **missing → unknown for that star (dropped), never a non-detection, zero flux or
  limit.** Outside the footprint or inside the 10″ margin → not counted (unchanged).
- **Comparison counts in fill exposures** use the same rules for the control stars: in play = in the footprint
  (10″) within the star's span; clean = clean detection from either source. A missing control row is "not a clean
  detection" in the count (no flux is assigned), so the frozen GOOD_FRAC exposure rule is unchanged.
- Every added observation is stored with expid, product URL, MAGZP, flags, separation and reason
  (`merge_log.json`, compact per combination; the full per-observation list stays in the scratchpad cache, not
  tracked, because of its size).

## 4. Calibration re-run with the same events
- **Same events:** the v6 draws (carrier key, t0, u0 for the 10 calibrated cells, seeds 1100 + cell index; control
  parameters, seeds 1200 + k) are replayed on the DR24-only series exactly as `stage_inject` drew them, and must
  reproduce the stored v6 flagged / recovered counts. The same (carrier, t0, u0, model) are then injected into the
  merged series. Carriers that newly qualify are not added; event numbers are kept.
- **Noise (fixed now):** the extra noise of an injected epoch is drawn from a generator seeded by (seed of the cell,
  event index, expid), so a shared exposure gets the same draw before and after the fill and the added epochs do not
  shift other draws. To separate data from noise effects, the DR24-only "before" is re-run with the same scheme and
  reported next to the stored v6 counts.
- **Same treatment for added points:** shifted by the model, noise from the exposure's maglimit, an epoch fainter
  than its maglimit becomes a non-detection for a DR24 exposure (unchanged) and **unknown (dropped) for a psfcat
  exposure**, the same as a real missing psfcat row.
- **Common-mode control:** the same 5 exposures (by expid) as the v6 draw get the 3-mag dimming and the 0.5
  comparison dim fraction.
- **Reported per cell:** v6 stored, before (new noise scheme) and after; per event: the bracket (≥ MIN_BRACKET
  non-dim epochs before and after the umbra), the first stage that stops it (sampling, pre-screen requirement or
  vetting test) before and after; and the count still failing for lack of observations (the official good
  exposures inside / around the umbra are fewer than the frozen requirement) versus other rules.
- **Burden:** the real calibration light curves (merged) are screened with the frozen rules.

## 5. Pass conditions (all required; computed by code)
- Each of the 10 cells recovers ≥ 0.5 (≥ 100 of 200) after the fill.
- Each control class (eclipse, fade, PSPL, common mode) passes vetting in ≤ 10 % of injections.
- The real calibration light curves leave ≤ 5 survivors after vetting.
- Data quality: every filled combination passed its overlap test; 0 search-member rows loaded during calibration.
- The tests for the new code pass and the review leaves no finding that changes these numbers.
- Time: at least 35 min left before the stop (search, image check ≤ 64 exposures, records and review).
- If any condition fails, the reason is recorded and the round stops. No threshold is changed and no sample is
  redrawn. If only time is short, the search is recorded as **not executed**, not as failed or empty.

## 6. Search (only if section 5 passes)
- One search of the 539 search groups with the merged data (psfcat measurements of search members are extracted
  from the cached files only at this step), no added sources, no rule changes.
- Image check as v6 PLAN section 6 (`../v6/image_script.py.txt` configured for v10): ≤ 5 flags, ≤ 12 exposures
  each, ≤ 64 in total and ≤ 50 MB; order: runs passing vetting first (by key), then flagged groups by the number of
  non-detections in their longest run; 2 exposures before, up to 8 inside (evenly spaced), 2 after. Classes: real
  sky change, data problem, ordinary explanation, undetermined. A missing record is never zero flux.

## 7. Scope of any result
- This is a method correction re-checked on the already-used calibration sample, not an independent validation.
  #339 (full W3 validation failed) stands. Efficiencies are catalogue-stage for the 10 cells only, not image-level
  completeness for full disappearance.
