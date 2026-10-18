# W3 ZTF v8 — data-entry check and minimal fill trial (Npix 7134, r band)

- **Scope:** checks only. No threshold changed, no new object searched, the search group not read, no injection
  re-run. No new photometry pipeline: the #338/#339 cutout method was reused.
- **Conclusions** are limited to this area (Npix 7134), DR24 and the frozen W3 rules.
- **Labels:** `observed` for HATS, IRSA LC API, TAP and IBE outputs; `derived` for counts and photometry.

## 1. Where the gap is: HATS vs another official entry vs the exposure log (`entry.json`)
- **Fixed sample:** the first 5 calibration groups (by key) whose r-band objectid lies in HATS row group 2, giving
  13 r-band objectids across fields 671, 672, 1669 and 1714.
- **HATS and the IRSA light-curve API agree exactly.** For every one of the 13 objectids, the number of detections in
  the gap (2022-10 to 2023-10) is the same in the HATS DR24 row group and in the API (e.g. 10/10, 22/22, 28/28, 1/1).
  The gap is not a reading or subset problem of this round.
- **The release lacks exposures in part of the area.** For the field 672 objectid 672208200001404, IRSA lists 110
  official exposures in the gap that cover it (64 with infobits 0), but DR24 has 10 detections. Four of the field
  1669 objectids have 6 official exposures and 1 detection each. The fifth, 1669214200003701, has 0 official
  exposures and 1 detection.
- Month labels in `entry.json` can be off by one at month edges (approximate binning). The gap totals use exact MJD
  limits.
- **Field 671 objectids have no such gap:** official 9–28 exposures against 10–28 detections.
- **The DR24 release notes** (IRSA, `ztf_release_notes_dr24.pdf`, read 2026-10-18) put the gap inside the release
  spans: public to MJD 60970 and private to MJD 60513. But "the single-exposures included in DR24 … do not cover all
  survey fields from both grids" (advisory 10). The missing exposures are therefore consistent with documented
  release coverage, not with a reading error. The reason for each individual exposure is not established.
- **No newer release on IRSA's S3** (`ipac-irsa-ztf/ztf/enhanced/` lists only `dr24`). The API serves the same data
  as HATS for these objects.
- **The 2,600 unmatched 2023 detections**, 10 of them sampled (`entry.json`, `unmatched_lookup`; exposure times in
  `unmatched_exptime.json`):
  - **5 are 299–300 s private exposures** whose detection lies 0.00155–0.00156 d (about 134 s) from the matched time.
    The code assumed 30-s exposures (mid-time = start + 15 s), so it was off by 134.5 s, more than the 86-s match
    tolerance. This is an **implementation error**, fixed here (section 4).
  - **1 is a 30-s exposure 0.153 d away.**
  - **4 are 30-s exposures whose nearest in-footprint exposure is 1.1–290 d away.** Their expid is in the area's
    official list, but not among the star's in-footprint exposures (it may fail the 10″ margin; not verified).

## 2. Ready-made products and services
- **IRSA LC API:** same release, same gaps (above).
- **Newer release:** none found on IRSA's S3 listing. A newer release outside S3 was not searched (unknown).
- **ZTF forced-photometry service:** needs a registered account, so this option was stopped under the round's
  rules. It was not tried.
- **Science images (IBE):** cutouts were obtained for all 36 trial exposures (section 3).

## 3. Limited image trial (`image_trial.json`, selection `image_trial_selection.json`, fixed before download)
- **Selection:** 3 calibration stars with an r-band objectid in field 672. Per star, 8 gap exposures plus 4 control
  exposures with a DR24 detection, 36 cutouts in total. A gap exposure is one with no DR24 detection with catflags 0
  for that objectid, so a flagged detection may still exist there. All 36 cutouts were obtained and all apertures
  are finite.
- **Gap flux** (`trial_ratios.json`; uncorrected, because no field stars were available): the target is present in
  all 24 gap exposures. The aperture flux is 0.62–1.16 × the median of the star's own 4 controls, so no gap
  exposure falls ≥ 1 mag below its controls. The field-corrected values in `image_trial.json` are NaN for that
  reason.
- **Control check, failed on two counts:**
  - At the 12 control exposures, the 2″-aperture flux is 0.54–0.75 × the DR24 flux of the same exposure. The
    random-aperture errors (≈ 0.4–2.7 % of the control flux) are far smaller than that spread.
  - The field-star search, run once per star on its first cutout (3 searches), found **0 comparison stars**.
    Same-exposure zero-point and background checks were therefore not possible.
- **Seeing dependence:** the 9 gap exposures with seeing > 3″ have ratios to the controls of 0.62–0.88. The lowest,
  0.62, is at seeing 3.96″.
- **Reading:** the method shows the star is present at roughly its usual level. It does **not** reproduce catalogue
  flux to catalogue precision, and its errors are not realistic. Under the round's rule, the failure cases are kept
  (the 9 seeing > 3″ gap exposures) and no PSF solver is built.

## 4. Implementation fixes (`scripts/w3_ztf.py`, tests)
- `exposure_times` now uses the `exptime` column (half the exposure) when present, and `stage_exposures` queries
  it. Without the column it falls back to 30 s.
- `match_exposures` returned −1 for a quadrant with exactly one exposure; this is fixed.
- One regression test was added, `test_v8_exposure_mid_time_uses_exptime`, which fails on the old code.
- **Effect on stored results:**
  - The `exptime` fix changes none. Every non-30-s r-band exposure in both areas has infobits ≠ 0 (`exptime.json`:
    542 in v6, 320 in v4; 0 usable), and the frozen infobits rule already excluded them.
  - The single-exposure fix was not re-run. It adds at most one epoch per such quadrant, against a frozen minimum of
    2 × MIN_BRACKET + MIN_RUN epochs per series.

## 5. Cost of filling the #350 baseline gaps (`cost_estimate.json`)
- **Scope:** the 49 ρ 0.01, blend 1, t_E 100 d cases whose post-event baseline is short, although IRSA lists ≥ 20
  official good exposures after the umbra.
- **Size:** 528–533 distinct star × exposure measurements across 47 carriers.
- **Unit costs measured in the trial:**
  - 1.69 s and 69.6 kB per cutout with photometry, so about 15 min and 37 MB;
  - 18.8 s of queries per star, so about 15 min.
- **Scope of that estimate:** one cell only. The other 9 cells were not estimated.
- **What it would not settle:** meeting the 20-epoch bracket does not mean the injections would be recovered. The
  measurements would also come from a method that failed its control check.
- **Reusable:** cutouts, TAP and API responses can be cached per star and exposure (not stored here).

## Downloads
- Metadata and light curves: about 69 MB (≤ 100 MB; byte counts in `entry.json`, `image_trial.json`,
  `exptime.json`):
  - HATS row group 2, objectid column 1.00 MB and light-curve column 65.85 MB;
  - API 0.99 + 0.33 MB;
  - TAP 0.50 + 0.31 MB plus aggregate queries.
- Cutouts: 2.51 MB (≤ 100 MB).
- Release notes PDF: 6.45 MB, documentation.

## Observation, not acted on (a frozen rule)
- The release notes call exposures with INFOBITS < 33554432 "likely … usable" (§12b). The frozen W3 rule keeps only
  infobits = 0.
- Example, field 671 objectids: 9–28 official exposures in the gap, of which only 2–9 have infobits 0, while DR24
  holds 10–28 detections (`entry.json`).
- Relaxing the rule would be a threshold change, so it was not made.

## Decision
- **The route stays paused** for this area and DR24.
- For field 672 in Npix 7134 (1 sampled objectid plus 3 trial stars), the gap exposures are absent from DR24 in both
  HATS and the API, so it is not a reading problem. The reason for each exposure is not established.
- The available images confirm the star is present. The reused cutout method fails its control: catalogue flux is
  not reproduced, no field stars are found, and the errors are not realistic. Filling about 530 star × exposures for
  one cell (about 30 min, about 37 MB, plus about 4 control exposures per star, about 188 more) would therefore rest
  on an unvalidated measurement.
- **Reopen if any of these is verified** (each check is cheap and uses this round's scripts):
  1. A newer ZTF release whose light curves contain the field 672 gap exposures. Check: `entry_script.py.txt` on the
     same 13 objectids shows release detections ≈ official infobits-0 exposures in the gap.
  2. Image photometry that passes the control. The cutout method must reproduce DR24 flux at control exposures to
     within the DR24 errors, with realistic errors, for example using catalogue stars as same-exposure comparisons.
     Running that check needs the owner's go-ahead, because it means building more than this round allowed.
  3. The owner decides to revisit the frozen infobits rule. That is a scientific threshold change, so it needs a new
     plan and a new validation.
