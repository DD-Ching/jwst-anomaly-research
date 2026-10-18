# W3 ZTF v7 — why recovery dropped from Npix 3838 to Npix 7134 (post hoc diagnosis)

- **Scope:** post hoc diagnosis only. Thresholds were not changed, the search group was not read, and no new area was
  added. It is not a new validation. #339 (failed) and #349 (failed; search not executed) stand.
- **Data:** the two partitions already used were re-read:
  - Npix 7134: 386.7 MB (`../v6/cost.json`);
  - Npix 3838: 7 light-curve columns of an 846.0 MB file (`../v4/cost.json`, upper bound) plus 5 object columns,
    79.9 MB (`../v5/cost.json`).
  - Total ≤ 1.31 GB, under the 2 GB cap.
- **Images:** none downloaded. 3 IBE metadata queries were made.
- **Search group:** `../v6/source_usage.jsonl` shows 0 search-member rows for every v7 script.

## 1. Stage trace reproduces the stored results (`trace_v4.json`, `trace_v6.json`)
- `trace_script.py.txt` replays `stage_inject`'s draws: same carriers, seeds 1100 + cell index, same calls.
- For each of the 200 injections per cell it records:
  - pre-injection sampling from the carrier's official good exposures, with no injected detection used;
  - whether the run was flagged;
  - the first pre-screen block and any further blocks;
  - the first failing vetting test;
  - whether the injection was recovered.
- The flagged and recovered counts equal the stored ones in all 10 cells, in both areas: `../v4/injections_v4.json`
  and `../v6/injections_v6.json`. Stage counts add up to 200 per cell.
- **Sampling-eligible group** (fixed from the frozen pre-screen requirements before the trace was read): ≥ MIN_RUN
  epochs on ≥ MIN_NIGHTS nights inside the umbra, and ≥ MIN_BRACKET epochs before and after it, all before
  injection. Recovery is reported for both groups and for all 200. It counts all carrier epochs, so it is a
  sampling indicator, not a strict necessary condition: for example, at blend 0.7, t_E 100 d in the old area 173 are
  flagged against 172 eligible.

## 2. How much of the gap is explained (`summary.json`)

| cell (ρ, blend, t_E) | old recovered | new recovered | old / new recovered if eligible | gap from sampling eligibility | gap from recovery given eligibility |
|---|---|---|---|---|---|
| 0.01, 1, 10 | 147 | 61 | 147/148, 61/65 | 77.9 | 8.1 |
| 0.01, 1, 30 | 170 | 80 | 170/175, 80/89 | 77.3 | 12.7 |
| 0.01, 1, 100 | 173 | 90 | 173/177, 90/93 | 81.3 | 1.7 |
| 0.1, 1, 10 | 124 | 80 | 124/138, 80/90 | 42.7 | 1.3 |
| 0.1, 1, 30 | 142 | 76 | 142/180, 76/92 | 72.7 | −6.7 |
| 0.01, 0.7, 10 | 142 | 80 | 142/146, 80/85 | 57.4 | 4.6 |
| 0.01, 0.7, 30 | 165 | 76 | 165/175, 76/91 | 70.2 | 18.8 |
| 0.01, 0.7, 100 | 148 | 61 | 147/172, 61/83 | 65.4 | 21.6 |
| 0.1, 0.7, 10 | 128 | 51 | 128/147, 50/62 | 68.5 | 8.5 |
| 0.1, 0.7, 30 | 120 | 48 | 120/172, 48/93 | 40.8 | 31.2 |
| **total (2,000 per area)** | **1,459** | **703** | — | **654.2** | **101.8** |

- **Sampling part** = conditional rate × (old eligible − new eligible). The rest is the difference in recovery among
  sampling-eligible injections (pre-screen and vetting). The split depends on the rate used:
  - with the new area's conditional rate: 654.2 sampling and 101.8 conditional (table);
  - with the old area's conditional rate: 704.2 and 51.8 (`summary.json`, `*_old_rate`).
  - The sampling part therefore explains about 654–704 of the 756 gap.
- **Sampling:** in the new area most injections lack the epochs the frozen pre-screen requires (umbra epochs and
  nights, and baseline before and after).
- **Recovery given eligibility** differs by more than 0.1 in 3 of 10 cells, all at blend 0.7: 0.94 vs 0.84, 0.85 vs
  0.73 and 0.70 vs 0.52 (old vs new).
- **Pre-screen and vetting given eligibility, about 52–102 of 756.** Vetting first-fail counts in the new area,
  summed over the 10 cells, are stored per cell:
  - **common-mode test: 20 injections** in the new area, 0 in the old. The v6 comparison pool is the control group
    only (sample members excluded), which is a v6 method change. The link from that change to the common-mode
    failures is not tested here.
  - **sharp-edge test:** for example 31 at ρ 0.1, blend 0.7, t_E 30 d (33 in the old area).
  - **achromatic test:** 12 and 16 at blend 0.7, t_E 30 and 100 d (7 and 21 in the old area).

## 3. The 72 ρ 0.01, blend 1, t_E 100 d injections (≥ 3 umbra nights, not flagged) (`cases_tE100.json`)
- **All 72 stop at MIN_BRACKET**, the frozen requirement of ≥ 20 non-dim epochs before and after the run. They are
  not blocked by run length or by nights.
- **59 have fewer than 20 series epochs after the umbra.** For 49 of them, IRSA lists ≥ 20 official good exposures
  after it (infobits 0, star ≥ 10″ inside the footprint), but the series lacks them. The median umbra end is
  MJD 60153 (2023-07-28).
- **13 have fewer than 20 series epochs before the umbra.** For 12 of them, the official exposures are also fewer
  than 20.

## 4. Why the series lacks those exposures (`drops_v6.json`, `month_v6.json`, `match_v6.json`; old area `*_v4.json`)
- **Monthly catalogue coverage, new area, official r-band quadrant exposures with infobits 0:**
  - **2022-10 to 2023-10:** 148 of 607 official quadrant exposures (24 %) have any catalogue detection of a loaded
    star (e.g. 2023-03: 2 of 41; 2023-10: 1 of 68; better in 2022-11, 30 of 43). The frozen good-exposure rule
    (≥ 50 % of in-play comparison stars detected cleanly) then drops them: 20,960 carrier-member × exposure pairs in
    2023 (`drops_v6.json`).
  - **2024:** the area was barely observed: 73 r-band quadrant exposures (`match_v6.json`), or 670 carrier-member ×
    exposure pairs (`yearly_v6.json`).
- **Old area:** every month until 2024-06 has detections in all exposures (`month_v4.json`).
- **Matching:** detection-to-exposure matching is complete in every year except 2023 in the new area, where 2,600 of
  29,164 detections do not match an exposure (`match_v6.json`). Another 529 detections have a non-finite time; all
  of them are unmatched. Whether the 2,600 would fill some empty exposures was not tested. They are kept as an open
  minimal case: Npix 7134, r band, 2023.
- **Images:** the check picked the same exposure twice, so only 2 distinct empty exposures were checked. IRSA's image
  search lists a science image for 1 of them (`ibe_check.json`; 2023-08-09). At least some of the images exist; how
  many is not established.

## 5. Classification of the cause
- **Implementation error:** none found that changes the outcome. The trace reproduces every stored count. The 2,600
  unmatched 2023 detections remain an open small item.
- **Missing sampling, the main cause:**
  - **in the catalogue:** DR24 light curves for these fields lack most exposures from 2022-10 to 2023-10; the image
    check is limited to 1 of 2 distinct exposures;
  - **in the survey:** 2024 was barely observed, and in some cases the event falls before the first observations.
- **Pre-screen and vetting given enough sampling, a minor cause:** about 52–102 injections depending on the rate
  used (common-mode, sharp-edge and achromatic tests; counts above). No rule was changed.

## 6. Recommendation (the route is paused)
- Another batch or another area in the DR24 catalogue would add more data of the same kind. It would not restore the
  missing post-event baselines, so it is not proposed.
- **Reopen if either holds:**
  1. a later ZTF data release whose light curves contain the 2022-10 to 2023-10 exposures for such fields. Check
     with `month_script.py.txt` on calibration data before any search data are read;
  2. the owner authorizes image-level forced photometry from the IRSA science images (an image is listed for 1 of 2
     distinct sampled exposures) to fill the post-event baseline of calibration carriers. The scale and cost of that are not
     estimated here (unknown).
- In both cases the frozen thresholds stay. Recovery is quoted over all injections, with the sampling groups shown
  alongside.
