# W3 ZTF v15 — Npix 3837 passes and is searched (4,778 groups, 0 vetted survivors); Npix 183 fails validation and is not searched

- **Plan:** `PLAN.md`, committed before any read. **Samples:** committed before any light curve.
- **Labels:** DR24 and TAP `observed`; counts `derived`; injections `simulated`.
- **Scope:**
  - the frozen v13 W3 rule (`--edge time`), the 10 calibrated cells, catalogue stage only;
  - no limit derived, and nothing here bears on wormholes in general.
  - Npix 3837 shares ZTF fields 826 and 1857 with v1–v5: **new sources, but partly shared exposures**, so it is not an
    independent-area validation.

## 1. Choice and samples (`sample_info.json`)
- **Npix 8170 was not used.** Only ≈ 738 unused eligible groups remain there. After a new 500-group calibration set
  that leaves ≤ 238 search groups for a 2.2 GB re-read.
- **Chosen by metadata and cost:** 3837 and 183, both qualified in #368.
  - **3837:** 5,278 eligible groups, split 500 calibration / 4,778 search. Refmag r 5/50/95 %: calibration 14.82 /
    16.46 / 17.43; search 14.80 / 16.43 / 17.42.
  - **183:** 2,303 eligible groups, split 500 / 1,803.
- **Dedupe** against every earlier calibration, development and search group (v1–v14): 0 shared objectids and none
  within 3″. Earlier reference stars lie in other, disjoint partitions.
- **Deviation:** `configure()` was called twice in one process, so the stored 183 `lc_partition` text and its cache
  folder carried 3837 values. Both were corrected, and the code was fixed with a unit test. The selection itself used
  partition 183: all 2,303 groups lie in it.

## 2. Validation (`3837/decision.json`, `183/decision.json`; frozen rule, new seeds, fade × 1,000)

| cell (ρ, blend, t_E) | 3837 | 183 |
|---|---|---|
| 0.01, 1, 10 | 139 | 100 |
| 0.01, 1, 30 | 166 | 136 |
| 0.01, 1, 100 | 159 | 139 |
| 0.1, 1, 10 | 130 | **96** |
| 0.1, 1, 30 | 142 | 149 |
| 0.01, 0.7, 10 | 135 | 114 |
| 0.01, 0.7, 30 | 143 | 132 |
| 0.01, 0.7, 100 | 124 | 129 |
| 0.1, 0.7, 10 | 130 | **91** |
| 0.1, 0.7, 30 | 141 | 127 |
| **fade passing vetting** | 66 / 1,000 (6.6 %, Wilson 5.2–8.3 %) | **154 / 1,000 (15.4 %, 13.3–17.8 %)** |
| eclipse, PSPL, common | 0, 0, 0 of 200 | 1, 0, 0 of 200 |
| real calibration survivors | 0 (30 flagged runs in 1 group) | 0 (no flags) |
| **verdict** | **passed** | **failed** (2 cells < 0.5; fade > 10 %) |

- **183 stopped.** Its 1,803 search groups were never loaded (`183/source_usage.jsonl`). Nothing was retried.
- **Why 183 failed** (post hoc trace on its calibration group with the same draws; `183/trace_posthoc.json`; the
  stored counts are reproduced):
  - In the two failing t_E 10 d cells, 50 and 51 injections have no dim epoch run inside the umbra, plus 14 and 20 at
    MIN_RUN, so the loss is mainly sampling. The metadata check had given 183 the lowest worst-cell margin (0.547).
  - The common-mode test removes 8–24 flagged injections per cell.
  - The fade-control excess (15.4 %) is not diagnosed further.
- **Fade control across validated areas, same rule:**
  - 2382: 17 / 200 (8.5 %);
  - 8170: 50 / 1,000 (5.0 %);
  - 3837: 66 / 1,000 (6.6 %);
  - 183: 154 / 1,000 (15.4 %, failed).
  - The rate depends on the area.

## 3. Search of Npix 3837 (`3837/search_holdout_v15_time.json`)
- **4,778 groups, 9,519 series:** 45 flagged runs in 4 groups, **0 pass vetting**. One group (826214400018444) holds
  37 of the runs.
- Every run fails the one-off test and/or the achromatic test (counts in the file). Since there is no final survivor,
  no image check was made (0 cutouts).

## 4. #380 pending items
- 656 rebuilt NDs and the undetermined 6th group stay pending. They come from #370 flagged runs that failed vetting
  (diagnostic sample; 0 survivors there). They are not re-checked here and not claimed excluded.

## 5. Cost and rate (`cost.json`)
- **Downloads ≤ 1,703.5 MB:**
  - LC 3837: 1,233.2 MB;
  - LC 183: 270.2 MB;
  - objects ≤ 186.2 MB;
  - TAP 13.9 MB;
  - cutouts 0 MB.
- **Per 1,000 searched groups in 3837** (session wall-clock, not stored timings):
  - about 1.9 min of compute;
  - about 280 MB downloaded;
  - 0 image checks.

## Deviations and notes
- The `configure()` bug and its fix (section 1).
- `183/trace_posthoc.json` is post hoc. It comes from the v7 trace script, run through `trace183_wrapper.py.txt` with
  the 8000 seed offset and the time edge rule.
- The 3837 sample is stored gzipped (`3837/sample.json.gz`; rebuild with `gunzip -k`), because the JSON is 1.4 MB,
  over the 1 MB rule.
