# W3 ZTF v13 — edge-duration rule passes validation; 2,000 groups searched; 0 vetted survivors; no image-confirmed dimming

- **Plan:** `PLAN.md`. The rule and the re-split were frozen in 7aae8c9, before the development comparison was looked
  at. The model-derived `edge_frac.json` was committed later (2818b97); its script was in 7aae8c9 and reproduces it
  exactly.
- **Labels:** model `model_prediction`; injections `simulated`; DR24, TAP and cutouts `observed`; counts `derived`.
- **Scope:**
  - one W3 vetting test was replaced; everything else is unchanged;
  - catalogue stage, the 10 calibrated cells, Npix 2382;
  - the validation set shares its area and exposures with the development set, so this is **not an independent-area
    validation**;
  - no abundance limit is derived.

## 1. Records checked first
- **#368 sample:** 3,000 groups (500 calibration and 2,500 search; 16,499 objectids). The "1,000" was only in the
  chat report.
- **The 80 cases:** confirmed. At ρ 0.1, blend 0.7, t_E 30 d, 80 of the 87 flagged-but-lost injections failed first on
  the sharp-edge test (`../v12/2382/trace_posthoc.json`). The cadence explanation was not tested.

## 2. The replacement rule (`edge_time()`, `--edge time`)
- **Edge extent in days:**
  - Before the run: from its first dim epoch back to the farthest significant intermediate detection.
  - After the run: from its last dim epoch forward to the farthest one.
  - Significant means 0.3 + 2σ ≤ Δm < 1, with σ from the epoch's maglimit; non-detections are never intermediate.
- **Pass** if both extents are ≤ 0.25 × the run span.
- **The 0.25** is the 99th percentile over u0 of the W3 model's own edge-to-span ratio at ρ 0.1, blend 0.7 (0.248;
  `edge_frac.json`, model only). The ratio is ≤ 0.008 at ρ 0.01.
- **No injection truth is used.**
- A unit test (`test_v13_edge_time_uses_time_not_count`) shows that dense edge epochs pass and a far one fails.

## 3. Development comparison (same events and noise; already-used v12 calibration group; `devcompare.json`)

| cell (ρ, blend, t_E) | old recovered | new recovered | old sharp-edge first failures | of these, recovered by the new rule |
|---|---|---|---|---|
| 0.01, 1, 10 | 150 | 150 | 0 | 0 |
| 0.01, 1, 30 | 173 | 176 | 4 | 3 |
| 0.01, 1, 100 | 150 | 157 | 9 | 7 |
| 0.1, 1, 10 | 130 | 142 | 15 | 12 |
| 0.1, 1, 30 | 122 | 174 | 57 | 52 |
| 0.01, 0.7, 10 | 138 | 139 | 1 | 1 |
| 0.01, 0.7, 30 | 155 | 156 | 1 | 1 |
| 0.01, 0.7, 100 | 136 | 145 | 12 | 9 |
| 0.1, 0.7, 10 | 117 | 130 | 23 | 17 |
| 0.1, 0.7, 30 | **84** | **148** | 80 | 64 |
| **total of 2,000** | **1,355** | **1,517** | 202 | 166 |

- **Recovered under the old rule but lost under the new one:** 4, all at ρ 0.1, blend 0.7, t_E 10 d.
- **Controls passing vetting (old → new):**
  - eclipse 0 → 0;
  - **fade 4 → 9** (of 200);
  - PSPL 0 → 0;
  - common mode 0 → 0.
- **Real calibration light curves:** 0 survivors under both rules (16 flagged runs in 3 groups).
- **Consistency check:** the count rule reproduces the v12 stored counts, and both rules flag the same events.

## 4. New validation (500 never-read groups, seeds 3100 + cell and 3200 + k; `2382val/decision.json`): **passed**

| cell | recovered / 200 |
|---|---|
| 0.01, 1, 10 | 138 |
| 0.01, 1, 30 | 174 |
| 0.01, 1, 100 | 164 |
| 0.1, 1, 10 | 130 |
| 0.1, 1, 30 | 166 |
| 0.01, 0.7, 10 | 143 |
| 0.01, 0.7, 30 | 168 |
| 0.01, 0.7, 100 | 134 |
| 0.1, 0.7, 10 | 128 (lowest, 0.64) |
| 0.1, 0.7, 30 | 147 |

- **Controls passing vetting:** eclipse 0, **fade 17 (8.5 %; limit 10 %)**, PSPL 0, common mode 0.
- **Real validation light curves:** 999 series, 6 flagged runs in 2 groups, 0 survivors.
- **The fade control rose from 4 (old rule, development set) to 9 (new rule, development set) and 17 (new rule,
  validation set).** The new rule admits more ordinary fades; the margin to the limit is 1.5 percentage points.
- **The validation set never touched the search group** (`2382val/source_usage.jsonl`: search_member_rows 0 for
  inject and search dev).
- **Coverage** (`coverage.json`, after the search): r-band series epochs 5 / 50 / 95 % were 1,611 / 2,021 / 3,822
  (validation; median 2,020.5) against 1,651 / 2,020 / 3,822 (search). Refmag r was 14.83 / 16.56 / 17.43 against 14.82 / 16.43 /
  17.42.

## 5. Search (`2382val/search_holdout_v13_time.json`)
- **One search of the 2,000 frozen groups:** 3,991 series, 96 flagged runs in 6 groups, **0 pass vetting**. All 96
  fail the one-off and/or the achromatic test:
  - one-off and achromatic: 90;
  - one-off only: 4;
  - achromatic only: 2.
- One group (724208400014500) holds 86 of the 96 runs.

## 6. Image check (frozen order; `2382val/image_check.json`, `image_check.png`, `column_check.json`)
- With no vetted survivor, the order took the flagged groups by non-detections in their longest run: 5 of 6 groups,
  50 cutouts (4.0 MB). The sixth (724208200014557) is **pending**, as are 26 non-detection exposures inside checked
  runs beyond the per-flag cap.
- **Automatic classes of the non-detection exposures:**
  - 17 "credible dimming" (class 4);
  - 4 "source present, catalogue missing it" (class 3);
  - 11 undetermined.
  - The automatic support rule called 3 groups a real sky change: 724208200007532, 723205400006654 and
    724208200006955.
- **Looking at the contact sheet, and in a post hoc column check** (`column_check_script.py.txt`), a dark detector
  column (column median > 5 robust σ below the cutout median) lies within the 2″ aperture in **17 of 17** class-4
  exposures (nearest column ≤ 0.55″ from the source). The same holds for **0 of 20** detected reference exposures
  (nearest 2.36″).
  - **Post hoc reading:** the non-detections coincide with a dark column at the source. This is consistent with a
    data problem in those exposures; it is **not image-confirmed dimming**.
  - In 4 of the 5 groups most cutouts show 2–3 dark columns. The catalogue appears to lose the source whenever one
    falls on it (not tested further).
  - The automatic "real sky change" labels are therefore not accepted (post hoc).
- **724208400014500** (86 runs): its non-detection exposures are class 3 (source present) or undetermined, so it is
  not a dimming.
- **Result:** no flag is an image-confirmed real dimming. No ordinary-versus-W3 comparison was needed. Nothing here is
  a W3 or wormhole candidate.

## 7. Cost (`cost.json`)
- **Downloads:**
  - Npix 2382 partition 1,846.5 MB (re-read; not stored);
  - objects ≤ 93 MB;
  - TAP 13.9 MB;
  - cutouts 4.0 MB.
  - Total ≤ 1,957.5 MB, under the 3 GB cap.

## Deviations and notes
- **Plan change (owner-authorized):** the re-split of the v12 search group into 500 validation and 2,000 search
  groups.
- **Post hoc:** `column_check` and `coverage`. The coverage script loaded search members after the search (one usage
  record).
- **Section numbers:** the headers of `split_script` and `decision_script` cite PLAN sections 4 and 5; the re-split
  and the pass rule are in PLAN section 3.
- **Lenient edge cases of the new rule,** applying equally to injections, controls and real data:
  - an epoch with a non-finite maglimit is never intermediate;
  - a run at the start or end of a series has an empty window on that side.
