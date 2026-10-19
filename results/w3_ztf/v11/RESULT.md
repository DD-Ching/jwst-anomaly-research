# W3 ZTF v11 — the area chosen by observing conditions passes validation; 488 search groups searched, 0 flags

- **Plan:** `PLAN.md`, committed in 09f2a12 before any metadata of this round was read.
- **Area choice:** committed in dbe3e68 before any light curve was read.
- **Validation:** committed in 2184f8b before the search group was read.
- **Labels:** objects, TAP, IBE and DR24 `observed`; counts `derived`; injections `simulated`.
- **Scope:**
  - No threshold, model or tool changed. No psfcat fill: DR24 sampling here is near-complete (section 2).
  - No abundance limit is derived.
  - The stars share fields, exposures and processing with the v1–v5 samples of the same area, so this is **not an
    independent-area validation**.
  - Efficiencies are catalogue-stage only, for the 10 cells.

## 1. Observing-condition check (`obs_conditions.json`; metadata only)
- **Candidates:** Npix 3837, 3838 and 3839 (order-4 pixel 959). Sources within 3″ of a v1–v3 or v4 sample group
  are excluded: 1,000 groups each in 3838 and 3839. That leaves 5,278, 3,633 and 3,600 eligible groups.
- **Positions:** 200 per area, drawn from the seed-7 draw of ≤ 1,000 groups (the would-be sample), not from all
  eligible groups (review note; the plan said "eligible groups"). For 3839 that draw is the v11 sample.
- **Exposures:** a median of 1,392–1,427 official infobits-0 r-band exposures per position inside the footprint.
  - 114–117 more per position (median) have infobits ≠ 0, so they have no usable product under the frozen rule.
  - **psfcat HEAD check:** 17 (3837), 19 (3838) and 15 (3839) of 20 sampled exposures returned 200. A failure and a
    timeout are not separated.
- **Necessary-sampling fraction** (2,000 draws per cell; noise-free model on the official exposures, frozen
  pre-screen; seasonal gaps kept):

| cell (ρ, blend, t_E) | 3837 | 3838 | 3839 |
|---|---|---|---|
| 0.01, 1, 10 | 0.737 | 0.734 | 0.7505 |
| 0.01, 1, 30 | 0.859 | 0.864 | 0.868 |
| 0.01, 1, 100 | 0.820 | 0.818 | 0.823 |
| 0.1, 1, 10 | 0.701 | 0.705 | 0.7145 |
| 0.1, 1, 30 | 0.821 | 0.820 | 0.833 |
| 0.01, 0.7, 10 | 0.724 | 0.719 | 0.7245 |
| 0.01, 0.7, 30 | 0.860 | 0.854 | 0.8625 |
| 0.01, 0.7, 100 | 0.824 | 0.819 | 0.825 |
| 0.1, 0.7, 10 | **0.687** | **0.6855** | **0.7025** |
| 0.1, 0.7, 30 | 0.816 | 0.8085 | 0.817 |

- **Choice:** the rule picks the highest worst-cell fraction, which is **Npix 3839 (0.7025)**. All three areas score
  ≥ 0.5. The tie-break was not needed. Its stored count is capped at 1,000 per area, so it could not have applied as
  written.
- **What is missing, by cell (3839, first missing requirement in `first_missing`):**
  - t_E 10 d: mostly no exposure during the dim part of the event (231–342 of 2,000) or fewer than 3 dim epochs
    (123–177).
  - t_E 100 d: mostly the 20-epoch bracket (330–337).

## 2. Data and sample
- **Sample** (`sample.json`): 1,000 groups, seed 7, split into 512 calibration and 488 search groups
  (crc32("v11:" + key)). There are 32,966 control (comparison) stars, none of them a v11 sample member
  (`sample_check.json`).
  - The nearest v1–v4 sample group is 5.4″ away. No groups lie within 3″, and no objectids are shared
    (`sample_check.json`).
- **Light curves:** 38,195 rows of the DR24 partition (`fetch_info.json`); TAP 38,026 exposure rows
  (`exposures_summary.json`).
- **Realised sampling** (`realised_sampling.json`), on the calibration carriers' actual series:
  - The series holds a median 98.2 % of each carrier's official infobits-0 in-footprint r-band exposures (10th–90th
    percentile 97.8–99.3 %).
  - The same indicator gives 0.721–0.889 per cell.
  - This is the near-complete DR24 coverage that Npix 7134 lacked (v7–v10).

## 3. Validation (frozen v6 rule, `decision.json` from `decision_script.py.txt`): **passed**

| cell | recovered / 200 | rate |
|---|---|---|
| 0.01, 1, 10 | 155 | 0.775 |
| 0.01, 1, 30 | 171 | 0.855 |
| 0.01, 1, 100 | 167 | 0.835 |
| 0.1, 1, 10 | 136 | 0.680 |
| 0.1, 1, 30 | 150 | 0.750 |
| 0.01, 0.7, 10 | 149 | 0.745 |
| 0.01, 0.7, 30 | 149 | 0.745 |
| 0.01, 0.7, 100 | 147 | 0.735 |
| 0.1, 0.7, 10 | 123 | 0.615 |
| 0.1, 0.7, 30 | 131 | 0.655 |

- **Controls passing vetting:** eclipse 0, fade 8 (4 %), PSPL 0 and common mode 0 of 200.
- **Real calibration light curves** (`search_dev_v11.json`): 1,022 series, 0 flagged runs, 0 survivors.
- **Carriers:** 512. The search group was never loaded during calibration (`source_usage.jsonl`,
  search_member_rows 0 in the three calibration records: inject, search dev, and `realised_script` logged as
  stage "?"). The fourth record is the holdout search, after the pass, with 2,752 search-member rows.
- **The necessary-sampling fraction is not a strict per-sample bound.** At ρ 0.01, blend 1, t_E 10 d, recovery is
  0.775 against an indicator of 0.7505 (3839) and 0.753 (realised). The 200 injections are a different draw from the
  2,000, so the binomial error of about ±0.03 (one σ; computed in the session, not stored) covers the difference. The indicator also ignores real non-detections
  already in the series. "Upper bound" holds on average, not draw by draw.
- **Pre-search review** (subagent): no finding affects the verdict. Recorded:
  - the position draw (section 1);
  - the tie-break;
  - the psfcat HEAD cause, not separated;
  - the comparison stars exclude only v11 members (as in v6).

## 4. Search (`search_holdout_v11.json`)
- **One search of the 488 search groups:** 969 series, **0 flagged runs, 0 flagged groups**. Nothing was added.
- **Image check:** not needed, because there were no flags. No cutout was downloaded (0 MB).
- **Reading:** under the frozen rules, no W3-like dimming run was flagged in these 488 groups. This is a
  catalogue-stage statement for the 10 calibrated cells. It is not an abundance limit (none derived) and says
  nothing about other model ranges. No source is a candidate of any kind.

## 5. Cost (`cost.json`)
- **Downloads:**
  - light-curve partition 944.9 MB (7 columns, not stored);
  - objects 92.6 MB (11 columns, read once and cached for the sample);
  - TAP 7.3 + 8.2 MB;
  - 60 HEAD requests;
  - cutouts 0 MB.
  - Total 1,053 MB, under the 1.5 GB cap. The plan's estimate (≈ 0.85 GB for the partition) was low.

## Deviations
- `obs_script.py.txt` line 31: the comment was corrected after the run ("all eligible groups" → the seed-7 draw).
  Code and outputs are unchanged.
- The `configure()` guard for Npix 3839 with vdir v11 was added after the plan, before the sample was written. It
  keeps the v1–v3 `sample.json` untouched; the review confirmed this.
