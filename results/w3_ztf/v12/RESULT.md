# W3 ZTF v12 — new-field area Npix 2382 fails validation (1 of 10 cells at 0.42); search not executed

- **Plan:** `PLAN.md`, committed in ad00408 with the dedupe list and the candidates, before any read.
- **Area check, budget decision and sample:** committed in 664cece, before any light curve.
- **Labels:** objects, TAP, IBE and DR24 `observed`; counts `derived`; injections `simulated`.
- **Scope:** the W3 negative-mass point-lens model only (10 calibrated cells, frozen thresholds, unchanged code
  apart from a `configure()` per-area folder option). No limit is derived. This round says nothing about wormholes
  in general.

## Verdict (`2382/decision.json`, frozen v6 rule by `decision_script.py.txt`)
- **Failed:** 9 of 10 cells reach 0.5. The cell ρ 0.1, blend 0.7, t_E 30 d recovers **84/200 (0.42)**.
- The other conditions hold:
  - controls passing vetting: eclipse 0, fade 4, PSPL 0, common mode 0 of 200;
  - real calibration light curves: 998 series, 16 flagged runs in 3 groups, 0 survivors.
- As planned, the area stopped. **Its 2,500 search groups were not searched**: they were never loaded
  (`2382/source_usage.jsonl`, search_member_rows 0). This is "not executed", not a search with zero candidates.
- No threshold was changed, no seed redrawn, and no second area read. Npix 8170 had already been dropped before any
  light curve by the budget rule (section 2).

## 1. Dedupe and candidates
- **Dedupe list** (`dedupe.json`): 4,072 earlier calibration or search groups (v1–v3, v4, v6/v10, v11), with
  positions and objectids. New groups share no objectid with them, and none lies within 3″.
- **Comparison stars of Npix 3838, 3839 and 7134** were read before but never searched. They are listed by rule.
  None of the chosen area's stars is among them, because it is a different sky area.
- **Candidates** (`candidates.json`, partition lists only): Npix 183, 2380, 2382 and 8170 (new parents), plus the
  fallback 3837 (old fields). That is 5 candidates.

## 2. Observing-condition and product check (`obs_conditions.json`; metadata only)

| Npix | worst-cell necessary fraction | median product ratio | eligible new groups | fields shared with earlier samples | qualifies |
|---|---|---|---|---|---|
| 183 | 0.547 | 1.005 | 2,303 | none | yes |
| 2380 | 0.4705 | 1.898 | 4,418 | none | no (< 0.5) |
| 2382 | **0.6965** | 1.038 | 4,028 | none | **yes, first** |
| 8170 | 0.6615 | 1.050 | 6,238 | none | yes, second |
| 3837 | 0.6975 | 1.025 | 5,278 | 826, 1857 | yes (fallback) |

- **Choice:** 2382 and 8170, the two best qualifying areas with no shared field.
- **Budget (`budget.json`):** their light-curve partitions are 1,846.5 and 2,207.9 MB. With 390 MB already spent on
  objects and TAP, both areas would need 4.4 GB, so only the first, **2382**, was used (PLAN section 7).
- **Product ratio:** the summed r-band `ngoodobs` per counted exposure. It can exceed 1 when a source has several r
  objectids. As a gap detector it is coarse.
- **psfcat HEAD:** 14–20 of 20 per area returned 200; the causes of failure are not separated.

## 3. Npix 2382 sample and realised sampling
- **Sample** (`2382/sample.json`, `samples_info.json`): 3,000 new groups, split 500 calibration / 2,500 search by
  crc32 rank. There are 17,853 control stars.
  - Refmag r percentiles (5/50/95): calibration 14.79 / 16.44 / 17.40; search 14.82 / 16.46 / 17.42.
  - Search-group epoch counts were not read, because its light curves were never loaded.
  - The `holdout: 1491` in the stored `info` is `select()`'s parity count from before the v12 rank split. The
    split used is in each group's `holdout` flag.
- **Realised sampling** (`2382/realised_sampling.json`):
  - calibration series hold a median 94.2 % of their official exposures;
  - r-band epochs per carrier, 5/50/95 %: 1,734 / 2,022 / 3,823;
  - realised indicator 0.717–0.892 per cell (0.8405 for the failing cell).

## 4. Injection results and the reason for the failure (`2382/injections_v12.json`; post hoc trace `2382/trace_posthoc.json`)

| cell (ρ, blend, t_E) | flagged | recovered / 200 | first vetting failure (trace) |
|---|---|---|---|
| 0.01, 1, 10 | 150 | 150 | — |
| 0.01, 1, 30 | 177 | 173 | sharp 4 |
| 0.01, 1, 100 | 160 | 150 | sharp 9, one-off 1 |
| 0.1, 1, 10 | 146 | 130 | sharp 15, one-off 1 |
| 0.1, 1, 30 | 180 | 122 | sharp 57, one-off 1 |
| 0.01, 0.7, 10 | 148 | 138 | achromatic 6, one-off 3, sharp 1 |
| 0.01, 0.7, 30 | 178 | 155 | achromatic 21, sharp 1, one-off 1 |
| 0.01, 0.7, 100 | 179 | 136 | achromatic 28, sharp 12, one-off 3 |
| 0.1, 0.7, 10 | 143 | 117 | sharp 23, achromatic 2, one-off 1 |
| **0.1, 0.7, 30** | **171** | **84** | **sharp 80, achromatic 7** |

- The trace (`../v7/trace_script.py.txt`, the same draws) reproduces the stored flagged and recovered counts.
- **The failing cell is not limited by sampling:**
  - 171 of 200 injections are flagged by the pre-screen;
  - the necessary fraction there is 0.828 (metadata) and 0.8405 (realised).
- **The loss is in vetting:** 87 flagged injections fail it, 80 of them first on the frozen sharp-edge test (≤ 1
  intermediate-depth epoch within one run length of the run).
- A finite source (ρ 0.1) dims and recovers gradually. The reading that this area's dense cadence (median 2,022
  r-band epochs per carrier) samples those edges more often is a hypothesis; it is not tested here.
- For comparison, in the v11 area (Npix 3839) the same cell recovered 131/200 (`../v11/injections_v11.json`); no
  trace was stored there.

## 5. Cost (`cost.json`)
- **Downloads:**
  - objects tables 372.4 MB (4 parents);
  - TAP 18.0 + 13.9 MB;
  - Npix 2382 light-curve partition 1,846.5 MB (not stored);
  - 100 HEAD requests;
  - cutouts 0 MB.
  - Total 2,250.8 MB, under the 3 GB cap.

## Deviations and notes
- **Post hoc:** `2382/trace_posthoc.json` is a failure diagnosis on the calibration group only.
- **Plan wording:** the candidate rule asked for "2, 2, 1 children", but only 4 partitions met it, because two
  parents had one eligible child each. That is recorded in `candidates.json`.
- **Code:** `configure()` accepts a per-area folder ("v12/2382"). This is plumbing only.
