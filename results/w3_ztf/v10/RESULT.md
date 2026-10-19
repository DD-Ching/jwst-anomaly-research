# W3 ZTF v10 — psfcat in the W3 chain: recovery rises by 161 of 2,000, but 8 of 10 cells stay below 0.5; search not executed

- Plan: `PLAN.md`, committed in 1b5b03f with the fill set (`fillset.json`) before any psfcat file was read.
- Scope: a method correction re-checked on the already-used v6 calibration sample (Npix 7134), not an independent
  validation. #339 (full W3 validation failed) and v6 stand. Efficiencies are catalogue-stage for the 10 cells
  only, not image-level completeness.
- Labels: psfcat, TAP and DR24 `observed`; counts `derived`; injections `simulated`.

## Verdict (`decision.json`, computed by `recal_script.py.txt`)
- **Failed: 2 of 10 cells reach 0.5** (ρ 0.01, blend 1, t_E 100 d: 124/200; ρ 0.01, blend 0.7, t_E 30 d: 102/200).
- The other conditions hold:
  - controls pass vetting in ≤ 10 %: eclipse 0, fade 16 (8 %), PSPL 0, common mode 0 of 200;
  - real calibration light curves (merged): 1,041 series, 1 flagged run, 0 survivors;
  - 0 search-member rows loaded (`source_usage.jsonl`).
- As the plan requires, the round stopped. **The 539-group search was not executed**: it is not a search with zero
  candidates. No psfcat value of a search member was extracted, and `search_script.py.txt` (written before the
  verdict) was not run. Nothing was re-tuned and no sample was redrawn.

## 1. Data entry
- **Fill set** (`fillset.json`): 724 exposures in 33 field / quadrant / filter combinations, where DR24 has no
  detection of any control star although the exposure is official, infobits 0 and covers a control star in play.
- **Downloads** (`download.json`, `download_pass2.json`, `cost.json`): 790 psfcat files, 183.6 MB, each
  downloaded once. One connection reset was retried on a second pass.
- **Overlap test per combination** (`overlap_v10.json`):
  - **32 of 33 passed.** Over the passed combinations, median Δ ranges from −0.019 to +0.006 mag, robust scatter
    from 0.002 to 0.013 mag, and pairs from 10 to 291.
  - **672_31_2 (field 672, CCD 8, q4, r) failed** on n ≥ 8: it has only 2 control stars, giving 3 pairs. Its 93
    fill exposures were not used.
- **Merge** (`merge_log.json`), over calibration and control rows:
  - **Used:** 631 files.
  - **Added:** 87,841 clean detections and 2,695 flagged or ambiguous (these count as unknown).
  - **Kept as DR24:** 261 star × exposure pairs where DR24 already had a detection.
  - **Missing:** 8,087 star × exposure pairs had no psfcat row. They are unknown, never a non-detection.
  - **Good quadrant exposures:** 14,483 → 15,111.
- **Calibration series** (`stats_v10.json`):
  - **r band:** 397,386 → 430,030 epochs (+32,644: 32,643 in fill exposures, 1 elsewhere);
  - **g band:** 512,761 → 525,866 (+13,105: 13,103 in fill exposures, 2 elsewhere);
  - The 3 epochs outside fill exposures are not separated. A likely mechanism is that psfcat points widen a star's
    span, which changes in-play counts and good status in a few DR24 exposures.
  - 0 epochs lost.
- Per-event stages and epoch counts before and after: `events_v10.json`.
- There are 509 carriers before and after, and no event's carrier was lost (no `carrier_lost` stage). No carrier
  was added.

## 2. Same events before and after (`injections_v10.json`)
- **Events are the same.** Replaying stage_inject's draws on the DR24-only series reproduces the stored v6 flagged and
  recovered counts in all 10 cells (`replay_mismatches_vs_v6: []`). Those (carrier, t0, u0) were injected into the
  merged series.
- **Noise scheme:** one generator per (cell seed, event, expid). "Before" is DR24-only with that scheme. The scheme
  alone moves the total from 703 (v6 stored) to 711.

| cell (ρ, blend, t_E) | v6 stored | before | after | bracket met before → after | sampling-eligible before → after |
|---|---|---|---|---|---|
| 0.01, 1, 10 | 61 | 61 | 64 | 104 → 145 | 65 → 77 |
| 0.01, 1, 30 | 80 | 80 | 93 | 95 → 131 | 89 → 113 |
| 0.01, 1, 100 | 90 | 90 | **124** | 93 → 138 | 93 → 134 |
| 0.1, 1, 10 | 80 | 80 | 81 | 122 → 156 | 90 → 101 |
| 0.1, 1, 30 | 76 | 76 | 88 | 102 → 146 | 92 → 116 |
| 0.01, 0.7, 10 | 80 | 82 | 91 | 121 → 155 | 85 → 95 |
| 0.01, 0.7, 30 | 76 | 79 | **102** | 96 → 140 | 91 → 115 |
| 0.01, 0.7, 100 | 61 | 59 | 95 | 84 → 136 | 83 → 130 |
| 0.1, 0.7, 10 | 51 | 52 | 62 | 113 → 149 | 62 → 74 |
| 0.1, 0.7, 30 | 48 | 52 | 72 | 105 → 148 | 93 → 114 |
| **total of 2,000** | **703** | **711** | **872** | **1,035 → 1,444** | **843 → 1,069** |

- Column definitions:
  - **bracket met:** the uninjected series has ≥ 20 epochs before and after the umbra.
  - **sampling-eligible:** the v7 definition, i.e. the bracket plus ≥ 3 umbra epochs on ≥ 3 nights.
- **Where the gain comes from** (per-cell stages and transitions in `injections_v10.json`): mainly the MIN_BRACKET
  pre-screen requirement. It stops 289 events before the fill and 171 after; for example 72 → 38 at ρ 0.01, blend 1,
  t_E 100 d. "No dim run in the umbra" goes from 647 to 602.
- **Vetting failures rise** where the added epochs make runs testable: 113 → 128 in total. Achromatic failures go
  from 35 to 46; the largest rise is at ρ 0.01, blend 0.7, t_E 100 d (20 → 29).
- **Fade control:** passing runs rise from 7 (before) to 16 of 200 after the fill, still under the 10 % limit.

## 3. What still fails after the fill (1,128 events; `stats_v10.json`)
- **No dim run in the umbra: 602.** For 527 of them the carrier has < 3 official infobits-0 exposures inside the
  umbra (unique r-band, in-footprint exposure times; g band not counted): too few official r-band exposures there.
- **MIN_RUN 197, MIN_NIGHTS 30.** Not checked against official exposures.
- **MIN_BRACKET: 171.**
  - For 86, the official exposures before or after the umbra are fewer than 20 (too few official r-band
    exposures on a side).
  - For 85, official r-band exposures are ≥ 20 on both sides, but the series has fewer than 20 non-dim epochs on a
    side. Possible reasons, not separated here: the frozen good-exposure rule, partly covered DR24 exposures (not in the fill rule), missing psfcat rows,
    the unvalidated combination, and dim epochs of the event itself.
- **Vetting: 128** (sharp-edge, achromatic, common-mode and one-off tests; per cell in `injections_v10.json`).
- **Implementation:** these checks found no error: the replay is exact, and two new unit tests cover the merge
  rules. They do not prove the absence of implementation errors.

## 4. Cost
- **Downloads:**
  - light-curve partition 299.9 MB, re-read and not stored;
  - objects columns ≤ 86.8 MB;
  - TAP 4.7 MB;
  - psfcat 183.6 MB.
  - Total ≤ 575 MB, under the 1.5 GB cap. No image cutouts.
- **Compute:** the re-run took 2.7 min (`decision.json`).

## Deviations and notes
- **Post hoc:** `stats_script.py.txt` (series epochs; official-exposure counts for the remaining blocks) was added
  after the verdict. It is diagnostic only.
- **Merge timing:** the mid-time is computed once per file at the quadrant stars' mean position, as `coverage()`
  does, instead of per star. This changed after a first, slow attempt was stopped during the merge, before any
  injection was run.
- **Colour term:** not applied, as in v9.
- **Deviation, too-faint injected epochs:** whether such an epoch is dropped or becomes a non-detection is decided
  by the exposure (fill or not), not by the point's source. The 261 DR24 detections kept inside fill exposures
  (`merge_log.json`) are therefore dropped instead of becoming non-detections, contrary to PLAN section 4.
  The direction of its effect on recovery is **not verified** (addendum below); it was not re-run.
- **Review fix after the first run:** a psfcat row with a non-finite magnitude now counts as flagged (it was
  unguarded). The ambiguous branch also no longer writes a placeholder magnitude. The whole re-run and the stats
  were repeated, and every committed number is from that repeat. Apart from the runtime field, its decision, injections, merge,
  overlap and stats outputs were identical to the first run (compared in the session, not stored).
- **Usage log:** every `source_usage.jsonl` record names `fillset.json` as the stage, because `context()` logs
  `sys.argv[1]`. The records are the recal runs (including the stopped one) and the stats runs, all with 0
  search-member rows. `stats_script.py.txt` re-runs the recal setup through runpy, so it rewrites
  `overlap_v10.json` and `merge_log.json` with identical content.
- **Not computed by code:** the time and the tests/review conditions of PLAN section 5. The verdict had already
  failed on the cells.

## Addendum 2026-10-19 (after #359; from stored files only)
- **Why "before" is 711 and v6 stored is 703.**
  - The two runs use the same events. The DR24 re-read has the same row count (below); byte identity was not
    checked.
  - The replay used the v6 sequential-noise scheme on the re-read inputs and reproduced the v6 per-cell flagged and
    recovered counts exactly (`decision.json` → `burden_calibration_merged.replay_mismatches_vs_v6: []`). The
    control passed counts were reproduced likewise (`injections_v10.json` `v6_replay`; v6 stored no flagged count
    for controls).
  - "Before" differs from that replay in **one thing only: the noise draws**. They come from a generator per
    (cell seed, event, expid) instead of one sequential generator per cell, so the realised noise differs.
  - The +8 (703 → 711) is therefore a noise-realisation difference of the same events. It is not a data or code
    change.
  - **What is known about the inputs:** the light-curve partition re-read gave 12,209 rows, the same as v6. The TAP
    re-query gave 21,632 rows, also the same, but its byte count was 4,660,671 against v6's 4,509,175. The v10-side
    values are from the session log, not stored, and the content was not compared. The code included the v8 fixes (exptime mid-time, single-exposure match). The exact replay shows
    that none of these changed the stored v6 counts. It does not show the inputs are byte-identical.
- **What +161 compares:** "before" (711) and "after" (872) use the same events, code, noise scheme and per-expid
  seeds. They differ only in the data: DR24 alone vs DR24 + validated psfcat fill (and the FILL rules that come with
  it). It is not a comparison with the v6 stored 703.
- **The 261 dropped points:** the earlier "conservative (it can only lower recovery)" is withdrawn.
  - Dropping a too-faint epoch instead of making it a non-detection removes a dim epoch. That can shorten a dim run
    (lower recovery). It can also remove a second dim run elsewhere, and so pass the one-off test (higher recovery).
  - The direction is **not verified**. The effect was not measured.
