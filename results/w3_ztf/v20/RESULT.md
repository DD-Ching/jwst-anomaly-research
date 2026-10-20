# W3 ZTF v20 — frozen queue of 9 from the v16 inventory: 7 started, 2 pass, 7,724 groups searched, 1 final survivor (undetermined, consistent with near-limit catalogue misses)

- **Plan:** `PLAN.md` and `queue.json` were committed before any light curve (2242ac9).
- **Labels:** DR24, TAP, IBE images and masks are `observed`; counts are `derived`; injections are `simulated`.
- **Scope:**
  - the frozen v13 W3 rule, the 10 cells, all v16 thresholds and cuts;
  - no limit derived; nothing here bears on wormholes;
  - the cumulative result is "no confirmed target event yet", not "no star dimmed".
- **Archived:** the two v16 sources (1624 / 708205400009336, 1625 / 709212400017577) are recorded as "no reliable
  evidence of a real dimming; some exposures still cannot be judged", per `../v19/RESULT.md`. They were not re-run.

## 1. Queue and verdicts (`summary.json`, each area's `decision.json`)

| pos | Npix | verdict | failing condition | fade passing / 1,000 | searched |
|---|---|---|---|---|---|
| 1 | 1356 | failed | cells: ρ0.01/0.7/100 d 97, ρ0.1/0.7/10 d 89 (of 200) | 66 | — |
| 2 | 2490 | failed | cell ρ0.1/0.7/10 d 95; fade 12.4 % | 124 | — |
| 3 | 2491 | failed | fade 13.9 % | 139 | — |
| 4 | 8126 | **passed** | — (lowest cell 110 / 200) | 70 | 4,608 groups (9,206 series), 35 runs in 8 groups, **1 survivor** |
| 5 | 2488 | failed | fade 12.8 % | 128 | — |
| 6 | 2489 | failed | cell ρ0.01/0.7/10 d 99; fade 10.6 % | 106 | — |
| 7 | 2900 | **passed** | — (lowest cell 137 / 200) | 54 | 3,116 groups (6,186 series), 11 runs in 2 groups, 0 survivors |
| 8 | 8124 | not started | download budget (projection 5,104.5 + 1,734.6 MB > 5,900) | — | — |
| 9 | 8127 | not started | download budget | — | — |

- **No failure was retried, tuned or diagnosed further.** The failed areas' search groups (6,243) were never read.
  Eclipse, PSPL and common-mode controls were 0 / 200 everywhere, and every area had 0 calibration survivors.
- **Descriptive only, not a rule:** the four failing 2488–2491 areas share new fields 760, 1759/1760 and 1796/1797,
  and all four fail the fade control (10.6–13.9 %).
- **Exposure independence** (`queue.json`):
  - 8126 shares fields 634 and 1679 with earlier samples; 1631 and 1678 are new.
  - 2900 shares 848 and 1856; 825 is new.
  - The new groups are deduplicated against all 36,062 earlier sample groups by position (3″) and objectid, so they
    are new physical sources. Their exposures in shared fields are not independent of earlier data.
- **Incident (deviation):**
  - The copied v16 driver deadlocked when its budget-stop rule fired: the stop event was logged while the lock was
    held. v16 never reached that rule.
  - No new area was started after 2900, as the rule intended.
  - The driver was stopped at 05:21. 2900 was resumed from its validation stage with the same commands and seed.
    Its data stages had already finished.
  - The exposures end and stop events were logged by hand, marked in `timeline.jsonl`. The driver is fixed.

## 2. Final survivor and image check (`8126/image_check.json`, `mask_check.json`, `image_check.png`; 8 science and 4 mask cutouts, 0.73 MB)
- **Source:** 8126 / 634210200007847, g band, 4 NDs (HMJD 58208.51–58210.37).
  - It is faint: the after-run detections give 66–72 µJy (2″ aperture), g ≈ 19.3 AB.
  - The ND image limits are 19.24–19.52.
- **Automated classes of the 4 NDs:** 2 × poor image and 2 × undetermined. There is no field-star reference: the
  reference is NaN.
- **Measured fluxes:** the 2″ fluxes in the NDs are 43–67 µJy (errors 7–8). The two before-run exposures give 44 and
  47 µJy.
- **Masks:** bits 1 and 7 (detection, pixel spike) in 3 NDs; no low-response, dead or saturated bits.
- **Visual reading** (`image_check.png`, separate from the classes): the source is faintly visible in every stamp,
  NDs included.
- **Reading:** consistent with near-limit catalogue misses; no image-confirmed drop. Formally **undetermined**, not
  a confirmed dimming.

## 3. Totals and cost (measured; `summary.json`, `timeline.jsonl`, `objects.json`)
- **Searched this round:** 7,724 new groups (15,392 series).
- **Cumulative with the v13 rule:** 25,497 groups. That is 17,773 from v13–v16 (CHANGELOG) plus 7,724; the groups
  are disjoint by the dedupe (not stored as one number).
- **Downloads:** 4,947.2 MB.
  - Objects tables are 344.4 MB. They were read again because the v16 files were deleted after #385.
  - Light curves come from parquet metadata (an estimate of bytes read); TAP 52.7 MB is measured.
  - Cutouts add 0.73 MB.
- **Time:**
  - the queue ran 05:00:44–05:29 UTC (29.1 min to the last logged event), plus 3.4 min for the objects tables;
  - stage-seconds of the two searched areas total 1,285.5 s;
  - **per 1,000 searched groups:** 166 stage-seconds, 3.8 queue minutes, 300.5 MB of light curves (searched areas
    only). Including the 5 failed areas, the round spent about 640 MB and 4.2 min per 1,000 searched groups (not
    stored).

## 4. What remains (`remaining.json`)
- **Not started (budget):** 8124 and 8127, with 10,978 unused groups (9,978 searchable) in 3,438.6 MB of light
  curves.
  - Both are sampling-qualified, with worst cells 0.636 and 0.657. That is no promise of a pass.
  - Both share fields 634, 635 and 1679 (8127 also 681) with earlier samples.
- **Failed areas:** 6,243 unread search groups. They are not searchable under the frozen rule.
- **Beyond these:** the v16 inventory of 25 partitions is now used up apart from 8124 and 8127.
