# W3 ZTF v22 — new inventory of 31 partitions (25 qualify, none sharing fields), frozen queue of 8: 1 passes, 1,882 groups searched, 0 survivors; the 8124 neighbour also sits on a low-responsivity column in the shared NDs

- **Pre-registration:** the plan and candidates were committed before any data read (78d4627, 06:54:10; the
  inventory started at about 06:54:20 by local file times, not stored). The queue and samples were committed before
  any light curve (d264f7e, 07:44:43; first light-curve fetch 07:45:02).
- **Labels:** DR24, TAP, IBE images, masks and psfcat are `observed`; counts are `derived`; injections are
  `simulated`.
- **Scope:**
  - the frozen v13 W3 rule, the 10 cells, the v16 thresholds and cuts;
  - no limit; nothing here bears on wormholes;
  - the cumulative result is "no confirmed target event", not "no star dimmed".
- **Undetermined W3 ZTF items** (this unit only; not a project total):
  - **tracked from v20–v21:** 8126 / 634210200007847 and 8124 / 634205400002525;
  - **archived in v19:** 1624 and 1625.

## A. 8124 neighbour check (`neighbour_check.json`; 80 s, 1.6 MB, 7 mask cutouts)
1. **Different sources:** yes. 634205400002482 and 634205400002525 are 9.85″ apart and have separate psfcat matches
   in the bracket detection 214109164.
2. **Shared NDs:** all 3 of the target's NDs (214310269, 214408748, 214508458) are also NDs of the neighbour. The
   neighbour has 4 more (7 in its longer window) and the target has none the neighbour lacks.
   - Both are r band, in the same field and quadrant (634, rcid 19).
   - NDs count only infobits-0 exposures (a clarification of PLAN A).
3. **The same detector problem?** In all 3 shared NDs both positions carry official mask bit 2 (LOW RESPONSIVITY)
   at the centre pixel, and each lies 0.04–0.50 px (x) from a bit-2 column in the cutout. Column identity and pixel
   positions were not stored.
   - The target's mask values come from this round's cutouts, not from v21 (deviation from PLAN A).
   - In 2 more neighbour NDs (214409141, 220251997) the neighbour also sits on bit 2.
   - In the other 2 neighbour NDs (214408356, 214409533) both positions show only bit 1, 0.66–1.25 px from the
     column.
- **What this shows:** two distinct sources missed in the same exposures, each on an official low-responsivity
  column. That is consistent with a detector problem causing both catalogue misses.
- **What it does not show:** it is not proof. No flux was measured. **8124 stays undetermined.** No veto rule was
  added.

## B. Inventory (`candidates.json`, `inventory.json`, `queue.json`)
- **Candidates:** 31 partitions under 8 parents never read before, chosen in farthest-point order from the used
  parents.
  - **Deviation (selection rule):** the rule was changed once before freezing, after reading only the HATS
    `partition_info.csv`. The first version gave 17 edge-parent candidates; the frozen one requires ≥ 3 eligible
    children per parent. v16 used even spacing in Npix.
- **Time (deviation):** 49.0 min, over the 30-min cap (319,648 TAP rows); the cause of the slowness was not
  measured.
- **Downloads:** 759.7 MB (objects 698.2 MB from parquet metadata, TAP 61.6 MB).
- **Results:**
  - 64,510 catalogue-eligible unused groups;
  - **25 partitions qualify**, holding 56,519 unused and 44,019 searchable groups. None shares a ZTF field number
    with v1–v21 sample objectids, judged on 200 sampled positions per area. A full-sample check for the 8 queued
    areas agrees (not stored); the other 17 were not checked in full;
  - 6 do not qualify.
- **Scope:** this round's inventory only, not the remaining supply of ZTF or the sky.

## C. Queue (`summary.json`, each area's `decision.json`)

| pos | Npix | verdict | failing (cells < 100/200; fade passing / 1,000) | searched |
|---|---|---|---|---|
| 1 | 1853 | failed | ρ0.1/1/10 d 89, ρ0.1/0.7/10 d 85; fade 117 | — |
| 2 | 1150 | failed | ρ0.1/0.7/10 d 96; fade 117 | — |
| 3 | 1148 | failed | ρ0.1/0.7/10 d 84; fade 102 | — |
| 4 | 1151 | failed | fade 114 | — |
| 5 | 5054 | failed | fade 104 | — |
| 6 | 4991 | failed | fade 138 | — |
| 7 | 2344 | **passed** | — | 1,882 groups (3,709 series), 2 runs in 2 groups, **0 survivors** |
| 8 | 4989 | failed | ρ0.1/1/10 d 93, ρ0.1/0.7/10 d 93; fade 112 | — |

- **Fade control:** every failing area exceeds it (10.2–13.8 %, limit 10 %). It was the only failing condition in 3
  (1151, 5054, 4991); the other 4 also failed ρ 0.1 / t_E 10 d cells. A similar fade excess occurred in v20
  (10.6–13.9 %). Eclipse, PSPL and common-mode controls were 0 / 200 everywhere, and every area had 0 calibration
  survivors.
- **No failure was retried, tuned or replaced.** The failed areas' 12,439 search groups were downloaded with the
  sample but never loaded by any stage (`source_usage.jsonl`). 1150 and 1148 each lack 1 sample member.
- **Deviation (execution setting only):** the driver ran with 1 worker after the v21 memory kill, set in a separate
  commit (3db7a4e) at the driver start (`timeline.jsonl`). It had no crash.
- **Exposure independence:** 2344's fields (632, 678, 679, 1676) share no ZTF field number with v1–v21 sample
  objectids (inventory sample and full sample; the full check is not stored). All 1,882 searched groups therefore
  come from fields with no earlier sample data.

## D. Survivors
- **0 final survivors**, so no image check was needed: 0 groups and 0 cutouts this round.
- Reliable real dimmings: 0.
- New undetermined: 0. The earlier undetermined items are unchanged (preamble).

## E. Totals and cost (measured or estimated, as marked)
- **Searched this round:** 1,882 new groups, all in fields with no earlier sample data (by field number).
- **Cumulative with the v13 rule:** 37,357 groups (35,475 + 1,882; not stored as one number). No confirmed target
  event.
- **Downloads:** 4,127.2 MB (`summary.json`; objects and light curves from parquet metadata, TAP measured), plus
  1.6 MB for part A.
- **Time:**
  - inventory 49.0 min;
  - queue 25.6 min (07:45:00–08:10:36);
  - part A 80 s.
- **Per 1,000 searched groups:**
  - in the passing area: 140 stage-seconds and 247 MB of light curves;
  - for the whole round (inventory, 7 failed areas, 1 passing): about 2.2 GB and 40 min (not stored);
  - review load: 0.

## F. What remains in this round's inventory (`remaining.json`)
- **Not queued:** 17 qualifying partitions with 38,198 unused groups (29,698 searchable) in 10,265.9 MB.
  None shares a field with earlier samples by the 200-position check; they were not checked in full.
- **Failed:** 12,439 search groups stay unsearched under the frozen rule.
- **Fade control:** 7 of 8 areas exceeded it (10.2–13.8 %). It was the only failing condition in 3; the other 4
  also failed short-t_E cells. Any change to the frozen calibration is a method decision for the owner; this round
  makes none.
