# W3 ZTF v6 — cross-sky validation failed (10 of 10 cells below 0.5); search group not read

- Plan: `PLAN.md`, frozen in commit 31baace before any light curve in the new area was read.
- Every number below is in a JSON file in this folder. Injections are `simulated`; everything else is `derived` from
  ZTF DR24 and IRSA metadata.

## Verdict (frozen rule applied by code: `decision.json`)
- **Failed.** Every one of the 10 calibrated cells recovered < 0.5 on the calibration group of the new area.
- The other two conditions hold:
  - controls: eclipse 0/200, fade 8/200, PSPL 0/200 and common mode 0/200 pass vetting;
  - the real calibration light curves give 1 flagged run and 0 survivors.
- As the plan requires, the round stopped. The **search group (539 groups) was not read**: `source_usage.jsonl`
  shows 0 search-member rows in every stage that loads light curves through `context()` (inject, calibration
  search, diagnostic). The fetch stage downloaded the whole partition, and the exposure stage read quadrant metadata
  only. Nothing was re-tuned.

## Recovered / injected per cell (`injections_v6.json`; 509 carriers from 533 calibration groups)

| ρ | blend | t_E 10 d | t_E 30 d | t_E 100 d | v4 (old area) for comparison |
|---|---|---|---|---|---|
| 0.01 | 1 | 61/200 (166 stars) | 80/200 (166) | 90/200 (163) | 147, 170, 173 |
| 0.01 | 0.7 | 80/200 (166) | 76/200 (161) | 61/200 (169) | 142, 165, 148 |
| 0.1 | 1 | 80/200 (159) | 76/200 (166) | — | 124, 142 |
| 0.1 | 0.7 | 51/200 (174) | 48/200 (165) | — | 128, 120 |

- In the table, "stars" is the number of distinct carrier stars used per cell. Injections on one star are not
  independent sky samples.
- These are catalogue-stage efficiencies, not image-level efficiencies for full disappearance.

## Why it failed: post hoc, calibration group only (`failure_posthoc.json`)
- **Cadence:**
  - r-band nights per carrier: median 411 (10th–90th percentile 262–544) over a median 2,751 d;
  - that is a median of 54.8 nights per year.
  - The old area's cadence was not stored on the same basis, so no ratio is quoted.
- **ρ 0.01, blend 1, t_E 30 d** (the only t_E 30 d cell diagnosed):
  - 112 of 200 injections were not flagged;
  - 82 of those 112 have fewer than 3 nights inside the umbra (median over the 112: 0), so most of the event fell
    in gaps between observations.
- **ρ 0.01, blend 1, t_E 100 d:**
  - 107 were not flagged; 35 of them have fewer than 3 umbra nights.
  - The other 72 have ≥ 3 umbra nights but still were not flagged (the median over all 107 is 6). **Cause not
    determined; not stored.**
- **Reading:** in the diagnosed ρ 0.01, blend 1, t_E 30 d cell, most unflagged injections have fewer than 3 umbra
  nights. At t_E 100 d a second loss inside the pre-screen remains unexplained. The v4 cadence was not measured on
  the same basis, so the v4 → v6 drop is not shown to come from cadence. The v4 capability **did not reproduce in
  Npix 7134**; the other eight cells were not diagnosed.

## What changed relative to the old area
- **Different:**
  - comparison stars now exclude every sample member (control group only; `scripts/w3_ztf.py`
    `save_comparisons`);
  - sky position (RA 184.3, Dec +32.8, b +80.7; `independence.json`) and ZTF fields (671, 672, 1669, 1713, 1714; none shared with
  v1–v5, `independence.json`). Therefore no exposure ID is shared.
- **Shared:** the camera, DR24 processing, the 2018–2024 period, the metadata and the rest of the code.

## Sources and cost
- **Sample:** 1,072 groups (all that pass the cuts; none added), of which 533 calibration and 539 search, plus a
  control group of 5,591 objects in the comparison range (count from the session cache, not stored). 5,212 of them
  have light curves in the partition (`source_usage.jsonl`, control rows).
- **Light curves:** 12,209 rows fetched (`fetch_info.json`). The fetch includes search-group members, but no stage
  loaded them.
- **Downloads:**
  - light-curve partition 299.9 MB and objects 86.8 MB compressed, total 386.7 MB (`cost.json`);
  - TAP 4.51 MB (`exposures_summary.json`);
  - cutouts 0 MB.
- **Deviation:** a bounded TAP retry was added after a connection reset. It is network handling only; the stage was
  re-run.

## Not done (by the frozen rule)
- No search of the search group and no image check. No abundance or limit. Nothing here bears on wormholes.
