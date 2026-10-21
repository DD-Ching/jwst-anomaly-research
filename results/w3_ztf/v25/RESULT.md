# W3 ZTF v25 — sibling prediction holds in 4 / 4 decided areas (partial; 2347, 4990 carry over); 1 area passes, 2,009 groups searched, 0 survivors

- **Pre-registration:** plan committed before any data read (edeb5d0); queue, objects, dedupe and the six samples
  committed before any light curve (first fetch 05:15:14 UTC, `timeline.jsonl`).
- **Labels:** DR24 objects and light curves, IRSA TAP `observed`; counts and rates `derived`; injections `simulated`.
- **Scope:** frozen v13 rule, v22 pass criteria, unchanged 10 % fade limit. No limit; nothing here bears on
  wormholes. The sibling tally is a queue-ordering test, not a physical result.

## Sibling prediction (PLAN section 1; `summary.json` → `sibling_tally`)
| npix | parent | predicted fade | fade passing / 1,000 (Wilson 95 %) | as predicted | calibration verdict |
|---|---|---|---|---|---|
| 2346 | 586 | ≤ 10 % | 72 (5.8–9.0 %) | yes | **passed** → searched |
| 1149 | 287 | > 10 % | 121 (10.2–14.3 %) | yes | failed (fade) |
| 2345 | 586 | ≤ 10 % | 67 (5.3–8.4 %) | yes | failed (cell ρ0.1 / blend 0.7 / t_E 10 d: 0.475 < 0.5) |
| 4988 | 1247 | > 10 % | 132 (11.2–15.4 %) | yes | failed (fade) |
| 2347 | 586 | ≤ 10 % | not run (time stop) | — | — |
| 4990 | 1247 | > 10 % | not run (time stop) | — | — |

- **Partial:** 4 / 4 decided as predicted; the frozen outcome needs all 6 (≥ 5 / 6 → H-sib usable as a queue rule).
  Already, at most 1 of the remaining 2 can go against the prediction for the outcome to hold.
- **Ordinary reading:** siblings share fields, CCD-quadrant seasons and crowding; this test cannot separate them.
- 2345 shows that a passing fade control does not guarantee a calibration pass: one injection cell missed by 5 of 200.
  No retry, seed change or replacement (PLAN section 4).

## Search (2346; `2346/search_holdout_v25_time.json`)
- 2,009 groups, 3,998 series, 7 flagged runs in 1 group, **0 survivors** after the frozen chain. No image check needed.
- Cumulative with the v13 rule: 39,366 groups (37,357 + 2,009; not stored as one number). No confirmed target event.

## Cost (`summary.json`)
- Downloads ≈ 2,704 MB (objects 266 MB + four areas' light curves from parquet metadata, TAP measured); caches deleted.
- Queue wall clock 19.5 min (05:15:14–05:34:46) for 4 areas; 1 worker.

## Carry-over (next run)
- Run `objects_script.py.txt` (parents 586 / 1247 suffice), then `area_sample_script.py.txt` for 2347 and 4990 to
  rebuild their control files in the data cache; the committed `sample.json.gz` must come out identical (same salt,
  same objects; check it). Then run the driver unchanged; it skips areas that already have `decision.json`.
- Then re-run `summary_script.py.txt` for the complete tally.
- 2026-10-21 06:16Z: carry-over run started (2347, 4990), cloud routine.
