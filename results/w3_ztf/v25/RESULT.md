# W3 ZTF v25 — sibling prediction holds in 6 / 6 areas (frozen outcome: H-sib usable as a queue rule); 2 areas pass, 3,926 groups searched, 0 survivors

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
| 2347 | 586 | ≤ 10 % | 73 (5.9–9.1 %) | yes | **passed** → searched (carry-over run) |
| 4990 | 1247 | > 10 % | 124 (10.5–14.6 %) | yes | failed (fade; also cell ρ0.1 / blend 0.7 / t_E 10 d: 93 / 200) |

- **Complete:** 6 / 6 as predicted ≥ 5 / 6 → **H-sib is usable as a queue rule** (PLAN section 1): future queues put
  untested parents and siblings of fade-passing parents first and do not fetch siblings of fade-failing parents.
  The six fade rates split cleanly by parent (586: 67–73; 287 / 1247: 121–132 per 1,000).
- **Ordinary reading:** siblings share fields, CCD-quadrant seasons and crowding; this test cannot separate them.
- 2345 shows that a passing fade control does not guarantee a calibration pass: one injection cell missed by 5 of 200.
  No retry, seed change or replacement (PLAN section 4).

## Search (`<npix>/search_holdout_v25_time.json`)
- 2346: 2,009 groups, 3,998 series, 7 flagged runs in 1 group, **0 survivors** after the frozen chain.
- 2347: 1,917 groups, 3,826 series, 8 flagged runs in 4 groups, **0 survivors**. No image check needed.
- Cumulative with the v13 rule: 41,283 groups (37,357 + 2,009 + 1,917; not stored as one number). No confirmed
  target event.

## Cost (`summary.json`)
- Downloads ≈ 2,704 MB (objects 266 MB + four areas' light curves from parquet metadata, TAP measured); caches deleted.
- Queue wall clock 19.5 min (05:15:14–05:34:46) for 4 areas, 1 worker; carry-over 5.7 min (06:17:30–06:23:12) for
  2 areas, 2 workers. `summary.json` `wall_clock_min` (68.0) spans both runs including the 42-min gap between them;
  use the two queue times. `downloaded_mb` (4,205) counts the objects tables once; the carry-over read them again
  (+266 MB).

## Carry-over run (2026-10-21, 06:16–06:25Z)
- `objects_script.py.txt` re-run: `objects.json` identical. `area_sample_script.py.txt` for 2347 and 4990: decompressed
  `sample.json.gz` and `sample_info.json` identical to the committed ones (committed gzip kept). Driver unchanged
  (`W3_WORKERS=2`); it skipped the four decided areas.
- Next: build the next queue from `../v22/inventory.json` under H-sib (parents untested or fade-passing first).
