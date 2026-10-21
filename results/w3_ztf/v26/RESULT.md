# W3 ZTF v26 — both parent heads fail the fade control (11.7 %, 14.3 %); H-sib skips their four siblings; the v22 inventory is exhausted; 0 groups searched

- **Pre-registration:** plan committed before any data read (318d678, 07:12:13 UTC); queue and scripts (07:13),
  objects, dedupe and the six samples committed before any light curve (07:15:14; first fetch 07:15:30,
  `timeline.jsonl`).
- **Labels:** DR24 objects and light curves, IRSA TAP `observed`; counts and rates `derived`; injections `simulated`.
- **Scope:** frozen v13 rule, v22 pass criteria, unchanged 10 % fade limit. No search, no limit; nothing here bears
  on wormholes. The head gate is a queue rule, not a physical result.

## Heads (PLAN sections 1, 3; `<npix>/decision.json`)
| npix | parent | fade passing / 1,000 (Wilson 95 %) | lowest injection cell | calibration verdict | siblings |
|---|---|---|---|---|---|
| 2878 | 719 | 117 (9.9–13.8 %) | ρ0.1 / blend 0.7 / t_E 10 d: 87 / 200 | failed (fade, cell) | 2876, 2879, 2877 skipped |
| 2813 | 703 | 143 (12.3–16.6 %) | ρ0.1 / blend 0.7 / t_E 10 d: 92 / 200 (all 10 cells 0.46–0.65) | failed (fade, cells) | 2812 skipped |

- PSPL, eclipse and common-epoch controls 0 / 200 in both; 0 calibration survivors. No retry, seed change or
  replacement (v22 PLAN section C).
- **Sibling agreement (PLAN section 2):** not tested; no parent reached two fade verdicts, because the head gate
  stopped both parents after one area. This is the rule working as written, so the check can only be made on a
  parent whose head passes.
- **Ordinary reading:** two more sky regions with fade-control rates above 10 %, as five of the 25 earlier areas'
  parents. The cause (fields, CCD-quadrant seasons, crowding, noise) stays undetermined (`../fade_diag/RESULT.md`).

## Supply
- With 1852 / 1854 / 1855 (parent 463) and 5052 / 5053 (parent 1263) excluded under H-sib, **no qualifying v22
  partition is left to run.** The next W3 ZTF batch needs a new-area inventory.
- Cumulative with the v13 rule unchanged: 41,283 groups searched, no confirmed target event.

## Cost (`summary.json`)
- Downloads ≈ 1,066 MB (objects 171.5 MB + two areas' light curves from parquet metadata + TAP); caches deleted.
  Skipped siblings would have cost 2,532 MB of light curves (inventory `lc_mb`).
- Queue wall clock 3.4 min (07:15:30–07:18:51), 2 workers.
- **Plan wording slip:** the PLAN section 3 column "searchable / GB" is in thousands of groups per GB
  (`queue.json` stores the per-GB values, e.g. 2,977 for 2878); the order is unaffected.
- **Driver fix after the run (review):** a crash in an area job now also marks its head as undecided, so a waiting
  sibling carries over instead of waiting forever. No stage crashed in this run, so the result is unaffected.
