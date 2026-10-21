# W3 ZTF v26 — H-sib queue on the v22 inventory: two untested parents (frozen v13 method, unchanged pass rule)

Frozen before any v26 data read (this commit). Labels as v22 / v25: DR24 / TAP products `observed`, counts
`derived`, injections `simulated`.

## 1. Queue rule (adopted by `../v25/RESULT.md`, 6 / 6 as predicted)
- Untested parents and siblings of fade-passing parents first; siblings of fade-failing parents are not fetched.
- **Remaining qualifying v22 partitions** (`../v22/remaining.json` minus the six v25 areas): 1852 / 1854 / 1855
  (parent 463; 1853 fade 11.7 %, `../fade_diag/parent.json`) and 5052 / 5053 (parent 1263; 5054 fade 10.4 %)
  are **not fetched** under H-sib. That leaves parent 703 (2812, 2813) and parent 719 (2876–2879), both untested.
- **Head-then-siblings:** per parent, one head area runs first (highest searchable groups per GB). Its siblings are
  fetched only if the head's fade control is ≤ 100 / 1,000 (the H-sib criterion; a head failing for another reason
  with fade ≤ 10 % still releases its siblings). A sibling whose head has no verdict in this run carries over.

## 2. Prediction (report only; no queue consequence beyond section 1)
- **H-sib sibling agreement:** in each parent with ≥ 2 fade verdicts, the fade rates span ≤ 4.0 percentage points
  (ASSUMPTION; earlier parents with ≥ 2 verdicts span 1.4–3.3 points: 586 1.4, 287 1.9, 1247 2.6, 622 3.3, 339 2.7,
  406 2.5, 2031 2.0, from `../fade_diag/parent.json` and `../v25/summary.json`). A wider span is
  reported as a counter-example to sibling agreement; it changes no frozen rule.
- **What changes the plan:** a head failing fade removes its siblings from this and every later queue (saves their
  light-curve budget). Two or more counter-examples to sibling agreement across v26 and later rounds → re-open the
  queue rule (H-sib) in a frozen plan.
- **Ordinary explanations** of a sibling pattern: shared ZTF fields / CCD-quadrant seasons, shared crowding and
  light-curve noise. Not separated here.

## 3. Queue (frozen order; seed offsets fixed)
| pos | npix | parent | role | lc_mb (inventory) | searchable / GB | seed_offset |
|---|---|---|---|---|---|---|
| 1 | 2878 | 719 | head | 441.0 | 2.98 | 56000 |
| 2 | 2813 | 703 | head | 440.2 | 2.17 | 57000 |
| 3 | 2876 | 719 | sibling of 2878 | 542.3 | 2.44 | 58000 |
| 4 | 2812 | 703 | sibling of 2813 | 442.6 | 2.13 | 59000 |
| 5 | 2879 | 719 | sibling of 2878 | 726.4 | 2.04 | 60000 |
| 6 | 2877 | 719 | sibling of 2878 | 820.4 | 1.86 | 61000 |

## 4. Inputs and batch (frozen v13 method; the v25 driver plus the head gate)
- Inventory metadata: `../v22/inventory.json` (no new inventory). Dedupe list: v1–v25 samples
  (`dedupe_script.py.txt`: the v25 list plus the six v25 area samples). Samples drawn as v25 (salt `v26:`),
  committed before any light curve is read.
- Validation, pass rule, on-pass search, on-failure handling and crash rule exactly as `../v22/PLAN.md` section C.
- **Stop rules:** no new area after the per-run deadline passed to the driver (cloud routine: ≈ 25 min after start)
  or once the run's projected download passes 3,000 MB; unreached areas carry over (same seeds).
- Light-curve caches are deleted after each area (cloud disk rule); survivors keep theirs for the image check of
  `../v22/PLAN.md` section D.
