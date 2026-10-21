# W3 ZTF v25 — sibling-prediction queue on the v22 inventory (frozen v13 method, unchanged pass rule)

Frozen before any v25 data read (this commit). Labels as v22: DR24 / TAP products `observed`, counts `derived`,
injections `simulated`.

## 1. Question and prediction (from `../fade_diag/RESULT.md`, post hoc there; tested here)
- **Hypothesis H-sib:** the fade-control pass rate of an order-5 area is set mostly by its order-4 parent, so
  siblings of a fade-passing parent pass the fade control (≤ 100 / 1,000) more often than siblings of failing parents.
- **Predictions (frozen):** 2345, 2346, 2347 (parent 586; sibling 2344 at 8.1 %) give fade ≤ 10 %;
  1149 (parent 287: 1148 / 1150 / 1151 at 10.2–11.7 %) and 4988, 4990 (parent 1247: 4989 / 4991 at 13.0–13.8 %)
  give fade > 10 %.
- **Ordinary explanations of a sibling pattern:** shared ZTF fields / CCD-quadrant seasons, shared stellar density
  and crowding, shared light-curve noise. A pass of H-sib does not separate these; it only decides the queue rule.
- **What changes the plan:**
  - ≥ 5 of the 6 areas as predicted → future queues are ordered "untested parents and passing-parent siblings
    first"; predicted-failing siblings are not fetched (saves their LC budget).
  - ≤ 4 of 6 as predicted → H-sib is not a usable queue rule; queues stay ordered by inventory metrics.
  - Areas not reached in this run are run by the next run unchanged (same seeds), and the tally is made only when
    all 6 have a verdict. A partial tally is reported as partial.
- **Scope:** the fade limit (10 %) is unchanged and stays an owner decision; nothing here changes any frozen rule.

## 2. Inputs
- Areas and inventory metadata: `../v22/inventory.json` (no new inventory). Dedupe list: v1–v22 samples
  (`dedupe_script.py.txt`: the v22 list plus the eight v22 area samples).
- Samples drawn as v22 (`area_sample_script.py.txt`: v22 script with paths and salt `v25:`), committed before
  any light curve is read.

## 3. Queue (frozen order; seed offsets fixed)
| pos | npix | parent | prediction | lc_mb (inventory) | seed_offset |
|---|---|---|---|---|---|
| 1 | 2346 | 586 | ≤ 10 % | 725.5 | 50000 |
| 2 | 1149 | 287 | > 10 % | 270.2 | 51000 |
| 3 | 2345 | 586 | ≤ 10 % | 880.6 | 52000 |
| 4 | 4988 | 1247 | > 10 % | 531.8 | 53000 |
| 5 | 2347 | 586 | ≤ 10 % | 947.8 | 54000 |
| 6 | 4990 | 1247 | > 10 % | 536.3 | 55000 |

Order alternates prediction classes, smallest download first within a class, so a partial run still tests both.

## 4. Batch (frozen v13 method; the v22 driver, paths changed)
- Validation, pass rule, on-pass search, on-failure handling and crash rule exactly as `../v22/PLAN.md` section C.
- **Stop rules:** no new area after the per-run deadline passed to the driver (cloud routine: ≈ 25 min after start)
  or once the run's projected download passes 3,000 MB; unreached areas carry over to the next run.
- Light-curve caches are deleted after each area (cloud disk rule); survivors keep theirs for the section D image
  check of `../v22/PLAN.md`.
