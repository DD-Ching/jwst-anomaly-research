# VENUS batch 4 with fixed `arc_diff`: result

Plan `PLAN.md` (before selection); `boxes.json` committed before any difference image. Owner instruction 2026-10-18,
branch B (the same-source test failed, `../../same_source/a2744_feasibility/RESULT.md`). Every number is `derived`
unless labelled; injected sources are simulations. An ordinary transient is not wormhole evidence; zero entries is not
an exclusion; calibration-failed boxes are not reliable zero-candidate regions.

## New coverage (`boxes.json`)
- 20 new 24″ boxes in 10 cluster fields (two each): MS0451, MACS0416, MACS1149, PLCK G165, AS1063, SMACS 0723, A370,
  Bullet, MACS1423, MACS0417; same ranking as batch 3 (`--prefer-extra-epochs`). Selection 4 min 43 s (session clock).
- Not reached again: RXCJ0018, RXJ2129, SPT0615, Sunrise (each pair still has free grid positions but ranks after the
  20th box under the unchanged rule).

## Calibration (`boxes/batch.json`, `totals.json`)
- **Passed: 13 boxes (2.08 arcmin²):** MS0451 ×2, MACS0416 2, MACS1149 ×2, SMACS ×2, A370 ×2, Bullet ×2, MACS1423 2,
  MACS0417 2.
- **Undeterminable: 7** (not counted as searched): vanish recovery below 0.70 in MACS0416 1 (9/14), PLCK G165 1
  (7/14), AS1063 1 (5/14), MACS1423 1 (7/14); too many control false hits in AS1063 2 (9/20); **timed out** (screen
  cap 900 s): PLCK G165 2 (screen and injection both at their caps) and MACS0417 1 (its injection finished with 0/14
  and 7/20 false hits). The time caps removed 2 boxes this time; they took 4,060 of 8,443 task-seconds.
- Pooled over the 13 passed boxes (review list with sign; 68 % Wilson): vanish 0.88 / 0.87 / 0.85 at 26 / 27 / 27.5 AB;
  partial fade 0.91 / 0.70 / 0.53; appear 0.86 / 0.80 / 0.69; unchanged controls with any hit 11/232 = 0.047. A370 1
  and Bullet 1 passed on small numbers (6/6 and 7/7 vanish injections). Same validity limits as batches 1–3.

## Screening and review (`boxes/review_queue.json`, `boxes/review/`)
- **204 hits → 70 not near a bright core → 18 review-list entries**: 7 in passed boxes, 11 in AS1063 1 (calibration
  failed). **10 reviewed** (capacity): all 7 passed-box entries and the first 3 of AS1063 1; 8 AS1063 1 entries are
  beyond capacity (`boxes/pending.json`), not judged.
- Each review: sheet and stamps looked at (`sheet_verdicts.json`), forced photometry in every covering mosaic, coded
  class, the stored #296 post-review checks and the #304 full-rule function (`final_classes.json`).

| Rank | Entry | Box | Final class |
|---|---|---|---|
| 1 | B4-A370-box1-h000 (39.968857, −1.571911) | passed | real change by rule → **not supported by matched controls** (below) |
| 2 | B4-SJ0723.3-7327-box1-h005 (110.809928, −73.463102) | passed | real change by rule → control check below |
| 3 | B4-SJ0723.3-7327-box1-h006 (110.815491, −73.462798) | passed | artifact (edge) |
| 4 | B4-A370-box2-h004 (39.980333, −1.583316) | passed | undetermined (residual chain along a bright arc edge) |
| 5 | B4-BULLET-box1-h008 (104.639893, −55.947082) | passed | no significant change |
| 6 | B4-A370-box1-h006 (39.966252, −1.567313) | passed | undetermined (≥ 3σ in other mosaics: F200W) |
| 7 | B4-SJ0723.3-7327-box1-h010 (110.802876, −73.464924) | passed | artifact (edge) |
| 8 | B4-AS1063-box1-h006 (342.174250, −44.527317) | failed | artifact (edge) |
| 9 | B4-AS1063-box1-h014 (342.171669, −44.521964) | failed | no significant change (ring residual pattern) |
| 10 | B4-AS1063-box1-h021 (342.176219, −44.523448) | failed | no significant change (ring residual pattern) |

## The two rule-level real changes: matched-control check (`confirm/`, plan `confirm/PLAN.md` fixed first)
- **B4-A370-box1-h000:** bright point-like source (F150W ~1,030 nJy in 0.15″) present in every programme. Its 0.15″
  F150W flux is unchanged between the search epochs (ratio 0.998); the screen hit is a core residual. Against 34
  matched compact controls: F150W tail 0.529 (decisive), F444W tail 0.382 → **consistent with the matched-control
  scatter; no evidence of a source-specific change** (`confirm/a370_h000/cz_*.ecsv`). No Gaia DR3 source within 2″,
  which cannot exclude a star at this brightness (AB ≈ 23.9 in F150W, below Gaia's limit).
- **B4-SJ0723.3-7327-box1-h005:** SMACS_H005_PLACEHOLDER

## Deviations (labelled)
- Control photometry: the first run used `--query-only` and so measured only the controls sharing the entry's
  position group (14 and 16 < 20); re-run without it before any tail was seen (tool usage error, no rule change).
- The control pick used the pre-declared 30″ fallback (16″ gave 8 and 4 controls).

## Cost (session clock; `boxes/*_task.json`, `totals.json`)
- Search + calibration: 03:15:43 → 03:58 wall on 4 workers; 8,443 task-seconds; passed boxes 173–272 s each (2,860 s
  in all); undeterminable 5,580 s, of which the two timed-out boxes 4,060 s; peak memory 1.40 GB per task.
- **Per reliable area:** 8,443 task-s / 2.08 arcmin² ≈ 4,060 task-s (1.1 CPU-hours) per calibrated arcmin²; 43 min wall
  on 4 cores for 2.08 arcmin².
- Review: 03:59 → 04:17 (photometry, sheets) plus sheet look and checks; control checks 04:20 →.

## Pending (no automatic follow-up)
- 8 AS1063 1 entries beyond capacity (calibration-failed box). Held undeterminable boxes now 33 (26 + 7).
- A370 2 h004 and A370 1 h006 undetermined.
