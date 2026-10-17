# VENUS batch 2 with fixed `arc_diff`: result

Plan: `PLAN.md`, committed before box selection; `boxes.json` committed before any difference image. Owner instruction
2026-10-17 (after #300); budget 10:55Z → 14:55Z. Every number is `derived` unless labelled; injected sources are
simulations. An ordinary transient is not wormhole evidence; zero entries is not an exclusion.

## New coverage (`boxes.json`)
- **20 new boxes** of 24″ × 24″ (3.2 arcmin²) in 11 cluster fields, each pairing an earlier public NIRCam epoch with VENUS
  (6882) in F150W + F444W; footprints and ≥ 0.9999 finite pixels confirmed in all four images (metadata only).
- **Fields and pairs:** the north and south neighbours of the batch-1 core boxes in Abell 370 (1208 o002, 2022-12 →
  2026-01), Abell S1063 (3293, 2024-09 → 2025-09), Bullet (4598, 2025-01 → 2026-01), MACS1423 (1208 o010, 2023-02 →
  2026-03), MACS0416 (1176 o211, 2022-10 → 2026-01), MACS0417 (1208 o006, 2022-10 → 2025-10), MACS1149 (1176 o231,
  2024-01 → 2026-05), MS0451 (5058, 2025-01 → 2025-09), PLCK G165 (1176 o251, 2023-03 → 2026-03); one El Gordo box
  (1176 o241, 2022-07 → 2026-08); the RXJ2129 core (2767, 2022-10 → 2025-10).
- **Same rule as batch 1** with the batch-1 boxes excluded. The same 10 pairs as batch 1 rank first (depth ties
  broken alphabetically, as in #296) and took 19 boxes; RXJ2129 (pair rank 10) took the last one. The pairs of
  RXCJ0018, SPT0615, SMACS and Sunrise rank after it and were not reached (`boxes.json`: `pair_rank`, `pairs_considered`).
- MACS0416 epoch 1 is now 1176 o211 (batch 1: o212): the lowest-obs_id tie-break that #296 put in the code.
- Pipeline versions: jwst 3.0.0 except epoch 1 of Bullet and MS0451 (2.0.1) and RXJ2129's epoch-1 F150W (2.0.1).

## Calibration (per box, fixed rule; `boxes/batch.json`, `totals.json`)
- **12 boxes passed** (1.92 arcmin²): A370 ×2, Bullet ×2, MACS1423 ×2, MACS0416 2, MACS0417 ×2, MACS1149 2, MS0451 1,
  RXJ2129 1.
- **8 boxes are undeterminable** (not reported as searched with zero candidates):
  - too many unchanged-control false hits: AS1063 1 (6/19) and 2 (5/20), MACS0416 1 (5/20);
  - too few vanish recoveries: El Gordo 1 (8/14), MACS1149 1 (7/14), MS0451 2 (6/14), PLCK G165 1 (5/8),
    PLCK G165 2 (7/14).
  - PLCK G165 2's screen also did not finish within the budget (see Deviations); that box was not searched.
- **Pooled over the 12 passing boxes** (review list with sign; 68 % Wilson intervals; `totals.json`):

  | Change | 26.0 AB | 27.0 AB | 27.5 AB |
  |---|---|---|---|
  | vanish (no flux left) | 59/70 = 0.84 [0.80, 0.88] | 59/68 = 0.87 [0.82, 0.90] | 49/61 = 0.80 [0.75, 0.85] |
  | partial fade (50 % drop) | 58/69 = 0.84 [0.79, 0.88] | 45/61 = 0.74 [0.68, 0.79] | 39/71 = 0.55 [0.49, 0.61] |
  | appear (brightening) | 51/63 = 0.81 [0.76, 0.85] | 59/73 = 0.81 [0.76, 0.85] | 40/69 = 0.58 [0.52, 0.64] |
  | unchanged controls, any hit | 9/205 = 0.04 [0.03, 0.06] | | |

- Valid only for point sources injected on galaxy and arc pixels, not within 0.8″ of bright cores, not extended
  changes, not fainter than 27.5 AB. Controls per box are 11–20, so single-box decisions rest on small numbers. The
  control fraction is not a real-sky false-alarm rate.

## Screening counts (`totals.json`, `boxes/review_queue.json`)
- 19 completed screens: **203 hits → 56 not near a bright core → 11 review-list entries** (morphology pass): 6 in passing
  boxes, 5 in undeterminable boxes.
- **10 reviewed in detail** (capacity 10): sheet looked at, forced photometry in every covering public mosaic (291
  measurements), the coded class, then the #296 full rule by function (`boxes/review/final_classes.json`).
- **1 pending:** AS1063 1 h002 (undeterminable box, max|S| 6.7; `boxes/pending.json`), not cleared.

| Rank | Entry | Box | z F150W / F444W (search epochs) | Final class | Note |
|---|---|---|---|---|---|
| 1 | MACS1423 2 h000 (215.949099, +24.086948) | passed | −29.1 / −12.5 | **real change** (fading) | point source on a z ≈ 0.53 galaxy's outskirts in 8 bands in 2023; not detected as a point source in 2026 |
| 2 | MACS0417 2 h002 (64.396407, −11.907040) | passed | +13.7 / +3.2 | **real change** (brightening) | point source in VENUS (seen in F150W, F210M, F300M; F444W z +3.2), not in 2022; on lensed image 1.1 of system 1 (z = 0.8718) |
| 3 | MS0451 1 h011 | passed | −0.4 / −5.4 | undetermined | LW-only source in two 5058 visits 4 d apart, not detected 4 d later; S peaks 0.19″ apart; persistence not tested |
| 4 | MS0451 1 h015 | passed | 0.0 / −2.7 | artifact (edge) | |
| 5 | Bullet 2 h019 (104.655304, −55.945271) | passed | +5.6 / +6.1 | **real change** (brightening) | faint point source in VENUS (4 bands, S/N 5–7), not in 2025 |
| 6 | A370 1 h005 | passed | −0.2 / −1.7 | undetermined | compact-core dipole; opposite-sign ≥ 3σ differences in other mosaics (F200W −3.7, F210M −3.0, F356W +3.2, F444W +3.4) |
| 7 | MACS1149 1 h012 | undeterminable | −37.8 / −11.7 | artifact (spike) | bright star's diffraction pattern |
| 8 | PLCK G165 1 h003 | undeterminable | −1.8 / −1.3 | undetermined | compact-core dipole; F090W z −9.0 between two 4744 mosaics |
| 9 | MACS1149 1 h051 | undeterminable | +10.2 / +2.7 | artifact (spike) | linear residual along a spike |
| 10 | AS1063 1 h001 | undeterminable | −2.3 / −2.6 | no significant change | |

## The three real changes (`docs/candidates/batch_venus2_real_changes.md`)
- Sky-fixed (S peaks 0.023–0.059″ apart), point-like in ≥ 3 bands of one visit, no Gaia star within 3″, no spike on the
  sheets. No reported transient found within 1″ in the sources searched; **TNS could not be queried** (HTTP 403).
- **Status of each:** real change; an ordinary transient (e.g. a supernova) is a reasonable hypothesis; identity and
  type are not confirmed. None is unexplained, none is an exotic candidate, and "not found" is not "first discovery".
- MACS0417 h002 lies on a multiply-imaged lensed galaxy; whether its counter-images show the transient with time
  delays was **not checked**.

## Revision 2026-10-17 (follow-up, `followup/RESULT.md`; the text above is kept)
- AS1063 1 h002 (the pending entry): reviewed; artifact (edge), and no significant change in any mosaic. Pending: 0.
- MS0451 1 h011: undetermined → **artifact**: a signal fixed on one NRCALONG detector pixel in all 18 exposures of
  visits 1–2, absent at the sky position in dithers that put it elsewhere.
- MACS0417 2 h002: still a real change; the counter-image test is feasible with public RELICS models but not testable
  with current epochs (counter-images lag by ≈ 650–830 d; re-checked 2026-10-17: 708–778 d sample range of 10 realisations, best fit 748 d, `followup/macs0417_h002/recheck/RESULT.md`).
- Final classes now: 3 real changes, 2 undetermined (A370 1 h005, PLCK G165 1 h003), 5 artifacts, 1 no significant
  change; 11 of 11 review entries reviewed.

## Deviations and corrections (labelled)
- **Selection bug fix** (before selection): `batch_pairs.overlaps()` now lets boxes touch a searched box (0.01″
  tolerance); the old test wrongly rejected 13 of 40 RA-neighbours of batch-1 boxes (`overlap_fix_check.json`). Tested; no threshold changed.
- **Post hoc cutoff:** PLCK G165 2's screen and injection runs took the per-stamp fallback path (stamp pixels above
  `max_pix`, so the D-094 shared design is skipped) and ran 38 min (injection, finished) and 83 min (screen, stopped).
  At 11:51Z, before any result of that box, a 12:45Z cutoff was set; the screen was stopped then (`screen_task.json`
  rc −15). The box is undeterminable by its calibration anyway (7/14).
- **Post hoc threshold in the full rule:** "sky-fixed" = F150W/F444W S peaks ≤ 0.1″ apart. #296 judged 0.03″ sky-fixed
  and 0.49″ not but did not state a number.
- **Review script:** the #296 post-review checks were a session script; it is stored unchanged except its paths
  (`boxes/review/post_review_checks_script.py.txt`).

## Cost (session clock; `boxes/*_task.json`, `totals.json`)
- Times are from `run_log.json` (session clock). Box selection: 4.8 min (metadata and finite-pixel reads).
- Search and calibration: 11:03Z → 12:45Z wall on 4 workers; 13,732 task-seconds, of which PLCK G165 2 took 7,308
  (2,294 s injection + 5,007 s stopped screen + 8 s sample and morphology); the other 19 boxes 6,424 task-seconds; peak memory 1.41 GB per task.
- Review: about 6 min computation (approximate) (sheets, 291 forced measurements, checks), 10 sheets and stamp sets looked at,
  about 20 min in all (not stored); catalogue and literature checks by a subagent in parallel (1,122 s wall).
- Records and PR: about 30 min (not stored).

## Pending
- AS1063 1 h002 (beyond capacity). PLCK G165 2 (screen not completed; undeterminable).
- The 8 undeterminable boxes of this batch and the 8 of batch 1 need a decision on re-calibration with more
  injections. Pairs not reached: RXCJ0018, SPT0615, SMACS, Sunrise.
- The three real changes: TNS when reachable (low cost, not polled); MACS0417 h002 counter-images (needs lens-model
  image positions); spectra are not available from public data now.
- MS0451 h011: detector-frame persistence check of the LW exposures; arXiv:2601.10097 (a transient search in MS0451)
  not read.
- No follow-on batch is started.

## Regenerate
- `python scripts/batch_pairs.py --pairs results/transients/sfft_trial/venus_pairs.json --exclude … --out boxes.json`
  (exclusions in `boxes.json` `rule`).
- `python scripts/batch_run.py run --boxes boxes.json --out boxes --scratch <dir> --workers 4`, then `summarize` and
  `totals`.
- `python scripts/batch_review.py --batch boxes --boxes boxes.json --scratch <dir> --capacity 10
  --sheet-verdicts boxes/review/sheet_verdicts.json`, then the stored post-review script and the full-rule function
  (`final_classes.json`).
