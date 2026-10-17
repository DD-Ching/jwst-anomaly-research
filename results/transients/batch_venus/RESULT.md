# Bounded VENUS batch search with fixed `arc_diff`: result

Plan: `PLAN.md`, committed before box selection. `boxes.json` was committed before any difference image. Owner
instruction 2026-10-17; budget 06:01Z → 10:01Z, used 06:01–06:45Z. Every number is `derived` unless labelled;
injected sources are simulations. An ordinary fading is not wormhole evidence; zero entries is not an exclusion.

## What was searched (new data)
- **Data:** 20 boxes of 24″ × 24″, 3.2 arcmin² in all, in 10 cluster fields. Each pairs an earlier public NIRCam
  epoch with VENUS (6882) in F150W + F444W (`boxes.json`).
- **Boxes per field:** two each, at the cluster core and its west neighbour (east for El Gordo):
  - Abell 370 (1208 o002, 2022-12 → 2026-01);
  - Abell S1063 (3293, 2024-09 → 2025-09);
  - Bullet (4598, 2025-01 → 2026-01);
  - El Gordo (1176 o241, 2022-07 → 2026-08; two new boxes beside B1–B3);
  - MACS1423 (1208 o010, 2023-02 → 2026-03);
  - MACS0416 (1176 o212, 2022-12 → 2026-01);
  - MACS0417 (1208 o006, 2022-10 → 2025-10);
  - MACS1149 (1176 o231, 2024-01 → 2026-05);
  - MS0451 (5058, 2025-02 → 2025-09);
  - PLCK G165 (1176 o251, 2023-03 → 2026-03).
- **Pipeline versions:** jwst 3.0.0 in every epoch, except Bullet and MS0451 epoch 1 (2.0.1).
- **Selection:** coverage, finite fraction ≥ 0.95 in all four images, and depth only.
  - Every pair's shallower F150W exposure is the same VENUS 1,245 s, so the depth ranking ties. The stable sort then
    took pairs in alphabetical order of target. This tie-break is labelled **post hoc**; it uses no difference
    information.
  - RXCJ0018, RXJ2129, SPT0615, SMACS and Sunrise fell beyond the 20-box cap and were not searched.
- **Already searched and excluded:** El Gordo B1–B3, the Sunrise box and the SMACS system-5 box (3 grid boxes
  dropped).

## Calibration (per box, fixed rule; `boxes/batch.json`, `boxes/*/inj_eval.json`)
- **12 boxes passed:** vanish 26.0 + 27.0 AB on the review list with sign ≥ 0.70, and ≤ 2 unchanged-control false
  hits. That is 1.92 arcmin².
  - A370 ×2, Bullet ×2, El Gordo 1, MACS1423 2, MACS0416 2, MACS0417 2, MACS1149 ×2, MS0451 2, PLCK G165 2.
- **8 boxes are undeterminable:**
  - too many control false hits: AS1063 ×2 (4/10 and 12/20), El Gordo 2 (3/20), MACS0416 1 (4/17);
  - too few vanish recoveries: MACS1423 1 (3/8), MACS0417 1 (3/7), MS0451 1 (4/10), PLCK G165 1 (9/14).
  - They are not reported as searched with zero candidates.
- **Sample size:** controls per box range 9–20 (fewer positions where the box has few bright pixels), so single-box
  pass decisions rest on small numbers.
- **Pooled over the 12 passing boxes** (review list with sign; 68 % Wilson intervals; `totals.json`):

  | Change | 26.0 AB | 27.0 AB | 27.5 AB |
  |---|---|---|---|
  | vanish (all flux gone) | 62/71 = 0.87 [0.83, 0.91] | 56/68 = 0.82 [0.77, 0.87] | 50/63 = 0.79 [0.74, 0.84] |
  | partial fade (50 % drop) | 54/68 = 0.79 [0.74, 0.84] | 47/63 = 0.75 [0.69, 0.80] | 37/71 = 0.52 [0.46, 0.58] |
  | appear (brightening) | 48/63 = 0.76 | 59/73 = 0.81 | 55/70 = 0.79 |
  | unchanged controls, any hit | 8/207 = 0.04 [0.03, 0.05] | | |

- **Where these efficiencies apply:**
  - point sources injected on galaxy and arc pixels;
  - not within 0.8″ of bright cores, which are flagged and not reviewed;
  - not extended changes;
  - not fainter than 27.5 AB.
  - The unchanged-control fraction is not a real-sky false-alarm rate.

## Screen and review (`boxes/review_queue.json`, `boxes/review/`)
- **Flags:** 113 hits in the 20 screens; 13 not near a bright core; **5 review-list entries.**
  - 3 are in passing boxes, 2 in undeterminable ones.
  - All 5 were reviewed in detail (capacity 10). Pending: 0.
- **Each detailed review:** the contact sheet, looked at; forced photometry in every covering public mosaic (180
  measurements); the coded class rule.

| Rank | Entry | Box status | z F150W / F444W (search epochs) | Class | Note |
|---|---|---|---|---|---|
| 1 | PLCK G165 2, h002 (171.805899, +42.481802) | passed | −27.1 / −13.4 | **real change** | point source in all 8 PEARLS bands on 2023-03-30 (S/N 14–41) and on 2023-04-22 (programme 4446, F090W 271 → 173 nJy, F150W 240 → 214 nJy); absent in 2025 (4744) and 2026 (VENUS) |
| 2 | MACS0416 2, h003 | passed | −0.35 / −2.87 | no significant change | F444W ridge along a bright elongated arc |
| 3 | MACS1423 2, h006 | passed | −0.89 / −1.03 | no significant change | weak S on a compact source |
| 4 | MACS0416 1, h004 | undeterminable | −1.10 / −2.82 | no significant change | F444W ridge along a thin arc |
| 5 | AS1063 1, h005 | undeterminable | −0.78 / −1.56 | no significant change | weak S on a compact source |

- **h002 has an ordinary explanation; it is not an unexplained change.**
  - Its multi-epoch light curve is that of a transient: present over at least 23 days in 2023 and fading in the blue
    bands, then gone by 2025. That is the light curve of a supernova or similar ordinary transient.
  - It lies about 1.2″ from a compact galaxy; no host is visible directly under it.
  - It is about 32″ from the published SN H0pe position in this cluster (11:27:15.60 +42:28:33.7). Whether it is one
    of SN H0pe's images or another transient is **not established**. The per-image coordinates were not found, and
    no transient-catalogue cross-match was done.
  - It is not sent to `/vet-candidate`: it is a catalogue-identification task, not an unexplained signal.
  - It is a fading recovered by the search, as the calibration says the search can. It is not wormhole evidence.
- **Not reviewed (by design):** 108 hits that are near bright cores or fail morphology. They are saved in each
  box's `screen/hits.ecsv` and `morph_all.json`, not deleted.

## Time used
- **Box selection** (MAST metadata, footprints, finite fractions, headers): about 5 min.
- **Computation:**
  - 40 `arc_diff` runs plus samples and morphology, 3,849 task-seconds;
  - 22 min of wall time on 4 workers (06:10–06:32Z);
  - peak memory 1.53 GB per task.
- **Review:** 2.4 min of computation (sheets and forced photometry), plus looking at 5 sheets and 1 stamp set and a
  literature check, about 10 min.
- **Records:** about 15 min.

## Pending
- Identify h002 in a transient catalogue (TNS) or against SN H0pe's published image positions.
- Boxes not searched in this batch: the undeterminable 8, which need a decision on whether to re-calibrate with more
  injections, and pairs beyond the cap.
- No follow-on batch is started.

## Regenerate
- `python scripts/batch_pairs.py --pairs ../sfft_trial/venus_pairs.json --exclude … --out boxes.json`
  (the exclusions are listed in `boxes.json` `rule`).
- `python scripts/batch_run.py run --boxes boxes.json --out boxes --scratch <dir> --workers 4`.
- `python scripts/batch_review.py --batch boxes --boxes boxes.json --scratch <dir> --capacity 10
  --sheet-verdicts boxes/review/sheet_verdicts.json`.
