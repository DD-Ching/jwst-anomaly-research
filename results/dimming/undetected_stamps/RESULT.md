# W3 dimming: the NEXUS `undetected` control pairs are visible sources lost to no-data aperture pixels (2026-10-20)

Plan frozen first (`PLAN.md`). `stamps.py` range-read 31 × 31 px SCI + DQ stamps of every NaN-flux exposure of
the NEXUS monitored pairs in e4–e7 (plus flag d30295) from S3 in 23 s (0 MB stored; `out/reasons.ecsv`,
contact sheet `out/stamps.png`). `score.py` scored the committed `persistence_run/nexus` tables
(`score.json`). Provenance: derived; thresholds are ASSUMPTIONs.

## Stamps (`score.json` `stamps`)
- The 13 monitored `undetected` pairs (11 held out, 2 train) have 24 covering exposures, all NaN flux: 22
  `aperture_no_data`, 2 `edge` (d45375 e5 and d45498 e6, centre within 3 px of the top row), 0 `not_reproduced`.
  In 21 of the 22 `aperture_no_data` exposures only 1–2 of the 13–29 aperture pixels are NaN; every NaN aperture
  pixel carries DQ DO_NOT_USE (`reasons.ecsv` `dq_bits_nan`). The median NaN fraction of these stamps is 0.020.
- By eye (contact sheet): every one of the 13 shows a compact source at the centre; the NaN pixels are scattered
  isolated pixels, some next to the core. They are on-sky sources that the D-039 test cannot measure, not
  non-detections.
- Over all monitored NaN-flux exposures in e4–e7 (76): 71 `aperture_no_data` (66 with ≤ 2 NaN pixels), 5 `edge`.
- Flag d30295 (e2) shows no source at the centre in its stamps, consistent with its `undetected` verdict; it
  stays vetoed as before.

## Exposure-count efficiency, measurable pairs only (`score.json`; ≥ 2 finite-flux exposures, mosaic S/N ≥ 10)
| split | class | all pairs `on_sky` / n (Wilson 95 %) | measurable `on_sky` / n (Wilson 95 %) | measurable fraction |
|---|---|---|---|---|
| train | ≤ 4 | 49 / 96 (0.412–0.608) | 49 / 62 (0.674–0.873) | 0.646 |
| train | > 4 | 94 / 94 (0.961–1.000) | 94 / 94 (0.961–1.000) | 1.000 |
| held out | ≤ 4 | 33 / 89 (0.278–0.475) | 33 / 45 (0.590–0.840) | 0.506 |

- **Frozen decision: R2.** The held-out ≤ 4 measurable lower bound (0.590) is below 0.80, so the deficit is not
  only measurability; the remaining 12 measurable held-out pairs are `inconclusive` (not split further here).
- Measurability alone moves the held-out ≤ 4 rate from 0.37 to 0.73: about half of the ≤ 4 deficit comes from
  NaN apertures.
- **Plan trigger (PLAN.md "What would change the plan"):** real, visible sources with only isolated flagged
  pixels in the aperture. By the plan, the fix is a NaN-tolerant aperture in `persistence_check`, not a new
  efficiency model. Not built here.

## What this changes
- Nothing in the current limits: `swap_inject` already counts these pairs as misses (ε = 0.412 for ≤ 4), which is
  conservative.
- A NaN-tolerant aperture (for example: accept ≤ 2 NaN aperture pixels, area-corrected flux, the same S/N rule)
  would raise the ≤ 4 efficiency and the swap gain. It changes the D-039 test, so it needs a frozen plan, a unit
  test, a re-run of the controls and the 37 targets in all three fields, and a re-scored efficiency before any
  limit uses it.

## Next
1. Frozen plan + NaN-tolerant aperture in `persistence_check.measure_pixel`; re-run NEXUS e1–e7 controls first
   (120-file epochs, about 2 min each) and re-score.
2. Split the 12 measurable held-out `inconclusive` pairs by cause (first exposure of a visit has `n_prior` = 0;
   one of two exposures below S/N 5).
3. NEXUS e2 controls (unchanged).
