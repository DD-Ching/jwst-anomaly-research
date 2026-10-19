# Result: forced-photometry vetting of review-list entries SMACS h039 and Sunrise o052 h005

Plan and rule: `PLAN.md` (frozen before the run). Every number is `derived` and sits in `outcome.json`,
`lightcurves.ecsv` or `ctrl_*/control_z.ecsv`; thresholds are ASSUMPTIONs.

## Frozen outcome (`outcome.json`, `scripts/review_pair_rule.py`)
| Entry | Pair | F150W Δ ± σ (nJy), z | F444W Δ ± σ (nJy), z | Outcome |
|---|---|---|---|---|
| h039 | 2736 o001 (2022-06) → 6882 o057 (2026-06) | +37.2 ± 14.7, 2.53 | +16.3 ± 8.4, 1.94 | **inconclusive** |
| h005 | 2282 o010 (2022-07-30) → 6882 o052 (2025-07-08) | −5.8 ± 9.8, −0.59 | −46.2 ± 6.6, −7.01 | **inconclusive** |

Each entry passes 2σ in one band only, so neither is `no_change` nor `change`.
**Deviation:** the rule was frozen in words in `PLAN.md`; the function and its unit test
(`tests/test_hit_lightcurve.py`) were written after the run, implementing that text unchanged.

## Stamps (looked at: `stamps_h039.png`, `stamps_h005.png`)
- h039: an isolated compact source in every band and epoch; a diagonal stripe (scattered light or a spike) crosses
  the lower part of the 2022 ERO stamps, away from the aperture. No visible change by eye.
- h005: a compact, slightly elongated galaxy core in every band; no visible change by eye.

## Post hoc: matched controls (`ctrl_h039/`, `ctrl_h005/`; `hit_lightcurve_controls.py pick` + `compare`)
Not in the plan. 40 compact peaks per entry from the epoch-A F150W mosaic (h039: 16″ box, 250–1000 nJy; h005:
30″ box, 400–2500 nJy, widened after a 16″ box gave 7), measured in the pair's mosaics with the same aperture.
The errors here are one robust value over the larger group box, so the raw z differ from the table above.

| Entry | Band | raw z | control median, robust std | z_cal | control tail fraction |
|---|---|---|---|---|---|
| h039 | F150W | 4.12 | −0.66, 8.82 | 0.54 | 0.575 |
| h039 | F444W | 2.94 | −0.74, 3.14 | 1.17 | 0.375 |
| h005 | F150W | −0.50 | 0.34, 1.67 | −0.50 | 0.600 |
| h005 | F444W | −7.11 | −2.01, 5.20 | −0.98 | 0.385 |

Against constant compact sources measured the same way, both entries sit inside the control spread
(|z_cal| ≤ 1.17; at least 37.5 % of controls are as far or farther from the median). Read narrowly: the
cross-program change of each entry is not larger than that of the matched controls. Control stamps were not kept.

Also post hoc (`lightcurves.ecsv`): within program 2282, h005 changes o010 → o120 by +108.8 nJy in F277W and
−123.8 nJy in F356W (opposite signs in adjacent bands; F090W/F115W within 6.5 / 23.5 nJy). Not interpreted.

## Conclusion
- Frozen outcome: **inconclusive** for both entries.
- Post hoc: neither changes beyond the matched-control scatter, so the review-list items are closed as
  "no change beyond control scatter (post hoc)", not as `no_change` under the frozen rule.
- Neither is a candidate. Both lie outside the arc mask (3.1″ / 9.8″) and away from catalogued images (3.71″ /
  10.46″), so they test no lensing prediction.
