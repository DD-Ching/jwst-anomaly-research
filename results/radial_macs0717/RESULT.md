# Radial (W1 exotic-lens) screen on MACS0717 with ASTRODEEP-FF photo-z: result (2026-10-17)

Plan: `PLAN.md` (committed before the run). Script: `scripts/radial_more.py --fields macs0717` (unchanged code;
`exotic_screens.cmd_radial` CLI defaults, `--max-radius` 100", Gaia DR3 `--spike-stars`, model `macs0717-cats`
= CATS v4.1). Labels: shapes `observed`; photo-z, arc classes, centres and null draws `derived`; predicted shear
axes `model_prediction`; every threshold an ASSUMPTION. Numbers: `summary.json`, `macs0717/macs0717-cats/radial*.{json,ecsv}`.

## Screen (flag rule: a centre with p_random < 0.05)
| | elongated (spike segs) | not behind lens | `anti` (not model-radial) | centres (dark) | null mean (p95) | max lines | p_random_max | flags |
|---|---|---|---|---|---|---|---|---|
| this run (photo-z, 100") | 411 (64) | 77 | 42 (37) | 4 (1) | 12.0 (21) | 5 | 0.69 | 0 |
| 2026-10-08 (no photo-z, 120"; docs/fields/macs0717.md) | 464 (65) | — | 74 (68) | 11 (8) | 16.7 (25) | 5 | >= 0.70 | 0 |

- **No centre is flagged.** The rule is not passed; nothing goes to `/vet-candidate`.
- The photo-z remove 77 elongated rows as not behind the lens; with the smaller radius (100" vs 120", a plan
  difference from the 2026-10-08 run, so the two rows are not a controlled comparison) the not-model-radial `anti`
  count falls from 68 to 37.
- `sheet_radial_top.png` (F200W, looked at): the 5-line centre lies between two bright cluster galaxies on a faint
  compact source (catalogue 764, 0.5" away); the contributing arcs are short faint segments, no dark centre. By
  eye, not measured.

## Limits of this result
- Photo-z cover only ~17 % of the shape rows (855 of 4,953, not stored; `results/w1_shear_macs0717/RESULT.md`), so
  most members outside the HST WFC3 footprint still enter.
- No injections for `radial` on MACS0717, so no detection limit is quoted.
