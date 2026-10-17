# Radial (W1 exotic-lens) discovery screen on the shear-only fields: result (2026-10-17)

Plan: `PLAN.md` (committed before any run). Script: `scripts/radial_more.py` (inputs `inject_shear.field_spec`,
screen `exotic_screens.cmd_radial` with the CLI defaults, `--max-radius` 100″, Gaia DR3 `--spike-stars`). Labels:
catalogue shapes `observed`; arc classes, centres and null draws `derived` (model_prediction for the predicted
shear axes); every threshold an ASSUMPTION. Per-field `<field>/radial.json`; table `summary.json`.

## Known case first
Abell S1063 re-run with the docs/fields/abells1063.md headline inputs (`abells1063_reproduce`): 436 elongated (26
spike segments), 51 `anti`, 46 not model-radial, 10 centres, null mean 8.2 (p95 13), max 4 lines, p_random_max
0.435: every value in that table's "no photo-z" row is reproduced.

## Screen (flag rule: a centre with p_random < 0.05, `radial_more.flagged`)
| Field | model | elongated (spike segs) | not behind lens | `anti` (not model-radial) | centres | null mean (p95) | max lines | p_random_max | flags |
|---|---|---|---|---|---|---|---|---|---|
| MACS0417 | `macs0417-canucs-maps` | 623 (58) | 88 | 71 (69) | 17 | 19.7 (29) | 5 | 0.17 | 0 |
| MACS1423 | `macs1423-canucs-maps` | 590 (54) | 112 | 41 (41) | 1 | 2.9 (7) | 3 | 0.915 | 0 |
| WHL0137 | `whl0137-relics-lenstool` | 260 (8) | 73 | 27 (26) | 2 | 1.4 (4) | 3 | 0.705 | 0 |
| RXJ2129 | `rxj2129-lenstool-maps` | 209 (14) | 14 | 21 (21) | 6 | 2.5 (7) | 3 | 0.825 | 0 |
| MACS0647 | `macs0647-clash-nfw` | 257 (29) | 63 | 36 (36) | 2 | 5.6 (11) | 6 | 0.075 | 0 |

- **MACS0417, MACS1423:** no centre is flagged (models passed the `arcs` alignment control).
- **WHL0137:** no flag. It repeats the 2026-10-08 Sunrise run (D-036 evidence: 29 not model-radial, max 3, p 0.885)
  with the Gaia spike veto added (9 stars, 8 segments dropped); the numbers stay close.
- **RXJ2129, MACS0647: no result.** Their models fail the `arcs` alignment control (PLAN), so the `anti` class is not
  trusted; the runs are kept for completeness. MACS0647's strongest centre (6 lines, p 0.075, not dark: 0.57″ from
  source 803) is not flagged either.

## Contact sheet (`sheet_radial_top.png`, looked at)
Each field's strongest centre in a 36″ F200W stamp, with the major axes of the arcs whose lines pass it. By eye:
- MACS0647: five of the six arcs lie along one straight faint linear feature running from the centre to the lower
  right, parallel to the bright streak in the stamp: segments of a streak, not galaxies (W1 noted streaks at this
  field's shear peaks too, `results/w1_streak_veto`). Not stored as a number; the field gives no result anyway.
- WHL0137: one of the three arcs lies on a diffraction spike of the bright star at the stamp's right edge.
- MACS0417, MACS1423, RXJ2129: the contributing arcs are galaxies scattered around the centre; nothing aligned.

## Reading
No exotic-lens candidate in the three control-passing fields; nothing for `/vet-candidate`. The `radial` screen has
now run on all twelve W1 shear fields (ten with a control-passing model). No injections were run on the new fields,
so the D-049 `radial` limits (eight fields) are unchanged; the shear screen's limits are 9–14× stronger
(docs/exotic_limits.md), so radial injections here would not move the headline.

## Limits of this result
- The spike veto does not remove streaks or every spike segment (WHL0137, MACS0647 above).
- Photo-z from CANUCS DR1 (MACS0417, MACS1423) and DJA (others); sources without a photo-z use the z grid.
- The MACS0417 map file failed its sha256 on the first fetch (truncated download) and passed on the retry.
- `radial.json` paths point to the session's scratch directory; inputs are recovered with
  `scripts/w1_stage_inputs.py --fields <field>`.
