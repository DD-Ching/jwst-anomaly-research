# SMACS ERO o001 → VENUS o057: tile-median registration check — result

Plan: `PLAN.md` (written before the run). Command: `summary.json` (the stored pair-1 command plus
`--reg-tile-arcsec 6`). Labels: shifts, S values and efficiencies are `derived`; the 0.5 px threshold is an
ASSUMPTION. The sheet (`sheet.jpg`: the four morph passes, then the old h063 and h039 positions as control rows)
was looked at before writing this.

## Registration
| Band | whole box (#237 run, `../o001_o057`) | tile median (6″ tiles) | difference |
|---|---|---|---|
| F150W | (2.11, −1.82) px | (0.20, −0.40) px, 15 tiles | 2.4 px ≈ 0.075″ |
| F444W | (1.04, −0.81) px | (0.28, −0.09) px, 16 tiles | 1.0 px ≈ 0.066″ |

The whole-box peak was pulled by a few tiles, as in El Gordo (#288), here without hitting the ±3 px window edge.
The difference exceeds the 0.5 px threshold, so the hit lists were compared (PLAN).

## Screen (frozen rule, `scripts/smacs_arcdiff_images.py`)
- **Positive control:** Atalanta A = h024, B = h025, both morph `pass` and `known` (≤ 0.021″); clean point
  brightenings in both bands on the sheet.
- Hits 73 → 76 (not `near_bright_core` 27 → 28). All-hits morph: pass 5 → 4, dipole 11 → 6, extended 12 → 23.
- Injection recovery (whole box) F150W 0.85 / 0.64 / 0.37 and F444W 0.82 / 0.65 / 0.29 at 27.5 / 28.0 / 28.5 AB
  (was 0.82 / 0.60 / 0.37 and 0.76 / 0.62 / 0.30).
- **Frozen rule: 0 `vet`** in both the default and the all-hits runs (`image_distances.json`,
  `image_distances_all.json`); the post hoc 2.0″ radius holds only Atalanta.
- New `off_image` passes: h073 (2.42″ from 6.2) and h015 (2.79″ from 18.2, all-hits only). On the sheet both are
  mixed-sign residuals on bright compact cores, the subtraction-residual pattern.
- The earlier passes h019, h039 and h063 are not passes now (dipole, rosette, and no hit within 4″ for h063); on
  the sheet h063's position shows no residual and h039 a weak dipole. h063 had no forced-photometry change ≥ 2σ
  (`../RESULT.md`), so nothing is lost.

## Outcome
The pending item from `../../fade_calib/RESULT.md` is closed: the registration changes by 1–2.4 px, the screen's
outcome does not (Atalanta recovered, 0 `vet`, no new pass within 2″ of a catalogued image). `../RESULT.md`
pair 1 stands; h061 (`edge` before) has no hit within 1″ now; it stays open as recorded there (its question is the 2024 F444W excess, outside this pair). Exploratory: one box, one pair.
