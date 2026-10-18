# W1 shear screen in PRIMER (UDS and COSMOS, program 1837, no lens model): result (2026-10-18)

Plan: `PLAN.md` (frozen before any screen). Code: the PRIMER tile entries in `scripts/inject_shear.py`
(`w1_footprints.tile_grid`, `footprints_{uds,cosmos}.json`), the catalogue and DJA zout pins in
`scripts/w1_stage_inputs.py`, tile-range streak masks in `scripts/w1_streak_mask.py`. Labels: shapes and pixels
`observed`; S, p, efficiencies and limits `derived`; injected lenses `simulated`; R and every threshold an
ASSUMPTION. Built over three cloud runs (2026-10-18, 07:10–09:50Z); inputs were re-staged and masks regenerated in
each, and every regenerated mask JSON (wall time aside) and every new summary's `screen` block equal the committed
ones (checked by `assert` in the session before copying; the checks are not stored separately).

## Gate P (`gate_p.json`)
All twenty tiles pass (photo-z share 0.592–0.873).

## Screens (`summaries/z0.3_<tile>.json`, `screen`; `screen_flags.json`)
Four tiles pass the E/B rule (p_random_max < 0.05 and S_max > max |S_×|): COSMOS x0y0 (S 4.44), UDS x0y2 (4.17),
x2y0 (4.85) and x3y0 (6.47). All four peaks lie within 20″ of a tile edge (`near_tile_edge`). The other sixteen do
not pass; they form the PRIMER part of the limit.

## The four flags (`sheet_smax_flags.png`, `vet/vet_<tile>.json`; frozen rule `w1_shear_more` addendum A)
| Tile | V1 rows | V2 rows | S without V1 ∪ V2 | p_rot | frozen outcome |
|---|---|---|---|---|---|
| primer-cosmos-x0y0 | 13 | 25 | 3.10 | 0.99 | explained |
| primer-uds-x0y2 | 0 | 27 | 0.14 | 1.0 | explained |
| primer-uds-x2y0 | 12 | 15 | 0.83 | 1.0 | explained |
| primer-uds-x3y0 | 0 | 31 | 5.11 | 0.0 | **not explained** |

**UDS x3y0** stays a flag, not a candidate (unvetted). On the sheet the peak sits on one sharp horizontal line lined
with catalogue detections. Removing its row band (V1r, post hoc, not part of the frozen outcome) gives S 2.99,
p_rot 0.99; 5 % of 200 rotated control bands drop S as far (`frac_le_no_v1r`). Visit check (post hoc,
`vet/x3y0_visit_check.json`, `vet/x3y0_visit_sheet.png`): within the 40″ cutout the combined mosaic c1002 differs
from visit o001 by at most 0.053 MJy/sr, so there it is effectively visit o001, and the visit test cannot separate
a satellite trail from sky. By eye the line ends about 15″ east of the peak with a sharp end (not a frozen test).
The decisive check needs the o001 L2 (`_cal`) exposures at the peak. Like every flagged tile, x3y0 is out of the
limit whatever its outcome (D-098).

## Limits (`w1_volume.py`; zero detections in the sixteen; Planck18)
Injections: Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5} × 3 masses × 100 lenses on each of the sixteen unflagged tiles
(`summaries/z<Z>_<tile>.json`). The twenty-four-field headline summaries reproduce
`results/w1_shear_ceers/limits_24fields_unflagged.json` exactly (checked in the session; not stored separately).

| Set | file | area (arcmin²) | n₉₅ at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (h³ Mpc⁻³) | × Takahashi & Asada |
|---|---|---|---|---|
| twenty-four fields (before) | `results/w1_shear_ceers/limits_24fields_unflagged.json` | 199.9 | 1.50 × 10⁻³ / 2.53 × 10⁻⁴ / 1.09 × 10⁻⁴ | 15.0 / 2.53 / 1.09 |
| **forty: + sixteen unflagged PRIMER tiles** | `limits_40fields_unflagged.json` | 446.2 | 9.32 × 10⁻⁴ / 1.50 × 10⁻⁴ / 6.38 × 10⁻⁵ | 9.32 / 1.50 / **0.638** |
| PRIMER only (sixteen) | `limits_primer_only.json` | 246.3 | 2.47 × 10⁻³ / 3.70 × 10⁻⁴ / 1.54 × 10⁻⁴ | 24.7 / 3.70 / 1.54 |

V_eff at 2 × 10¹³ M☉ grows from 88,969 to 151,703 Mpc³ (+71 %). Area check (`area_check.json`): the sixteen tiles'
screened grid area (246.3 arcmin²) is 0.937 of the tile area their mosaic covers (262.9 arcmin², share × tile), so
PRIMER does not overstate its area on that measure. Tiles do not overlap each other; PRIMER does not overlap the
earlier fields.

Plan check ("what would change the plan"):
- **n₉₅ at 2 × 10¹³ M☉ is below 10⁻⁴ h³ Mpc⁻³** (6.38 × 10⁻⁵, 0.638×): under the stated assumptions W1 now reaches
  the Takahashi & Asada (2013) value at that mass; the plan's next step is lower masses, not more area;
- the PRIMER tiles add 71 % to V_eff at 2 × 10¹³ M☉ (more than 20 %).

## Responsivity sensitivity (`r_sensitivity.json`, Z = 0.3, sixteen tiles)
Injections at Z = 0.3 repeated with `--fixed-r 0.231` (`summaries_r0.231/`):

| Mass (M☉) | R = 0.414 | R = 0.231 | ratio |
|---|---|---|---|
| 2 × 10¹² | 68 | 19 | 0.279 |
| 8 × 10¹² | 320 | 49 | 0.153 |
| 2 × 10¹³ | 638 | 116 | 0.182 |

Per tile and mass the ratio spans 0.0–0.667 (47 cells with a non-zero adopted count). **As for GOODS-N and CEERS,
the PRIMER contribution depends strongly on R.** A forty-field limit at R = 0.231 was not computed (only Z = 0.3 was
re-run, as the plan says); with these ratios the PRIMER V_eff would fall well below its R = 0.414 value, so the
sub-Takahashi & Asada reading holds only if R ≳ 0.414 (D-101 supports that for blank fields; not shown for PRIMER
itself).

## Limitations
- The comparison depends on R = 0.414, borrowed from a cluster fit in another program (D-099). The image
  simulations of #329 give R_sim ≥ 0.94 for blank fields, so the borrowed value is conservative under that
  simulation; the factor-two gap between cluster-fit and simulated R is unexplained.
- Takahashi & Asada (2013) is a different method and data set (a published comparison value), and both limits
  assume W1 lenses uniform in comoving volume; the ratio is not a like-for-like test.
- Efficiencies are linear between the five grid redshifts and zero outside them (conservative).
- UDS x3y0 is unresolved (above); flagged tiles (four of twenty) are excluded from the limit, which loses area but
  does not bias it.
- **Deviation:** stacking both DJA zouts for every tile of a survey instead of only straddling tiles (run 2).
