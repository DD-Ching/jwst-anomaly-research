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
| primer-uds-x3y0 | 0 | 31 | 5.11 | 0.0 | **not explained** (post hoc: bright-star artefact, below) |

**UDS x3y0** stayed a flag after the frozen rule; the `_cal` check below explains it post hoc. On the sheet the peak sits on one sharp horizontal line lined
with catalogue detections. Removing its row band (V1r, post hoc, not part of the frozen outcome) gives S 2.99,
p_rot 0.99; 5 % of 200 rotated control bands drop S as far (`frac_le_no_v1r`). Visit check (post hoc,
`vet/x3y0_visit_check.json`, `vet/x3y0_visit_sheet.png`): within the 40″ cutout the combined mosaic c1002 differs
from visit o001 by at most 0.053 MJy/sr, so there it is effectively visit o001, and the visit test cannot separate
a satellite trail from sky. By eye the line ends about 15″ east of the peak with a sharp end (not a frozen test).
The decisive check needs the o001 L2 (`_cal`) exposures at the peak. Like every flagged tile, x3y0 is out of the
limit whatever its outcome (D-098).

**`_cal` exposure check (2026-10-18, `vet/x3y0_cal_plan.md` frozen before reading pixels; `scripts/w1_cal_trail.py`).**
The script searched the 22 visits of observation 001 (`--search-dirs`) and found four F200W `_cal` exposures that cover
the peak (`exposures_found` in the JSONs): visit 001 dithers 1–2 (nrca4, MJD 59960.255 / .262) and visit 011
dithers 1–2 (nrcb3, MJD 59965.999 / 59966.005). Each is resampled onto the c1002 grid.
- **Frozen outcome (`vet/x3y0_cal_check.json`, the plan's ±10″ stamps, band on the peak row):** contrasts 2.47 /
  8.05 / 4.50 / 15.72 against the 5σ call. The two "in" exposures are both dither 2, so the frozen rule returns
  **"trail: only in exposures of one dither position": x3y0 explained.**
- **The sheets do not support the trail mechanism** (`vet/x3y0_cal_sheet.png`, `vet/x3y0_wide_cal_sheet.png`).
  Visit 011 dither 1 shows the line clearly but scores 4.50, because stray-light streaks enter its control band.
  Visit 001 dither 2 shows no coherent line: its track has residual rms 2.26 px (40″ run) and 11.8 px (wide run),
  against 0.30–0.31 px for visit 011. The band-median contrast is therefore a weak per-exposure statistic, and the
  frozen "explained" holds for a different reason than the rule states.
- **Post hoc (`vet/x3y0_wide_cal_check.json`, ±70″ stamps, `--star`):**
  - The line is radial from Gaia DR3 2489648647202619264 (G 10.07, 57.6″ west of the peak, on the line's track;
    queried from the Gaia archive 2026-10-18). The line fits of both visit-011 dithers pass 6.0 / 5.3 px
    (0.19 / 0.17″) from the star, extrapolated over 1845 px.
  - In visit 011 both dithers put the line on the same sky row: fit rows at the stamp centre 0.77 / 1.26 px,
    slopes 0.0130 / 0.0129 ± 0.0003. The detector offset between the two dithers is (53.9, 17.5) px
    (`x_linear_wcs`, `y_linear_wcs`), so a detector-fixed line would move by about 17 rows (17.5 − 0.013 × 53.9).
  - At the same sky positions (stamp x 2000–2800), the visit-011 chunk peaks are 0.40–0.89 MJy/sr above the
    chunk median. Visit 001's are 0.009–0.028, with no coherent track (rms 11.8–12.1 px).
  - By eye, visit 011 also shows a fan of fainter streaks converging toward the star. The star lies 9.5–11.1″
    outside the first track chunk of the two visit-011 stamps ((x₀ − 25 − 395.6) × 0.0312″).
- **Reading (post hoc):** the line is a linear optical artefact of the G 10.07 star. It is present when the star
  sits just off nrcb3 (visit 011) and absent at the same sky position in visit 001. Whether it is a diffraction-spike
  wing or stray light is not determined. This is the row band whose removal (V1r) drops S to 2.99 (p_rot 0.99).
  **x3y0 is explained as bright-star contamination.** The frozen `_cal` rule also returns "explained", but it names
  a trail, which the sheets do not support. The earlier frozen V1 ∪ V2 outcome stays "not explained". It is not a
  candidate.
- **Why the screen's Gaia spike veto missed it:** the star lies 152″ from the tile centre (not stored), inside the
  tile's 4′ Gaia cone. Its uncapped spike law gives 291″ (`spike_law_uncapped_arcsec`). External stars are capped
  at `EXTERNAL_SPIKE_CAP` = 60″ (`scripts/exotic_screens.py`), and the peak sits at 57.6″, so most of the arm's
  sources within the aperture lie beyond the cap. The star's presence in the run's star table is inferred from the
  cone, not checked. In the mosaic the arm is traced to stamp x ≈ 2675, about 71″ from the star (chunk peak ≥ 0.1
  MJy/sr; chunk lists in the wide JSON; the distance is not stored). Visit-011 coverage ends there too, so the arm's
  full length is not measured.
- **Deviation:** the plan's third branch ("fixed in detector pixels") had no computable threshold. The detector-fixed
  test above is post hoc. A first run used 40″ stamps and a band centred one row above the peak (contrasts 1.21 /
  5.09 / 2.87 / 7.81). It was superseded by the plan's settings; its outputs were overwritten and are not stored.

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
- UDS x3y0 is explained post hoc as a G 10.07 star's optical artefact (above); flagged tiles (four of twenty) are excluded from the limit, which loses area but
  does not bias it.
- **Deviation:** stacking both DJA zouts for every tile of a survey instead of only straddling tiles (run 2).
