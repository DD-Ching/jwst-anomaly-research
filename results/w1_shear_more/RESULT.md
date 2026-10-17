# W1 shear screen (D-053) on MACS0417, MACS1423 and WHL0137: result (2026-10-17)

Plan: `PLAN.md` (committed before the run; addendum A committed before the vetting computation). Labels: shapes
`observed`; ε, S maps, R `derived`; cluster shear `model_prediction` (CANUCS DR1 Lenstool maps for MACS0417 and
MACS1423, RELICS Lenstool v1 for WHL0137); injections `simulated`; limits `derived` from `simulated`; every
threshold is an ASSUMPTION (D-053 defaults, unchanged).

## Inputs (sha256 match `results/arcs/*/RESULT.md`; not stored)
F200W pipeline catalogues `jw01208-o006_t003` (MACS0417), `jw01208-o010_t005` (MACS1423), `jw02282-o010_t001`
(WHL0137); CANUCS DR1 `macs0417-clu` / `macs1423-clu` photometry (converted with `inject_radial.canucs_zout`); DJA
v7.5 sunrise zout (streamed out of its tarball). Gaia DR3 G < 17 within 4′ of each model reference point: 13 / 9 / 9
stars. Fields enter `inject_shear.py` as shear-only entries (`SHEAR_FIELD_OVERRIDES`; the radial screen's FIELDS is
unchanged). Commands: `exotic_screens.py --model <model> shear ... --max-radius 100` (`screen_shear_*.json`) and
`inject_shear.py --fields <field>` (seed 50, 200 per mass; `inject_summary_*.json`).

## Screen
| Field | sources used | R | S_max (p_rot) | max S_× / max −S_× (p_rot) | D-053 rule |
|---|---|---|---|---|---|
| MACS0417 | 2,279 | 0.488 ± 0.034 | 3.56 (0.065) | 3.24 (0.34) / 3.20 (0.415) | not passed |
| MACS1423 | 2,268 | 0.604 ± 0.033 | 3.88 (0.01) | 3.64 (0.065) / 3.23 (0.59) | **passed** |
| WHL0137 | 340 | 0.533 ± 0.057 | 2.75 (0.72) | 3.02 (0.215) / 3.18 (0.03) | not passed |

- Known-case check: R is positive and well above zero in all three (the cluster shear is seen). MACS1423's 0.604 lies
  above the 0.41–0.49 of the seven earlier fields; PLAN's "far outside" was not a computable condition.
  **Deviation:** MACS1423 is left out of the headline limits (below; the plan said null fields replace the
  seven-field limits) and shown only as an alternative, because of this R and because its hit is removed only by a
  test not run through injections.
- WHL0137 keeps 340 of 1,340 resolved rows after the |g| < 0.5 and κ < 1 cuts; the cause was not measured.
- Sheets looked at (F200W 30″; `sheet_macs1423_smax.png`, `sheet_macs0417_whl0137_smax.png`): the MACS0417 and
  WHL0137 maxima sit at a mosaic edge (WHL0137's aperture is about half off the image). By eye, not measured.

## MACS1423 hit: vetting (addendum A, `scripts/w1_peak_vet.py`, `vet_macs1423.json`)
The sheet shows a vertical stripe through the aperture centre and an m = 18.2 source (not a Gaia G < 17 star)
8.9″ from the peak with diffraction spikes. **Post hoc** test (designed after seeing the sheet, frozen before computing): V1 removes the 15 used sources
within 1″ of the stripe columns (61 columns, 442–502, about 1.9″ wide, found by the 5-MAD column-median rule;
zero-valued fill is ignored and a zero MAD flags nothing, guards added after code review); V2 removes the 5 used sources within
15″ of the bright source whose axis points at it (±15°).

| Removed | S at the peak | p_rot (same 200 null draws) |
|---|---|---|
| none | 3.88 | 0.01 |
| V1 (stripe) | −0.07 | 1.0 |
| V2 (star) | 3.90 | 0.01 |
| V1 ∪ V2 | 0.14 | 1.0 |

**Frozen outcome: explained by streak/spike contamination (screen-level).** The 15 sources on the stripe carry
the signal (S 3.88 → −0.07 without them); the star's spikes do not (3.88 → 3.90). The stripe's origin (detector
artefact, trail, scattered light or mosaic seam) was not determined. Not a candidate; nothing for `/vet-candidate`. This is the second
field (after Abell S1063) where a linear streak mimics the W1 sign.

## Injections (`inject_summary_*.json`; recovered / 200)
| \|M\| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
|---|---|---|---|---|---|
| MACS0417 | 0 | 0 | 5 | 51 | 101 |
| MACS1423 | 0 | 0 | 9 | 72 | 140 |
| WHL0137 | 0 | 0 | 1 | 12 | 44 |

## Limits (`combined_9fields.json`, `combined_10fields.json`)
- Check: `--combine-only` on the seven earlier summaries reproduces `results/w1_shear_smacs_elgordo/combined_7fields.json`
  exactly (not stored; equality checked in the session).
- **Headline, nine fields (MACS1423 left out), 55.1 arcmin²:** Σ(W1) < 4.9 × 10³ / 8.2 × 10² / 4.1 × 10² deg⁻² at
  2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (border-corrected 5.3 × 10³ / 9.0 × 10² / 4.6 × 10²); was 5.3 × 10³ / 9.6 × 10² /
  5.0 × 10² from seven fields.
- Alternative, ten fields (MACS1423 counted as null after the streak test), 61.6 arcmin²: < 4.3 × 10³ / 7.0 × 10² /
  3.5 × 10² deg⁻² (border-corrected 4.7 × 10³ / 7.7 × 10² / 3.9 × 10²). Using it needs a streak veto run through
  the injections too (scripts/CLAUDE.md: a vetting test must be calibrated).
- Still blind at ≤ 2 × 10¹¹ M☉.

## Limits of this result
- Same caveats as D-053 (catalogue isophotal moments, no PSF-anisotropy model, one B map per field).
- The CANUCS maps are preliminary (MACS1423: 5 spectroscopic systems); no frame offset to the JWST catalogue was
  measured for any of the three models.
- The streak test is a post hoc, screen-level vetting step (its rule was frozen before computing); it was not run
  through injections. A first version compared column medians with the median of column medians instead of the
  frozen "cutout median" (19 sources, S 0.14 without them, same outcome); the committed JSON uses the frozen rule.
