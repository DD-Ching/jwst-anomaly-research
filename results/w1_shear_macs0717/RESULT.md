# W1 shear screen (D-053) on MACS0717 with ASTRODEEP-FF photo-z: result (2026-10-17)

Plan: `PLAN.md` (committed before the runs). Labels: shapes and pixels `observed`; photo-z, ε, S maps, R, masks
`derived`; cluster shear `model_prediction` (CATS v4.1); injections `simulated`; limits `derived` from `simulated`;
every threshold an ASSUMPTION (D-053 defaults and the streak-veto settings, unchanged).

## Inputs
- F200W VENUS catalogue `jw06882-o029_t063` (sha256 `099c3c38…`, data/manifests/macs0717.ecsv), 4,953 shape rows.
- ASTRODEEP-FF (VizieR J/A+A/607/A30 `m0717cla` + `m0717clz`, accessed 2026-10-17): 4,082 rows; kept: `relFlag` = 1 or a
  spectroscopic z (`f_zbest` = 1), 2,928 rows (2,866 + 62; `zout_counts.json`, live test in `tests/test_astrodeep_zout.py`) (`scripts/astrodeep_zout.py`;
  file sha256 `b072ff2e…`, written twice, identical; pinned in `w1_stage_inputs.PHOTOZ` by its column-value hash
  `astrodeep_zout.content_sha256` `ef1a24d3…`, re-measured 2026-10-17). 855 shape rows (of 4,953) get
  a photo-z within 0.3″ on the final zout (808 before review 2; not stored); JWST − ASTRODEEP median offset (−0.03″, +0.06″) over 1,498 pairs within 1″ (not stored).
- Gaia DR3 G < 17 within 4′ of the VENUS pointing: 34 stars. Streak mask (`masks/macs0717_streak_mask.json`, 2 trails,
  19,621 binned pixels; the FITS mask is gitignored, re-make it with `w1_streak_mask.py --fields macs0717`).

## Screen (`inject_summary_macs0717.json` → `screen`)
- 2,251 weak-shear rows (`n_weak_shear`), `n_veto_used` 92, `n_streak_vetoed` 132; R = 0.229 ± 0.026 (other fields 0.41–0.48
  except MACS0647 0.23).
- S_max 3.14 at (109.36939, +37.72272), p_rot 0.28; max |S_×| 3.28.
- **D-053 rule (p_rot < 0.05 and S > max |S_×|): not passed.** Nothing for `/vet-candidate`.
- `sheet_smax.png` (F200W, looked at): the S_max aperture holds two bright cluster galaxies with diffuse/tidal light,
  a vertical image stripe through it, and two spiked stars at its rim (by eye, not measured; no causal reading).

## Injections (seed 50, 200 per mass; recovered / 200)
| |M| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
|---|---|---|---|---|---|
| MACS0717 | 0 | 0 | 0 | 7 | 9 |

Screened area 22,207 arcsec² (6.17 arcmin²). MACS0717 is nearly blind (efficiency ≤ 0.045).

## Limits: the headline stays at twelve fields
83 % of the MACS0717 shape rows have no photo-z (855 of 4,953 matched, not stored). `inject_radial.combine` puts a
field in the headline (photo-z) set when it has any photo-z file, but rows without photo-z are lensable at z_s = 2,
which biases a field's efficiency high (D-049). So MACS0717 is **not** added to the headline. Counted anyway
(`combined_13fields_veto.json`, 76.6 arcmin², optimistic for MACS0717): 4.4 × 10³ / 7.1 × 10² / 3.4 × 10² deg⁻² at
2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉, against the twelve-field headline 4.4 × 10³ / 7.2 × 10² / 3.5 × 10²
(`results/w1_streak_veto/combined_12fields_veto.json`): changes ≤ 1.5 %. Still blind at ≤ 2 × 10¹¹ M☉.
Fields of the thirteen-field file: the twelve of D-097 plus macs0717 (`--combine-only --fields` with all thirteen;
`DEFAULT_FIELDS` unchanged).

## Code-review fixes (post hoc, before the PR)
- Review 1: the first run kept `relFlag` ≠ 1 rows as NaN and floored z160 at 0 (`attach_photoz` turns it into NaN);
  the adapter dropped those rows and floored z160 at 0.01. Re-run, same seed: every number unchanged.
- Review 2: spectroscopic rows with `relFlag` = 0 (62) were dropped; they are now kept. Re-run, same seed: p_rot
  0.23 → 0.28, R 0.228 → 0.229, 8 × 10¹² M☉ recoveries 6 → 7; the numbers above are from this run. Review 2 also
  moved MACS0717 out of the headline (above).
- **Deviations from PLAN.md:** (1) spectroscopic rows with `relFlag` = 0 are kept (the plan set every `relFlag` ≠ 1
  row to NaN); (2) the planned thirteen-field combination is not the headline (post hoc criterion, no threshold;
  the twelve headline fields were not checked against it); (3) the critical-curve check of the peak was not run
  (the screen did not pass). No final zout row has zbest < 0.01, so the z160 floor never inverts a range.

## Limits of this result
- Photo-z cover only the HST WFC3 footprint; the other rows enter as in the unchanged code. The low R (0.23) is
  **not explained here**; a hypothesis, untested: unflagged cluster members without photo-z dilute the calibration.
- HST-based photo-z on JWST shapes; SIQR → 16/84 % conversion and the `relFlag` cut are ASSUMPTIONs.
- Sheet (by eye): the trail finder misses a bundle of bright parallel diagonal streaks at the mosaic's north-east
  corner (`masks/macs0717_streak_mask.png`, top left), and the vertical masked trail runs along a mosaic gap.
- p_rot rests on 200 rotation draws.
