# W1 shear screen (D-053) on MACS0717 with ASTRODEEP-FF photo-z: result (2026-10-17)

Plan: `PLAN.md` (committed before the runs). Labels: shapes and pixels `observed`; photo-z, ε, S maps, R, masks
`derived`; cluster shear `model_prediction` (CATS v4.1); injections `simulated`; limits `derived` from `simulated`;
every threshold an ASSUMPTION (D-053 defaults and the streak-veto settings, unchanged).

## Inputs
- F200W VENUS catalogue `jw06882-o029_t063` (sha256 `099c3c38…`, data/manifests/macs0717.ecsv), 4,953 shape rows.
- ASTRODEEP-FF (VizieR J/A+A/607/A30 `m0717cla` + `m0717clz`, accessed 2026-10-17): 4,082 rows, 2,866 with
  `relFlag` = 1 (printed by the script, not stored) (`scripts/astrodeep_zout.py`; column-bytes sha256 `a3a9f660…`, not stored as a file). 808 shape rows get
  a photo-z within 0.3″ (not stored); JWST − ASTRODEEP median offset (−0.03″, +0.06″) over 1,498 pairs within 1″ (not stored).
- Gaia DR3 G < 17 within 4′ of the VENUS pointing: 34 stars. Streak mask (`masks/macs0717_streak_mask.json`, 2 trails,
  19,621 binned pixels; the FITS mask is gitignored, re-make it with `w1_streak_mask.py --fields macs0717`).

## Screen (`inject_summary_macs0717.json` → `screen`)
- 2,251 weak-shear rows (`n_weak_shear`), `n_veto_used` 92, `n_streak_vetoed` 132; R = 0.228 ± 0.026 (other fields 0.41–0.48
  except MACS0647 0.23).
- S_max 3.14 at (109.36939, +37.72272), p_rot 0.23; max |S_×| 3.28.
- **D-053 rule (p_rot < 0.05 and S > max |S_×|): not passed.** Nothing for `/vet-candidate`.
- `sheet_smax.png` (F200W, looked at): the S_max aperture holds two bright cluster galaxies with diffuse/tidal light,
  a vertical image stripe through it, and two spiked stars at its rim. Read as ordinary contamination (by eye,
  not measured).

## Injections (seed 50, 200 per mass; recovered / 200)
| |M| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
|---|---|---|---|---|---|
| MACS0717 | 0 | 0 | 0 | 6 | 9 |

Screened area 22,207 arcsec² (6.17 arcmin²). MACS0717 is nearly blind (efficiency ≤ 0.045).

## Thirteen-field limits (`combined_13fields_veto.json`; 76.6 arcmin²)
Σ(W1) < 4.4 × 10³ / 7.1 × 10² / 3.4 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (border-corrected 4.9 × 10³ /
8.0 × 10² / 3.9 × 10²), against twelve fields 4.4 × 10³ / 7.2 × 10² / 3.5 × 10² (`results/w1_streak_veto/
combined_12fields_veto.json`): changes ≤ 1.2 %. Still blind at ≤ 2 × 10¹¹ M☉.

## Limits of this result
- Photo-z cover only the HST WFC3 footprint; the other rows enter as in the unchanged code. The low R (0.23) is
  **not explained here**; a hypothesis, untested: unflagged cluster members without photo-z dilute the calibration.
- HST-based photo-z on JWST shapes; SIQR → 16/84 % conversion and the `relFlag` cut are ASSUMPTIONs.
- Sheet (by eye): the trail finder misses a bundle of bright parallel diagonal streaks at the mosaic's north-east
  corner (`masks/macs0717_streak_mask.png`, top left), and the vertical masked trail runs along a mosaic gap.
- p_rot rests on 200 rotation draws.
