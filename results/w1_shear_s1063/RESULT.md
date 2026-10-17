# W1 shear screen (D-053) on Abell S1063: result (2026-10-17)

Plan: `PLAN.md` (committed before the run). Labels: shapes `observed`; ε, S maps, R `derived`; cluster shear
`model_prediction` (CATS v4.1); injections `simulated`; limits `derived` from `simulated`; every threshold is an
ASSUMPTION (D-053 defaults, unchanged).

## Inputs (all sha256-checked against SOURCES.md / `data/manifests/abells1063.ecsv`)
- F200W GLIMPSE pipeline catalogue `jw03293-o001_t001` (4,238 rows).
- DJA v7.5 eazy zout (`abells1063-grizli-v7.5-fix.eazypy.zout.fits`, from the 75 MB tarball; tarball deleted).
- Gaia DR3 G < 17 within 4′ of (342.1832, −44.5309): 11 stars (`scripts/gaia_stars.py`), as in
  docs/fields/abells1063.md; sha256 `79330a36…` (`screen_shear.json`). The injection run used a copy of the same file
  (sha256 checked equal in the session; not stored).
- Command: `exotic_screens.py --model abells1063-cats shear --catalog <F200W cat> --photoz <zout> --spike-stars
  <gaia> --max-radius 100`; injections `inject_shear.py --fields abells1063` (seed 50, 200 per mass), with the
  shear-only input override `SHEAR_FIELD_OVERRIDES` (the radial FIELDS entry has no photo-z and stays unchanged).

## Screen (`screen_shear.json`, `shear_peaks.ecsv`)
- 80 spike segments dropped; 1,440 weak-shear sources (0.07 per arcsec², not stored); R = 0.472 ± 0.034 (in the 0.41–0.48 range
  of the other fields).
- S_max 3.61 at (342.19413, −44.53478), p_rot 0.07; max |S_×| 3.33 (p 0.29); second peak 3.52 (p 0.105).
- **D-053 rule (p_rot < 0.05 and S > max |S_×|): not passed. Null.** Nothing for `/vet-candidate`.
- Post hoc (not in the plan): `sheet_top_peaks.png` (F200W, S3 byte-range read, looked at): both top-peak apertures are crossed by long, straight
  streaks (several parallel ones at the maximum, one at the second peak) and hold cluster members. Linear structures
  pointing through an aperture centre are the spike-type mimic of the W1 sign; the Gaia veto does not cover streaks
  from stars outside 4′ or G < 17 cuts (by-eye reading, not measured).

## Injections (`inject_summary_abells1063.json`; recovered / 200)
| |M| (M☉) | 2 × 10¹⁰ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ |
|---|---|---|---|---|---|
| Abell S1063 | 0 | 0 | 17 | 69 | 91 |

Screened area 5.62 arcmin².

## Five-field limits (`combined_5fields.json`; 35.8 arcmin²)
Post hoc (not in the plan): the four earlier fields were re-run in this session with the same code and seeds
  (`rerun_4fields/`): recoveries (14/7/3/9, 72/50/30/45,
146/112/74/73) and four-field limits reproduce docs/exotic_limits.md exactly.
- Σ(W1) < 5.8 × 10³ / 1.1 × 10³ / 5.7 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (was 7.7 × 10³ / 1.3 × 10³ /
  6.5 × 10²); border-corrected 6.3 × 10³ / 1.2 × 10³ / 6.2 × 10². Still blind at ≤ 2 × 10¹¹ M☉.

## Limits of this result
- Same caveats as D-053 (catalogue isophotal moments, no PSF-anisotropy model, one B map per field).
- The DJA photo-z puts 13 of 46 matched `arcs.txt` images at z < 0.6 (docs/fields/abells1063.md), so some background
  sources are dropped as foreground: conservative for the limit.
