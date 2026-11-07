# SMACS 0723: 2736 NIRISS / MIRI coverage of the open records and Noirot+2023 grism redshifts

Script: `scripts/smacs_niriss_miri_coverage.py` (reuses the polygon geometry of `scripts/grism_4043_coverage.py`;
matching unit-tested in `tests/test_smacs_niriss_miri_coverage.py`). Output: `summary.json`. No plan file: this
is a coverage and catalogue lookup with no outcome rule; there are no images, so no contact sheet.

## Questions
1. TASKS Later "NIRISS imaging of program 2736 (F115W, F200W)": which open SMACS records does it cover?
2. Could a second instrument at the 2022 epoch test the 2736 NIRCam position of `f200w_384` (the remaining
   "2736 systematic" explanation in `results/smacs_jwst3/epoch3_384/RESULT.md`)?
3. Open question of `docs/candidates/jw02736-o001_t001_nircam_f200w_2925.md`: is its NED redshift spectroscopic?

## Data
- Observed: MAST level-3 `s_region` footprints of 2736 o003 NIRISS (F115W, F200W direct images; GR150R/C WFSS)
  and o002 MIRI (F770W, F1000W, F1500W, F1800W); one `astroquery.mast` query, UTC in `summary.json`.
- Observed: Noirot+2023 NIRISS WFSS redshift catalogue (190 rows, sha256 pinned in the script; SOURCES.md).
  `Q_FLAG` 1 = multiple features, 2 = single feature (paper Table B1; Fig. 5 caption: "single spectral features
  aided by photometry").
- Record positions: `results/grism_4043/coverage.json` (9 candidate records plus `1816`, `384`).

## Result (derived)
- 5 of 11 records are inside the NIRISS footprints (`1032`, `1096`, `2559`, `2925`, `829`); 4 of them also inside
  MIRI (`2925` is not). `829` is 0.3″ inside the MIRI edge.
- **`f200w_384` and `f200w_1816` are outside every 2736 NIRISS and MIRI footprint.** Together with the 2″ MAST
  cone (only NIRCam 2736 / 6882 / 12774 and one 9276 NIRSpec image), no second 2022 instrument covers `384`;
  the "2736 systematic" explanation cannot be tested from 2022 data in MAST (question 2: no).
- Noirot+2023 entries within 0.5″ (ASSUMPTION: match radius): `1096` (ID 908, 0.025″, z_grism 1.3618 ± 0.0013,
  Q_FLAG 2) and `2925` (ID 1100, 0.292″, z_grism 1.9807 ± 0.0017, Q_FLAG 2). Both NED redshifts equal these
  values, so both are **single-feature NIRISS grism redshifts**, not multi-line or slit redshifts (question 3).
- `1032`, `2559` and `829` are inside the footprint but have no Noirot entry within 0.5″ (nearest 8.64″, 4.47″,
  1.81″). Being inside a footprint does not imply an extracted spectrum.

## Reading
- `2925` (verdict "known background galaxy"): the verdict stands; its redshift is a single-line grism value.
- `1096` (shear check `results/smacs_1096_shear/RESULT.md` used z 1.3618): the redshift is single-feature; the
  shear result depends on it being correct (not re-tested here).
- `384` stays "needs spectroscopy" (D-017 addendum 3); no archival route to a 2736-independent 2022 position.
- A NIRISS re-extraction of `1032`, `2559` or `829` is possible from the 2736 WFSS data but is not attempted
  (no open question on those records needs a redshift: `829` already has a NIRSpec z).

## Limitations
- Footprint membership uses the MAST polygons only (no detector gaps or trace geometry).
- Q_FLAG meaning is from the paper's table description; the single feature's identity per source is not checked.
