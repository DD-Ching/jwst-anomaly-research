# Vetting 1316 and 2184 (Abell S1063 `arcs` anti): result (2026-10-16)

Plan: `PLAN.md` (committed before the run of `scripts/arcs_vet_noise.py`). Labels: injected sources are
simulations; re-measured axes and rates `derived`; predicted axes `model_prediction`; HST shapes `observed`.

## Measurement-noise injections (frozen plan; `noise_injections.json`, `.ecsv`)
- Gate passed: the m = 24 copies come back aligned in ≥ 0.9667 of valid injections for all 24 rows.
- **1316: P_anti = 0 / 60; 2184: P_anti = 0 / 60** → both `not noise-explicable` under the frozen rule (< 0.05).
- Sample: 3 bgsub anti rows observed, 0.017 expected from measurement noise alone → `excess remains`.
- Limitation: the injected ellipticity is each row's catalogue ellipticity (0.58 and 0.55 for the targets); a
  Gaussian profile without PSF convolution; only the bgsub (moment) measurement is calibrated, not the pipeline
  segmentation axis that defined `anti`.

## Sensitivity (post hoc; deviation from the plan)
- Injected e = 0.2 (`noise_injections_e0.2.json`): the all-row gate **fails** (min 0.45; the failing rows are
  large galaxies, `semimajor_px` 14–24 px against the 3 px weight, `.ecsv`). For the two targets the bright copies
  are aligned in 60 / 60 each, and P_anti = 2 / 60 (1316) and 1 / 60 (2184), both below 0.05 (per-target gate,
  post hoc).
- Injected e = 0 (`noise_injections_e0.0.json`): no true axis, so the gate is meaningless (P_anti 26 / 60 and
  25 / 60, near a random axis's 1/3), so this test cannot exclude a near-round true source.
- Reading: if the true sources are as elongated as observed (or e ≥ 0.2), measurement noise at their flux does not
  turn an aligned axis anti. A near-round true source remains a possible explanation; the intrinsic-shape null
  (`../arcs.json`, P = 0.0021 for the sample) is the test for that.

## Independent HST axis (post hoc; `hst_shapes.json`)
- HFF-DeepSpace SExtractor shapes (VizieR `J/MNRAS/506/928`), separate data and noise: 1316 PA 49.4°, 65.9° from
  the prediction (`anti`); 2184 PA 65.3°, 86.1° from it (`anti`). Both are at S/N ≈ 1 in that catalogue (σ_mag 1.07
  and 1.26), so this agreement is weak. ASSUMPTION: the HFF mosaics are north-up.

## Looked at (F200W 6″ raw, 1″-filtered and 0.75″ zoom; figure not stored)
- 1316: two compact knots about 5 px (0.16″) apart along the measured axis (by eye): a blend of two faint sources
  or a two-clump galaxy is a leading ordinary explanation; a deblend has not been run.
- 2184: a faint irregular blob about 2″ from the giant tangential arc, whose direction matches the prediction.

## Status
Both stay **open, unvetted, not candidates**. Tested and passed: measurement noise (observed elongation), ICL,
model error (four models). F444W axis: untestable. Untested: a blend or clumps (1316 by eye), photo-z failure (an all-VizieR
query within 1.5″ found no spectroscopic redshift for either; the nearest is the giant-arc image at z 1.2282, 1.34″
from 2184; not stored), fragment of the giant arc (2184).

## Next
Deblend 1316 into its two knots and measure each knot's colour and axis; a spectroscopic redshift for 2184 (MUSE
cube of the field, if public through the ESO archive) or a photo-z from a second code.
