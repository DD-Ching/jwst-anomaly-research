# W1 real-image metacal R, v3 (empirical PSF_d) — result (2026-10-18)

Plan (frozen first): `../PLAN_v3.md`. Script: `scripts/w1_r_realmetacal.py` (`epsf3`, `validate3`, `real3`);
tests: `tests/test_w1_r_realmetacal.py`. Numbers: `epsf.json` (gate E, star list, ePSF image), `posthoc_validate.json`,
`posthoc_r_real.json` (both **post hoc**). Sheet: `epsf_sheet.png` (looked at).
Data: `observed` MAST jw01181-o001_t001 F200W i2d (S3, streamed; nothing stored); Gaia DR3 `gaia_source` (ESA TAP).

## Gate E (frozen): **failed** → v3 stops; no R from this plan
- 6 Gaia DR3 sources fall inside the mosaic's pixel bounds (≥ 20 px from the edge). 4 fail the selection: G 15.97
  (RUWE 6.45), RUWE 1.97, RUWE 1.66, and a 2-parameter solution. 2 pass (G 18.53 and 20.53); both stamps pass every
  stamp test and the size test. 2 < 5 stars → gate E fails (`epsf.json` `gate_e_pass` false).
- The ePSF FWHM part would pass: half-maximum area 0.0782″ (in [0.060″, 0.090″]); GalSim FWHM 0.0764″ vs the
  analytic PSF_d ⊗ pixel 0.0740″ (reported).
- Sheet: both stars are clean and unsaturated and show the six-armed JWST diffraction spikes out to the stamp
  edge (16 px); the analytic PSF_d has none. The median ePSF keeps the spikes.
- **Code fix after the frozen run (deviation, outcome unchanged):** the flat-core test first took its 3 × 3 core
  at the brightest pixel, which on a saturated plateau is the plateau corner (a unit test caught it); it now uses
  the stamp centre, as `w12_niq_epsf.build_epsf` does. Re-run: identical star decisions and ePSF image.

## Post hoc diagnostic with the two-star ePSF (not the frozen rule; not adoptable)
Run although gate E failed (`--posthoc`), to size the PSF-model term.
- **Simulated gate (ePSF as scene PSF and PSF_d): V1 fails.** Metacal 0.969 ± 0.235 vs direct 1.634 ± 0.055;
  |Δ| 0.665 > tolerance 0.647. V2 sign passes. So the chain is not validated with this PSF either.
- **Real tiles** (22 tiles; gate decided on the ePSF plan variant: median 3,505.5 detections; excluded (0, 1024),
  (1024, 8192), (1024, 9216)):

| Variant | tiles | n | R | ±1σ tiles | R_11 / R_22 | rule function |
|---|---|---|---|---|---|---|
| ePSF, fixnoise | 19 | 1,032 | 0.707 | 0.169 | 0.466 / 0.948 | conservative |
| ePSF, no fixnoise | 19 | 4,194 | 0.906 | 0.043 | 0.784 / 1.029 | (adoptable; not the estimator) |
| analytic PSF_d (v2 estimator) | 19 | 850 | 1.006 | 0.134 | 0.736 / 1.276 | conservative |

- The analytic row does not reproduce v2 exactly (1.056 ± 0.156, 18 tiles, 804 sources) because the tile gate is
  set by the plan variant, which is now the ePSF (other detection counts; (1024, 8192) out, (3072, 1024) and
  (3072, 2048) in). The two agree within 1σ.
- In the same tiles the ePSF lowers R by 0.30 relative to the analytic PSF_d, the same sign as PSF_d × 1.15
  (−0.18 in v2). Hypothesis (not tested): the spike wings widen the deconvolution kernel, so the reconvolved
  galaxies are rounder.

## Reading
The frozen v3 gives no R: GOODS-N F200W o001 has too few usable Gaia stars (2) for the plan's ePSF. The post hoc
two-star diagnostic fails the simulated gate and moves R down to 0.71 ± 0.17. Every estimate (v1, v2, v3 post hoc)
stays above 0.414 − 2σ, so the borrowed **R = 0.414 stays** and the blank-field limits are unchanged. The
real-image R is model-dependent at the ±0.3 level, so this metacal route does not reach the 0.8 adoption bar.

## Limitations
One field and filter; two stars; one PSF for the whole mosaic; spikes truncated at 16 px; the cubic recentring
broadens the ePSF slightly; white-noise fixnoise; selection response not separated; metacal tile images of v3 not
looked at (v2's `sheet.png` covers the same tiles with the analytic PSF).
