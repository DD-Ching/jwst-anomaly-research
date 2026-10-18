# W1 blank-field responsivity R by metacalibration of real F200W images — result (2026-10-18)

Plan (frozen first): `PLAN.md`. Script: `scripts/w1_r_realmetacal.py`; tests: `tests/test_w1_r_realmetacal.py`.
Numbers: `validate.json` (frozen gate), `validate_pooled_posthoc.json` (post hoc), `r_real.json` (the plan's run
and its two reported variants), `r_real_posthoc_excl_tile.json` (post hoc). Contact sheet: `sheet.png` (`sheet` subcommand).
GalSim 2.8.5, photutils 3.1.1 (not stored in the JSON).

## Validation gate (simulated)
- **V1** (one 1024² scene, seed 7): metacal R = 1.442 ± 0.226 vs direct shearing 1.732 ± 0.064; tolerance 0.642;
  **pass**. **V2** sign: a true g1 = +0.03 raises mean Re ε by 0.036; **pass**. Gate passed.
- The gate is weak: its tolerance (0.64) is set mostly by the metacal error. Post hoc, four scenes pooled
  (`validate_pooled_posthoc.json`): metacal 1.519 ± 0.118 vs direct 1.844 ± 0.039 (medians 1.292 vs 1.402).
  In this simulation metacal reads lower than direct shearing (by 0.33 in the mean, 0.11 in the median), so a
  metacal R on real data is more likely to understate the response than to overstate it (simulated; one profile
  type, isolated galaxies).

## Real data (`derived` from observed MAST jw01181-o001_t001 F200W i2d)
Six 1024² tiles passed the ≥ 98 % coverage rule (`r_real.json` `tiles`); 323 sources pass every plan cut.

| Variant | n | R | ±1σ (bootstrap) | median | R_11 / R_22 | frozen rule |
|---|---|---|---|---|---|---|
| plan (fixnoise, PSF_d 0.071″) | 323 | 0.953 | 0.239 | 0.835 | 0.951 / 0.954 | conservative |
| no fixnoise (reported) | 881 | 0.979 | 0.067 | 1.068 | 0.866 / 1.093 | conservative |
| PSF_d FWHM × 1.15 (reported) | 323 | 0.886 | 0.237 | 0.852 | 0.736 / 1.036 | conservative |

- **Rule outcome (frozen, `rule_outcome`): conservative.** The plan's R − 2σ is 0.475, above 0.414 but not by
  much; the measurement does not exclude R ≈ 0.5.
- The no-fixnoise variant is more precise (fixnoise doubles the noise, so fewer sources reach S/N 10) but its
  R_11 ≠ R_22 (0.866 vs 1.093) shows the expected noise-shear anisotropy; it is not the plan's estimator.
- Per tile (plan): R 0.604–1.851 with errors 0.133–0.811 (`per_tile`); by S/N bin 0.783–1.311 (`by_snr`).

## Contact sheet and a post hoc exclusion (deviation)
`sheet.png` (two tiles: original, metacal unsheared, 1p − 1m). Tile y0 = 0, x0 = 1024 shows a radial ringing
pattern after deconvolution near its lower-right corner (a feature sharper than PSF_d, by eye; not identified), and
its clipped noise is inflated: 94 detections vs 1,755–2,370 in the other tiles, PSF σ estimate 0.29 px vs
1.30–2.27 px. Excluding it (post hoc, `r_real_posthoc_excl_tile.json`; error approximate): plan R = 1.004 ± 0.270
(n 276), no fixnoise 1.002 ± 0.072, PSF_d × 1.15 0.940 ± 0.255; all "conservative". The other tile shown has no
such pattern. The tool has no bad-pixel or sharp-feature mask before deconvolution (a gap).

## Reading
On real GOODS-N F200W images, with the pipeline-like shapes and an analytic PSF, the metacalibration response is
0.95 ± 0.24 (plan), consistent with the simulated R_sim ≥ 0.94 (D-101) and above the borrowed 0.414 by the frozen
rule. The GOODS-N, CEERS and PRIMER limits stand as computed with R = 0.414 (no limit changed). The error is too
large to adopt R ≈ 1 for the limits; the cluster-fit gap (R 0.41–0.60, D-053 / D-098) remains unexplained.

## Limitations
One field, one filter, 323 sources; analytic PSF_d (Airy ⊗ Gaussian), not an empirical drizzled PSF; white-noise
fixnoise on drizzle-correlated noise; no bad-pixel mask before deconvolution; selection response not separated; R
belongs to images reconvolved with a 0.10″ Gaussian, broader than the catalogue's PSF; the gate's simulation is
isolated n = 1 galaxies only.
