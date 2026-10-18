# W1 blank-field responsivity R from image simulations — result (2026-10-18)

Plan (frozen first): `PLAN.md`. Script: `scripts/w1_r_sim.py`; tests: `tests/test_w1_r_sim.py`. Numbers:
`r_sim.json` (the plan's settings) and `r_sim_rawthr.json` (post hoc threshold variant). All R values here are
`simulated`; the GOODS-N row distribution used as weights is `derived` from the MAST pipeline catalogue
(jw01181-o001_t001 F200W) and the DJA v7.3 photo-z. GalSim 2.8.5, photutils 3.1.1 (not stored in the JSON).

## Checks
- **C1** (noiseless, no PSF, no kernel, n = 0.5, half-light radius 8 px): R = 0.998 (`check_c1` in both JSON
  files; `tests/test_w1_r_sim.py` asserts |R − 1| < 0.03).
- **C2** a noiseless round n = 1 galaxy through PSF and kernel, with the stamp PSF σ subtracted: |ε| < 0.02 (unit
  test; the value is not stored).
- **C3** sign: g1 > 0 raises Re ε, g2 > 0 raises Im ε (unit test); every cell R is positive (`r_sim.json`).
- The galsim/photutils tests are skipped in CI (its install has neither); they passed locally (4 passed).

## Deviations from the plan (labelled)
- The PSF model measures FWHM 0.071″ (`psf_fwhm_arcsec_measured`), not the plan's 0.066″.
- Post hoc, added after the first run: (a) the real field's PSF σ (1.327 px, what the screen subtracts) beside
  the plan's noiseless-stamp σ (1.767 px; that stamp is thresholded at 1e-3 of its peak, a deeper isophote than the
  galaxies', so it over-subtracts); (b) `_resolved90` weights, which keep only rows whose nearest cell has ≥ 90 %
  of its pairs used; (c) the jwst-style threshold run. After code review a pair counts only if both 90°-rotated
  copies are detected and resolved in all sheared versions (the first run kept them separately, which breaks the
  shape-noise cancellation); all numbers below are from the re-run.
- The weight rows pass the screen's S/N ≥ 10, resolved and spike-veto cuts; the GOODS-N streak mask and PSF χ are
  not applied, so they approximate the screen's used rows.

## Result (rule: R_sim ≥ 0.414 − 2σ → "conservative", `w1_r_sim.rule_outcome`)
Weighted by the 2,259 GOODS-N o001 rows (`r_sim.json`, plan threshold):

| PSF σ subtracted | n = 1 | n = 4 | n = 1, resolved90 | n = 4, resolved90 | rule |
|---|---|---|---|---|---|
| stamp (plan, 1.767 px) | 1.369 ± 0.028 | 1.343 ± 0.018 | 1.182 ± 0.007 | 1.137 ± 0.009 | conservative (all) |
| field (post hoc, 1.327 px) | 1.017 ± 0.006 | 0.965 ± 0.007 | 1.006 ± 0.005 | 0.940 ± 0.007 | conservative (all) |

- Rows without a usable cell (`n_rows_unrepresented`): 0–238 of 2,259 in the plain columns, 43–523 under
  resolved90.
- Per cell (field PSF σ): 36 of 40 cells have ≥ 2 used pairs; R 0.836–2.390. The 27 cells with ≥ 90 % of pairs
  used span 0.836–1.404. Cells where few pairs stay resolved exceed 1 (a hypothesis, not tested: the pairs that
  stay resolved after the PSF subtraction are a selected subset, and the subtraction amplifies their ε).
- **Post hoc threshold variant** (`r_sim_rawthr.json`; 3σ of the unsmoothed noise on the smoothed image, the jwst
  `source_catalog` convention as we understand it, ASSUMPTION, not re-read this run; a higher isophote): weighted R
  1.034–2.057 in all eight columns, every one "conservative"; 18 cells usable with the field PSF σ (0.926–4.13).
- **Reading:** in every variant the weighted simulated response of the pipeline-like shapes is ≥ 0.940, more than
  twice the borrowed R = 0.414. In this simulation, injections that paint W1 shear with R = 0.414 understate the
  shape change a real W1 lens would cause, so the GOODS-N and CEERS limits (D-099) are conservative with respect
  to R. No limit is changed in this PR (the plan's rule).
- **The gap to the cluster fits** (R 0.41–0.60, D-053/D-098) is about a factor two. The fitted slope also absorbs
  lens-model, photo-z (z = 2 default without one) and member/foreground dilution; which of these makes the gap is
  a hypothesis, not tested here.

## Limitations
Isolated stamps (no blending or neighbours inside the isophote), Airy ⊗ Gaussian PSF (FWHM 0.071″, no spikes, not
STPSF), white noise (drizzled noise is correlated), Sérsic profiles only, a fixed σ_e = 0.25. The weighted R depends
strongly on the PSF σ subtracted (0.94–1.02 vs 1.14–1.37) and on the threshold, but none of the variants comes near 0.414.

## Next
A real-image check: shear GOODS-N / CEERS F200W i2d cutouts (S3 byte ranges) with metacalibration (deconvolve an
STPSF model, shear, reconvolve) and re-measure them like the pipeline; and run the simulation on a cluster field's
rows to see whether the fitted R gap comes from the measurement or from the lens model.
