# W1 blank-field responsivity R from image simulations — result (2026-10-18)

Plan (frozen first): `PLAN.md`. Script: `scripts/w1_r_sim.py`; tests: `tests/test_w1_r_sim.py`. Numbers:
`r_sim.json` (the plan's settings) and `r_sim_rawthr.json` (post hoc threshold variant). All R values here are
`simulated`; the GOODS-N row distribution used as weights is `derived` from the MAST pipeline catalogue
(jw01181-o001_t001 F200W) and the DJA v7.3 photo-z. GalSim 2.8.5, photutils 3.1.1 (not stored in the JSON).

## Checks
- **C1** (noiseless, no PSF, no kernel, n = 0.5, half-light radius 8 px): R = 0.998 (`tests/test_w1_r_sim.py`
  asserts |R − 1| < 0.03).
- **C2** a noiseless round n = 1 galaxy through PSF and kernel, with the stamp PSF σ subtracted: |ε| < 0.02 (test).
- **C3** sign: g1 > 0 raises Re ε, g2 > 0 raises Im ε (test); every cell R is positive (`r_sim.json`).

## Result (rule: R_sim ≥ 0.414 − 2σ → "conservative", `w1_r_sim.rule_outcome`)
Weighted by the 2,259 used rows of GOODS-N o001 (S/N ≥ 10, resolved, spike veto; real-field PSF σ 1.327 px):

| PSF σ subtracted | n = 1 | n = 4 | rule |
|---|---|---|---|
| stamp (plan: noiseless PSF-only stamp, 1.767 px) | 1.621 ± 0.014 | 1.522 ± 0.023 | conservative |
| field (deviation, post hoc: the screen's own 1.327 px) | 1.024 ± 0.006 | 0.963 ± 0.007 | conservative |

- Per cell (field PSF σ): 38 of 40 cells have resolved pairs; R ranges 0.684–2.781. The S/N 200 cells and the r_h ≥ 4 px, S/N 50 cells sit
  at 0.84–0.97; small, faint cells exceed 1 (a hypothesis, not tested: only the pairs still resolved after the PSF
  subtraction are used, and the subtraction amplifies their ε). The lowest cell is n = 1, r_h = 10 px, S/N 10 (0.684 ± 0.162, 42 pairs).
- **Post hoc threshold variant** (`--raw-threshold`: 3σ of the unsmoothed noise on the smoothed image, the jwst
  `source_catalog` convention as we understand it, ASSUMPTION, not re-read this run; a higher isophote): weighted R 1.54–2.37 in all four columns, rule
  "conservative"; only 19 cells keep resolved pairs (field PSF σ), cell R 0.926–11.8.
- **Reading:** under every variant the simulated response of the pipeline-like shapes is ≥ 0.96, more than twice
  the borrowed R = 0.414. The injections that paint W1 shear with R = 0.414 therefore understate the shape change
  a real W1 lens would cause in this simulation, so the GOODS-N and CEERS limits (D-099) are conservative with
  respect to R. No limit is changed in this PR (the plan's rule).
- **The gap to the cluster fits** (R 0.41–0.60, D-053/D-098) is about a factor two. The fitted slope also absorbs
  lens-model, photo-z (z = 2 default without one) and member/foreground dilution; which of these makes the gap is
  a hypothesis, not tested here.

## Limitations
Isolated stamps (no blending or neighbours inside the isophote), Airy ⊗ Gaussian PSF (FWHM 0.071″, no spikes, not
STPSF), white noise (drizzled noise is correlated), Sérsic profiles only, a fixed σ_e = 0.25. The weighted R depends
strongly on the PSF σ subtracted (1.0 vs 1.6) and on the threshold, but none of the variants comes near 0.414.

## Next
A real-image check: shear GOODS-N / CEERS F200W i2d cutouts (S3 byte ranges) with metacalibration (deconvolve an
STPSF model, shear, reconvolve) and re-measure them like the pipeline; and run the simulation on a cluster field's
rows to see whether the fitted R gap comes from the measurement or from the lens model.
