# W1 blank-field responsivity R by metacalibration of real F200W images — plan

Frozen 2026-10-18 before any number. Script: `scripts/w1_r_realmetacal.py`. Trigger: D-101 handoff and "Revisit
if" (a real-image calibration); the forty-field W1 limit (0.638× Takahashi & Asada, #331) needs R ≳ 0.414.

## Question and hypothesis
- **H:** the measurement response R of the pipeline-like shapes (the `w1_r_sim.measure` emulation: kernel 2 px,
  3σ smoothed-noise threshold, 5 px, no deblending; PSF σ subtracted by `exotic_screens.intrinsic_ellipticity`)
  on real GOODS-N JADES 1181 o001 F200W i2d data is ≥ 0.414.
- **Ordinary reasons it could be lower than the simulation (R_sim ≥ 0.94, D-101):** blending and neighbours inside
  isophotes, real (clumpy) profiles, drizzle-correlated noise, the real PSF.

## Method (metacalibration, Huff & Mandelbaum 2017 / Sheldon & Huff 2017; GalSim implementation)
- Data (`observed`): 1024 × 1024 px SCI tiles of `jw01181-o001_t001_nircam_clear-f200w_i2d.fits` (MAST S3, streamed,
  not stored), on a fixed 1024-px grid; the first six tiles (row-major from the lower-left) with ≥ 98 % finite,
  WHT > 0 pixels. Sky per tile = median of sigma-clipped pixels, subtracted.
- Metacal per tile (`derived`): GalSim `InterpolatedImage` (x_interpolant Lanczos(15)) ⊗ Deconvolve(PSF_d), sheared
  by g = ±0.01 in each component (and 0), ⊗ PSF_r. PSF_d = `w1_r_sim.psf()` (Airy 2 µm ⊗ Gaussian, FWHM 0.071″;
  ASSUMPTION for the drizzled PSF). PSF_r = round Gaussian, FWHM 0.10″ (ASSUMPTION; broader than PSF_d ×(1+2g), so
  the measured R belongs to slightly broader-PSF images, which lowers R for small sources: conservative for H).
- Noise symmetrisation ("fixnoise", ASSUMPTION white noise): a white-noise image of the tile's clipped-MAD σ, rotated
  90°, put through the same operation, rotated back and added to each output.
- Measurement on each of the five images separately with `w1_r_sim.measure`'s settings (noise σ = that image's
  clipped MAD). Sources are matched to the unsheared image's by the sheared position (first order in g) within
  1.5 px. PSF σ subtracted = `exotic_screens.psf_sigma_px` on the unsheared catalogue (the screen's rule).
- Per source R_ii = (ε_i(+g) − ε_i(−g)) / 2g; R = mean of (R_11 + R_22)/2 over sources with S/N ≥ 10 (unsheared),
  resolved in all five images, not touching the tile edge (8 px) and matched in all four sheared images. Error:
  bootstrap over sources (1000). Median also reported.

## Validation gate (simulated, before any real-image number is read)
- **V1** a 1024² GalSim scene, Sérsic n = 1 galaxies on a 48-px grid with `w1_r_sim` sizes {2.5, 4, 6} px HLR and
  S/N {20, 50, 200}, intrinsic ε ~ σ 0.25, PSF_d, white noise. R_direct: the same galaxies drawn at g = ±0.01 each
  with the same noise (`w1_r_sim` style). R_mc: the full metacal chain above on the g = 0 scene.
  Pass: |R_mc − R_direct| ≤ 0.10 R_direct + 2σ (ASSUMPTION). R_direct uses PSF_r for its images, so both measure the
  response of the same PSF.
- **V2** sign: a scene with true g1 = +0.03 gives mean Re ε above the g = 0 scene.
- Fail → the tool is not trusted; no real-image number is used for a decision and the PR says so.

## Rule (frozen)
- R_real ≥ 0.414 − 2σ: the borrowed R is conservative or matched; GOODS-N / CEERS / PRIMER limits stand (no change).
- R_real < 0.414 − 2σ: the blank-field limits are optimistic; re-run them with R_real (a new PR).
- Reported, not used to tune: per-tile R, the gap to R_sim, the no-fixnoise variant, a PSF_d FWHM ×1.15 variant.

## Limitations stated in advance
One field and filter; PSF_d analytic, not the drizzled empirical PSF; white-noise fixnoise vs correlated noise;
selection response (detection in sheared images) not separated; R is for PSF_r images, not the catalogue's own.
