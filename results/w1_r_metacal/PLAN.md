# W1 blank-field responsivity R from image simulations of the catalogue shapes — plan

Frozen 2026-10-18 before any simulated number. Script: `scripts/w1_r_sim.py`. Revisit trigger: D-099 "Revisit
if" (a blank-field R calibration); handoff of #328 (the borrowed R dominates the W1 systematics).

## Question and hypothesis
The blank-field W1 screens (GOODS-N, CEERS) use R = 0.414, borrowed from the Abell 370 F200W cluster fit (D-099);
0.231 (MACS0647) is the sensitivity case. A cluster fit measures the moment response *and* lens-model, photo-z and
member dilution, so it need not equal the pure measurement response that the injections should use.
- **Hypothesis H:** the measurement response of the pipeline-like shapes (jwst `source_catalog`-style isophotal
  moments on kernel-smoothed data, Gaussian-PSF subtraction of `exotic_screens.intrinsic_ellipticity`) at F200W
  sampling, over the size / S/N range of the used rows, is ≥ 0.414 (then the borrowed R is conservative).
- **Ordinary reasons it could be lower:** isophotal truncation at faint S/N, the detection kernel rounding small
  sources, and the PSF-σ estimate (1st percentile of minor-axis σ) mis-subtracting.

## Simulation (simulated; every setting an ASSUMPTION)
- GalSim exponential (n = 1) and de Vaucouleurs-like (n = 4) Sérsic galaxies, isolated 64 × 64 px stamps,
  0.031″ px, F200W-like Gaussian⊗Airy PSF (FWHM 0.066″), white noise.
- Grid: half-light radius {1.5, 2.5, 4, 6, 10} px × isophotal S/N {10, 20, 50, 200}.
- Intrinsic ellipticity |e| ~ Rayleigh(σ = 0.25 per component), truncated < 0.8; each galaxy drawn at 90°-rotated
  pairs and at true shear g = ±0.02 (both components), with the same noise for ±g (noise cancels in the
  difference).
- Measurement mimics the pipeline: Gaussian kernel FWHM 2 px (`w1_psf_smooth.KERNEL_SIGMA_PX`), detection at
  3σ of the smoothed noise, npixels 5, no deblending; photutils `SourceCatalog(data, segm, convolved_data=…)`
  `semimajor_sigma`, `ellipticity`, `orientation`; the source nearest the stamp centre. PSF σ for the subtraction
  as the pipeline does it: minor-axis σ of a noiseless, PSF-only stamp measured the same way (the 1st-percentile
  star rule in the limit). Shapes → `exotic_screens.intrinsic_ellipticity`.
- R per cell = (⟨ε(+g)⟩ − ⟨ε(−g)⟩) / 2g, both components averaged; error from the per-pair scatter.

## Checks (before any conclusion)
- **C1** noiseless, PSF-free, kernel-free, large round Gaussian: R within 0.03 of 1.
- **C2** the PSF-σ subtraction: a noiseless round galaxy convolved with a Gaussian PSF of known σ gives |ε| < 0.02.
- **C3** shear sign: ε(+g) − ε(−g) has the sign of g in both components.

## Rule (frozen)
Weight the cell R by the joint (semimajor_px, S/N) distribution of the GOODS-N o001 used rows (S/N ≥ 10, resolved;
`inject_shear` shapes) mapped to the grid by nearest cell. Call R_sim that weighted mean.
- R_sim ≥ 0.414 − 2σ: the borrowed R is conservative or matched; GOODS-N/CEERS limits stand (no change made here).
- R_sim < 0.414 − 2σ: the blank-field limits are optimistic; re-run them with R_sim (a new PR) and record it.
- Either way the cluster-fit vs simulation gap is reported, not used to tune anything in this PR.

## Limitations stated in advance
Isolated stamps (no blending, no neighbours in the isophote), Gaussian-like PSF (no spikes, no STPSF), white noise
(drizzle correlates noise), Sérsic profiles (no clumps). A real-image metacalibration on the mosaics is the next
stage if this one is inconclusive.
