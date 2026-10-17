# W1: does detection smoothing, jitter and drizzle rounding explain STPSF's over-predicted star χ? (plan, 2026-10-17)

Claimed by a cloud routine session. TASKS "JWST focus" 1 (W1); follows `results/w1_stpsf/RESULT.md` "Limitations"
and "Next". Written before any smoothed model PSF was computed (frozen at this commit; deviations get a labelled
addendum).

## Hypothesis and ordinary alternatives
- **H (processing):** the catalogue's morphology is measured on the detection image smoothed by the source-catalog
  kernel, from drizzled data with telescope jitter. Applying these roundings to the STPSF F200W model brings the
  area-matched model χ of the three Abell 2744 PSF stars (now 0.05–0.07 too large in |χ|) into agreement.
- **Alternatives:** the STPSF wavefront/spider model itself differs from the 2022 PSF; the three objects are not
  clean stars; isophote-level ↔ area matching is too coarse.
- **What changes the plan:** if S3 passes, the smoothed model is a validated PSF for a per-source forward model of
  χ_PSF (next cycle). If it fails for every tried smoothing, the PSF-model route to the floor stops and the next
  floor step is PSF-weighted moments on pixels.

## Method
Re-use `scripts/w1_stpsf_floor.py` PSFs (same 24 positions, 10″, `DET_DIST`). Convolve each with (a) a Gaussian
jitter of σ_j, (b) a box of one mosaic pixel (drizzle footprint, approximate), and (c) the catalogue's detection
kernel (Gaussian, FWHM read from the catalogue metadata or the jwst default 2.0 px if absent; ASSUMPTION). Then
re-run the area-matched star check.

## Rule (ASSUMPTION, computable)
- **S3:** a smoothing variant is consistent if, for each of the 3 stars, |χ_model − χ★| ≤ 0.05 per component at
  the star's own `isophotal_area`. Variants tried (all reported): kernel only; kernel + box; kernel + box + σ_j ∈
  {0.005″, 0.010″}. No variant is tuned beyond this list.
