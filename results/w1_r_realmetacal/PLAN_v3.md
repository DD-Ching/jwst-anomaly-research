# W1 real-image metacal R, v3 (empirical drizzled PSF_d) — plan

Frozen 2026-10-18 14:12 UTC (commit time; first draft said 14:25) before any v3 number. Script: `scripts/w1_r_realmetacal.py` (`epsf3`, `validate3`,
`real3`). v2 (`PLAN_v2.md`, `v2/RESULT.md`, #340): R = 1.056 ± 0.156 (18 tiles), R − 2σ 0.744 < 0.8; PSF_d × 1.15
moved R by −0.18, so the analytic PSF_d is the limiting term. v1/v2 numbers and code paths are unchanged.

## Question
With PSF_d taken from stars of the same mosaic (JADES 1181 o001 F200W i2d) instead of the analytic Airy ⊗
Gaussian, is R − 2σ ≥ 0.80, so that the borrowed R = 0.414 of the blank-field W1 limits can be raised?

## Empirical PSF (frozen; ASSUMPTIONs)
1. Stars: Gaia DR3 `gaia_source` within the mosaic's footprint, 17 ≤ G ≤ 21, `ruwe` < 1.4 and a 5- or 6-parameter
   solution (`astrometric_params_solved` > 3; rejects most compact galaxies), at least 20 px from the mosaic edge.
2. Stamp: 33 × 33 px (half 16) cut around the Gaia position mapped by the i2d WCS; dropped if any pixel is
   non-finite or WHT ≤ 0, if the brightest pixel is > 3 px from the stamp centre, or the core is flat-topped
   (3 × 3 core min > 0.85 × max: saturation). Background = median of the stamp border; recentred by the
   thresholded flux-weighted centroid (cubic spline shift), normalised to unit sum.
3. Size outliers: stamps whose second-moment σ differs from the median of the kept stamps by > 20 % are dropped
   (neighbours, galaxies).
4. ePSF = pixel-wise median of the kept stamps, clipped at 0, unit sum. It is a drizzled effective PSF (pixel
   included), so it is deconvolved directly (GalSim `InterpolatedImage`, no extra `Pixel`).
5. **Gate E (computable):** ≥ 5 kept stars **and** the ePSF FWHM (half-maximum area, √(4A/π) × 0.031″) in
   [0.060″, 0.090″]. Fail → v3 stops; no R is read.

## Simulated gate (before the real run)
V1 and V2 of `PLAN.md` with the v2 changes, where the scene PSF and PSF_d are both the ePSF image (scene drawn
without an extra pixel, since the ePSF already includes it); same tolerance. Fail → no v3 real number is used.

## Real run
Exactly `real2` (tiles, mask, floor, unique matching, tile gate decided once on the plan variant, tile bootstrap),
with variants: **plan = ePSF, fixnoise**; reported: ePSF without fixnoise; analytic PSF_d (the v2 plan, as a
same-run check that v2 reproduces).

## Rule (frozen, `rule_outcome_v2` on the plan variant with the tile-bootstrap σ)
- R − 2σ ≥ 0.80: **adoptable** — a new frozen plan re-runs the blank-field limits with R = R − 2σ.
- 0.414 − 2σ ≤ R and R − 2σ < 0.80: **conservative** — borrowed 0.414 stays.
- R < 0.414 − 2σ: **optimistic** — re-run the limits lower (new PR).

## Limitations stated in advance
One field and filter; one PSF for the whole mosaic (no spatial variation, mixed detector/dither PSFs); a
median stack of few stars carries noise into the deconvolution; the cubic recentring broadens the ePSF slightly;
white-noise fixnoise on correlated noise; selection response not separated.
