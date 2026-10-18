# W1 real-image metacal R, v2 (mask, PSF σ floor, unique matching, all tiles) — plan

Frozen 2026-10-18 13:15 UTC before any v2 number. Script: `scripts/w1_r_realmetacal.py` (`validate2`, `real2`).
v1 (`PLAN.md`, `RESULT.md`, #337): R = 0.953 ± 0.239 (323 sources, 6 tiles); handoff named the fixes below.
v1 numbers and code paths are unchanged.

## Question
Does a cleaner, larger measurement on the same mosaic (JADES 1181 o001 F200W i2d, MAST S3, streamed) give a
2σ lower bound on R high enough to raise the borrowed R = 0.414 of the GOODS-N / CEERS / PRIMER W1 limits?

## Changes from v1 (all fixed before running)
1. **Tiles:** every 1024² tile of the fixed grid with ≥ 98 % finite, WHT > 0 pixels (v1: the first six).
2. **Mask before deconvolution:** non-finite or WHT ≤ 0 pixels, dilated by 2 px, are filled with white noise of the
   tile's clipped σ (seeded) instead of zeros; sources whose unsheared centroid is within 16 px of a masked pixel
   are dropped (ASSUMPTION: 16 px ≈ 8 PSF_d FWHM).
3. **PSF σ floor:** the subtracted PSF σ is max(`es.psf_sigma_px` estimate, 0.97 px), 0.97 px = the 0.071″
   measured FWHM of PSF_d / 2.355 (ASSUMPTION). The raw estimate is stored.
4. **Unique matching:** a sheared-image source matched by more than one unsheared source is dropped for all of
   them.
5. **Tile quality gate (frozen):** a tile is excluded from R if its raw PSF σ estimate is below the floor or its
   unsheared detection count is below half the median over accepted tiles (the v1 ringing tile had 94 vs
   1,755–2,370 and σ 0.29 px). Excluded tiles are listed.
6. **Error:** primary 1σ = bootstrap over tiles (2000 draws, tiles weighted by their sources); the source
   bootstrap is reported.

Unchanged: PSF_d analytic (`w1_r_sim.psf`, parameter 0.066″), PSF_r 0.10″ Gaussian, g = ±0.01, white-noise
fixnoise, `w1_r_sim.measure` settings, S/N ≥ 10, 8-px edge, 1.5-px match radius. The empirical drizzled PSF is
**not** in this step (deferred; the PSF_d × 1.15 variant is the reported sensitivity).

## Validation gate (simulated, before the real run)
V1 and V2 of `PLAN.md` re-run with changes 3 and 4 (seed 7; same tolerance). Fail → no v2 real number is used.

## Rule (frozen), on the plan variant with the tile-bootstrap σ
- R − 2σ ≥ 0.80: **adoptable** — open a new frozen plan to re-run the blank-field limits with R = R − 2σ
  (ASSUMPTION: the lower bound, not the mean, is adopted).
- 0.414 − 2σ ≤ R and R − 2σ < 0.80: **conservative** — borrowed 0.414 stays; limits unchanged.
- R < 0.414 − 2σ: **optimistic** — re-run the limits lower (new PR).

Reported, not used for the rule: no-fixnoise variant, PSF_d × 1.15 variant, per-tile R, by-S/N R.

## Limitations stated in advance
One field and filter; analytic PSF_d; white-noise fixnoise on correlated noise; selection response not separated.
