# W1 real-image metacal R, v2 — result (2026-10-18)

Plan (frozen first): `../PLAN_v2.md`. Script: `scripts/w1_r_realmetacal.py` (`validate2`, `real2`, `sheet --v2`);
tests: `tests/test_w1_r_realmetacal.py`. Numbers: `validate.json`, `r_real.json`. Contact sheet: `sheet.png`.
Data: `observed` MAST jw01181-o001_t001 F200W i2d, streamed from S3 (22 tiles read in 203 s; nothing stored).

## Gate (simulated)
V1 metacal 1.442 ± 0.226 vs direct 1.732 ± 0.064, tolerance 0.642: **pass**; V2 sign +0.036: **pass**. Identical to
v1 to the printed precision: in the isolated-galaxy scene the floor never binds (estimate 2.19 px) and no match is
shared. The gate stays weak (its tolerance is set by the metacal error).

## Real data (`derived`)
22 tiles pass the coverage rule; the frozen tile gate excludes 4: (0, 1024) 96 detections, raw PSF σ 0.29 px
(the v1 ringing tile); (1024, 9216) 32; (3072, 2048) 139; (3072, 1024) 679 (median 1,860).

| Variant | tiles | n | R | ±1σ tiles | ±1σ sources | median | R_11 / R_22 | frozen rule |
|---|---|---|---|---|---|---|---|---|
| **plan** (fixnoise, PSF_d 0.066″) | 18 | 804 | **1.056** | **0.156** | 0.162 | 0.916 | 0.773 / 1.338 | **conservative** |
| no fixnoise (reported) | 18 | 2,647 | 1.075 | 0.052 | 0.039 | 1.087 | 0.970 / 1.181 | (adoptable) |
| PSF_d × 1.15 (reported) | 18 | 799 | 0.875 | 0.162 | 0.148 | 0.924 | 0.456 / 1.294 | conservative |

- **Rule outcome (frozen): conservative.** R − 2σ = 0.744 < 0.80, so the borrowed R = 0.414 stays and the
  GOODS-N / CEERS / PRIMER limits are unchanged. v2 and v1 (0.953 ± 0.239) agree; the error fell by 1.5×.
- The no-fixnoise variant would pass the adoption bar but is not the plan's estimator (noise-shear bias not
  removed); it is reported only.
- By S/N (plan): 10–20 2.90 ± 1.01 (n 65), 20–50 0.81 ± 0.32, 50–200 1.01 ± 0.14, ≥ 200 0.81 ± 0.17.
- **R_11 < R_22** in every variant (plan 0.77 vs 1.34). The simulated gate shows the same sign (1.18 vs 1.71), so
  it is a property of the chain (hypothesis: the 90°-rotated fixnoise or the pixel-grid interpolation), not of
  the data; not investigated. The screen uses the mean of both, as here.

## Contact sheet (looked at)
`sheet.png`: three excluded tiles and one kept. All three excluded tiles show a radial starburst of ringing after
deconvolution (lower right of (0, 1024); right edge of (1024, 9216); left edge of (3072, 2048)) and a suppressed
background in the metacal image, so the inflated clipped σ explains their low counts. (3072, 2048) has no masked
pixel, so the noise fill of change 2 does not remove the ringing: its source is a sharp feature in the data
(not identified). The kept tile (1024, 1024) shows no ringing. (3072, 1024) was not looked at.

## Reading
With a mask, a PSF σ floor, unique matching and 18 tiles, the real-image metacal response of the pipeline-like
shapes is 1.06 ± 0.16, consistent with R_sim ≥ 0.94 (D-101) and v1. The bar for adopting it (R − 2σ ≥ 0.8) is not
reached; the PSF_d × 1.15 shift (−0.18) is as large as the statistical error, so the PSF model, not the tile
count, is now the limiting term.

## Limitations
One field and filter; analytic PSF_d; white-noise fixnoise; selection response not separated; the ringing source
is not removed, only gated; the tile gate thresholds are ASSUMPTIONs.
