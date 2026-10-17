# W1: detection-kernel, jitter and drizzle rounding vs STPSF's over-predicted star χ (2026-10-17)

Plan and rule S3: `PLAN.md` (frozen before any number). Script: `scripts/w1_psf_smooth.py`; numbers:
`smooth_chi.json`. Labels: model PSFs are `model_prediction` (STPSF 2.2.0, data 2.2.0, F200W, `DET_DIST`, 10″, the
same 24 grid positions as `results/w1_stpsf`); star χ and areas are `derived` (`results/w1_stpsf/stars_check.json`).
Smoothing widths and the S3 tolerance are ASSUMPTIONs.

## Checks
- **Reproduction:** the unsmoothed variant (`none`) gives model |χ| 0.106 / 0.229 / 0.251 at the three stars'
  areas, as `results/w1_stpsf/stars_check.json`.
- **Smoothing code:** Gaussian σ add in quadrature on an elongated test Gaussian (`tests/test_w1_psf_smooth.py`).
- **Kernel:** current jwst `tweakreg_catalog` measures morphology on data convolved with a Gaussian of FWHM
  `kernel_fwhm` (default 2.0 px; source code read 2026-10-17, SOURCES.md). That the Abell 2744 catalogue used the
  default is an ASSUMPTION (the catalogue's metadata were not checked). STPSF already includes 0.007″ jitter; the
  jitter variants add to it. The drizzle box is approximated by its second-moment Gaussian (σ = 1/√12 px).

## Result (S3: each star within 0.05 per component; `smooth_chi.json`)
| variant | total σ (px) | star 623 Δ (re, im) | star 7644 Δ | star 8291 Δ | S3 |
|---|---|---|---|---|---|
| none | 0 | −0.053, +0.004 | −0.071, −0.006 | −0.085, +0.109 | fail |
| kernel | 0.85 | −0.037, +0.009 | −0.051, +0.003 | −0.021, +0.103 | fail |
| kernel+box | 0.90 | −0.033, +0.006 | −0.051, −0.000 | −0.012, +0.105 | fail |
| kernel+box+jit5 | 0.91 | −0.032, +0.008 | −0.051, −0.000 | −0.012, +0.105 | fail |
| kernel+box+jit10 | 0.95 | −0.034, +0.008 | −0.048, −0.003 | −0.012, +0.105 | fail |

(Δ = star − model, pixel frame.) S3 fails for every variant. Per star (`s3_star`): 623 passes in every smoothed
variant; 7644 passes only in kernel+box+jit10 (−0.048); 8291 fails in all, so it alone decides the jit10 outcome. The smoothing lowers the over-prediction of the real
component (the axis both share) by 0.016–0.023 for stars 623 and 7644 and by 0.06–0.07 for star 8291. Star 8291's
imaginary component stays 0.10 off in every variant: its orientation differs from the model's (15.5° vs ≈ 179°,
`results/w1_stpsf`), and an isotropic smoothing cannot rotate χ.

## Conclusion (narrow)
The detection kernel (if the default) accounts for part of STPSF's |χ| excess at the stars' areas, not all of it;
with the third star's orientation unexplained, the smoothed model is not a validated PSF for a per-source χ_PSF
forward model. Per the plan, the PSF-model route to Abell 2744's B-mode floor stops here. W1 limits unchanged.

## Limitations
- Three stars; the area-matched level comes from a 36-level grid, so neighbouring variants can pick the same level.
- The box and jitter terms are Gaussian approximations; the PSF is not drizzled from dithered exposures.
- Whether star 8291 is a clean, unsaturated star was not re-checked here.

## Next
PSF-weighted moments on pixels (e.g. GalSim HSM with an empirical PSF; `/reuse-check` first) instead of catalogue
isophotal moments, or accept the 4.21 floor and spend W1 time on new fields.

## Regenerate
```
curl -sSL <stpsf-data-LATEST.tar.gz URL, SOURCES.md> | tar xz -C <dir>   # 88 MB, streamed, deleted after
STPSF_PATH=<dir>/stpsf-data python scripts/w1_psf_smooth.py --out results/w1_psf_smooth/smooth_chi.json
```
Time: 1.7 min on 4 processes.
