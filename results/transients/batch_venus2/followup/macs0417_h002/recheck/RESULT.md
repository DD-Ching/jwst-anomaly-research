# MACS0417 h002 time-delay re-check: result

Owner instruction 2026-10-17 (after #308); plan `PLAN.md`. Model quantities are `model_prediction`; the realisation ranges
are sample ranges of the first 10 glafic realisations (map000–map009), not confidence intervals. Files:
`recheck_delays.json` (per model, interpolation order 1 and 3, old and new values), `summary.json`.

## Checks
- **Units and axes:** ∂ψ/∂θ from the glafic ψ map equals the deflection maps (median ratio 1.0000 ± 0.0002 on both pixel
  axes, best fit and map000; not stored beyond this text). ψ is in arcsec² on the same axes, normalised like α to
  D_LS/D_S = 1, so scaling both by f = 0.4299 is consistent.
- **z_l = 0.443:** the arXiv:1810.13439 abstract (Mahler et al., RELICS) states "MACS J0417.5-1154 at z=0.443".
- **Reviewer's claim (verified, not assumed):** #308 used φ = ½|θ − β(θ)|² − ψ with each pixel's own β. At an exact
  image this equals φ(θ; β0) = ½|θ − β0|² − ψ, but at a grid point with source-plane residual δβ it is off to first
  order (≈ α·δβ). The re-check fixes β0 = β(h002) and root-finds β(θ) = β0 from each grid candidate.
- **Self-reference:** the delay of h002 relative to itself is exactly 0 in all 11 models (the #308 grid point nearest
  h002 showed −14 to −65 d, an artifact of the per-pixel β).
- **Convergence:** converged images have source-plane residuals ≤ 3e-14″; interpolation order 1 vs 3 differs by
  ≤ 3.6 d. One 7.3–7.6″ grid candidate per model does **not** converge (residual 0.06–0.39″): a local minimum of
  |β − β0|, not an image. It includes the best fit's 7.34″ entry (826 d) that set the #308 upper bound.

## Old (#308) vs new delays relative to h002 (days)
| Image | #308 best fit | #308 realisations | New best fit | New realisations (sample range) |
|---|---|---|---|---|
| ≈ 7.6″ (one image, not two) | 826.4 (spurious), 751.1 | 649.1–817.0 (mixed with spurious) | 748.0 | 708.1–777.7 |
| ≈ 9.3″ | 1,025.5 | 1,067.9–1,183.2 | 1,021.5 | 1,056.3–1,136.3 |
| ≈ 25.4″ | −4,814.8 | −5,258.6 to −4,504.5 | −4,824.2 | −5,318.7 to −4,515.9 |

## What changes
- **7.6″ image:** in VENUS (2025-10-06) it shows image 1.1's state of **2023-08-20 to 2023-10-28** (#308: about
  2023-07 to 2023-12). No observation covers image 1.1 then, so the VENUS non-detection stays **uninformative**; the
  conclusion "not testable with current epochs" is unchanged.
- **Next useful epoch of the 7.6″ image:** **2027-09-14 to 2027-11-22** (#308: about 2027-07 to 2028-01).
- **Position:** the converged glafic 7.6″ image lies 0.82–1.10″ from catalogued 1c.2 (#308's 0.41″ was the spurious
  candidate); lenstool positions were not re-checked here. Position uncertainty is of order 1″.
- At the converged best-fit position (0.047″ from the #308 `g2_7.6` aperture) the 2022 → 2025 photometry shows no brightening
  (F150W 94.5 → 93.4 nJy; `../counter_image_lightcurves.ecsv`), consistent with the delays.
- The real change at h002 itself is not affected.
