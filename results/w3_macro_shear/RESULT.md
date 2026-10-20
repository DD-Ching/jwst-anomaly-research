# W3 negative-mass microlens inside a macro-lensed arc: result (2026-10-20)

All numbers are **simulated / model_prediction** (inverse ray shooting); no data were used. Plan frozen first:
`PLAN.md` (gates G1–G3; addendum 1 replaced the mis-specified G3 before G3′ was computed). Code:
`exotic_sim.macro_magnification_map`, `exotic_sim.macro_hole_halfwidths`, `scripts/w3_macro_shear.py`;
numbers in `result.json`, maps in `maps.png` (looked at).

## Gates
- **G1** (isolated lens vs `total_magnification`, radial profile, μ < 20): max deviation 0.08 % (+) / 0.05 % (−).
  Pass. Inside r < 1.95 the negative-mass map is exactly 0.
- **G2** (positive mass, κ 0, γ 0.3; analytic det J = 0 caustic vs the ray-shooting ridge): 4,000 / 4,000 caustic
  points have a local μ > 5 μ_macro; 100 % of pixels with μ > 20 μ_macro lie within 2 pixels of the caustic. Pass.
- **G3** as written (window mean = μ_macro within 2 %): **failed** (+14.3 % / −17.3 %). The gate was wrong, not
  the tool: a point lens adds a fixed source-plane excess π(R√(R² + 4s) − R²) (addendum 1).
- **G3′** (that excess on the isolated maps, R = 2.5 and 3): agrees to < 0.01 %, both signs. Pass.

## Measurements (negative mass; λ1 = 1 − κ − γ along the macro stretch, λ2 = 1 − κ + γ = 1)
| Case | μ_macro | umbra (μ / μ_macro < 0.1) area, θ_E² | half-width along stretch / across | predicted 2√λ (on axis) | min μ / μ_macro, ρ = 0.1, u₀ = 0 track |
|---|---|---|---|---|---|
| isolated | 1 | 12.49 | 1.99 / 1.99 | 2 / 2 | 0 |
| minimum | 10 | 4.47 | 0.625 / 2.00 | 0.632 / 2 | 0 |
| minimum | 100 | 1.55 (2.23 incl. side regions) | 0.20 / 2.00 | 0.20 / 2 | 0 |
| saddle | −10 | 0 | — | 0 (no on-axis umbra) | 0.91 (point source 0.88 / 0.69 off-axis) |
| saddle | −100 | 0 | — | 0 | 0.98 (point source ≥ 0.86) |

## What this shows
- **Minimum-parity macro images:** the W3 "vanishing" survives macro magnification. The umbra keeps zero flux,
  but along the macro stretch it narrows to the predicted half-width 2√λ1 (≈ 2 / √μ_macro when λ2 = 1);
  across it stays 2 θ_E. A source crossing along the stretch axis vanishes for a time shorter by √λ1 than the
  isolated formula 2 t_E √(4 − u₀²) gives.
- **Saddle-parity macro images:** there is **no umbra**. A negative-mass microlens dims a saddle image by at
  most ~30 % (point source) or ~10 % (ρ = 0.1) in these cases, with caustic spikes: it looks like ordinary
  positive-mass microlensing near a saddle (the known mimic in docs/exotic_lensing.md), not a vanishing.
- Hence the isolated-lens W3 injections (`inject_light_curve`) are a fair approximation only for minimum
  images and only for motion across the stretch. The plan's decision rule: hole area at μ_macro = 100 is 12 %
  of the isolated value in the minimum (between the 10 % and 50 % thresholds) and 0 % in the saddle, so the
  doc's W3 amplitudes for arcs are revised (docs/exotic_lensing.md notes).

## Not shown / limits
- One microlens, no microlens population (κ_* = 0), λ2 = 1 only, two |μ_macro| values; point and ρ = 0.1 disk
  sources; geometric optics. Light curves were not timed against any observed caustic transient.
- The 2√λ1 on-axis rule is a short derivation (minimum of λ t + 1/t), confirmed here numerically; the literature
  was not searched for a published version in this run, so no novelty is claimed.
- No search changed: existing W3 limits (dimming screens, ZTF, MOA) are for unlensed or macro-unmodelled
  sources and keep their wording.

## Next
- If a W3 search ever targets caustic-transient stars in arcs (e.g. the Dragon arc), inject through this map
  with the image's (κ, γ) and parity; saddle images cannot test W3 by vanishing.
