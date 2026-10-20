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
| Case | μ_macro | umbra (μ / μ_macro < 0.1) area, θ_E² | half-width along stretch / across | predicted 2√λ (on axis) | min μ / μ_macro over u₀ = 0, 0.5, 1 tracks, ρ = 0.1 (interior; ±100: 5 of 6 tracks) | no-lens control, ρ = 0.1 interior (point, whole map) |
|---|---|---|---|---|---|---|
| isolated | 1 | 12.49 | 1.99 / 1.99 | 2 / 2 | 0 | 0.995–1.005 (0.79–1.23) |
| minimum | 10 | 4.47 | 0.625 / 2.00 | 0.632 / 2 | 0 | 1.000–1.008 (0.89–1.23) |
| minimum | 100 | 1.55 (2.23 incl. side regions) | 0.20 / 2.00 | 0.20 / 2 | 0 | 0.998–1.002 (0.89–1.33) |
| saddle | −10 | 0 | — | — (formula covers minima only) | 0.85 | 1.000–1.008 (0.89–1.23) |
| saddle | −100 | 0 | — | — (formula covers minima only) | 0.96 | 0.998–1.002 (0.89–1.33) |

Sampling noise (post hoc, after review): the ray grid is regular, so single-pixel μ aliases. The no-lens
control (`sign = 0`, same grids) gives the noise range in the last column; the ρ = 0.1 range excludes the
convolution's padded borders (second review: with them it was 0.89–1.12); the ρ = 0.1 track minima use the same
interior (third review; values unchanged). Zero-flux umbrae and their widths, and
the ρ = 0.1 saddle minima, are far outside it; point-source saddle minima (0.69 / 0.86) are not resolved from the
point-source control and are not quoted.

## What this shows
- **Minimum-parity macro images:** the W3 "vanishing" survives macro magnification. The umbra keeps zero flux,
  but along the macro stretch it narrows to the predicted half-width 2√λ1 (≈ 2 / √μ_macro when λ2 = 1);
  across it stays 2 θ_E. A source crossing along the stretch axis vanishes for a time shorter by √λ1 than the
  isolated formula 2 t_E √(4 − u₀²) gives.
- **Saddle-parity macro images:** there is **no umbra** (no pixel below 0.1 μ_macro in either case). For
  ρ = 0.1 the deepest dimming is 15 % at μ_macro = −10 and 4 % at −100 (control interior within ±0.8 %). That looks like ordinary microlensing near a saddle (the known mimic in docs/exotic_lensing.md),
  not a vanishing.
- Hence the isolated-lens W3 injections (`inject_light_curve`) are a fair approximation only for minimum
  images and only for motion across the stretch. The plan's decision rule: hole area at μ_macro = 100 is 12 %
  of the isolated value in the minimum (between the 10 % and 50 % thresholds) and 0 % in the saddle, so the
  doc's W3 amplitudes for arcs are revised (docs/exotic_lensing.md notes).

## Deviations from the plan (labelled)
- Only the negative-mass sign was measured (positive mass appears only in gates G2 / G3′), and only ρ = 0.1 (no
  ρ = 0.01 disk: the 0.005–0.01 θ_E pixels cannot resolve it); track minima are stored, not full light curves.
- The u₀ = 1 track across the stretch is outside the ±1 θ_E window interior at μ_macro = ±100 (`null` in
  `result.json`), so those rows take the minimum over 5 of the 6 tracks.
- G2's thresholds (local μ > 5 μ_macro in a 5 × 5 pixel maximum, μ > 20 μ_macro within 2 pixels, both ≥ 90 %)
  were set with the code, after the plan's "within two source pixels" wording; both fractions came out 100 %.
- G1 compares radial-bin means, as the plan's "radial profile" says; single pixels scatter by up to ±9 % near
  r = 2.5–2.9 (aliasing, see the control).
- The no-lens control and the alignment of both ray grids to the fine box (no seam) were added after review; the
  re-runs left every gate and umbra number unchanged. `macro_hole_halfwidths` now raises outside minimum images
  (the saddle "no umbra" is a simulated result, not this formula).

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
