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
| Case | μ_macro | window \|β_x\| × \|β_y\|, θ_E | umbra (μ / μ_macro < 0.1) area, θ_E² | half-width along stretch / across | predicted 2√λ (on axis) | umbra half-width along stretch at β_y = 0.5 / 1.0 / 1.5 (ellipse from on-axis) | min μ / μ_macro, ρ = 0.1: six tracks (±100: 5) / whole window interior | no-lens control, ρ = 0.1 interior (point, whole map) |
|---|---|---|---|---|---|---|---|---|
| isolated | 1 | 6 × 3 | 12.49 | 1.99 / 1.99 | 2 / 2 | 1.93 / 1.73 / 1.32 (1.93 / 1.72 / 1.31) | 0 / 0 | 0.995–1.005 (0.79–1.23) |
| minimum | 10 | 1.9 × 3 | 4.48 | 0.63 / 2.00 | 0.632 / 2 | 0.62 / 0.60 / 0.55 (0.61 / 0.55 / 0.42) | 0 / 0 | 0.999–1.001 (0.89–1.11) |
| minimum | 100 | 1 × 3 | 1.55 (2.23 incl. side regions) | 0.20 / 2.00 | 0.20 / 2 | 0.195 / 0.195 / 0.195 (0.19 / 0.17 / 0.13) | 0 / 0 | 0.998–1.002 (0.89–1.33) |
| saddle | −10 | 1.9 × 3 | 0 | — | — (formula covers minima only) | — | 0.85 / **0.79** | 0.999–1.001 (0.89–1.11) |
| saddle | −100 | 1 × 3 | 0 | — | — (formula covers minima only) | — | 0.96 / **0.89** | 0.998–1.002 (0.89–1.33) |

Sampling noise (post hoc, after review): the ray grid is regular, so single-pixel μ aliases. The no-lens
control (`sign = 0`, same grids) gives the noise range in the last column; the ρ = 0.1 range excludes the
convolution's padded borders (second review: with them it was 0.89–1.12); the ρ = 0.1 track minima use the same
interior (third review; values unchanged). Fifth review: the windows now span whole bins (μ_macro = 10 window
1.897 → 1.9 θ_E; before, the real bin width was 0.25 % off `pixel` and μ was biased by that much; the area
4.47 → 4.48 and the ±10 control ranges moved, nothing else), the minimum over the whole window interior is
stored next to the six tracks (for the saddles it is deeper: 0.79 / 0.89 against 0.85 / 0.96), and the umbra
half-width along the stretch is stored at β_y = 0.5, 1.0 and 1.5. Zero-flux umbrae and their widths, and the
ρ = 0.1 saddle minima, are far outside the control; point-source saddle minima are not quoted.

## What this shows
- **Minimum-parity macro images:** the W3 "vanishing" survives macro magnification. The umbra keeps zero flux,
  but along the macro stretch it narrows to the predicted on-axis half-width 2√λ1 (≈ 2 / √μ_macro when λ2 = 1);
  across it stays 2 θ_E. It is **not** an ellipse: off axis its half-width along the stretch stays near 2√λ1
  (0.62–0.55 at μ_macro = 10 for β_y = 0.5–1.5; a constant 0.195 at 100), wider than an ellipse with the
  on-axis half-widths. So a source moving along the stretch at |u₀| ≲ 1.5 vanishes for about
  2 × 2√λ1 t_E, much shorter than the isolated 2 t_E √(4 − u₀²); a source moving across the stretch vanishes
  only if |u₀| < 2√λ1 (0.63 / 0.20; the u₀ = 1 track at μ_macro = 10 and the u₀ = 0.5 track at 100 do not dim).
- **Saddle-parity macro images:** there is **no umbra** inside the windows (|β_x| ≤ 1.9 / 1 θ_E, |β_y| ≤ 3 θ_E;
  no pixel below 0.1 μ_macro). For ρ = 0.1 the deepest dimming in the window interior is 21 % at μ_macro = −10
  and 11 % at −100 (control interior within ±0.2 %). That looks like ordinary microlensing near a saddle (the
  known mimic in docs/exotic_lensing.md), not a vanishing. Outside the windows was not sampled.
- Hence the isolated-lens W3 injections (`inject_light_curve`) are not a fair approximation inside arcs: for
  minima they overstate the cross-section for motion across the stretch (by 1 / √λ1) and the duration for motion
  along it, and in saddles they predict a vanishing that does not occur. The frozen rule covers < 10 % (revise) and ≥ 50 % (keep); the measured hole
  area at μ_macro = 100 is 12 % (central) or 18 % (all regions) of the isolated value, between the two, so the
  rule gives no verdict. The doc is revised on the simulation itself (deviation, labelled: the rule did not
  decide it).

## Deviations from the plan (labelled)
- Only the negative-mass sign was measured (positive mass appears only in gates G2 / G3′), and only ρ = 0.1 (no
  ρ = 0.01 disk: the 0.005–0.01 θ_E pixels cannot resolve it); track minima are stored, not full light curves.
- The rule's outcome falls between its two thresholds (12–18 %); the revision of docs/exotic_lensing.md is a
  judgement on the maps, not the rule's verdict (fifth review).
- The u₀ = 1 track across the stretch is outside the ±1 θ_E window interior at μ_macro = ±100 (`null` in
  `result.json`), so those rows take the minimum over 5 of the 6 tracks.
- G2's thresholds (local μ > 5 μ_macro in a 5 × 5 pixel maximum, μ > 20 μ_macro within 2 pixels, both ≥ 90 %)
  were set with the code, after the plan's "within two source pixels" wording; both fractions came out 100 %.
- G1 compares radial-bin means, as the plan's "radial profile" says; single pixels scatter by up to ±9 % near
  r = 2.5–2.9 (aliasing, see the control; not stored).
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
