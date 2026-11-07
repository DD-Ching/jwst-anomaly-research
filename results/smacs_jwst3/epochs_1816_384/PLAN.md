# Two-epoch check of SMACS `f200w_1816` and `f200w_384` (2736 → VENUS 6882 o057)

Frozen before any catalog is read (claim commit). TASKS Next M2 "Second epochs"; the open records left by
D-022 addenda 6–7 (point-source screen not validated) and D-110 (no grism coverage).

## Hypotheses
- H-star: both compact no-fit sources (DJA `flux_radius` 2.276 / 1.906 px, ref_mag 25.90 / 26.31) are faint
  Galactic stars. A halo or thick-disc star at 26 mag (several kpc to tens of kpc) would move ≳ 1–10 mas/yr,
  i.e. ≳ 4–40 mas over the 4.0 yr baseline (model_prediction, order of magnitude only).
- H-gal: compact background galaxies (no motion; constant flux).
- H-var: AGN or transient (flux change, no motion).

## Ordinary explanations
Frame-tie error, centroid noise at S/N ≈ 10–30, deblending differences between pipeline versions, and
zero-point differences (cancel in the reference median, D-017).

## Method (D-017 tool, unchanged)
`scripts/epoch_compare.py`, epoch 1 = program 2736 o001 level-3 CAT (jwst 3.0.0, as #766), epoch 2 =
`jw06882-o057_t057_nircam_clear-<band>` for every band shared with 2736 (F150W, F444W; F200W/F277W/F356W
only if VENUS has them). Defaults: references S/N ≥ 30 within 60″, match 0.3″.

## Frozen outcome rules (ASSUMPTIONs)
- Motion: residual ≥ 3σ (tie + centroid error) in **both** F150W and F444W, with the same direction within
  45°, → "moves: star" (H-star supported). Residual < 2σ in both with a 3σ bound ≤ 15 mas → "no motion
  detected" (bound stated; H-star at ≳ 4 mas/yr disfavoured, not excluded at smaller μ).
  Anything else → "undetermined".
- Variability: |Δm − reference median| ≥ 5 × reference scatter in both bands → "variable".
- Target rows must be mutual nearest neighbours in both epochs; otherwise "not measured" for that band.

## What would change the plan
A "moves" result closes the record as a star (no lens interest). "No motion" + no variability leaves a
compact galaxy (closes the star route; record stays open only for spectroscopy). "Variable" → `/vet-candidate`.
