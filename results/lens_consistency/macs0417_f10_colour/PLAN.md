# MACS0417 family 10: is CANUCS 1115938 the third image? Catalogue colour test (plan, frozen 2026-10-22 02:27Z before any measurement)

**Question.** `../macs0417_f10_zscan` (outcome A′) puts the third family-10 image on CANUCS DR1 1115938 only at
z 1.975–2.100, above 1115938's photo-z interval 1.72–1.88. Is 1115938 the third image at a z the photo-z
underestimates (H1), or an unrelated z ≈ 1.8 galaxy (H2)?

**Prediction.** Lensing is achromatic: under H1, 1115938's spectral energy distribution equals 10.1's and 10.2's up to
a flux ratio. Under H2 the colours may differ. Equal colours are necessary for H1, not sufficient.

**Method (observed fluxes; derived χ²).**
- Catalogue: CANUCS DR1 MACS0417 photometry (sha256 ca6b48aa…90dd). Fluxes `FLUX_COLOR03_TOTAL_<band>` and errors,
  bands F090W, F115W, F150W, F200W, F277W, F356W, F410M, F444W (JWST wide + F410M); a band is used for a pair only
  when both fluxes and errors are finite and both errors > 0.
- Counterparts: each Jauzac+19 image (JWST frame, `../macs0417_frame` offset) → the nearest CANUCS source within
  0.3″ (ASSUMPTION); 1115938 by id.
- Error floor: σ² = err² + (0.05 f)² per band (ASSUMPTION, a common photo-z floor).
- Pair χ²: min over r of Σ_band (f_a − r f_b)² / (σ_a² + r² σ_b²), dof = n_bands − 1; p from χ². "Differ" if
  p < 0.01 (ASSUMPTION); needs n_bands ≥ 5, else "untestable".
- Known cases: every Jauzac family other than 10 with ≥ 2 counterparts, first two images by id, both with F200W
  S/N ≥ 10. Gate (ASSUMPTION): ≥ 8 testable pairs and ≥ 80 % "consistent". If the gate fails, the test decides nothing.
- Target pairs: (10.1, 1115938) and (10.2, 1115938); (10.1, 10.2) reported as a reference.
- Decision function `decide(gate, target_classes)`: gate fails → "not validated"; any target pair "differ" →
  "H2 favoured: colours differ" (the #453 flag then stands as `absent` at the photo-z and the A′ landing is a chance
  alignment); both "consistent" → "H1 not excluded"; else "untestable".

**What would change the plan.** "H2 favoured" sends family 10's third image to a second-model test (RELICS Lenstool)
as a model-residual question; "H1 not excluded" closes the #453 flag as ordinary pending a z_spec.
