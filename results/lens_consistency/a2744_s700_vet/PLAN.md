# Abell 2744 system 700.1 (UNCOVER 78): vetting plan

Frozen before any measurement (2026-10-20, cloud routine). Follows `../a2744_s700_uncover/RESULT.md` ("Reading").

**Question.** Are 700.1a and 700.1b (UNCOVER 78.x1 / 78.x2, three nested knot pairs) two lensed images of one
source across a critical curve that both models miss, or two parts of one image?

**Hypotheses and ordinary explanations.**
- H1 (ordinary, one image): a/b are opposite sides of one galaxy. Prediction: knot colours may differ between a
  and b.
- H2 (ordinary, fold pair from an unmodelled perturber): a cluster member near the arc adds a local critical
  curve. Prediction: equal colours of mirrored knots; a z ≈ 0.31 galaxy within a few arcsec.
- H3 (model residual, not yet explained): equal colours, fold geometry and no plausible perturber.

**Tests.**
1. Colours. Forced photometry (1″ median high-pass, r = 0.15″ apertures, local MAD noise) of the six knot positions
   in the UNCOVER `jw02561-o001_t003` mosaics, every available wide band (S3 byte ranges). For each knot pair,
   χ² of the colour difference (all bands, flux ratio a/b free) against the noise. ASSUMPTION: p < 0.01 for any
   pair means "colours differ" (supports H1); otherwise "colours consistent" (H1 not excluded: equal colours are
   necessary for H2/H3, not sufficient).
   Known case first: the same procedure on a catalogued UNCOVER multiple-image pair of similar brightness must
   give "consistent".
2. Perturber. The nearest galaxies within 5″ of the a/b midpoint in the DJA (or UNCOVER) photo-z / spec-z
   catalogue: a source with z within 0.1 of 0.308 (ASSUMPTION) and F200W brighter than 22 AB (ASSUMPTION) is a
   plausible perturber (supports H2 as possible, not proved).
3. Look at the multi-band stamps before writing anything.

**What changes the plan.** "Colours differ" closes 700.1 as one image (not a candidate). "Consistent" plus a
member perturber keeps it ordinary (H2) and records a perturber-model test as the next step. "Consistent" with no
perturber keeps it undetermined and records it in the candidates list for a lens-model expert.

**Specification (added 22:20Z, before any photometry was run).**
- Bands: F115W, F150W, F200W, F277W, F356W, F410M, F444W of `jw02561-o001_t003`. Stamp half-width 2″.
- Position: the UNCOVER catalogue position, re-centred once on the F200W high-pass peak within 0.1″; the same
  position in every band.
- Noise: per stamp, 200 random r = 0.15″ apertures with centres 0.6–1.8″ from the target on the high-pass image;
  σ = 1.4826 × MAD of their sums.
- Pair χ²: min over r of Σ_band (f_a − r f_b)² / (σ_a² + r² σ_b²), dof = n_bands − 1 (bands with finite values
  only); p from χ². "Differ" if p < 0.01 (ASSUMPTION).
- Known cases: every UNCOVER system other than 78 with a spectroscopic z, taking its first two images, when both
  have F200W S/N ≥ 10 (at most the first 20 such pairs by system order). Gate (ASSUMPTION): ≥ 80 % of them give
  "consistent". If the gate fails, the colour test is "not validated" and decides nothing.
- Perturber: DJA v7.2 eazy zout (`z_spec` if > 0, else `z_phot`) within 5″ of the 78.2 a/b midpoint; F200W
  magnitude from an r = 0.5″ aperture on the stamp minus its median beyond 1.5″ (ASSUMPTION; a crude total).
- Decision function `decide(colour_gate_passed, any_pair_differs, perturber_found)`:
  not validated → "undetermined (colour test not validated)" plus the perturber result; differs → "one image (H1)";
  consistent and perturber → "ordinary perturber possible (H2)"; consistent, no perturber → "undetermined (H3)".
