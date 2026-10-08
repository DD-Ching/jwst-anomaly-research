# Abell 2744 family 4, c images (4.1c underluminous, 4.2c undetected) — vetting record

## 2026-10-08 (cloud run, session 01F4SgDLLX7dpbxi9Kj6gEmC)

- **Source:** multiple-image family 4 of Abell 2744 (Bergamini+2023b, `abell2744-bergamini23`; D-036). The knots
  4.1 and 4.2 are two clumps of one galaxy (z 3.577, spectroscopic). The flag arose in the exotic-lens screens
  (D-031/D-034; docs/fields/abell2744.md):
  - predicted 4.2c (RA 3.57972, Dec −30.40835, μ 8.7) is undetected;
  - catalogued 4.1c is underluminous.
  This is the pattern a demagnifying lens would produce, so it was vetted before any interpretation.
- **Score:** none (not a ranked pipeline source). Its forced-photometry and flux-ratio classes are `derived`.
- **Data:** UNCOVER `jw02561-o001_t003` F150W, F200W, F277W, F444W `_i2d`, read by S3 byte ranges.

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Image edge / low weight | A | does not explain | all positions well inside the mosaic; finite ERR | observed |
| Diffraction spikes, bright-star neighbour | A | does not explain | none within 5″ (cutouts) | observed |
| Detected in one band only (cosmic ray, snowball) | A | does not explain | 4.1c detected in all four bands (S/N 16–77 after BCG subtraction) | observed |
| Light of the brightest cluster galaxy (BCG) | B | does not explain on its own | BCG modelled with elliptical isophotes (photutils 2.3.0, sigma-clipped) and subtracted. 0.5″ apertures give L(4.1c)/L(4.1a) = 0.12 / 0.24 / 0.21 (F150W/F277W/F444W); L(4.1b)/L(4.1a) = 0.88–0.94. A 0.6″ high-pass gives 0.13 / 0.12 / 0.18 | derived |
| Aperture or knot mismatch | B | does not explain | 4.1a/4.2a/4.1b/4.2b are compact knots, 0.49–0.53″ apart. The deficit holds at r = 0.25″ and 0.5″ | derived |
| 4.2c hidden near 4.1c (model offset) | B | consistent | 4.2c shifted by 4.1c's own model error (0.59″): high-pass S/N < 2 in all bands. The model puts the knots' c images 2.5″ apart, against 0.5″ in a/b | derived |
| Dust | C | partially | 4.1c is about 0.6 mag redder (F150W vs F444W ratio), but the deficit is ≥ 4× in F444W too, so dust alone does not explain it | derived |
| **Magnification error near a member's critical curve** | C | **explains** | μ(4.1c) changes with the nearby member potential 34423 (σ 199 km/s, 2.9″ from 4.2c): σ × 0.7–1.4 and r_cut × 0.5–3 give μ = 3.9, 5.6, 7.6, 14.8 (best), 28.7, and parity flips (−12, −4.7, −3.1). μ(4.1a) and μ(4.1b) change by < 10 %. The independent CATS v4.1 model gives μ(4.1c) = 7.3. The observed deficit needs \|μ\| ≈ 2–4, which is inside this range. The same perturbations move 4.2c by 0.2–1.95″ | model_prediction |
| Microlensing by intracluster stars | C | not needed | the knots are extended (~100 pc); not tested further | assumption |

**Verdict:** `known population: lensed images near a cluster member's critical curve`. The family-4 c-image flux
deficit and the undetected 4.2c are explained by the published model's magnification and position systematics
next to member 34423. That is ordinary galaxy-scale lensing that the best-fit model does not constrain. Status:
**explained; not a candidate.** Any exotic interpretation would be a `hypothesis` with no supporting evidence
after these tests.

**Rule for later fields (folded into the screens' vetting order):** before reporting an under- or overluminous
image, re-evaluate its μ under ±30 % changes of the nearest member potential, and under an independent model where
one exists. Treat a μ that moves by more than 2× as untestable.

**Regeneration:**
- BCG subtraction: scratch script, photutils isophotes on an 8″ stamp around potential 34423.
- μ variations: `dataclasses.replace` on `LensModel.components`.
- CATS: `MapLensModel.from_fits` on the `lensmodel.HFF_CATS["abell2744"]` maps.
- All inputs are pinned by sha256 (DECISIONS D-030, D-035).

**Open questions:** posterior (`bayes.dat`) spread of μ(4.1c); a model with 34423 free.
