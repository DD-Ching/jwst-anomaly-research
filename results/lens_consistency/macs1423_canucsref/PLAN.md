# MACS1423: extra-image test with CANUCS DR1 reference fluxes and a local-noise rule (plan, frozen)

Frozen 2026-10-22 (commit time of this file), before any CANUCS DR1 MACS1423 flux or pixel was read. Follows
`../macs0417_canucsref/RESULT.md` (gate C2 passed 37 / 40 there) and its handoff: the same reference on MACS1423 with
a local-noise `undetectable` rule frozen first (on MACS0417 the field-median depth made a BCG-core prediction look
testable when the local noise was 1.4× its flux). Script: `scripts/m1423_canucsref.py`.

## Hypothesis, ordinary explanations, expected value
- Question 1 (method): does the CANUCS DR1 reference pass gate C2 on a second field? A pass makes it a validated
  reference for the remaining CANUCS fields; a fail means MACS0417's pass was field-specific.
- Question 2 (field): are the CANUCS v1 extra predictions of MACS1423 present at the predicted flux?
- Expected before running: MACS1423 has 2 testable extra predictions (system 1 `missing` |μ| 17.3, `candidate`
  |μ| 8.8; `../macs0417_extra/predicted_macs1423.ecsv`), both with F200W flux in the earlier forced test
  (`recovered` / `confused`). The 2 `demagnified` central predictions (μ 0.39 / 0.43) are excluded by rule 8, as
  on MACS0417. So a flag is unlikely; the main value is question 1.
- Ordinary explanations for an `absent`: model μ or position error near critical curves, low surface brightness,
  deblending or BCG-model residuals, a reference with unrelated light. Exotic readings are not considered.
- What would change the plan: C2 fails → do not use CANUCS DR1 as the reference on other fields without a new
  gate. C2 passes and a flag → contact sheet, then `/vet-candidate`.

## Inputs
- observed: CANUCS DR1 MACS1423 cluster-field photometry catalogue
  `hlsp_canucs_jwst-hst_multi_macs1423-clu_multi_v1_photometry-cat.fits.gz` and the F200W 40 mas BCG-subtracted
  mosaic `hlsp_canucs_jwst_nircam_macs1423-clu-40mas_f200w_v1_bgsub-sci.fits.gz` (DOI 10.17909/18nv-np70; sha256 in
  `summary.json`; streamed to the scratch disk, not committed).
- model_prediction: `../macs0417_extra/predicted_macs1423.ecsv`, `control_macs1423.ecsv` and the `unpredicted`
  list of its `summary.json` (CANUCS v1 maps; unchanged). The published image list is used with no frame offset
  (ASSUMPTION: the F200W peaks lie a median 0.06″ from it, `../macs0417_extra/RESULT.md`).

## Frozen rules
1. Rules 1–8 of `../macs0417_arcref/PLAN.md` with the column mapping of `../macs0417_canucsref/PLAN.md`, both
   verbatim (same code: `m0417_arcref.run`, `m0417_canucsref.to_cat`).
2. **Local noise (new; ASSUMPTIONs):** at every control and target position, σ_loc = σ-clipped standard deviation
   of 200 random 0.3″-radius aperture sums on the BCG-subtracted F200W mosaic, centres 1.5–4″ away, redrawn when
   within 0.5″ of a CANUCS source or touching a NaN pixel (method of `../macs0417_canucsref/posthoc.json`; seed
   20261022). A test classed `recovered`, `confused`, `absent` or `ambiguous` by rule 1 becomes
   `undetectable_local` when F_pred < 5 σ_loc, and an `absent` becomes `ambiguous` when F_pred < 10 σ_loc.
   F_pred is a total flux and σ_loc an aperture noise, so the rule is lenient for extended images (stated, not
   corrected).
3. Gate C2 and the flag (rules 7–8) are evaluated on the classes after rule 2.
4. Reported for comparison, not a rule: gate C2 and flags without rule 2.
