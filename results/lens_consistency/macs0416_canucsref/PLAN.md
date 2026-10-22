# MACS0416: CANUCS DR1 reference extra-image test with the containing-source rule (plan, frozen)

Frozen 2026-10-22 (commit time of this file), before any MACS0416 counter-image prediction, CANUCS DR1 flux or
pixel was read. Follows `../macs1149_canucsref/RESULT.md` (gate C2 passed on a third field; handoff: Abell 370 and
MACS0416 next with this reference and rule). Script: `scripts/m0416_canucsref.py`; tests
`tests/test_m0416_canucsref.py`.

## Why MACS0416 before Abell 370
The pinned CANUCS Abell 370 model is a source-plane fit with image-plane rms 2.32″ on its 115 images (D-044), so
gate K (≤ 1.0″) is expected to fail there. The CANUCS MACS0416 model is an image-plane fit, rms 0.512″ on 303 images
(D-044), so gate K is expected to pass and the test can return a reading. Abell 370 stays next, with gate K run first.

## Hypothesis, ordinary explanations, expected value
- Question 1 (method): does the CANUCS DR1 reference with the containing-source rule pass gate C2 on a fourth field,
  with many more controls (111 families) than MACS1149 (23 gold images)?
- Question 2 (field): are the extra images that the CANUCS MACS0416 model predicts present at the predicted flux?
- Ordinary explanations for an `absent`: model μ or position error near critical curves (MACS0416 has long merging
  arcs at fold caustics, e.g. CATS system 26 in `docs/fields/macs0416.md`), low surface brightness, deblending or
  BCG-model residuals, ICL, a reference with unrelated light, faint images below the local noise. Exotic readings are
  not considered.
- What would change the plan: gate K fails → the model is not used, stop. C2 fails → this reference route is not
  validated on MACS0416; no flag is read. C2 passes and a flag → contact sheet, then `/vet-candidate`.

## Inputs
- model_prediction: `macs0416-canucs`, the CANUCS DR1 Lenstool best-fit parameter file of Rihtarsic et al. 2025
  (pinned with its multiple-image file and input file, D-044). **Deviation from the three earlier fields:** a
  parameter-file model, not the DR1 deflection maps; D-044 pinned and validated it (χ²pos within 4 % of Lenstool's).
- observed: the 303 images of the model's multiple-image file (`lensmodel.MACS0416_CANUCS["arcs.dat"]`), copied to
  `canucs_m0416_images.ecsv`. **ASSUMPTION:** every image is treated as `z_spec` (the file's z is used as fixed; the
  model treats these z as fixed constraints). The images are the model's own constraints, so gate K checks the
  reproduction, not an independent prediction (as for MACS0417's Jauzac+19 list).
- observed: CANUCS DR1 MACS0416 cluster-field photometry catalogue and F200W 40 mas BCG-subtracted mosaic (DOI
  10.17909/18nv-np70; sha256 pinned in the script and written to `summary.json`; not committed). JWST pipeline
  level-3 F200W catalogue and mosaic of 1208 `o004_t002` (S3; the predicted / control tables and stage-2 forced test
  only).
- ASSUMPTION: frame offset 0 between the image list and the JWST/CANUCS frame (D-044 measured (−0.008″, +0.052″),
  not pinned). Reported, not a gate: the median offset of each image from its rule-1 CANUCS counterpart centroid.

## Frozen rules (same code as `../macs1149_canucsref/PLAN.md`)
1. **Gate K** (`m0417_images.gate_k`, 0.25″ deflection grid as on the earlier fields): image-plane rms of the 303
   images ≤ 1.0″ and ≥ 80 % reproduced within 3.0″. Fails → stop.
2. Predicted and control tables from `m0417_extra.run_field` unchanged (grid 0.1″, pipeline catalogue depth); its
   pipeline forced test and gate C are reported, not used.
3. Rules 1–8 of `../macs0417_arcref/PLAN.md` with the column mapping of `../macs0417_canucsref/PLAN.md`, rule 2
   (local noise) of `../macs1423_canucsref/PLAN.md` and rule 4 (containing source) of `../macs1149_canucsref/PLAN.md`,
   all verbatim (same functions and constants).
4. Gate C2 (rule 7) and the flag (rule 8) on the classes after the local-noise and containing-source rules.
5. Reported for comparison, not rules: gate C2 and flags without the containing-source rule, and without the
   local-noise rule.
