# Abell 370: BUFFALO model A extra-image test with CANUCS DR1 reference fluxes (plan, frozen)

Frozen 2026-10-22, before any counter-image prediction, CANUCS flux or JWST pixel of Abell 370 was read for this test.
Follows `../abell370_buffalo/RESULT.md` (gate L passes; gate K2 on the model's own gold list passes, rms 0.910″).
Script: `scripts/a370_buffalo_extra.py`; tests `tests/test_a370_buffalo_extra.py`.

## Hypothesis, ordinary explanations, expected value
- Question 1 (method): does the CANUCS DR1 reference route with the local-noise and containing-source rules pass gate
  C2 on Abell 370 with a map model fitted in the image plane?
- Question 2 (science, only if C2 passes): does model A predict a magnified counter-image of a gold system, in the
  footprint and unmatched to any gold image, where CANUCS F200W shows no source (flag, rule 4)?
- Ordinary explanations for a flag: the image is predicted only by this model (model systematics; other published
  models disagree), it is below the local depth, blended with a cluster member or the BCG residual, or the system's z
  is wrong. Each is tested by `/vet-candidate` before any reading.
- Expected: C2 passes (as on MACS1149 and MACS0416) and 0–2 flags, each with an ordinary reading.
- What would change the plan: C2 fails → the route is not validated on Abell 370 and no flag is read. A flag →
  contact sheet, then `/vet-candidate`, before any wording beyond "flag".

## Inputs
- model_prediction: `abell370-buffalo-a-maps` (pinned in `lensmodel`; as `../abell370_buffalo/PLAN.md`).
- observed: the 104 images of the BUFFALO gold list (`a370_buffalo.BUFFALO_GOLD`, sha256 pinned), copied to
  `buffalo_sl_gold.ecsv`. **ASSUMPTION:** z fixed. These are the model's own constraints (as for MACS0416).
- observed: CANUCS DR1 Abell 370 cluster-field photometry catalogue and F200W 40 mas BCG-subtracted mosaic; file names
  and URLs in the script (`BASE`, `CAT_NAME`, `IMAGE_NAME`; DOI 10.17909/18nv-np70). Their sha256 are pinned in the
  script (`CAT_SHA`, `IMAGE_SHA`) at first download, from the file bytes only, before the script runs; the script
  refuses to run while they are unset. Not committed; in cloud sessions deleted after the run.
- observed: JWST pipeline level-3 F200W catalogue and mosaic of 1208 `o002_t001` (S3; `TARGPROP` ABELL370; byte-range
  cutouts only).
- ASSUMPTION: frame offset 0 between the HST-frame image list and the JWST/CANUCS frame (not measured for Abell 370).

## Frozen rules (same code as `../macs0416_canucsref/PLAN.md`)
1. **Gate K** (`m0417_images.gate_k`, 0.25″): rms ≤ 1.0″ and ≥ 80 % within 3.0″ (repeats gate K2). Fails → stop.
2. Predicted and control tables from `m0417_extra.run_field` unchanged; its pipeline forced test and gate C are
   reported, not used.
3. Stage 3 = `m0416_canucsref.stage3` (rules 3–5 of `../macs0416_canucsref/PLAN.md`, same functions and constants).
4. Gate C2 and the flag on the classes after the local-noise and containing-source rules.
5. Reported, not rules: the same without the containing-source rule (`without_contain/`).
