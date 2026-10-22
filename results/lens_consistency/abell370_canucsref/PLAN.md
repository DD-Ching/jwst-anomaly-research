# Abell 370: CANUCS DR1 reference extra-image test, gate K first (plan, frozen)

Frozen 2026-10-22 (commit time of this file), before any Abell 370 counter-image prediction, CANUCS DR1 flux or pixel
was read. Follows `../macs0416_canucsref/RESULT.md` (handoff: "then Abell 370 (gate K first; expected to fail)").
Script: `scripts/a370_canucsref.py`; tests `tests/test_a370_canucsref.py`.

## Hypothesis, ordinary explanations, expected value
- Question (method): can the pinned CANUCS Abell 370 Lenstool model (`abell370-canucs`, Gledhill et al. 2025, D-044)
  serve the extra-image test, i.e. does it reproduce its own 115 images in the image plane (gate K)?
- Expected: **fail**. D-044 measured `validate` rms 2.32″ (max 17.9″); the model was fitted in the source plane.
  `m0417_images.gate_k` differs from `validate` (0.25″ deflection grid, best of the solved images per family, 3.0″
  match), so the number can differ, but not by a factor > 2 without a bug.
- Ordinary explanations for a fail: source-plane fitting (image-plane positions not constrained), our solver near
  critical curves, the fixed-z ASSUMPTION below. None of them is a lensing anomaly; a fail is a model-route result.
- What would change the plan: gate K passes → stages 2–3 of `../macs0416_canucsref/PLAN.md` run in the next cycle
  with the CANUCS DR1 Abell 370 catalogue and F200W mosaic pinned by sha256 first. Gate K fails → the model is not
  used for the extra-image test, no prediction or flux is read, and the Abell 370 line of this route stops until an
  image-plane-validated model exists (the HFF CATS v4 list already fails D-035: rms 10.7″, map only).

## Inputs
- model_prediction: `abell370-canucs`, the CANUCS DR1 Lenstool best-fit parameter file (pinned with its
  multiple-image and input files, D-044).
- observed: its 115-image multiple-image file (`lensmodel.ABELL370_CANUCS["arcs.dat"]`), copied to
  `canucs_a370_images.ecsv`. **ASSUMPTION:** every image is `z_spec` (the file's z are used as fixed), as for
  MACS0416.

## Frozen rule
1. **Gate K** (`m0417_images.gate_k` via `m0417_images.run_field`, 0.25″ grid as on the earlier fields): image-plane
   rms of the 115 images ≤ 1.0″ and ≥ 80 % reproduced within 3.0″. Fails → stop.
2. Reported, not rules: the per-family rms and the `post_hoc` drop-unreproduced check that `run_field` writes.
