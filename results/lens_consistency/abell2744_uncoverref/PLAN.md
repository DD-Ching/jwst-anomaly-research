# Abell 2744: UNCOVER-reference extra-image test, stage 1 (gate K) (plan, frozen)

Frozen 2026-10-23 (commit time of this file), before this stage's gate K run. Follows the handoff of D-034 addendum i
(#475): the five CANUCS-reference fields are used, so the next field needs another reference catalogue; Abell 2744
with UNCOVER is the first to check, gate K first. Script: `scripts/a2744_uncoverref.py`; tests
`tests/test_a2744_uncoverref.py`.

## Hypothesis, ordinary explanations, what changes the plan
- Question (method): does the UNCOVER v2.0 map model (`abell2744-uncover2`, D-041 / `../abell2744-uncover2/`)
  reproduce its spectroscopic multiple images under the same gate K as the CANUCS fields? A pass lets a later cycle run
  the extra-image stages (predicted / control tables, reference forced test, gate C2 with the frozen shared-pair rule,
  flag) on Abell 2744 with an UNCOVER photometric reference. A fail stops the Abell 2744 extra-image route with this model.
- Ordinary explanations for a fail: map resolution (100 mas) and the grid search near critical curves, the source
  estimator (`m0417_images.Solver`, not the WLS source of `a2744_uncover.py`), knot-level ids. Exotic readings are not
  considered at this stage.

## Inputs
- model_prediction: UNCOVER DR4 lensing maps v2.0 (Furtak+2023 model update), `alpha_x` / `alpha_y`
  low-resolution 100 mas, sha256 pinned in `scripts/a2744_uncover.py` (`UNCOVER_V2`), z_lens 0.308, flat H0 70,
  Om 0.3 (ASSUMPTION as before; maps at D_LS/D_S = 1). Registered here as `abell2744-uncover2-maps`.
- observed: the model's own multiple-image table (same release, `UNCOVER_V2["images"]`, 186 rows). Rows with
  `Redshift_type` = spec are `z_spec` = true (172); `phot` rows (9) get a fitted z as in `m0417_images`; rows typed
  `None` (5) are left out. Ids are rewritten so `lensmodel.image_family` gives the knot-level source of
  `a2744_uncover.uncover_system` ("S.KI" -> "S.K" + letter I; "S.I" -> "S" + letter I).
- ASSUMPTION: frame offset 0 (the table is in the UNCOVER JWST frame).

## Known before freezing (disclosed)
- `../abell2744-uncover2/summary.json`: own-list image-plane rms 0.84″ on 181 images (WLS source, every system at its
  table z, 1.5″ match). Gate K here uses a different estimator and only the 172 spectroscopic images, so that number
  does not decide it, but the outcome is not blind. These images are the model's own constraints: a pass validates
  the model on its constraints only (as BUFFALO K2, D-034 addendum h), not on an independent list.

## Frozen rules
1. **Gate K** (`m0417_images.gate_k` via `m0417_images.run_field`, 0.25″ grid): image-plane rms of the spectroscopic
   images ≤ 1.0″ and ≥ 80 % reproduced within 3.0″. Fails → stop (no extra-image stage with this model).
2. Reported, not rules: `extra_predicted` / `fewer_predicted` families, the fitted z of the phot systems, and
   `drop_unreproduced` (post hoc in `run_field`). No forced photometry or flag in this stage.
3. Reference catalogue (availability only, no flux read): whether an UNCOVER photometric catalogue with per-band
   fluxes, errors and segmentation, and a BCG-subtracted mosaic, is reachable from the cloud session and pinnable by
   sha256. Recorded in `RESULT.md`; the stage-2 plan freezes its column mapping.
