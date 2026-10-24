# MACS0416: forced photometry at the MARS-only extra images (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any MARS image position of this stage or any pixel / flux at
it was read. Follows D-034 addendum ar (`../m0416_mars_extra/RESULT.md`): over the 40 families with a CANUCS
extra prediction, the MARS maps make one image matched to neither an observed image nor a CANUCS extra
prediction, in D115a (|μ_MARS| 29.2; position not stored). Script `scripts/m0416_mars_only_phot.py`; test
`tests/test_m0416_mars_only_phot.py` (both before the run).

## Inputs
- model_prediction: MARS HFF free-form maps of MACS0416 (Cha & Jee 2023; the zip pinned in
  `scripts/mars_hff_inputs.py`; 1.515″ pixels; gate H 0.514″ on 77 held-out images).
- observed: the CANUCS multiple-image list (`../macs0416_canucsref/canucs_m0416_images.ecsv`); the CANUCS extra
  predictions (`../macs0416_canucsref/predicted_macs0416.ecsv`, model_prediction, used only to define "MARS-only").
- observed: CANUCS DR1 MACS0416 photometry catalogue and F200W 40 mas bgsub image, the sha256-pinned files of
  `scripts/m0416_canucsref.py`.

## Question, ordinary explanations, what changes the plan
- Question: is there flux at the position where only the MARS model predicts an image?
- Ordinary explanations for "no flux": the MARS map is coarse (1.515″) near a galaxy-scale critical curve, so its
  |μ| and position there are unreliable (CANUCS does not predict the image); BCG / member light residuals.
- Ordinary explanations for "flux": an unrelated galaxy or cluster member within the radius (rule 4 calls that
  "confused" when brighter than 3 × F_pred).
- What changes the plan: target "absent" with gate C2 passing → a MARS-side `lensing_residual_input`, vetted with
  `/vet-candidate`, first test the map resolution; target "recovered" → a possible uncatalogued counter-image,
  reported only (needs a colour / morphology check against D115a before any reading); anything else → no reading,
  the D115 two-model route stops here.

## Frozen rules (functions with a unit test before the run)
1. Solve as addendum ar rule 2 (frame offset of addendum ap applied to the CANUCS positions, mean back-traced source,
   `find_images` on a 0.25″ grid), for the same 40 families. Pairs MARS image ↔ observed image: one-to-one, greedy
   by distance, within 3.0″ (`match_pairs`; same matching as addendum ar's `match_observed`).
2. Target = a MARS image paired with no observed image and farther than 3.0″ from every CANUCS extra prediction of
   its family (`mars_only`). Positions go back to the CANUCS frame by subtracting the addendum ap offset.
3. Photometry is `m0417_arcref` with the CANUCS catalogue (`m0417_canucsref.to_cat`), as in
   `../macs0416_canucsref/` (counterpart, units, reference with its S/N > 5, |μ| ≤ 50 and 3× spread rules,
   `contain=True`, classify with the catalogue depth), with these changes: each observed image's |μ| is the MARS
   |μ| of its paired MARS image (unpaired observed images are left out); every tested position is the MARS image
   position (controls too, so the controls test MARS's positions and |μ|); the search radius is 1.5″ for every
   position (ASSUMPTION: one MARS map pixel).
4. Controls: every unit of the 40 families at its first image's paired MARS position, |μ| = sum of the unit's MARS
   |μ|, reference from the family's other units. Gate C2: `m0417_arcref.gate_c2` with its frozen thresholds.
5. Then the local-noise rule (`m1423_canucsref.apply_local`, seed 20261022) on controls and targets.
6. Reading (`reading`): target `class_local` "absent" and C2 passes → `lensing_residual_input`; "recovered" →
   `seen_where_only_mars_predicts`; otherwise `no_reading`.

## Stop rules
- No rule change after the first MARS position of this stage; departures labelled deviation or post hoc.
- A cutout of every target (and of its family's observed images) is looked at before RESULT.md is written.
- Data: the MARS zip and the CANUCS files are deleted after the run; only `summary.json`, `controls.ecsv`,
  `targets.ecsv` and a small PNG are stored.
