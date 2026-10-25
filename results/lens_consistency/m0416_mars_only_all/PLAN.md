# MACS0416: MARS-only images in every CANUCS family (plan, frozen)

Frozen 2026-10-25 (commit time of this file), before any MARS image position of the families new to this stage, or
any pixel / flux at one, was read. Extends D-034 addendum as (`../m0416_mars_only_phot/`), which searched only the 40
families with a CANUCS extra prediction (its RESULT "Limits": "MARS-only images of families without a CANUCS extra
prediction were not searched"). Script `scripts/m0416_mars_only_all.py`; test `tests/test_m0416_mars_only_all.py`
(both before the run).

## Inputs (as addendum as)
- model_prediction: MARS HFF free-form maps of MACS0416 (Cha & Jee 2023; zip pinned in `scripts/mars_hff_inputs.py`;
  1.515″ pixels; gate H passes on 77 held-out images, D-034 addendum ap).
- observed: CANUCS multiple-image list (`../macs0416_canucsref/canucs_m0416_images.ecsv`); CANUCS extra predictions
  (`../macs0416_canucsref/predicted_macs0416.ecsv`, model_prediction, used only to define "MARS-only").
- observed: CANUCS DR1 MACS0416 photometry catalogue and F200W 40 mas bgsub image (sha256 pinned in
  `scripts/m0416_canucsref.py`).

## Question, ordinary explanations, what changes the plan
- Question: in the CANUCS families without a CANUCS extra prediction, does MARS predict images that no observed image
  accounts for, and is there flux at them?
- Ordinary explanations for "no flux": the 1.515″ MARS map is unreliable near critical curves and in the BCG cores
  (positions and |μ| there are model geometry); a demagnified or central image below the depth ("undetectable"); BCG /
  member light residuals.
- Ordinary explanations for "flux": an unrelated galaxy or member within 1.5″ ("confused" when brighter than
  3 × F_pred); flux of a catalogued image of the same source (checked post hoc as in addendum as).
- What changes the plan: a target `absent` with gate C2 passing → `lensing_residual_input`, vetted with
  `/vet-candidate` (first test: map resolution, i.e. distance to the MARS critical curve, and whether CANUCS makes a
  counter-image there); a target `recovered` → possible uncatalogued image, reported only (post hoc distance to the
  family's observed images first); otherwise no reading. Gate R failing → stop, no reading.

## Frozen rules (functions with unit tests before the run)
1. Families: every system of the CANUCS image list whose images all back-trace with finite source positions.
   `split_families` partitions them into the addendum-as set (systems with a CANUCS extra prediction) and the new set
   (all others).
2. Solve, pair and select targets exactly as addendum as rules 1–2 (`m0416_mars_only_phot.solve`: frame offset of
   addendum ap, mean back-traced source, `find_images` on the 0.25″ grid, `match_pairs` within 3.0″, `mars_only`
   farther than 3.0″ from every CANUCS extra prediction; for new families the extra list is empty).
3. **Gate R** (known case, reproduction): the addendum-as set must give the addendum-as target list, i.e. one target in
   D115a within 0.1″ of (64.044742, −24.060901) and no other (`gate_r`). Fails → stop.
4. Photometry, controls, gate C2 and the local-noise rule: addendum as rules 3–5 unchanged (`m0417_arcref` with the
   CANUCS catalogue, MARS positions and |μ|, 1.5″ search radius, ASSUMPTION; `gate_c2` frozen thresholds;
   `m1423_canucsref.apply_local`, seed 20261022), over the controls of **all** solved families. Gate C2 of the new
   families' controls alone is reported, not a rule.
5. Reading (`m0416_mars_only_phot.reading`) per new-family target: `class_local` "absent" and C2 passes →
   `lensing_residual_input`; "recovered" → `seen_where_only_mars_predicts`; otherwise `no_reading`.

## Stop rules
- No rule change after the first MARS position of the new families; departures labelled deviation or post hoc.
- Every target with a reading other than `no_reading` gets a cutout (and its family's observed images), looked at
  before RESULT.md is written; if there are more than 12, the sheet shows the 12 with the largest F_pred / σ_loc.
- Data: the MARS zip and the CANUCS files are deleted after the run; only `summary.json`, `controls.ecsv`,
  `targets.ecsv` and a small PNG are stored.
