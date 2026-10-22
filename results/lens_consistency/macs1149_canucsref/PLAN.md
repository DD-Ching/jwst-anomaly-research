# MACS1149: CANUCS DR1 reference extra-image test with a containing-source rule, third field (plan, frozen)

Frozen 2026-10-22 (commit time of this file), before any MACS1149 model prediction, CANUCS DR1 flux or pixel was
read. Follows `../macs1423_canucsref/RESULT.md` (gate C2 failed on arc-referenced controls, 2 / 4) and its handoff:
a rule that also accepts the source whose ellipse contains the position, frozen and gated on a field other than
MACS0417 and MACS1423. Script: `scripts/m1149_canucsref.py`; tests `tests/test_m1149_canucsref.py`.

## Hypothesis, ordinary explanations, expected value
- Question 1 (method): with the containing-source rule, does the CANUCS DR1 reference pass gate C2 on a field not
  used to design it? A pass makes reference + rule usable on the remaining CANUCS fields (Abell 370, MACS0416); a
  fail stops this reference route.
- Question 2 (field): are the CANUCS v1 MACS1149 extra predictions present at the predicted flux?
- Ordinary explanations for an `absent`: model μ or position error near critical curves, the frame offset between
  the HST-frame image list and the JWST/CANUCS frame, low surface brightness, deblending or BCG-model residuals, a
  reference with unrelated light. Exotic readings are not considered.
- What would change the plan: gate K fails → the model is not used, stop. C2 fails → stop this reference route.
  C2 passes and a flag → contact sheet, then `/vet-candidate`.

## Inputs
- model_prediction: CANUCS DR1 MACS1149 Lenstool v1 maps (best fit, 100 mas, D_LS/D_S = 1; registered as
  `macs1149-canucs-maps`, sha256 in `src/jwst_anomaly/lensmodel.py`), z_lens 0.543 (as `macs1149-cats`).
- observed: Sharon & Johnson HFF v4 gold image list (`../macs1149_sharon/sharon_gold_images.ecsv`, 23 images, all
  with spectroscopic z; copied here with `z_spec` = true). CANUCS DR1 MACS1149 cluster-field photometry catalogue
  and F200W 40 mas BCG-subtracted mosaic (DOI 10.17909/18nv-np70; sha256 in `summary.json`; not committed). JWST
  pipeline level-3 F200W catalogue and mosaic of 1208 `o008_t004` (S3, for the predicted / control tables only).
- ASSUMPTION: frame offset 0 between the image list and the JWST/CANUCS frame (as MACS1423). Reported, not a gate:
  the median offset of each gold image from its rule-1 CANUCS counterpart centroid.

## Frozen rules
1. **Gate K** (model validation, `m0417_images.gate_k`): image-plane rms of the 23 gold images ≤ 1.0″ and ≥ 80 %
   reproduced within 3.0″. Fails → stop.
2. Predicted and control tables from `m0417_extra.run_field` unchanged (CANUCS maps, pipeline catalogue depth, grid
   0.1″); its pipeline forced test and gate C are reported, not used.
3. Rules 1–8 of `../macs0417_arcref/PLAN.md` with the column mapping of `../macs0417_canucsref/PLAN.md`, and rule 2
   (local noise) of `../macs1423_canucsref/PLAN.md`, all verbatim (same code), except:
4. **Containing-source rule (new; ASSUMPTION):** per test (control or target), the sources whose flux ratio
   F / F_pred enters the class are those whose centroid lies within the search radius **or** whose 3σ ellipse
   (rule 1's ellipse) contains the tested position; the reference's own counterpart stays excluded. Function
   `m0417_arcref.compared_sources`, tested before the run.
5. Gate C2 (rule 7) and the flag (rule 8) on the classes after rules 3–4 and the local-noise rule.
6. Reported for comparison, not rules: gate C2 and flags on MACS1149 without rule 4; and, post hoc, MACS0417 and
   MACS1423 re-run with rule 4 (both fields were seen while designing it, so they are not gates).
