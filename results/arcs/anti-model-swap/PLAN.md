# `arcs`: do the pooled anti rows persist under a second published model? (plan, 2026-11-08)

Frozen before any second-model class of a primary anti row was read (only the primary anti labels and the table
sizes were listed). Script: `scripts/arcs_model_swap.py`; tests: `tests/test_arcs_model_swap.py`. No network, no
new data: the stored `arcs.ecsv` of each field's deciding run (as in `../pooled-tally/`) and of its second run.

**Question.** The pooled anti excess (`../pooled-tally/`, 10 vs 3.22; variant B 6 vs 2.27) has untested ordinary
explanations; one is primary-model error in the predicted stretch direction. Abell S1063's 4 anti rows persist
in three other HFF models (`../abells1063-cats/second_models.json`); this checks the replication fields.

**Hypothesis H_model (ordinary).** The pooled anti rows are anti because of errors of the deciding model's shear
direction. Prediction: under an independent model of the same field the anti rows lose the anti class far more
often than the aligned rows lose the aligned class.

**Inputs (deciding run → second run).** MACS0417 CANUCS v1 → RELICS Lenstool v1; WHL0137 RELICS Lenstool → RELICS
GLAFIC; RXC J0600 RELICS GLAFIC → RELICS Lenstool; MACS0308 RELICS GLAFIC → RELICS Zitrin LTM. MACS1423 (4 anti)
and RXC J0018 (0 anti) have no second run and are reported untested.

**Matching.** A primary row matches a second-run row with the same `label` whose position agrees within 0.1″
(ASSUMPTION); otherwise `unmatched`. The second-run class is its `orientation_class` if that row is a strong-shear
background row there, else `not_strong_background`.

**Rules (functions in the script, unit-tested before the run).**
1. Control C: primary aligned strong-shear background rows of the four fields, matched; retention r = share whose
   second-run class is `aligned`. Gate: r ≥ 0.70 (ASSUMPTION, the `arcs` control threshold); else outcome
   `undetermined` (models disagree in general).
2. Tested anti rows: m = matched primary anti rows; k = those `anti` in the second run.
3. Outcome (pooled over the four fields, combined): gate passes and m ≥ 4 and k ≥ m − 1 → `persists` (H_model
   disfavoured for these rows); gate passes and k ≤ m / 2 → `model-dependent` (H_model supported); otherwise `mixed`;
   m < 4 → `undetermined`.
Secondary (reported, not rules): pooled variant-B count recomputed with the second runs in those four fields
(`lens_consistency` shape-noise expectation, σ_ε 0.26) and the aligned → anti flip rate of the control.

**What would change the plan.** `model-dependent`: the excess is a model-error effect; S1063 then needs its own
explanation (it persists in three models). `persists`: model error is disfavoured; the next test is blends / clump
pairs or photo-z contamination of these rows.
