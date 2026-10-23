# Abell 2744: UNCOVER-reference extra-image test, stage 1b (gate K with the WLS source) (plan, frozen)

Frozen 2026-10-23 (commit time of this file), before the stage-1b run. **Deviation, chosen after a failure:** stage 1
(`PLAN.md`) failed gate K (rms 4.90″, `summary.json` `stage1_gate_K`) with the mean back-traced family source of
`lensmodel.imageplane_residuals`. Stage 1's RESULT named the estimator as the likely ordinary cause; this stage swaps
it for the magnification-weighted least-squares (WLS) source of D-041 (`a2744_uncover.imageplane`, unchanged code),
as BUFFALO K2 (D-034 addendum h) added a gate after a failure. A pass is labelled as such everywhere it is quoted.

## Hypothesis, ordinary explanations, what changes the plan
- Question (method): with the WLS source, do the UNCOVER v2.0 maps reproduce their 172 spectroscopic images under
  the unchanged gate K? A pass lets a later cycle freeze stage 2 (UNCOVER DR3 reference, extra-image stages,
  gate C2 with the frozen shared-pair rule) using the WLS source throughout. A fail stops the Abell 2744
  extra-image route with this model; the remaining option is Bergamini+23 (stage-1 RESULT option b).
- Ordinary explanations of a fail: map resolution (100 mas), the 0.25″ grid search near critical curves,
  knot-level ids, the spec images' own z. Exotic readings are not considered at this stage.

## Known before freezing (disclosed)
- D-041 (`../abell2744-uncover2/summary.json`): WLS source, rms 0.84″ on 181 own-list images (1.5″ match, every
  system at its table z, the grid of `a2744_uncover.py`). Stage 1 post hoc: median offset 0.66″, three single-image
  families (34.1, 80.1, 74). So a pass is expected and the outcome is not blind; it validates the model on its own
  constraints only.

## Frozen rules
1. **Gate K-WLS** (function `a2744_uncoverref.gate_k_wls`, unit-tested before the run): `a2744_uncover.imageplane`
   on the 172 spectroscopic images (same image list `uncover2_images.ecsv`, same `m0417_images.Solver` 0.25″ grid,
   frame offset 0), then the unchanged `m0417_images.gate_k`: rms ≤ 1.0″ **and** ≥ 80 % within 3.0″ (ASSUMPTIONs,
   combined over all 172 images). Fails → stop.
2. Reported, not rules: per-family predicted-vs-catalogued counts (`extra_predicted` / `fewer_predicted`), images not
   reproduced within 3″, `shared_match` count. No forced photometry, no flag in this stage.
3. Stage 2 (only after a pass) uses the WLS source for predicted images, and needs its own frozen plan.
