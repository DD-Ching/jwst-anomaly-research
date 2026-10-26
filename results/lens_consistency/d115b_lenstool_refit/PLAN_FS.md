# D115b stage F-S: local re-fit on Lenstool's source-plane χ² (plan, frozen 2026-10-26, before any F-S run)

## Question and why this stage
Stage F asks whether lowering member 8515's σ, which removes the D115b extra image, can be absorbed by re-fitting
the cluster-scale halos O1–O6. The image-plane χ² is discontinuous at `best.par` (RESULT_V3P.md), so every local
re-fit on it was uninformative. Stage VS found the source-plane χ² (`forme 0`) without jumps above 2 at `best.par`
(RESULT_VS.md), and noted that `best.par` is not its stationary point. F-S re-fits O1–O6 on that objective and reads
the image plane (lost images, proximity penalty, D115b flag) only at each end point.

## Hypothesis, ordinary expectations and what would change the plan
- H (ordinary): a halo re-fit absorbs the σ(8515) change; some tested f has a clean end point (no lost image, no
  proximity penalty, no D115b flag image) at Δχ²_S ≤ 6. D115b's extra image is then a model-degeneracy effect.
- Against it: the D14a proximity penalty and lost images of stage F may persist after a halo re-fit, or the re-fit
  may restore the extra image.
- Plan change: `cheap` → D115b is explained by a model degeneracy (σ of one member), vetting closes `explained`
  after review. `costly` → the extra image is a stable model prediction on the tested grid; next is an independent
  second-model check before any further reading. `undetermined` → no reading; the cause is listed.

## Runs (tool `scripts/d115b_lenstool_refit.py --stage fs`)
- Objective file: the stage VS run file (`build_vs`: stage V′ file with `forme 0`), O1–O6 substituted (28
  parameters, `free_params`, u in [0, 1] as PLAN_F2L.md), 8515 at f × 120.600977 km/s. Parameters are written with 9
  decimals (RESULT_VS.md review). Objective = the last printed `INFO: chi2`; a run that fails, exceeds 60 s or prints
  none scores 10⁶ (ASSUMPTION).
- Optimiser: bounded Nelder–Mead as PLAN_F2L.md (initial simplex 0.05 in u, `xatol` 1e-4, `fatol` 0.01), `maxfev`
  2,000, wall cap 900 s per run (ASSUMPTIONs; cloud slot).
- Four runs in parallel, one Lenstool thread each (4 cores):
  - **A** f 1.00 from `best.par` (reference optimum χ²_A);
  - **B** f 1.00 from the PLAN_F2L.md perturbed start (u_best ± 0.03 alternating, clipped);
  - **C** f 0.80 and **D** f 0.75 from `best.par`.
- f 0.70 (PLAN_F2L.md grid) is deferred to a later run (cloud slot); rule FS below reads the tested grid only.
- End-point image-plane check (C and D, and A for reference): stage V′ file (`forme -1`) at the end point's
  parameters and f, `image 1` mode; recorded: `chimul`, summed `nwarn`, `chipenalty`, the D115b flag images
  (`flag_images`). Reported χ² differences in the image plane are not gated.

## Gate VFS (function `gate_vfs`, unit-tested before the run)
- tol = max(2, 0.01 χ²_A) (ASSUMPTION: VS's jump floor, PLAN_F2L.md's 1 %).
- `uninformative` iff B's start has no value or its χ² ≤ χ²_A + tol (start already at the optimum).
- `pass` iff |χ²_B − χ²_A| ≤ tol (end values; both must have values).
- `undetermined` iff not pass and either A or B was ended by budget or wall cap; `fail` otherwise.

## Rule FS (function `rule_fs`, unit-tested before the run; applied only if VFS passes)
- An f is **clean** iff its end point's image-plane check has `nwarn` 0, `chipenalty` 0 and no flag image.
- Δχ²_S(f) = χ²_S(f) − χ²_A (ASSUMPTION: the image-plane threshold 6 of rule F carried to the source plane).
- `cheap` iff some tested f is clean with Δχ²_S ≤ 6; `costly` iff every tested f is clean with Δχ²_S > 6; else
  `undetermined`.

## Not part of this stage
No change to D115b's vetting status by this run alone; no BayeSys; no other potentials freed.

## Provenance
All χ² and parameters are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044);
the image list is `observed`; tolerances, budgets and thresholds are ASSUMPTIONs.
