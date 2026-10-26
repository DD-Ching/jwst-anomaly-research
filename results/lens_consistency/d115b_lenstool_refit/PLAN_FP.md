# D115b stage F-P: local re-fit on the source-plane χ² plus an image-plane feasibility penalty (plan, frozen 2026-10-26 19:35Z, before any F-P run)

## Question and why this stage
Stage F asks whether lowering member 8515's σ, which removes the D115b extra image, can be absorbed by re-fitting
the cluster-scale halos O1–O6. Stage F-S (RESULT_FS.md, D-109 addendum 8) minimised the unconstrained source-plane
χ² (`forme 0`) and every end point lost images or picked up proximity penalties in the image plane. Its handoff asks
for an objective that keeps image-plane feasibility. F-P adds a penalty per lost or merged image to the
source-plane χ².

## Hypothesis, ordinary expectations and what would change the plan
- H (ordinary): a halo re-fit absorbs the σ(8515) change; at some tested f the end point has no lost image, no
  proximity penalty and no D115b flag image, at Δχ²_S ≤ 6. D115b's extra image is then a model-degeneracy effect.
- Against it: the D14a merge (from f 0.90, RESULT_D14A.md) and the lost D115a images may persist after a constrained
  halo re-fit, or the re-fit may restore the extra image.
- Plan change: `cheap` → D115b is explained by a model degeneracy (σ of one member); vetting closes `explained` after
  review. `persists` or `costly` → the extra image is a stable prediction on the tested grid; next is an independent
  second-model check. `undetermined` → no reading; the cause is listed. A VFP `fail` or `undetermined` with B never
  reaching L = 0 closes the local derivative-free route for this model in cloud slots (record it in D-109).

## Objective (tool `scripts/d115b_lenstool_refit.py --stage fp`)
- Parameters: O1–O6 (28, `free_params`, u in [0, 1] as PLAN_F2L.md), 8515 at f × 120.600977 km/s, 9 decimals.
- Per evaluation, two Lenstool runs: (1) the stage VS file (`forme 0`) → χ²_S = last `INFO: chi2`; (2) the stage V′
  file (`forme -1`) → `nwarn` (summed warnings in `chires.dat`) and `chipenalty`.
- L = nwarn + round(chipenalty / 3748.4382) (one penalty unit = (30″ / 0.49″)², RESULT_D14A.md).
- g = χ²_S + P · L with P = 100 (ASSUMPTION: F-S traded about 12 lost / merged images for −53 in χ²_S, under 5 per
  image, so 100 per image forbids that trade). A failed, stalled (> 60 s) or value-less run gives g = 10⁶.
- Optimiser: bounded Nelder–Mead as PLAN_FS.md except initial simplex 0.003 in u (ASSUMPTION: V3′ found 33 / 56
  single ±0.1 % offsets feasible against 7 / 56 at ±3 % (RESULT_F2L.md), so the simplex starts near the feasible
  scale), `xatol` 1e-4, `fatol` 0.01, `maxfev` 1,000, wall cap 1,080 s per run (ASSUMPTIONs; cloud slot).

## Runs (four in parallel, one Lenstool thread each, 4 cores)
- **A** f 1.00 from `best.par` (reference g_A, L_A).
- **B** f 1.00 from u_best + 0.003 × (+1, −1, +1, …), clipped to [0.001, 0.999] (ASSUMPTION: V3′ found every
  all-parameter start ≥ 0.1 % infeasible, so B tests whether the objective recovers feasibility).
- **C** f 0.80 and **D** f 0.75 from `best.par` (both start with the D14a merge, RESULT_D14A.md).
- End-point image check (all runs): stage V′ file in `image 1` mode at the end point; recorded: `chimul`, `nwarn`,
  `chipenalty`, the D115b flag images (`flag_images`). Image-plane χ² differences are reported, not gated.

## Gate VFP (function `gate_vfp`, unit-tested before the run)
- tol = max(2, 0.01 g_A,end) (ASSUMPTION, as gate VFS).
- `undetermined` iff A has no end value.
- `uninformative` iff B's start has no value or g_B,start ≤ g_A,end + tol.
- `pass` iff L_A,end = 0, L_B,end = 0 and |g_B,end − g_A,end| ≤ tol.
- otherwise `undetermined` iff A or B was ended by budget or wall cap; `fail` otherwise.

## Rule FP (function `rule_fp`, unit-tested before the run; applied only if VFP passes)
For each tested f (C, D): `clean` iff L_end = 0 and no flag image at the end point; `flagged` iff L_end = 0 and a
flag image present. Δχ²_S(f) = χ²_S,end(f) − χ²_S,A,end (ASSUMPTION: threshold 6 of rule F carried to the source
plane, as rule FS).
- `cheap` iff some f is clean with Δχ²_S ≤ 6;
- else `persists` iff every f is flagged;
- else `costly` iff every f is clean with Δχ²_S > 6;
- else `undetermined`.

## Not part of this stage
No change to D115b's vetting status by this run alone; no BayeSys; no other potentials freed; f 0.70 deferred.

## Provenance
All χ² and parameters are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044);
the image list is `observed`; P, tolerances, budgets and thresholds are ASSUMPTIONs.
