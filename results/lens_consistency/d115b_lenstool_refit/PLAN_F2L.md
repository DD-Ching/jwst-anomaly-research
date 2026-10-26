# D115b stage F2-L: local optimiser from best.par (plan, frozen 2026-10-26 15:25Z, before any optimiser run)

## Question and why this stage
Stage F2 asks whether a re-fit with the cluster-scale potentials free keeps every catalogued image while member
8515's σ is low enough to remove the D115b extra image (RESULT_D14A.md "Reading"). Lenstool's own sampler cannot
start from `best.par` and wrote 0 samples in 900 s (RESULT_F2.md, gate V2 `undetermined`). This stage replaces it
with a local optimiser started from `best.par`, driving Lenstool `inverse 0` + `chi2 1` (the stage V′ run file, one
χ² in about 0.5 s). Before any D115b reading, the optimiser must show on the fiducial model (8515 at f 1.00) that it
returns to the published optimum from a perturbed start (gate V3). This run does **only V3**; F2-L proper runs in a
later cycle and only if V3 passes.

## Hypothesis and what would change the plan
- H: bounded Nelder–Mead in the 28 normalised O1–O6 parameters, started from a perturbed `best.par`, returns to
  χ² ≤ 344.42 × 1.01 within the evaluation budget.
- Ordinary expectations against it: 28 dimensions are many for Nelder–Mead (it can stall); the χ² surface is cut by
  lost-image and proximity-penalty cliffs.
- Plan change: V3 `pass` → F2-L proper (below) runs next. V3 `fail` or `undetermined` → no D115b reading with this
  optimiser; the next step is a different local method (e.g. Powell, or a gradient-free method with restarts),
  planned separately. V3 `uninformative` (the perturbed start already within the tolerance) → the perturbation is
  too small and is enlarged in a new frozen plan.

## Free parameters
- O1–O6 = `best.par`'s first six `potentiel` blocks = `input.par`'s `limit 1`–`6` (as stage F2-B). Every parameter
  in those limit blocks is free with its own uniform bounds (all are type 1 in `input.par`): 28 parameters.
  Lenstool keyword `ellipticity` is `ellipticite` in `best.par`. When `core_radius` or `cut_radius` is set, its
  `_kpc` companion line is removed so the arcsec value is the only one Lenstool reads.
- The potfile, the other potentials (including 8515) and the source redshifts stay at `best.par`. Deviation from
  the CANUCS free set, stated: this is a local re-fit of the cluster-scale halos only.
- Coordinates: each parameter mapped to u = (p − lo) / (hi − lo) in [0, 1].

## Objective (lost and merged images counted, as RESULT_D14A.md requires)
- One evaluation = one Lenstool `inverse 0` run of the stage V′ file with the parameters substituted.
- **Feasible** iff Lenstool exits 0, `chires.dat` has a `chimul` row, the summed `nwarn` is 0 (a lost image raises
  `nwarn`, RESULT_VPRIME_F.md) and `chipenalty` is 0 (a merged image adds the proximity penalty, RESULT_D14A.md).
- Objective = `chimul` when feasible, else 10⁶ (ASSUMPTION: an infeasible point is never accepted).

## Optimiser (ASSUMPTIONs)
- `scipy.optimize.minimize(method="Nelder-Mead")` with bounds [0, 1] per coordinate, initial simplex steps 0.05 in
  u, `maxfev` 1,500, `xatol` 1e-4, `fatol` 0.01; wall cap 1,500 s (checked between evaluations).
- Lenstool on 4 OpenMP threads (the session's cores); evaluations are sequential.

## Gate V3 (function `gate_v3`, unit-tested before the run)
- Perturbed start: u₀ = u_best + 0.03 × s_k, s_k = +1 for even k and −1 for odd k (k = parameter order: O1…O6, each in
  its `limit` block order), clipped to [0.001, 0.999] (ASSUMPTION: 3 % of each range).
- `uninformative` if the start is feasible with χ² ≤ 344.42 × 1.01 (= 347.86) or the start is infeasible.
- Otherwise `pass` iff the best feasible χ² found ≤ 347.86; `undetermined` if the budget or wall cap ended the run
  above it while χ² was still falling by > 1 per 100 evaluations over the last 300 evaluations; else `fail`.
- Recorded: start χ², best χ², evaluations, wall time, infeasible-evaluation count, max |u_best_found − u_best| and
  the parameter with it.

## F2-L proper (frozen now, runs only after V3 passes; not in this run)
- For f in 0.80, 0.75, 0.70 (stage F grid points where the extra image is absent without a re-fit), 8515 at
  f × 120.600977 km/s, the same optimiser from unperturbed `best.par`.
- At the optimum: Δχ² = χ² − 344.42, feasibility as above, and the D115b flag test (`flag_images`, image mode).
- Rule F2-L: `cheap` iff some f has a feasible optimum with no flag image and Δχ² ≤ 6 (ASSUMPTION, as rule F);
  `costly` iff every f has a feasible optimum with no flag image and Δχ² > 6; else `undetermined`.

## Provenance
All χ² and parameters are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044); the
image list is `observed`; thresholds are ASSUMPTIONs.
