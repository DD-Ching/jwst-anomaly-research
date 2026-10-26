# D115b stage V3′: gate V3 with small steps (plan, frozen 2026-10-26 16:25Z, before any V3′ run)

## Question and why this stage
Gate V3 (PLAN_F2L.md) was `uninformative`: the 3 % start was infeasible, and post hoc only 7 / 56 single ±3 %
offsets keep every image (RESULT_F2L.md). V3′ asks the same question with steps that fit the narrow feasible set:
does bounded Nelder–Mead, started from a feasible perturbed `best.par`, return to χ² ≤ 344.42 × 1.01 (= 347.86)?
Only V3′ runs here; F2-L proper (frozen in PLAN_F2L.md, unchanged) runs in a later cycle and only if V3′ passes.

## Hypothesis and what would change the plan
- H: with a perturbation and initial simplex of the same small size, Nelder–Mead returns within tolerance.
- Ordinary expectations against it: 28 dimensions are many for Nelder–Mead; the feasible set is cut by lost-image and
  proximity-penalty cliffs, so many simplex moves are rejected.
- Plan change: `pass` → F2-L proper next. `fail` / `undetermined` → no D115b reading with this optimiser; a different
  method is planned separately. `uninformative` (no ladder step gives a feasible start above tolerance) → this
  hard-constraint local route is closed; the next step is a different handling of the constraint.

## What changes from PLAN_F2L.md (everything else is as frozen there: free parameters, objective, feasibility,
## bounds, `maxfev` 1,500, `xatol` 1e-4, `fatol` 0.01, wall cap 1,500 s, 60 s run timeout, gate tolerance)
- **Start ladder (ASSUMPTION):** fractions 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02 of each range, in
  this order; start u₀ = u_best + frac × (+1, −1, +1, …), clipped to [0.001, 0.999]. The chosen frac is the
  **first** one whose start is feasible with χ² > 347.86. The ladder evaluations are recorded and are not part of
  the optimiser budget.
- **Initial simplex step = the chosen frac** (was 0.05), same sign rule as PLAN_F2L.md.
- Gate V3′ (function `gate_v3p`, unit-tested before the run): `uninformative` if no ladder step qualifies; otherwise
  the PLAN_F2L.md rule on the optimiser run (`pass` iff best feasible χ² ≤ 347.86; `undetermined` if the budget or wall
  cap ended the run while χ² still fell by > 3 over the last 300 evaluations; else `fail`).

## Provenance
All χ² and parameters are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044); the
image list is `observed`; thresholds and the ladder are ASSUMPTIONs.
