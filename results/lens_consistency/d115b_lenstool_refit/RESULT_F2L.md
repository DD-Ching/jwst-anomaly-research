# D115b stage F2-L: gate V3 of a local optimiser from best.par (result, 2026-10-26)

Plan `PLAN_F2L.md` (committed 15:25Z with the tool and gate tests, bf5d7211, before any optimiser run). Tool
`scripts/d115b_lenstool_refit.py --stage f2l-v3`; gate `gate_v3`, the parameter readers and the feasibility rule are
tested in `tests/test_d115b_lenstool_refit.py`. Numbers: `result_f2l_v3.json` (gate), `posthoc_f2l_v3.json` (post
hoc). All χ² are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044); member 8515
stays at its `best.par` σ (f 1.00), so this stage gives no D115b number.

## Result
- **Gate V3: `uninformative`** (`gate_V3`). The perturbed start (each of the 28 O1–O6 parameters moved by 3 % of its
  `input.par` range, alternating sign) is **infeasible** (`start_feasible` false), and so were all 299 distinct Lenstool
  evaluations (`n_infeasible` 299 of `n_evaluations` 299). Nelder–Mead saw a flat objective (10⁶) and stopped in
  104.6 s (`runtime_s`; its message "Optimization terminated successfully." is not a convergence,
  `optimizer_message`).
- The optimiser was therefore not tested; no D115b reading follows.
- Tool check before the run (not gated, not stored): `best.par`'s own values written back through the tool (with
  the `_kpc` companion lines removed) give χ² 344.42, feasible, in 0.29 s.

## Post hoc (after the gate; `posthoc_f2l_v3.json`)
- **Start:** `chimul` 73229.24, of which `chipenalty` 56226.57 (merged images; the `penalty_images` list names
  images in D1, D6e, D8, D14a, D17, D19, D115a, D115b and D34), and family D33 loses images (`warn_families`).
- **Single offsets of ±3 % from `best.par`, one parameter at a time:** 7 of 56 stay feasible (`n_single_feasible`).
  The feasible ones are O4 core_radius −, O5 v_disp +, O5 cut_radius +, O6 angle_pos −, O6 cut_radius +, O6
  ellipticity ±; their χ² runs from 344.57 (O6 ellipticity −) to 422.31 (O4 core_radius −). Every O1, O2 and O3
  offset is infeasible.
- Reading: under the frozen feasibility rule (no lost image, no proximity penalty), the feasible set around
  `best.par` is narrower than 3 % of the `input.par` ranges along most O1–O4 directions. A local optimiser with this
  hard rule needs steps far below 3 % of the range, or a different handling of the constraint.

## What this changes
- The plan listed `uninformative` only for a start already within tolerance (perturbation too small). Here the start
  is infeasible, so the perturbation is too large for the hard constraint. The next step is a new frozen plan
  (V3′): perturbation and simplex steps of order 0.3 % of the range (ASSUMPTION to be frozen), with the same gate,
  before any F2-L run at f < 1.
- D115b stays `inconclusive`.

## Deviations
- The post hoc writer (`--stage f2l-v3-posthoc`, `posthoc_f2l_v3`, `_warn_families`) and its test were added after
  gate V3 was read.
- `/code-review` fixes after the first run, then both stages re-run (Lenstool is deterministic here): the start is
  no longer evaluated twice (first run: 300 evaluations, 125.5 s, not stored), each Lenstool run has a 60 s timeout
  counted as infeasible (`run_timeout_s`, deviation from the plan's between-evaluation wall check), the best-point
  deviation is reported only for a feasible point (the first run stored the start's 0.03), post hoc Lenstool failures
  are recorded instead of aborting. Gate outcome, start and all 56 single-offset rows are unchanged.
- The plan's "Plan change" bullet did not cover an infeasible start; the follow-up above is chosen post hoc.
