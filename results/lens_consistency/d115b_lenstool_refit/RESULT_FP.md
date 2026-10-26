# D115b stage F-P: local re-fit on the source-plane χ² plus an image-plane feasibility penalty (result, 2026-10-26)

Plan `PLAN_FP.md` (frozen 19:11Z, 6680a728), tool `scripts/d115b_lenstool_refit.py --stage fp` with `fp_lost`,
`fp_objective`, `gate_vfp` and `rule_fp` tested in `tests/test_d115b_lenstool_refit.py` (f92d5f37), all committed
before the run. Numbers: `result_fp.json` (gate) and `posthoc_fp.json` (post hoc); end-point `chires.dat` files
`chires_fp_<run>_end.dat.gz`. All χ² and parameters are `model_prediction` of re-fitted variants of the pinned
`macs0416-canucs` model (D-044); the image list is `observed`; P = 100, tolerances and budgets are ASSUMPTIONs.

## Result
- **Gate VFP: `undetermined`.** Run B (f 1.00, start 0.3 % off `best.par`) began with L = 3 (g 671.33) and ended
  at L = 2 (g 466.26, χ²_S 266.26) when its 1,000-evaluation budget ran out; it never reached L = 0
  (`first_eval_with_lost_0` null). Run A (f 1.00 from `best.par`) stayed feasible and ended at g = χ²_S 225.85
  (start 227.92), also budget-ended. Not pass, and a run was budget-ended, so `undetermined`. Rule FP is not applied.
- Runs C (f 0.80) and D (f 0.75) started at L = 4 (χ²_S 220.94 / 220.68) and ended at L = 2 (χ²_S 218.29 / 218.25).
  C was budget-ended; D stopped by Nelder–Mead's own tolerances after 703 evaluations. Every evaluation returned a
  value (`n_no_value` 0). Largest end-point move from `best.par`: 0.0014 (A), 0.0056 (B), 0.0025 (C), 0.0020 (D) in u.
- End-point image check: A has both D115b flag images (expected at f 1.00), `nwarn` 0, no penalty, `chimul` 356.97.
  C and D have no flag image, `nwarn` 2, no penalty, `chimul` 367.23 / 371.43 (lost images count 0 there, D-109
  addendum 2). B has penalties on D1d.1 and D34.2 (`chipenalty` 7496.88).

## Post hoc (after the gate was read; `posthoc_fp.json`)
- At the C and D end points both warnings sit in family D115a (`warn_families` {D115a: 2}) and the D14a proximity
  penalty of the starts (RESULT_D14A.md) is gone. So within about 0.25 % of each parameter range of `best.par`, the
  penalised halo re-fit removes the D14a merge at f 0.80 / 0.75 but keeps two D115a images lost. D115a is the
  companion knot whose catalogued images stage F and the vetting record's σ × 0.7 `validate` run already lost.
- Run D ended there by its own convergence test, so under these settings L = 2 is a local stopping point, not a
  measured minimum over the halo parameters.

## Reading
The local re-fit is not validated (VFP `undetermined`: the penalised objective did not bring a 0.3 % start back to
a feasible point in 1,000 evaluations). By the plan, a VFP `undetermined` with B never reaching L = 0 closes the
local derivative-free route for this model in cloud slots. No D115b number follows; D115b stays `inconclusive`.
Post hoc, the lost D115a images, not D14a, are what the local re-fit could not repair at f 0.80 / 0.75.

## Deviations
- The plan header's time stamp was first written as 19:35Z (typo) and corrected to the commit time 19:11Z in the
  tool commit, before the run.
- The post hoc writer was an inline script over the stored end-point `chires.dat` files (`posthoc_fp.json`), using
  the tool's `read_chires`, `read_penalty` and `_warn_families`.
