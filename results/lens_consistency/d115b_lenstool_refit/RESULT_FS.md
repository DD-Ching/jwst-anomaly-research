# D115b stage F-S: local re-fit on Lenstool's source-plane χ² (result, 2026-10-26)

Plan `PLAN_FS.md`, tool `scripts/d115b_lenstool_refit.py --stage fs`; `gate_vfs`, `fs_clean` and `rule_fs` were
tested in `tests/test_d115b_lenstool_refit.py`. All of it was committed in 34ab3934 before the run. Numbers are in
`result_fs.json`. All χ² are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044).
The image list is `observed`. No D115b number follows from this stage.

## Result
- **Gate VFS: `undetermined`.** Runs A and B (f 1.00, from `best.par` and from the 3 % perturbed start) both used
  the full 2,000-evaluation budget (`stopped_by_budget`; 660 s each, under the 900 s wall cap). They ended at
  source-plane χ² 174.900 (A) and 193.797 (B). The difference is above the tolerance of 2, so they do not agree. The
  frozen gate reads a budget-ended disagreement as `undetermined`, so rule FS is not applied.
- Starts: A 227.918 (= VS0), B 11727.27; C (f 0.80) 220.94 → 172.37; D (f 0.75) 220.68 → 172.58. Every evaluation
  returned a value (`n_no_value` 0).

## Post hoc (not gated; the plan reports end-point image-plane values without a rule)
- **None of the four end points is clean in the image plane.** Each one has lost images (`nwarn` 6–9) and a
  proximity penalty (`chipenalty` 7496.88–14993.75, so image-plane `chimul` is 7831–15335 against 344.42 at
  `best.par`).
  - A: penalty images D1d.1, D14a.2, D77.2, D33.5.
  - C and D: D14a.1, D33.5, D34.2.
- So lowering the source-plane χ² by 53.02 (A: 227.92 → 174.90; difference not stored) moved the model to configurations that do not
  reproduce the catalogued images. Run A's end point lies far from `best.par` (max |Δu| 0.4152, O6 cut_radius).
  Read this way, the unconstrained source-plane χ² is not a usable proxy for this model's image-plane fit at these
  distances. VS tested only its continuity at offsets ≤ 1 %.
- The D115b flag images are absent at the A, C and D end points and present (D115b.1) at B's. With images lost
  elsewhere, this says nothing about D115b.

## Reading
- With the frozen budget, the local source-plane re-fit is not validated (VFS `undetermined`). Its end points also
  break the image-plane solution (post hoc). The source-plane route, as frozen, gives no D115b reading. D115b stays
  `inconclusive`.
- Not tested: a larger budget; an objective that keeps image-plane feasibility (for example the source-plane χ²
  with a lost-image or penalty term, or the image-plane χ² with D52b's matching fixed, PLAN_VS.md's `fail` branch);
  f 0.70 (deferred in the plan).

## Deviations
None from the plan. Notes: the four runs were concurrent with one Lenstool thread each, as planned. The image-plane
reading of runs A and B is reported for reference, as planned.
