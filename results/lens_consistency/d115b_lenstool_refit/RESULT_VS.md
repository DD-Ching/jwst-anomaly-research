# D115b stage VS: continuity of Lenstool's source-plane χ² at `best.par` (result, 2026-10-26)

Plan `PLAN_VS.md`, tool `scripts/d115b_lenstool_refit.py --stage vs`, `build_vs`, `series_vs` and `gate_vs` tested
in `tests/test_d115b_lenstool_refit.py`; all committed in 7e95ce42 (17:13Z) before the run. Numbers: `result_vs.json`.
All χ² are `model_prediction` of the pinned `macs0416-canucs` model (D-044) and its single-parameter variants, with
8515 at its `best.par` σ; the image list is `observed`. No D115b number follows from this stage.

## Result
- **Gate VS: `pass`.** With `forme 0` (source-plane χ², amplification averaged per system) Lenstool 8.15.5 prints
  χ² 227.918425 at `best.par` (VS0).
- All 56 (parameter, sign) series were evaluated: 28 O1–O6 parameters × 2 signs, none skipped and no value missing
  (56 × 6 offsets; `runtime_s` 97.8 s). Each of the 56 is `smooth`: no jump under the frozen rule.
- The changes are real and scale with the step. At 3 × 10⁻⁵ of a range they run from ±0.0001 (O5 cut_radius,
  O6 ellipticity) to ±0.058 (O1 and O2 v_disp). Most series grow by about 3.3× per step at the small steps (linear).
- For comparison, in the image plane (RESULT_V3P.md) the same 3 × 10⁻⁵ offset of O2 angle_pos or O3 v_disp changed
  χ² by +12.1. Here those offsets change it by +0.0087 and +0.0452.

## Reading
- In this neighbourhood (offsets up to 1 % of each range, one parameter at a time), the source-plane χ² is a
  continuous objective. A local optimiser on it can be validated with a χ² tolerance, which the image-plane χ² did
  not allow.
- Post hoc (not part of the gate): the opposite signs give changes of opposite sign at the smallest steps (for
  example O1 v_disp +0.0579 / −0.0569). So `best.par` is not a stationary point of the source-plane χ². That is
  expected for an image-plane optimum. A source-plane re-fit therefore moves away from `best.par`, so its validation
  must use the source-plane optimum, not `best.par`.
- Untested here: joint (multi-parameter) moves, offsets above 1 %, and whether a source-plane optimum keeps the
  D115b extra image. The source-plane χ² is not compared to the image-plane χ².

## Deviations
- None from the plan. The micromamba binary came from its conda-forge package (`micromamba-2.10.0-0.tar.bz2`,
  sha256 checked against conda-forge's repodata; D-109 addendum 2). Piping micro.mamba.pm into tar was refused
  this run by the session's permission check.
