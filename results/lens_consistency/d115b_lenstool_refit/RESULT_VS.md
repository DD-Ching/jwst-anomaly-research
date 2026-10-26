# D115b stage VS: continuity of Lenstool's source-plane χ² at `best.par` (result, 2026-10-26)

Plan `PLAN_VS.md`, tool `scripts/d115b_lenstool_refit.py --stage vs`, `build_vs`, `series_vs` and `gate_vs` tested
in `tests/test_d115b_lenstool_refit.py`; all committed in 7e95ce42 (17:13Z) before the run. Numbers: `result_vs.json`.
All χ² are `model_prediction` of the pinned `macs0416-canucs` model (D-044) and its single-parameter variants, with
8515 at its `best.par` σ; the image list is `observed`. No D115b number follows from this stage.

## Result
- **Gate VS: `pass`.** With `forme 0` (source-plane χ², amplification averaged per system) Lenstool 8.15.5 prints
  χ² 227.918425 at `best.par` (VS0).
- All 56 (parameter, sign) series were evaluated: 28 O1–O6 parameters × 2 signs, none skipped and no value missing
  (56 × 6 offsets; `runtime_s` 97.8 s). Each of the 56 is `smooth` under the frozen rule: no jump larger than
  |Δχ²| 2 (the frozen floor) at offsets up to 10⁻³. Smaller jumps are not tested: only 3 series (O1, O2 and O4 v_disp,
  + sign) exceed 2 at all by 10⁻³.
- The changes are real and scale with the step. At 3 × 10⁻⁵ of a range they run from ±0.0001 (O5 cut_radius,
  O6 ellipticity) to ±0.058 (O1 and O2 v_disp). Most series grow by about 3.3× per step at the small steps (linear).
- For comparison, in the image plane (RESULT_V3P.md) the same 3 × 10⁻⁵ offset of O2 angle_pos or O3 v_disp changed
  χ² by +12.1. Here those offsets change it by +0.0087 and +0.0452.

## Reading
- In this neighbourhood (offsets up to 10⁻³ of each range, one parameter at a time), the source-plane χ² shows no
  jump above 2, against +12.1 in the image plane. So a local optimiser on it can be validated with a χ² tolerance
  of a few units, which the image-plane χ² did not allow. Jumps below 2 are not excluded.
- Post hoc (not part of the gate): the opposite signs give changes of opposite sign at the smallest steps (for
  example O1 v_disp +0.0579 / −0.0569). So `best.par` is not a stationary point of the source-plane χ². That is
  expected for an image-plane optimum. A source-plane re-fit therefore moves away from `best.par`, so its validation
  must use the source-plane optimum, not `best.par`.
- Untested here: jumps below 2, joint (multi-parameter) moves, offsets above 1 %, and whether a source-plane optimum keeps the
  D115b extra image. The source-plane χ² is not compared to the image-plane χ².

## Limitations (from the review, after the run)
- `result_vs.json` stores Δχ² rounded to 4 decimals and not the raw χ². The 3 × 10⁻⁵ values of O5 cut_radius and
  O6 ellipticity (±0.0001) are one rounding unit.
- `set_param` writes values to 6 decimals, so a 3 × 10⁻⁵ step of a parameter with a range below about 1 is quantised
  by a few per cent. Neither changes the gate (no Δ is near 2), but F-S should write more digits and keep raw χ².
- The objective is the printed `INFO: chi2`. For VS0 it matches `chires.dat`'s `chitot` 227.92 (`chires_vs0.dat.gz`,
  added after the run).
- Post-run code fixes (tested; they do not change this run's outcome): `series_vs` treats a short series as
  `incomplete`; `build_vs` requires the `forme -1` line to be in the `image` section.

## Deviations
- None from the plan. The micromamba binary came from its conda-forge package (`micromamba-2.10.0-0.tar.bz2`,
  sha256 checked against conda-forge's repodata; D-109 addendum 2). Piping micro.mamba.pm into tar was refused
  this run by the session's permission check.
