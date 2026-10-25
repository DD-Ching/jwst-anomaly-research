# SN Eos 101.2 / 101.1 ratio from HST F110W 2024 − 2012 difference images: result (2026-10-25)

Plan `PLAN.md` (frozen before any pixel was read). Tool `scripts/sneos_hst.py`, tests `tests/test_sneos_hst.py`,
both committed before the run. Numbers from `summary.json`; stamps in `contact_sheet.png`. Fluxes are observed
(difference images); R_model is a model_prediction; the Table S1 ratio is derived from published photometry.

## Result (frozen rules)
- **Q1: "agrees"** with R_model 1.162. R_HST = 1.655, σ_ln 0.152 (scatter term), from 4 of 6 PID 17504 visits.
  |ln R_HST − ln R_model| = 0.354 (not stored) is below max(3σ_ln = 0.456, floor 0.30): the 3σ bound decides, so
  the precision is too low for this test to separate 1.16 from 1.8.
- **Q2: "both".** R_HST is within 3σ of the Table S1 value 1.805 (−0.55σ, not stored) and of our JWST 1.382
  (+1.04σ, not stored). This run cannot tell whether the HST 2024 ratio differs from the JWST one.
- Registration: offsets before the pre-shift 0.002–0.010″, after it 0.002–0.006″ (rule ≤ 0.03″): all 6 registered.
  Gate I: 6 / 6 per visit, all visits. Gate D: `ifaa51010` and `ifaa54010` fail (101.1 S/N 3.1 and 3.0; the two
  306 s visits); the other four pass (101.1 S/N 5.1–6.5).
- Per used visit R: 2.16 ± 0.44, 1.10 ± 0.23, 1.99 ± 0.42, 1.56 ± 0.37 (`ifaa53/52/50/49`, MJD 60380–60384).
  The two unused visits give 1.71 ± 0.64 and 2.48 ± 0.90. Second-template check on `ifaa53010`: R 1.99 with
  `ibt2a9030` vs 2.16 with `ibt2a2020`.
- PSF star: Gaia DR3 6766145724110284544 (the same star as #541), inside every visit.

## Reading
- Our own HST pixels give a ratio above the model on average, as Table S1 does, but with σ_ln 0.15 the result is
  equally consistent with the JWST 1.38. The 101.1 image is near the detection limit in single 300–600 s visits
  (S/N 3–6.5), which sets the precision.
- Contact sheet viewed: both images are point sources in every difference panel, 101.2 the brighter in all six (R > 1 in every visit);
  no other comparable residual in the 6″ stamps; no visible host residual at the SN positions.
- No microlensing, host or PSF test beyond the frozen gates; nothing here is evidence for anything beyond an
  ordinary lens model with microlensing.

## Post hoc and deviations
- **Deviation (implementation bug, fixed before the result was read as final):** the first run shifted the template
  with `scipy.ndimage.shift(cval=nan)`; the NaN edge row, zeroed inside `register_offset`, biased the re-measured
  offset to about 0.1″, so 5 / 6 visits failed registration (first-run output, not stored: 1 visit used, R 2.16). Shown on
  `ifaa52010` in a session check (not stored): after a 0.04 px pre-shift the re-measured offset was −0.71 px with a
  NaN edge and 0.006 px with `mode="nearest"`. The run uses `mode="nearest"`. A synthetic unit test did not
  reproduce the edge bias; it checks sub-pixel registration only.
- **Post hoc (code review):** with the frozen 0.6″ exclusion a random aperture's 0.40–0.60″ annulus can include an
  SN image. Re-run with 0.85″ (`posthoc_minsep085.json`; not the frozen result): 3 visits used (`ifaa53010` 101.1
  S/N 4.68 fails gate D), R 1.478, σ_ln 0.173, Q1 "agrees", Q2 "both". The decisions do not change; R moves by about
  1σ, so the central value depends on this choice.
- **Deviation:** rule 3 names `sneos_ratio.aperture_pair`; the run uses `sneos_hst.aperture_fluxes`, the same
  algorithm with this plan's radii and no error map (σ comes from rule 3's empirical term).
- Template and science pixel scales are equal (0.12825″, `summary.json` "templates"); the run now raises otherwise.
- **Display only:** the contact-sheet stretch removes each panel's median; `summary.json` is identical with and
  without it.
- Rule 3's error has no pipeline ERR (drz files carry none); the empirical random-aperture σ replaces it, as planned.

## Limits
One template epoch per visit (CLASH 2012), single-visit depth, aperture photometry only, one lens model, no
microlensing model. Data: SCI sections by S3 byte range from 8 `_drz.fits` (about 12.7 MB each, not stored); nothing stored.
