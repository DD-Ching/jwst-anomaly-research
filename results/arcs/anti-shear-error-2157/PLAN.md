# `arcs`: MACS0308 2157 under the Zitrin LTM-Gauss realisations (plan, 2026-11-08)

Frozen before any realisation deflection or κ value was read (only the RELICS `range/` directory listing, file
sizes and the FITS headers of map000 and the best fit were looked at). Script: `scripts/arcs_shear_error_2157.py`;
tests: `tests/test_arcs_shear_error_2157.py` (strip reader and conversion tested offline before the run). Follow-up
of D-024 addendum 13, which left this row `untested`.

**Question.** Is the `anti` class of MACS0308 2157 (deciding run GLAFIC v2; also `anti` in the Zitrin best fit,
addendum 8) robust to the Zitrin LTM-Gauss model's own statistical error (100 RELICS `range/` realisations)?

**Hypothesis H_err (ordinary).** As in addendum 13: model error in the predicted stretch direction or |g| makes the
row anti. Prediction: the row often loses the `anti` class across realisations.

**Inputs.** RELICS Zitrin LTM-Gauss v1 `range/` maps `zitrin-ltm-gauss-mapNNN_v1_{x,y}-pixels-deflect.fits`
(3000 × 3000 float32, 0.06″ px, 36 MB each, DOI 10.17909/T9SP45) and `_kappa.fits` for maps 000 / 049 / 099; the best
fit `zitrin-ltm-gauss_v1_{x,y}-{pixels,arcsec}-deflect.fits` for validation. **Byte-range strip read** (no file is
stored): rows j covering the tested row and the control rows plus a 16 px margin, full width, then sliced in
columns; the WCS is the header WCS sliced to the strip. Deflection in arcsec = pixel deflection × the header
pixel scale (0.06″), D_LS/D_S = 1 (RELICS readme v3, ASSUMPTION as for the other RELICS maps). Rows: the deciding
GLAFIC run's stored `arcs.ecsv` (observed axis, photo-z range); the tested row 2157 and the run's `aligned`
strong-shear background control rows.

**Validation (route gates; failing either → the whole test is `undetermined`).**
- V1 unit conversion: on the best-fit strip, max |pixel × scale − arcsec| ≤ 1 × 10⁻³″ in both axes (ASSUMPTION).
- V2 known case: a model built from best-fit arcsec strips, run through `orientation_table` on the stored Zitrin
  run's strong-shear background rows, reproduces every stored `orientation_class` and each `offset_min` within
  0.5° (ASSUMPTION; the stored run used the full maps).

**Rules (addendum 13's frozen functions, unchanged, imported from `scripts/arcs_shear_error.py`).**
1. Map gate on maps 000 / 049 / 099: κ from our strip model against the published realisation κ strip,
   `lens_consistency.map_check`, median relative difference ≤ 0.05.
2. Control gate: median over realisations of the share of control rows still `aligned` ≥ 0.70.
3. f = share of realisations in which 2157 is `anti`: `robust` f ≥ 0.68; `model_sensitive` f ≤ 0.32; else `partly`;
   any failed gate (V1, V2, rule 1, rule 2) → `undetermined`.
4. Outcome: `range_label` result mapped by addendum 13's `outcome` over this one row.
Secondary (reported, not rules): 16 / 50 / 84 % of offset_min and g_min; the share of realisations in which 2157 is
`aligned`.

**What would change the plan.** `model_sensitive`: 2157 leaves the open set and only RXC J0600 581 remains.
`robust`: the remaining ordinary routes for 581 and 2157 are a spectrum (membership) or a further independent model.
