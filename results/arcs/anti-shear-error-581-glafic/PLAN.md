# `arcs`: RXC J0600 581 under the RELICS GLAFIC v3 realisations (plan, 2026-11-08)

Frozen before any realisation deflection or κ value was read (only the RELICS `glafic/v3/range/` directory
listing, the file sizes and the FITS headers of map000, map098 and the best fit were looked at, plus the stored
GLAFIC run table). Script: `scripts/arcs_shear_error_581g.py`; tests: `tests/test_arcs_shear_error_581g.py` (unit
conversion, resampling and gate functions tested offline before the run). Follow-up of D-024 addendum 13, which
tested 581 only under the Lenstool v1 realisations (`robust`) and named the GLAFIC model, in which 581 was first
classed `anti`, as having no realisations. The listing shows that it has them.

**Question.** Is the `anti` class of RXC J0600 581 (deciding run GLAFIC v3, offset_min 63.6° against the 60° cut,
stored `rxc0600-relics/rxc0600-relics-glafic/arcs.ecsv`) robust to the GLAFIC model's own statistical error?

**Hypothesis H_err (ordinary).** The row sits 3.6° past the anti cut in the best fit; GLAFIC's statistical error in
the predicted stretch direction makes it anti. Prediction: the row often falls below 60° across realisations.

**Inputs (observed products of a model, labelled `model_prediction`).** RELICS GLAFIC v3 `range/` maps
`glafic-mapNNN_v3_{x,y}-pixels-deflect.fits`, NNN = 000–098 (99 maps; 600 × 600 float32, 0.4″ px, 1.4 MB each,
DOI 10.17909/T9SP45), and `_kappa.fits` for maps 000 / 049 / 098; the best fit `glafic_v3_{x,y}-{pixels,arcsec}-
deflect.fits` (2400 × 2400, 0.1″ px). Each file is streamed into memory, hashed (sha256) and dropped; only the
three κ gate maps touch the disk, briefly, for `map_check`.

**Deviations from addendum 13's frozen constants (labelled, decided before the run).** 99 realisations instead of
100 and κ gate maps 000 / 049 / 098 instead of 000 / 049 / 099, because only map000–map098 are published.

**Unit and resolution route (gates; failing any → `undetermined`).**
- V1 unit convention on the best fit: max |pixel × own header pixel scale − arcsec| ≤ 1 × 10⁻³″ in both axes
  over the full map (ASSUMPTION). Realisation deflection in arcsec = pixel deflection × its own header pixel scale
  (0.4″), D_LS/D_S = 1 (RELICS readme, ASSUMPTION). The κ gate (rule 1) checks this convention on the
  realisations themselves: a wrong scale changes κ by a constant factor.
- V2 resolution: the best-fit arcsec maps resampled (bilinear) onto the 0.4″ realisation grid of map000, run
  through `orientation_table` on the 12 stored strong-shear background rows, reproduce every stored
  `orientation_class`, with max |Δ offset_min| ≤ 2° (ASSUMPTION; looser than addendum 14's 0.5° because the grid
  is 4× coarser).
- V3 centring: at the 12 rows, the median over rows of |median-over-realisations α − best-fit α| / |best-fit α| ≤
  0.10 (ASSUMPTION): the realisations scatter about the best fit.

**Rules (addendum 13's frozen functions, unchanged, imported from `scripts/arcs_shear_error.py`).**
1. Map gate on maps 000 / 049 / 098: κ from our realisation model against the published realisation κ,
   `lens_consistency.map_check`, median relative difference ≤ 0.05.
2. Control gate: median over realisations of the share of the 11 `aligned` control rows still `aligned` ≥ 0.70.
3. f = share of realisations in which 581 is `anti`: `robust` f ≥ 0.68; `model_sensitive` f ≤ 0.32; else `partly`;
   any failed gate (V1, V2, V3, rule 1, rule 2) → `undetermined`.
4. Outcome: `range_label` mapped by addendum 13's `outcome` over this one row.
Secondary (reported, not rules): 16 / 50 / 84 % of offset_min and g_min; share of realisations in which 581 is
`aligned`; 581 under the resampled best fit (V2 model).

**What would change the plan.** `model_sensitive`: 581 leaves the open set and only MACS0308 2157 (membership
`undetermined`, addendum 17) remains. `robust`: 581 is anti in both published RXC J0600 models and under both
models' statistical errors; the remaining ordinary routes are spectroscopy (membership) or an independent model.
Neither outcome is a lensing anomaly claim; the row is not a candidate.
