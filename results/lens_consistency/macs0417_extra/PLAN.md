# MACS0417 / MACS1423: are the images CANUCS v1 predicts beyond the catalogued ones there? (plan, frozen)

Frozen 2026-10-21 before any forced photometry. Follows `results/lens_consistency/macs0417_images/RESULT.md`
(`extra_predicted`: MACS0417 families 8, 9, 10, 11, 17; MACS1423 system 2). Script: `scripts/m0417_extra.py`.

## Hypothesis and ordinary explanations
- Question: for each catalogued family, does the CANUCS v1 map model predict an image (|μ| ≥ 0.5) where the
  F200W data show nothing, although the model predicts it ≥ 10σ? Such an image is a `lensing_residual_input`
  (input for vetting, not a candidate).
- Ordinary explanations, in order: the published family z is wrong (fitted z for non-spectroscopic families);
  the extra image sits on a cluster member or in the BCG/ICL light; the reference image's flux is not a clean
  compact flux; the model is poorly constrained there (CANUCS v1 image catalogue not public); a chip gap / edge.
- What would change the plan: a flag that survives the halo filter and both models' reach would go to
  `/vet-candidate`; if gate C fails, the forced-photometry route is not valid on these fields and no reading is made.

## Inputs (observed / model_prediction)
- Images: `../macs0417_images/jauzac19_images.ecsv` (MACS0417), `../macs0417_images/kaleidoscope_m1423_images.ecsv`
  (MACS1423). Family z (ASSUMPTION): `z_used` of `../macs0417_images/imageplane_<field>_<field>-canucs-maps.ecsv`
  (z_spec, else the CANUCS-fitted z of the earlier run).
- Model: `macs0417-canucs-maps`, `macs1423-canucs-maps` (pinned, unchanged).
- Catalogue and image: JWST 1208 F200W pipeline `_cat.ecsv` (sha256 in SOURCES.md) and `_i2d.fits` by S3 byte
  range: o006_t003 (MACS0417), o010_t005 (MACS1423).

## Method (unchanged code: `lens_consistency.predict_counter_images`, `forced_check`, D-029 / D-034 rules)
- Every family with ≥ 2 traced images is back-traced, its mean source re-imaged; match 1.5″, footprint 5″,
  |μ| < 0.5 = `demagnified` (not tested).
- Forced photometry of every predicted image not matched to a catalogued one, as D-029 (reference = least
  magnified compact catalogued image, search radius from the family's residuals, ≤ 3″).

## Frozen rules (computable; `m0417_extra.flag`, `m0417_extra.gate_c`)
- **Gate C (positive control, per field):** the catalogued images themselves, run through the same
  `forced_check` as targets; of those with `forced_class` in {recovered, absent, confused, ambiguous} (i.e. a
  testable prediction), the fraction `recovered` or `confused` must be ≥ 0.8. Fails → no reading for that field.
- **Flag (per predicted image):** `image_class` not in {observed, demagnified, outside} **and**
  `forced_class` == "absent" → `lensing_residual_input`. Every other image: `no_flag`.
- No flag is a candidate before `/vet-candidate`. Thresholds are the D-029 ones (ASSUMPTIONs).
