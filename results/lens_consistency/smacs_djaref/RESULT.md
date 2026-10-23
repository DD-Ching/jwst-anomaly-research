# SMACS 0723 reference-catalogue extra-image test

## Stage 1 (2026-10-23): gate K passes; the DJA v7.4 reference carries every stage-2b input

Frozen plan: `PLAN_1.md` (committed before the run). Script `scripts/smacs_djaref.py`; numbers from `summary.json`.

- **Rule 0 (ASSUMPTION):** systems without `z_m_limit` in `input.par`: 1, 2, 3, 5, 7 (16 images), as `z_spec`.
- **Rule 1, gate K** (`smacs0723-iclv2`, mean back-traced source, 0.25″ grid): **passes**, rms 0.321″, 16 / 16
  within 3″ (thresholds ≤ 1.0″ and ≥ 80 %). Own-constraint validation only (the model was fitted to these
  positions); the fixed-z set is a proxy for spectroscopic, not a cross-checked spec-z list.
- **Rule 2, reference inputs:** the DJA v7.4 `smacs0723-grizli-v7.4-fix_phot.fits` (38,626,560 bytes, 13,554 rows,
  407 columns; headers only, HTTP range) has every column and card El Gordo stage 2b read; `APER_1` 12.5 px.
  **Passes.**
- **Stage-2 inputs, not read** (`extra_predicted`, catalogued / predicted images per family, fixed-z systems at their
  fixed z, the rest at `m0417_images.fitted_z`): 3 (4 / 5), 8 (2 / 3), 9 (2 / 3), 11 (2 / 3), 13 (4 / 5), 16 (2 / 3),
  17 (3 / 5), 25 (3 / 5), 26 (2 / 7). No family has fewer predicted than catalogued images. System 8's fitted z sits
  at the grid's upper edge (10.0), as its free z in `best.par` (11.76, D-024) is unconstrained.
- **Implementation note (not a rule change):** `fitted_z` needs `z_lens`; `lc.MODELS` has none for this model, so
  the script takes `best.par`'s 0.39. It affects only the recorded fitted z of free systems, not gate K.

**Next:** stage 2 under a frozen plan: El Gordo stage 2b's rules on this reference (`flux_auto × tot_corr` unit
check, grid check, frame tables with `m0417_extra.run_field` on a GO 2736 F200W tile, gate C2 with shared pairs,
local noise).
