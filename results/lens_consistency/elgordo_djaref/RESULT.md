# El Gordo reference-catalogue extra-image test

## Stage 1 (2026-10-23; `PLAN_1.md`, `scripts/elgordo_djaref.py`, `summary.json`)
- **Rule 1, gate K: passes** (`stage1_gate_K.elgordo-caminha23.gate_K`). The 56 spectroscopic Caminha+2023 images,
  model and list in the JWST frame (D-034 rule 1): image-plane rms 0.754″ (≤ 1.0″), matched fraction 1.0 within 3″
  (≥ 0.8). Own-constraint validation only (derived; the model was fitted to these positions). No family has fewer
  predicted than catalogued images; 15 families have more (`extra_predicted`: 12 with 2 → 3, families 6, 12 and 17
  with 2 → 5 / 3 → 5 / 3 → 5). These are stage-2 inputs (model_prediction); none was read against the data here.
- **Rule 2, stage-2 reference inputs: fails as written** (`stage1_reference_inputs`). The DJA v7.0
  `elgordo-grizli-v7.0-fix_phot.fits` is reachable (17,170,560 bytes, 6,775 rows, 359 columns; headers only, HTTP
  range) and has the segment shape columns, but no per-band total-flux column. **Post hoc** (read after the failure):
  each band has aperture fluxes (`<band>_flux_aper_0..2`) and an aperture-to-total correction `<band>_tot_corr`
  (18 bands); the detection apertures are 9 / 12.5 / 17.5 / 25 px in diameter (`APER_0..3`).
- **Data:** only the header bytes (≤ 2 MB) were read, in memory; nothing kept.
- **Next:** a frozen stage-2 plan that defines the reference flux from these columns (e.g. `flux_aper_k ×
  tot_corr` of the containing segment, with its own unit test) and checks it on the catalogued controls through gate
  C2 before any extra-image class is read; or another reference catalogue with segment totals.
