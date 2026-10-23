# El Gordo reference-catalogue extra-image test

## Stage 2b (2026-10-23; `PLAN_2b.md`, `--stage 2b`, `summary.json` `stage2b`)
- Rule 2a (with `tot_corr`): 23.9 mag → µJy (outcome known before freezing; disclosed). Rule 2b passes as in stage 2.
- **Gate C2 (shared pairs): passes** (`stage2b.frozen.gate_C2_shared`): 14 / 15 testable controls recovered or
  confused, arc-referenced 7 / 8. **No flag.** DJA rows kept: 6,654.
- Targets (`target_classes`, 16 rows): recovered 1, undetectable 7, no_reference 4, undetectable_local 3,
  inconsistent_reference 1; so 1 of 16 predicted extra images was testable at all (low testability, not a null).
- Reported, not rules: without the local-noise rule C2 passes, no flag; without the containing-source rule C2 fails
  (`without_contain`), no flag; plain `gate_c2` passes.
- **Contact sheet looked at** (`sheet_2b.png`: F200W, 6″ stamps, red circle = `search_arcsec`; 15 controls whose
  `class_local` is recovered / confused / absent and the 4 targets that are recovered or `undetectable_local`; made
  by a session script, not committed, by S3 range reads). Controls: each recovered control has its arc or knot in
  the circle; 10b (`confused`) sits on a diffraction spike; 12a (`absent`) shows only faint knots. No control
  stamp contradicts its class.
- **Target 23 (`missing`, the only testable target), post hoc:** it is the model's image of catalogued 23c, not an
  additional image: `frame/summary.json` lists 23c as `unpredicted`, and `readings_elgordo.ecsv` puts the model's
  23c 1.538″ from the catalogued one. `sheet_2b_t23.png` (10″ stamp, median / MAD stretch) shows no compact
  m ≈ 23 source (`m_pred` 22.99) at the predicted position; the frozen `recovered` comes from 3 sources up to
  2.25″ away (`max_ratio` 0.58). Reading: a model position error of about 1.5″ for 23c, as for Abell 2744 78.3c
  (D-034 addendum n); not a candidate, nothing to vet.
- **Outcome:** El Gordo by this route gives no flag and tests 1 of 16 predicted extra images; it is not a null over
  the other 15 (low testability: 7 undetectable, 4 no_reference, 3 undetectable_local, 1 inconsistent_reference).

## Stage 2 (2026-10-23; `PLAN_2.md`, `scripts/elgordo_djaref.py --stage 2`, `summary.json` `stage2`)
- **Rule 2b, grid: passes** (`rule_2b`): rotation 0.0002°, determinant negative (east left), pixel scale
  0.0400″ on 6,775 rows.
- **Rule 2a, flux unit: fails as written** (`rule_2a`): median `mag_auto + 2.5 log10(flux_auto)` = 23.69 mag over
  6,751 rows, outside 23.9 ± 0.05 and 28.9 ± 0.05. The run stopped there: no F200W flux row, mosaic pixel, frame
  table or prediction was read; gate C2 not run.
- **Post hoc** (`posthoc_rule_2a`, read after the failure): `mag_auto + 2.5 log10(flux_auto × tot_corr)` = 23.9 mag
  (5–95 %: 23.9–23.9; median `tot_corr` 1.21), so the flux columns are µJy (their `TUNIT`) and `mag_auto` includes
  the Kron-to-total factor. The header `ZP` = 28.9 does not apply to these columns.
- **Data:** DJA table (sha256 `ac1cf706…`, 17,170,560 bytes) read into memory twice (run and post hoc check); nothing kept.
- **Next:** stage 2b, a frozen plan with rule 2a in that form (disclosed as seen), otherwise PLAN_2.md unchanged.

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
