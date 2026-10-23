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

## Stage 2 (2026-10-23): gate C2 fails (16 / 23, arc-referenced 7 / 12); no flag is read

Frozen plan: `PLAN_2.md` (committed before the run; El Gordo stage 2b's rules plus deviations D1–D3, made before
any read). Run `scripts/smacs_djaref.py --stage 2`; numbers from `summary.json` `stage2`.

- **Rule 2a** (unit): 23.9 mag over 13,433 rows → µJy. **Rule 2b** (grid): rotation 0.028°, east left,
  0.0400″ per pixel. Both pass; 8,592 of 13,554 DJA rows kept after the flux and mask rules.
- **Frame tables** (`frame/`, GO 2736 o001_t001 F200W, reported, not used): pipeline gate C 37 / 39.
- **Gate C2 (shared pairs, frozen): fails.** All testable controls 16 / 23 (0.70), arc-referenced 7 / 12 (0.58);
  both need ≥ 0.8. Control classes: 14 recovered, 2 confused, 1 ambiguous, 6 absent, 8 undetectable_local,
  19 inconsistent_reference, 10 no_reference. Reported, not rules: without the containing-source rule 14 / 23;
  without the local-noise rule 23 / 31 (0.74), arc 11 / 17; plain `gate_c2` as the frozen result. All fail.
- **Outcome:** the DJA v7.4 reference is not validated on SMACS 0723, so no target is read as a flag (`flags` and
  `flags_z_edge` empty). Of 11 predicted extra-image targets none is testable with this reference anyway
  (4 inconsistent_reference, 6 no_reference, 1 undetectable_local).
- **Contact sheet looked at** (post hoc, `sheet_2.png`, `--stage sheet`: F200W, 6″ stamps, red circle =
  `search_arcsec`, the 9 absent / ambiguous / confused controls). By eye, a source lies inside the circle at 5 of
  the 6 `absent` controls (systems 12, 13, 18, 19, 21; 25 shows little above the bCG halo); `posthoc_absent_controls`
  counts all 6 as `absent_sources_too_faint` (a source is there, fainter than the predicted range). So the failing
  controls are flux-ratio misses of this reference route (the reference image's flux × μ ratio over-predicts),
  not missing images. Not diagnosed further (cause not isolated: μ ratios of the fitted-z systems 12, 13, 18, 19,
  21, 25, the 0.5″ aperture total on extended arcs, or the ICLv2 μ).

**Next:** the reference route is closed for SMACS 0723 with this model and reference. Next field: WHL0137 (gate K
first).
