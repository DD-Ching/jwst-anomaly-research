# `arcs`: photo-z contamination of the 3 single-component pooled anti rows (plan, frozen before the run)

Continues `results/arcs/pooled-census/RESULT.md` (D-024 addendum 9). Exploratory: the rows were chosen after the
pooled counts were seen.

## Question and hypotheses
- Rows: WHL0137 2194 (z_cl 0.566), RXC J0600 581 (z_cl 0.43), MACS0308 2157 (z_cl 0.356); z_cl are the `z_lens`
  values of `scripts/lens_consistency.py` (assumption, as in every `arcs` run).
- H_contam (ordinary): a row is a cluster member or foreground galaxy that the photo-z put behind the cluster; it is
  then unlensed and anti with P 1/3, so it is not evidence against the model shear.
- H_bg: the row is background; the anti orientation stays unexplained by this test.
- Stored photo-z (observed catalogue values): 2194 DJA eazy z160 1.14; 581 DJA eazy z160 0.675; 2157 RELICS BPZ
  zbmin 0.671. All three lie above z_cl, so H_contam needs a photo-z failure.

## Inputs (all public; none stored in git)
- RELICS HST HLSP multi-band catalogues `hlsp_relics_hst_acs-wfc3ir_<cluster>_multi_v1_cat.txt` for `whl0137m08`,
  `rxc0600m20`, `macs0308p26` (DOI 10.17909/T9SP45; archive.stsci.edu; read once, deleted after the run).
- NED cone (1.0″) through `jwst_anomaly.crossmatch.query_matches(services=("ned",))`; DJA v4.6 NIRSpec table rows
  within 0.5″ (`macs0647_specz.TABLE`, the bytes used by #781).
- Row positions: the stored `arcs_anti.ecsv` of each field's primary-model run. Matching to RELICS: nearest within
  1.0″ after the stored HST-to-JWST offset where one exists (`photoz_offset.json`), else none (ASSUMPTION).

## Tests (in this order; first decisive one wins per row)
1. **Spectroscopic redshift** (observed): a NED redshift within 1.0″ whose NED redshift flag is not photometric, or
   a DJA grade-3 z_spec within 0.5″. `member_or_foreground` if z_spec ≤ z_cl + 0.03 (1 + z_cl); `background`
   otherwise (ASSUMPTION: 0.03 (1 + z) is ≫ a cluster velocity dispersion).
2. **Red sequence** (derived, photometric): per cluster, fit F606W − F814W vs F814W (straight line, 3σ-clipped, 3
   iterations) to RELICS sources with stel < 0.8, F814W < 23, and BPZ |zb − z_cl| ≤ 0.05 (1 + z_cl). σ_RS = clipped
   scatter. A row is `red_sequence_like` if |Δcolour| ≤ 2 σ_RS (ASSUMPTION), else `off_sequence`.
   Validation (computed, reported, gate): false-member rate f_bg = fraction of RELICS sources with zbmin > z_cl +
   0.1, F814W < 23, stel < 0.8 that are `red_sequence_like`; member recall = fraction of the fitting set within
   2 σ_RS. Gate: recall ≥ 0.8 and f_bg ≤ 0.3, else the red-sequence test is `undetermined` for that cluster.
3. **Second photo-z** (observed catalogue): RELICS BPZ zbmin for the DJA-based rows (2194, 581); none exists for 2157
   besides the stored BPZ. `photoz_disagrees` if zbmin ≤ z_cl + 0.03 (1 + z_cl).

## Per-row label
- `contaminated` if test 1 says member_or_foreground, or (no spec-z and red_sequence_like and (photoz_disagrees or
  no second photo-z)).
- `background_spec` if test 1 says background.
- `background_photometric` if no spec-z, off_sequence and not photoz_disagrees.
- `undetermined` otherwise (incl. gate failure or no RELICS match).

## Frozen outcome
- ≥ 2 of 3 `contaminated` → "open part of the pooled anti excess explained by photo-z contamination (exploratory)".
- 0 `contaminated` and ≥ 2 `background_*` → "photo-z contamination not supported; the rows stay open".
- otherwise → "mixed / undetermined".

What would change the plan: `contaminated` rows drop out of the `arcs` anti count (next: re-run the pooled tally
without them); if the rows stay background, next is a second-band shape (F150W / F444W) and the models' local shear
error, not more redshift work.
