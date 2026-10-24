# HST-era MACS1149 models (CATS v4.1, Sharon v4cor) vs the Schuldt et al. 2024 new images (held-out; plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any residual is computed. Follows D-034 addendum at
(`../mars_m1149_heldout/`): MARS failed gate H on this list (2.51″) while CANUCS v1, which saw the list, gave 0.849″.
Question: is that miss MARS-specific, or do HST-era parametric models fitted before the list also miss it?
Script `scripts/m1149_hst_heldout.py`; tests `tests/test_m1149_hst_heldout.py` (written before the run).

## Inputs
- model_prediction: `macs1149-cats` (CATS v4.1, pinned, D-035) and `macs1149-sharon` (Sharon v4cor maps; use limited
  to its 200″ map footprint, `../macs1149_sharon/RESULT.md`), loaded with `lens_consistency.load_model`.
- observed: the "seen" list of each model: CATS `arcs.dat` (`lens_consistency.image_list`) and Sharon's gold list
  (`../macs1149_sharon/sharon_gold_images.ecsv`).
- observed: `../mars_m1149_heldout/schuldt24_new_images.ecsv` (18 rows, Schuldt+2024 Table 2), unchanged.
  ASSUMPTION: listed z are fixed; positions are in the HST frame both models use (the CATS frame offset is < 0.02″,
  D-034, so no offset is applied unless rule 1 finds one).

## Frozen rules
Rules 1–5 of `../mars_hff_heldout/PLAN.md`, the same functions (`frame_offset`, `split_seen`, `heldout_residuals`,
`gate_h`, `outcome` in `scripts/mars_hff_heldout.py`) and thresholds (pair 1.0″, ≥ 5 pairs for an offset, seen 0.5″,
gate H rms ≤ 1.0″ and ≥ 80 % within 3.0″, ≥ 5 tested images), applied to each model separately, with that model's
seen list in place of the MARS catalogue. Grid step 0.25″ (as `mars_hff_inputs.GRID_STEP`).

## Readings fixed in advance
- Both pass: the MARS miss is MARS-specific; CATS / Sharon are validated second models on MACS1149 for a later
  two-model extra-image stage against CANUCS v1.
- Both fail with the same images > 3″: the families themselves (z, identification, MUSE positions) or a shared
  HST-era constraint gap are the likelier cause; no flag reading.
- Split result: report only; the passing model is the candidate second model.
- A large residual is **not** a flag in this stage. Ordinary causes first: model incompleteness away from the
  HST-era constraints, the 4–6 z families outside the core, MUSE-only positions (203, 204, 401), the mean-source
  estimator, Sharon's map edge.

## Stop rules
- No rule change after the first residual; any departure is labelled deviation or post hoc.
- No flag, flux or extra-image reading; no image pixels read.
- Data: model files are fetched by the existing pinned loaders into the data root (gitignored), not committed.
