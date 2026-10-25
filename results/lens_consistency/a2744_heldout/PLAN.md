# Abell 2744: Bergamini+2023 and MARS vs UNCOVER v2.0 images absent from their catalogues (held-out; plan, frozen)

Frozen 2026-10-25 (commit time of this file), before any residual is computed. Follows D-034 addenda ap / au / aw:
held-out validation (gate H) passes so far only for MARS on MACS0416. Abell 2744 has a JWST-era reference route
(UNCOVER v2.0, `../abell2744_uncoverref/`) and an early two-model comparison against Bergamini+2023
(`../abell2744-uncover2/crossmodel.json`), but neither second model was ever tested on images it was not fitted to.
Script `scripts/a2744_heldout.py`; tests `tests/test_a2744_heldout.py` (written before the run).

## Inputs
- model_prediction: `abell2744-bergamini23` (Bergamini et al. 2023b Lenstool, JWST GLASS-era, pinned), loaded with
  `lens_consistency.load_model`; its "seen" list is its own `arcs.dat` (`lens_consistency.image_list`).
- model_prediction: MARS Abell 2744 free-form map (Cha & Jee 2023, Zenodo 10.5281/zenodo.7804575, md5 checked,
  arcsec reading), loaded as in `mars_hff_heldout.field_block`; "seen" = MARS catalogue rows of Abell 2744.
  Its training gate K failed on rms (1.112″, addendum ao), so its gate H is **reported only**: a MARS pass cannot
  make it a validated model here.
- observed: `../abell2744_uncoverref/uncover2_images_spec.ecsv` (UNCOVER v2.0 images with spectroscopic z, unchanged).
  ASSUMPTION: listed z are fixed.

## Question, ordinary explanations, what changes the plan
- Question: does the Bergamini+2023 model reproduce the UNCOVER v2.0 spectroscopic images it did not use?
- Bergamini passes gate H: a held-out-validated second model on Abell 2744, independent of UNCOVER v2.0 on those
  images; the next stage is a frozen two-model extra-image check (UNCOVER-only vs both) re-using the stored
  crossmodel tables. Fails: no validated second model on Abell 2744; the route waits for another model.
- Ordinary causes of a fail, first: model incompleteness away from its constraints, UNCOVER family / z assignments,
  knot ids that sit < 0.5″ apart (a held-out knot may be called seen), the mean-source estimator.
- A large residual is **not** a flag in this stage; no flag, flux or extra-image reading; no image pixels read.

## Frozen rules
Rules 1–5 of `../mars_hff_heldout/PLAN.md`, the same functions (`frame_offset`, `shift`, `split_seen`,
`heldout_residuals`, `gate_h`, `outcome` in `scripts/mars_hff_heldout.py`) and thresholds (pair 1.0″, ≥ 5 pairs
for an offset, seen 0.5″, gate H rms ≤ 1.0″ and ≥ 80 % within 3.0″, ≥ 5 tested images), applied to each model
separately with its own seen list. Grid step 0.25″. Family = `lensmodel.image_family` of the UNCOVER id (knots
"1.1a" / "1.2a" are separate families, as in the reference route).

## Readings fixed in advance
- Bergamini passes: candidate second model (above). Bergamini fails: report; no second model.
- Images > 3″ in both models: family-level causes (z, identification) are likelier; no flag reading.
- MARS: reported; compared image by image with Bergamini (post hoc only if not listed here: the shared > 3″ list
  above is in the plan).

## Stop rules
- No rule change after the first residual; any departure is labelled deviation or post hoc.
- Data: model files go to the data root (gitignored) / the scratchpad (MARS zip, deleted after the run).
