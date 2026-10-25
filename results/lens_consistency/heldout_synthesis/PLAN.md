# Held-out residuals across fields and models: synthesis (plan, frozen)

Frozen 2026-10-25 (commit time of this file), before any pooled, cross-model or within-family number is computed.
Inputs are the stored per-image held-out tables of D-034 addenda ap, at, au, aw, ax, ay; no new data, no network.
Script `scripts/heldout_synthesis.py`; tests `tests/test_heldout_synthesis.py` (written before the run).

## Question, ordinary explanations, what changes the plan
- Question: across the 11 stored model × field held-out tests, is any JWST-era image missed in a way that model
  imprecision does not explain, i.e. missed by two independent published models **and** isolated within its family?
  And how large must an image-position offset be to stand out from published-model error on each field?
- Ordinary explanations of a miss: model incompleteness away from its constraints, coarse maps, wrong listed z,
  a misidentified image (wrong family assignment), the mean-source estimator. A shared, isolated miss points first
  to misidentification; it is a vetting target, not a candidate.
- A flag sends the image to `/vet-candidate`. No flag closes the position-residual reading of these tables; the
  per-field floors (rule 4) become the stated limit of the lens-model-residual line until finer models appear.

## Inputs (all model_prediction residuals of observed positions; thresholds ASSUMPTION)
`results/lens_consistency/<dir>/heldout_<model>.ecsv`, columns `image_id system type z dtheta_arcsec n_predicted
inside`. Primary models (passed their own training gate and read in a verified map unit): MARS on Abell 370,
MACS0416, MACS1149; CATS v4.1 and Sharon v4cor on MACS1149; RELICS glafic v2 and zitrin-ltm-gauss v2 on SMACS 0723;
Bergamini+23 on Abell 2744; Caminha+23 on El Gordo. Reported only: MARS on Abell 2744 (gate K failed), RELICS
lenstool v1 on SMACS 0723 (map normalisation unverified, addendum aw).

## Frozen rules (functions with unit tests before the run)
Only rows with `inside` true enter. A row with no predicted image (`n_predicted` 0 or `dtheta_arcsec` not finite)
is a miss with offset +inf.
1. **Per-model statistics** `model_stats`: n, rms (finite rows), median, 68th and 95th percentile (inf counted),
   fraction > 3.0″.
2. **Isolated miss** `isolated_misses`: an image with offset > 3.0″ whose family (`system`) has ≥ 2 other tested
   images in the same model, all ≤ 1.0″ (ASSUMPTION).
3. **Shared miss** `shared_misses`: per field, an `image_id` with offset > 3.0″ in ≥ 2 primary models of that field.
4. **Flag** `flags`: a shared miss that is an isolated miss in at least one of those models.
5. **Floor** `field_floor`: per field, the 95th-percentile offset of the primary model with the lowest rms. An
   image-position offset below it does not stand out from published-model error at about 5 % per image
   (ASSUMPTION: held-out images represent future unknown images).

## Stop rules
- No rule change after the first number is computed; any departure is labelled deviation or post hoc.
- Fields with one primary model (Abell 370, MACS0416, Abell 2744, El Gordo) cannot give a shared miss; their
  isolated misses are reported, not flagged.
