# W3 dimming: NEXUS `o014` depth/footprint test and forced photometry on the 12 `on_sky` persistence targets

Frozen 2026-10-19 17:25 UTC, before any forced flux was read (handoff of #377,
`results/dimming/persistence_run/RESULT.md` Next 1). Screen code, thresholds and `forced` unchanged.

- **Targets.** The 12 NEXUS compact single-epoch vanish flags whose D-039 persistence verdict was `on_sky`
  (`persistence_run/nexus/verdicts.ecsv`). All are catalogue-detected only in epoch 2 (`jw05105-o014`).
- **Hypotheses.** (H0) Each is a real source present in every epoch that the other epochs' pipeline catalogues
  miss (catalogue depth, deblending or footprint), i.e. not a vanish. (H1) It is fainter by the screen's vanish
  criterion in the other epochs in forced photometry at fixed position.
- **Ordinary explanations tested here:** (a) footprint: the other epochs' F200W mosaics do not cover the position
  (`forced` cut flags / NaN); (b) depth: `o014` is deeper than the other epochs, so a source at its catalogue limit
  is not catalogued elsewhere. Not tested here: SN / AGN / variable star, which go to `/vet-candidate` for survivors.
- **Measurement.** Rebuild the NEXUS screen (catalogues from S3 by sha256, Gaia mask from the ESA TAP as in #376/#377),
  then `dimming_screen.py forced` with the 12 targets as the only flags (every other flag marked
  `catalogue_veto` in a copy of `lightcurves.ecsv`; the 12 targets' `catalogue_veto` cleared because D-039 `on_sky`
  replaces `single_epoch` for them) and the default 300 matched controls.
- **Depth (derived, descriptive).** Per epoch, the median forced F200W 1σ error of the controls (µJy); and per
  target, its `o014` forced S/N and its S/N expected in each other epoch at constant flux.
- **Decision rule (computed by the pipeline, not by eye).** A target survives iff `survives_artefact_tests`
  (`forced_confirmed` with the vanish flag reappearing in forced photometry, cutout-clean, no chromatic veto;
  `catalogue_veto` cleared as above). Survivors go to `/vet-candidate`. A non-survivor is explained by the first of:
  `forced_vanish` False with finite other-epoch fluxes → (b) catalogue depth / deblend (source present in forced
  photometry); other-epoch fluxes NaN or cut-flagged → (a) footprint / image quality.
- **Calibration requirement.** The run counts only if `forced_summary.json` `calibrated` is true and
  `controls_flagged` ≤ 1 % of controls (ASSUMPTION, as in the earlier fields).
- **What would change the plan:** any survivor (vetting before the D-039 swap); a calibration failure (stop and
  record).
