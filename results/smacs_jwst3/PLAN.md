# SMACS 0723 (program 2736) re-run on MAST's jwst 3.0.0 reprocessing — plan (written before any fetch)

TASKS "Next (M2)": "Re-fetch program 2736 after MAST reprocesses it with jwst >= 3.0". MAST serves jw02736
level-3 products built with jwst 3.0.0 since at least 2026-10-26 (CHANGELOG, #374).

- **Hypothesis:** the nine SMACS candidate records in `docs/candidates/jw02736-*` (all from jwst 2.0.1 catalogs)
  describe real sources, not products of the 2.0.1 source catalog step.
- **Ordinary explanation tested:** processing (pipeline-version) artefacts: deblending, segmentation and aperture
  changes between jwst 2.0.1 and 3.0.0 (D-010).
- **Measurement:** fetch the 3.0.0 `_cat.ecsv` catalogs only (no i2d), run `configs/reference_sample.yaml`, and for
  each candidate position find the nearest 3.0.0 source in its band (match radius 0.15″, ASSUMPTION) and its
  rank in the new run.
- **What changes the plan:** a candidate with no 3.0.0 counterpart, or one whose recorded anomaly feature changes
  by more than a factor 2 (ASSUMPTION), gets a "processing-dependent" note in its record; new top-20 entries that
  were absent before are listed for vetting. Science and control (CEERS, 3.0.0) become comparable (D-022).
