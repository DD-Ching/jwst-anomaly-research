# Halo-subtracted flux ratios on Abell 370 (CANUCS model): one system testable, no anomaly

Exploratory (gates set on other clusters; nothing here is held out). Plan and the pre-run deviation: `PLAN.md`.

- **Inputs:** `abell370-canucs` (D-044); CANUCS 1208 o002 t001 F150W + F444W `_i2d` by S3 byte range; unchanged
  `halo_phot.py` gates; ±1 mag end-to-end injections (simulated).
- **Per-system gate (deviation, frozen before photometry):** `lens_consistency.systems_passing_rms` on this run's
  `validate/imageplane.ecsv` (global image-plane rms 2.32″, `validate/validate.json`): 27 of 45 systems, 64 of 115
  images, pass rms ≤ 0.5″ with every image solved.
- **Tested (derived, `halo_phot.json`):** 2 of 64 images, both of system 9 (9.1, 9.3), in both bands; 62 `untested`.
  Residuals 0.09–0.23 mag (F444W medfilt 0.094, annulus 0.228; F150W 0.202 / 0.232): `consistent`, no
  `flux_outlier`. 9.2 is untested; its stamps show several bright neighbours within 1.5″ (looked at).
- **`fluxratios --halo-phot` (F444W):** 2 images per route, 0 `flux_outlier`, 0 colour outlier, `flux_tension`
  untested on all 64 (no σ_sys fit from 2 images): `fluxratios_f444w_{medfilt,annulus}.json`.
- **Reach (simulated, `halo_phot.json` "injections"):** the −1 mag (brighter) injection is recovered as
  `flux_outlier` for 1 of the 2 tested targets per band; +1 mag (fainter) makes both targets `untested` (no recovery).
- **Why so few:** of the 64 F444W images, sequential gate losses are |μ| > 20: 6, concentration: 38, axis ratio:
  14, S/N: 4 (computed in the session from `halo_phot.json`'s stored fields, not stored as a table).
- **Stamps looked at** (`system9_stamps.png`, 3″): 9.1 and 9.3 are compact sources; 9.1 has a second faint knot
  ≈ 0.4″ away in F150W.
- **Outcome:** no anomaly among the 2 images tested; the test has almost no power in Abell 370 (one system, a
  ±1 mag change recovered in at most 1 of 2 cases). Not a limit on the cluster's flux ratios.
