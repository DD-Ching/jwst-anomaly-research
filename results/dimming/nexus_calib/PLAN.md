# W3 dimming screen: calibrate the NEXUS forced stage (D-052 "Revisit if"; handoff of #367)

Frozen 2026-10-19 10:19 UTC (commit 73785280), before any run.

- **Question.** NEXUS (8 epochs) is the one D-052 field whose forced stage was never calibrated; its uncalibrated
  W3 efficiency (0.07–0.24) is what could move the headline. Does its forced stage pass the D-052 gate, and does it
  stay null?
- **Hypothesis.** 0 surviving events (the 2026-10-08 forced run measured 100 flags and confirmed 0).
- **Ordinary explanations for any confirmed flag:** saturated-star wings and spikes (`bright_neighbour`), mosaic
  edges and low weight, deblending, persistence (single epoch), SNe / AGN (SIMBAD/NED, TNS).
- **Speed (code change, no change to any measurement).** The sequential forced stage needs ≈ 16 images × 1,351
  targets of S3 range reads (> 1 h at Abell 2744's rate). `forced --workers N` reads the (band, epoch) images
  concurrently; everything after the reads is unchanged. **Validation gate:** `forced --workers 6` on Abell 2744
  must reproduce `results/dimming/abell2744_calib/forced_summary.json` exactly (except `wall_s`); if not, the
  option is not used.
- **Run.** Unchanged thresholds (`Params` defaults): `fetch`, `screen`, `forced --n-controls 300` (all flags,
  `--max-flags 0`), then `inject` with the measured forced inflation.
- **Calibration gate (D-052 rule, unchanged):** ≥ 50 controls in some epoch pair and control zero points in
  every epoch. Passed → NEXUS joins the headline limits; failed → reported separately.
- **What changes the plan:** a confirmed flag that survives the artefact tests goes to `/vet-candidate`; a failed
  gate with ≥ 300 controls means the field stays indicative.
