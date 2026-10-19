# W3 dimming screen: calibrate the Abell 2744 forced stage (D-052 "Revisit if": S3 cutout jobs work again)

Frozen 2026-10-19 09:11 UTC (commit 7953081), before any run.

- **Question.** D-052's headline W3 limits use MACS0416 only (34 sources); Abell 2744's forced stage was never
  calibrated because the S3 multi-target jobs failed on 2026-10-08. Can the Abell 2744 forced stage now be
  calibrated, and does the field stay null?
- **Hypothesis.** 0 surviving events (the 2026-10-08 forced run measured 32 flags and confirmed 0).
- **Ordinary explanations for any confirmed flag:** saturated-star wings and spikes (`bright_neighbour`), mosaic
  edges and low weight, deblending, persistence (single epoch), SNe / AGN (SIMBAD/NED, TNS).
- **Run.** Unchanged code and thresholds (`Params` defaults): `fetch`, `screen`, `forced --n-controls 300`
  (sequential, one S3 job), then `inject` with the measured forced inflation.
- **Calibration gate (D-052 rule, unchanged):** ≥ 50 controls in some epoch pair and control zero points in
  every epoch. Gate passed → Abell 2744 joins the headline limits; failed → reported separately, as now.
- **What changes the plan:** a confirmed flag that survives the artefact tests goes to `/vet-candidate`; a failed
  gate with ≥ 300 controls means the field's epoch overlap is too small and the field stays indicative.
