# MACS1149: D-041 cross-model rule with the pinned `macs1149-sharon`: gate passes; one frozen flag (system 16 at its fitted z), explained by D-037

Plan: `PLAN.md` (frozen and committed before the run). Script: `scripts/m1149_crossmodel.py`. Numbers: `summary.json`.
Positions and μ are `model_prediction` (CATS v4.1, Sharon v4cor); catalogued positions `observed`; forced classes
`derived` (CATS `images` run of 2026-10-14). No photometry was re-measured. Thresholds are ASSUMPTIONs.

**Gate G (known cases): passes.** Sharon reproduces 38 / 39 catalogued images of the 14 tested systems within 1.5″
(rule ≥ 0.8). The one miss is 16.1 at the CATS-fitted z 4.419. Every catalogued image and every CATS row lies inside
the Sharon maps.

**Frozen outcome (16 rows):** `no_flag` 15, `lensing_residual_input` 1, `untestable` 0.
- D-041 readings: `both models` 11, `model-dependent (CATS only)` 4, `catalogued image misplaced by CATS` 1 (13).
- The 11 `no_reference` rows read the same as in the scratch-Sharon run of 2026-10-14 (`no_reference_refs.json`;
  0 disagreements in `summary.json`), so the pinned entry (z_lens 0.543) changes none of those D-041 readings.
- Knot 1192 (`confused`): Sharon predicts an extra image 0.44″ away (D-037 quoted 0.56″ from the scratch run);
  not absent, no flag.

**The flag: system 16's CATS image at (177.41115, 22.389165), `absent`, μ_CATS 15.8.** At the CATS-fitted z 4.419
Sharon also places an extra image 0.90″ away (μ 30.2), so the frozen rule returns `lensing_residual_input`.
Ordinary explanation, from committed numbers and the existing record:
- At z 4.419 Sharon does not reproduce the catalogued pair: it reproduces 16.2 only (`summary.json`, `reproduced`).
  So both models' extra image at this z comes from a source position that fails the system's own constraints
  (D-037 addendum 2026-10-14: pair rms 4.3″ in CATS and 3.9″ in Sharon at 4.419).
- Pre-declared secondary at z 2.7 (D-037 addendum: z ≈ 2.5–2.9 from both models and the CANUCS photo-z of 16.2):
  Sharon's nearest extra image is 20.2″ away and the reading is `model-dependent (CATS only)`.
- D-037's "fitted redshift" rule already classes this row untestable; the cutouts were looked at then (empty sky,
  −0.4σ after high-pass). The frozen outcome function here did not include that rule (a plan gap, stated as a
  deviation below).

**Reading.** Not a candidate: the flag rests on a model-fitted redshift at which neither model fits the system.
MACS1149's forced-photometry verdicts now rest on pinned inputs only; the D-037 tally (0 surviving) is unchanged.

**Deviations.** The frozen `outcome` omitted D-037's fitted-redshift rule; the system 16 reading above applies it and
the pre-declared secondary after the run (post hoc reading, not a new rule). Diagnostic only, post hoc: at z 2.7
Sharon reproduces only 16.1 under the nearest-image assignment, while D-037's own z scan reproduced both 16.1 and
16.2 within 1.0″ at z 2.5–3.3 (different assignment: two images per pair, `m1149_sys16.py`); not used.

**Limits.** No new photometry, no injection: this is a model-geometry cross-check of existing forced classes. The gate
uses the CATS image list, which both models were partly fitted to, so it checks that the solve works, not that the
models predict independent data.
