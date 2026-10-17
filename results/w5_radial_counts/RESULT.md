# W5 count deficit at the `radial` centres: result (2026-10-17)

Plan: `PLAN.md` (committed before any run). Code: `scripts/w5_radial_counts.py`. Labels: catalogue rows `observed`;
counts, ratios and nulls `derived`; Q_H and θ_E `model_prediction`; thresholds ASSUMPTIONs. Numbers below are in
`summary.json` (planned run) and `vet_near2/summary.json` (post hoc diagnostic).

## Run as planned (`summary.json`)
- 30 `radial` centres in four fields, all inside the footprint proxy (S1063 10, MACS0417 17, MACS1423 1, WHL0137 2);
  radial arc rows removed (32, 80, 41, 26).
- Stacked counts N_obs vs the control mean (5–95 %): R = 3″: 66 vs 34.5 (22–50), Q 1.91; 6″: 226 vs 128.7 (101–159),
  Q 1.76; 9″: 485 vs 285.7 (240–333), Q 1.70. Deficit p_low ≥ 0.9995 at every R → frozen outcome **null** (no
  deficit).
- **Power is 0–0.0035 at every M and R, and the plan's deficit direction was wrong (failed approach).** At these depths
  (m_lim 27.3–28.2) the fields' own background counts are shallow, so the n = 1 negative-mass model predicts an
  **excess** inside R, not a deficit: Q_H = 1.6–2.9 at 3″ for 1 × 10¹² M☉, up to 8.4–17 for 1 × 10¹³ M☉
  (`fields[].radii[].q_hypothesis`; demagnified images crowd into a smaller area faster than they drop below the
  flux limit). The "~10 % of galaxies inside θ_E are missing" figure of the W5 recipe holds for shallow-survey
  (DR10) slopes only. A deficit search at JWST cluster depth cannot test H.

## Post hoc: the observed excess (not in the frozen rule)
The excess has the sign H predicts, so it was looked at further; nothing below is a test under the plan.
1. **Arc neighbours.** Dropping also the rows within 2″ of a radial arc (`vet_near2/`, the radius of the W1 stack
   vet) lowers Q to 1.34 (3″, p_high 0.104), 1.23 (6″, 0.0575), 1.26 (9″, 0.0075). Most of the excess sits next to
   the arcs that define the centres (selection; fragments or companions of those arcs).
2. **Profile shape.** The remaining Q is nearly flat in R (1.34 → 1.26 from 3″ to 9″). H at 1 × 10¹² M☉ falls from
   1.6–2.9 at 3″ to 1.05–1.17 at 9″ (`q_hypothesis`); a flat excess out to 9″ reads as a local overdensity around
   the centres (environment or centre selection), not a θ_E-scale profile. This comparison is by eye on the stored
   numbers, not a fitted test.
3. Per field (`vet_near2`): S1063 9″ 162 vs 99.6, MACS0417 167 vs 151.6, MACS1423 13 vs 13.4, WHL0137 5 vs 11.6;
   the residual excess comes mostly from Abell S1063.

No contact sheet was made: under the frozen rule there are no flags.

## Conclusion
No deficit at the `radial` centres (frozen rule: null), but the test has no power because H predicts an excess at
JWST depth. The W5 count item at these centres is closed as specified. An excess test would need a new frozen plan
with a control matched to local source density (the centres sit in overdense spots by construction) and
injections; given the shear stack already excludes ≥ 1 × 10¹² M☉ at every centre, it is low priority.
