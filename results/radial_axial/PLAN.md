# Radial (W1) screen: centre-level axial (collinear-streak) rule — plan (2026-10-17, written before any run)

Handoff of `results/radial_streak/RESULT.md`: `findsat_mrt` masks miss the faint MACS0647 streak, and an arc-level
veto thins the null. This unit tests a centre-level rule instead (D-TBD).

- **Hypothesis (to falsify):** the MACS0647 top `radial` centre is a streak, so its contributing arcs lie on one axis
  through the centre; a dark radial lens instead puts radial images at spread azimuths.
- **Ordinary explanations kept in view:** satellite/scattered-light streaks, diffraction spikes (already vetoed by
  `--spike-stars`), chance alignments of field galaxies.
- **Rule (ASSUMPTION, fixed now):** for a centre, φ_i = azimuth of each contributing arc around the centre; axial
  concentration A = |mean exp(2iφ_i)|. A centre with ≥ 3 lines and A ≥ 0.9 is `axial` and does not count, in the
  data and in every null draw alike (so the null is not thinned).
- **Gate K (known case):** MACS0647's 6-line centre must be `axial`. If not, the rule is not adopted.
- **Gate I (injections):** injected radial patterns recovered by the existing radial injections must lose ≤ 5 % of
  their recoveries to the rule. If more, not adopted.
- **What changes the plan:** a control-passing field whose flag status changes under the rule is reported and looked
  at before anything else.
