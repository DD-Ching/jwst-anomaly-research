# W3 dimming screen: regenerate MACS0416 forced + inject, then a joint three-field `combine`

Frozen 2026-10-19 13:25 UTC, before any run (handoff of #371 / #369).

- **Question.** The three-field headline (docs/exotic_limits.md, D-052 addendum c) adds MACS0416 exposures as
  rounded in docs/exotic_limits.md, because MACS0416's `limits.ecsv` is not committed. Regenerating it lets
  `combine` produce the joint per-source and per-deg² limits from committed per-field tables.
- **Hypothesis.** MACS0416 stays null (D-052: 0 surviving events) and its regenerated exposures reproduce the
  rounded values used so far.
- **Run.** Unchanged code and thresholds (`Params` defaults): `fetch`, `screen`, `forced --n-controls 300
  --workers 16`, `inject`; then `combine --field macs0416,abell2744,nexus` over the three committed `limits.ecsv`.
- **Reproduction gate (function `repro_gate`, decided by code, not by eye):** per t_E at ρ = 0.1, the regenerated
  MACS0416 `exposure_source_yr` is within 10 % (ASSUMPTION) of the D-052 value used in the headline
  (1.2, 4.7, 9.9, 18.7, 53.7, 62.6 at t_E 0.01, 0.03, 0.1, 0.3, 1, 3 yr), and the forced stage is calibrated.
  Pass → the joint `combine` table replaces the summed headline. Fail → report both; the headline keeps the old
  sum and the difference is investigated (the 2026-10-08 forced run measured 375 targets with 300 controls but its
  fluxes were re-calibrated offline, so a difference is possible and is itself a result).
- **Ordinary explanations for any new forced-confirmed survivor:** saturated-star wings / spikes, edges, low weight,
  persistence (single-epoch), cosmic rays (JUMP, chromatic veto of #371, post hoc there), SNe / AGN.
  A survivor goes to the exposure-level check (`scripts/dimming_exposure_check.py`) before any limit is quoted.
- **What changes the plan:** a survivor that passes the exposure-level check → `/vet-candidate` and owner
  notification; the headline then becomes conditional.
