# W3 dimming: cost and value of the D-039 persistence test in place of the `single_epoch` veto

Frozen 2026-10-19 15:20 UTC, before any run (handoff of #375; D-052 "Revisit if", second bullet).

- **Question.** The `single_epoch` veto removes every vanish of a source caught in only one detection-band epoch.
  Replacing it with the D-039 level-2 persistence test would restore full-vanish sensitivity, but costs one
  `_cal` exposure set per tested source. Is the swap worth building, field by field?
- **Measurement (scoping only; no code change to the screen, thresholds unchanged).**
  1. Workload: per field (MACS0416, Abell 2744, NEXUS), the number of `screen` flags that are vetoed by
     `single_epoch` alone, split into compact (point-like) and all, and how many have the vanish flag.
  2. Value: re-run `inject` with the per-copy `single_epoch` veto off (what-if, `simulated`): the gain in
     100 % / 50 % dimming efficiency and in W3 efficiency per t_E at ρ = 0.1 over the committed limits.
- **Ordinary explanations the swapped-in test must still cover:** persistence afterimages (D-039), cosmic rays
  (JUMP, chromatic veto), SNe / AGN (vet-candidate), edges.
- **Decision rule (ASSUMPTION).** Build the swap if the what-if tightens the joint headline (any t_E, ρ = 0.1) by
  ≥ 20 % and the compact single-epoch-only workload is ≤ 50 sources over the three fields. Otherwise record the
  numbers and leave the veto in place.
