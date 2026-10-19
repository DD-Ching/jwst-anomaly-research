# W3 dimming: weight-map coverage term in the screen and injection chain (gates passed; corrected headline)

Run 2026-10-19 22:13–22:18 UTC under the frozen `PLAN.md` (`run.py`, `summarise.py`). Catalogues re-fetched from
MAST, Gaia from the ESA TAP (the `persistence_scope` route), F200W `WHT` grids from S3 row reads. Numbers:
`summary.json`, `runs/<field>/summary.json`, `limits_combined.ecsv`, `runs/<field>/limits_*.ecsv`.

## Gates (frozen, computed by `run.py`)
- **R0 passed** (all three fields): with the term off, the veto-on exposures reproduce the committed `limits.ecsv`.
- **R1 passed** (all three fields): with the term on, the production `_screen` reproduces the validated gate V2
  what-if (`edge_coverage/runs/<field>/limits_veto_on_masked.ecsv`).

## Effect on the chain (derived)
| Field | weight-edge (source, epoch) pairs | injection sources off → on | area deg² off → on | compact unvetoed flags off → on |
|---|---|---|---|---|
| MACS0416 | 3,035 | 34 → 24 | 0.00299 → 0.00238 | 10 → 9 |
| Abell 2744 | 3,199 | 62 → 53 | 0.00389 → 0.00317 | 1 → 0 |
| NEXUS | 18,212 | 139 → 125 | 0.0221 → 0.0211 | 12 → 4 |

- The term removes 10 compact unvetoed screen flags (macs0416-d04053, abell2744-d05190 and 8 NEXUS uids listed in
  `runs/nexus/summary.json`) and adds none, so no new stamps needed looking at. None of the removed uids is one of
  the 5 `nexus_calib` artefact-test survivors or the 12 `nexus_o014` targets.

## Corrected headline (`limits_combined.ecsv`; derived, injections simulated; zero surviving events assumed as before)
Ratio of the corrected to the committed (`macs0416_regen/limits_combined.ecsv`) per-source limit:

| t_E (yr) | 0.01 | 0.03 | 0.1 | 0.3 | 1 | 3 |
|---|---|---|---|---|---|---|
| ρ = 0.1 | 0.82 | 1.25 | 1.30 | 1.18 | 1.19 | 1.33 |
| ρ = 0.01 | 0.92 | 1.28 | 1.20 | 1.19 | 1.07 | 1.29 |

- At t_E ≥ 0.03 yr the committed limits were optimistic by 7–33 % (per source) and the area limits by
  −3 to +23 % (`area_limit_ratio_corrected_over_committed`, range 0.78–1.23 over all cells). At t_E = 0.01 yr the
  corrected limit is tighter (×0.82–0.92); why is not tested.
- Corrected per-source limits at ρ = 0.1: 0.359 (0.01 yr), 0.140, 0.0565, 0.0341, 0.0246, 0.0209 per source per yr
  at t_E 0.03 … 3 yr.
- D-039 swap gain with the term on (veto-off / veto-on exposure): ×1.80–3.52 (`swap_gain_veto_off_over_on_range`).

## Reading (narrowest statement)
- The production screen and `inject` now uncover F200W weight-map edge epochs; the corrected three-field headline
  replaces `macs0416_regen/limits_combined.ecsv` and is up to 33 % weaker per source at t_E ≥ 0.03 yr.
- ASSUMPTIONs: weight threshold 0.5, 32-pixel grid, F444W footprint follows F200W, committed forced inflations.
- Limitation: the condition of `nexus_calib/headline.json` still holds (valid only if none of the 5 NEXUS
  artefact-test survivors is a W3 event; #371 vetted them, 0 W3 events).
