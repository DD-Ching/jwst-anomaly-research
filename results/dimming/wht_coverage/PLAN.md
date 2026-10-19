# W3 dimming: weight-map coverage term in the screen and injection chain (frozen 2026-10-19 22:15 UTC)

Follow-up of `results/dimming/edge_coverage/RESULT.md` "Gate V2" (passed; decision `over_counts`).

## Hypothesis and what would change the plan
- Hypothesis: the committed W3 dimming headline limits (`results/dimming/macs0416_regen/limits_combined.ecsv`,
  MACS0416 + Abell 2744 + NEXUS) are optimistic because the catalogue-count coverage proxy counts mosaic-edge and
  low-weight epochs as covered (gate V2 what-if: veto-on exposure ×0.75–1.22 of committed; not a result).
- Ordinary explanations of a change: none needed; this is a bookkeeping correction. The term can only remove
  epoch coverage, so it can lower or raise (via fewer flagged baselines) the injected-source exposure.
- Plan change: if a reproduction gate fails, no corrected headline is quoted and the mismatch is diagnosed.

## Change (production, `scripts/dimming_screen.py`)
`_screen` applies the gate-V2-validated grid rule (`edge.edge_pair`): F200W i2d `WHT` every 32 px, 3 × 3 node
minimum / median positive WHT < 0.5 in a proxy-covered epoch -> that (source, epoch) is uncovered in both bands
(`weight_edge`, `mask_weight_edges`). Epochs without an F200W image are not tested. The mask is applied after the
noise scale and the ordinary (catalogue) columns and before the flags, the order of `edge_coverage/run.py`.
`screen` and `inject` use it by default; `--no-weight-coverage` restores the old proxy. Unit test:
`tests/test_dimming_screen.py::test_wht_rel_weight_edge_and_mask`.

## Gates (computable; `run.py` decides, all three fields each)
- **R0:** with the term off, the veto-on `exposure_source_yr` reproduces each field's committed `limits.ecsv`
  (`np.allclose`, rtol 1e-9).
- **R1:** with the term on, the veto-on `exposure_source_yr` reproduces
  `results/dimming/edge_coverage/runs/<field>/limits_veto_on_masked.ecsv` (rtol 1e-9), i.e. the production code is
  the validated what-if.
- Both pass for all three fields -> the corrected per-field `limits.ecsv` and `combine_limits` headline are quoted
  (derived; injections simulated), with the ratio corrected / committed per (t_E, rho). Otherwise no corrected
  headline.

## Also recorded (descriptive, no gate)
- Screen flags with the term on vs off: `flagged & ~catalogue_veto & monitored` per field, uids added / removed.
  Any added uid gets its stamps looked at before any statement.
- Veto-off exposure with the term on (the D-039 swap gain, already measured as a what-if in `edge_coverage`).

ASSUMPTIONs: thresholds 0.5 (weight), 32 px (grid), F444W footprint follows F200W; forced inflations are the
committed ones (`inject_summary.json` of each field).
