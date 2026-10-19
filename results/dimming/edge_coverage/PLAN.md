# W3 dimming: does the injection chain count mosaic-edge positions as covered? (plan, frozen before data)

Frozen 2026-10-19 19:22 UTC (cloud routine), before any weight map was read. From
`results/dimming/nexus_o014/RESULT.md` Next 2.

- **Question.** The screen calls an epoch "covered" for a source when ≥ 3 catalogue rows lie within 5″
  (`Params.footprint_min`, `footprint_radius_arcsec`). `forced` found 11 / 12 NEXUS `on_sky` vanish targets at or
  beyond the other epochs' mosaic edges, so the proxy counts edge positions as covered. Does it also do so for the
  injection set, and by how much does that inflate the veto-off what-if exposure of `persistence_scope` (the value
  case for the D-039 swap)?
- **Hypothesis H1:** a material share of injection pairs sit at mosaic edges / low weight and the what-if gain is
  inflated. **H0 (ordinary):** the proxy is good away from the few edge sources; the gain stands.
- **Measurement.** Per field (NEXUS, MACS0416, Abell 2744), rebuild the screen exactly as `persistence_scope`
  (S3 catalogues, ESA Gaia mask). For each F200W epoch read the i2d `WHT` extension on a 32-pixel grid (S3 row
  reads). For each master source and epoch, `wht_rel` = the smallest WHT of the 3 × 3 grid nodes around the source
  divided by the median positive WHT of that epoch's grid. **Edge pair** (rule `edge.edge_pair`, unit-tested):
  proxy-covered (finite F200W flux) and `wht_rel < 0.5` (the `cut_low_weight` threshold of `forced`), NaN counting
  as 0.
  1. `f_edge`: share of edge pairs among the proxy-covered pairs of the injection set (`monitored` and unflagged).
  2. What-if (`simulated`): set edge pairs uncovered (NaN flux and error, not detected) in both bands (ASSUMPTION:
     the F444W footprint follows F200W; both channels see the same field), re-derive flags and `monitored`, and
     re-run `efficiency_table` with the per-copy `single_epoch` veto on and off (`persistence_scope/whatif.py`
     logic, same seed, same forced inflation).
- **Validation gate V (before any reading).** The rule must call an edge pair, in at least one epoch where the
  target is not catalogue-detected, ≥ 9 of the 11 NEXUS targets that `nexus_o014` placed at other epochs' mosaic
  edges (`cut_edge` or `cut_low_weight`). Fails → the weight-map tool is not trusted; record and stop.
- **Decision rule (ASSUMPTION; function `edge.decision`).** With V passed: if the edge-masked veto-off joint
  exposure (three fields, ρ = 0.1) is below the unmasked veto-off exposure by > 10 % at any t_E, the chain
  over-counts edges, and the swap needs a weight-map coverage term in screen and injections before its gain is
  quoted (the corrected gain is reported). Otherwise the scoped gain stands (still an upper bound, as scoped).
