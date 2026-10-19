# W3 dimming: does the injection chain count mosaic-edge positions as covered? (gate V failed)

Run 2026-10-19 19:15–19:22 UTC under the frozen `PLAN.md`. Screen rebuilt from the S3 catalogues and the ESA Gaia
mask (the `persistence_scope` route). The `veto_on` exposures reproduce the committed `limits.ecsv` of all three
fields exactly (`summary.json` `reproduces_committed_veto_on`). F200W `WHT` read on a 32-pixel grid, 15 images,
4–158 s per field (`runs/<field>/edge_summary.json`). Numbers: `summary.json`, `runs/`, `gate_v_diagnostic.json`.

## Gate V: failed (frozen outcome)
The rule `edge.edge_pair` calls an edge pair in an undetected epoch for **5 of the 11** NEXUS targets that
`nexus_o014` placed at mosaic edges; the plan required ≥ 9. **By the plan, the weight-map tool is not trusted and
no decision is read.** `summary.json` `decision` = `no_reading_gate_v_failed`.

## Why the gate failed (post hoc diagnostic, `gate_v_diagnostic.json`)
- The gate set was mis-specified. `forced`'s `cut_edge` / `cut_low_weight` are OR'ed over **all** epochs, including
  epochs the screen's proxy already calls uncovered, which do not matter here. Three of the six uncalled targets
  (d13073, d20291, d22021, all in the earlier `catalogue_depth_or_deblend` class, i.e. detected by `forced` in the
  covering epochs) have `wht_rel` 0.807–0.998 in every proxy-covered epoch other than `o014` (e2).
- All 4 targets that `nexus_o014` found with "no other F200W epoch measurable" (d10423, d18183, d29073, d30717)
  are called: their only other proxy-covered epoch has `wht_rel` 0.
- Not called and not explained by the gate set: d29638 (`wht_rel` 0.97–0.99 in its covered e4 / e7, where `forced`
  saw a low-weight stamp: a feature finer than the 32-pixel grid is one possible reading, not tested), d15524 and
  d24470 (covered e0 `wht_rel` 0.516 / 0.752, just above the 0.5 cut).
- Sheet (`sheet.png`, MACS0416 and Abell 2744, looked at): the red edge pairs lie on outer mosaic borders, chip-gap
  borders, and some fully off the data.

## Measured anyway (post hoc; the tool did not pass its gate, so these are not results)
| Field | injection sources | proxy-covered pairs | edge pairs | f_edge | sources with an edge pair |
|---|---|---|---|---|---|
| MACS0416 | 34 | 90 | 16 | 0.178 | 15 |
| Abell 2744 | 62 | 125 | 11 | 0.088 | 11 |
| NEXUS | 139 | 393 | 52 | 0.132 | 45 |

What-if with edge pairs set uncovered (`simulated`, ρ = 0.1, joint over the three fields): veto-off exposure ×0.81–0.87
of the unmasked veto-off exposure (by t_E); veto-on exposure ×0.75–1.22 of the committed one. The veto-off / veto-on
gain is ×1.80–3.37 masked against ×1.99–3.36 unmasked. If a re-gated tool confirms this, the frozen rule would read
`over_counts`, and the committed headline limits would also be optimistic by up to about 25 % at t_E = 3 yr (not
established).

## Reading (narrowest statement)
- No decision. The question stays open: the plan's gate was mis-specified, not shown to be failed by the tool on
  the cases it should catch.
- The proxy counts some positions as covered that the `WHT` maps place outside the data (sheet; the 4 "no other
  epoch measurable" NEXUS targets).

## Next
1. Re-freeze gate V on cases that test the rule: the per-epoch `forced` stamp flags (`cut_low_weight`, `cut_edge`)
   in proxy-covered epochs only, for the 37 targets and the 300 controls of each field, and require agreement
   with `edge_pair` both ways (false-edge rate on interior controls too). A finer grid (8 px) where needed.
2. If it passes, read the decision; if `over_counts`, add a weight-map coverage term to `build_light_curves`
   coverage before the D-039 swap and re-run `inject` for the headline.
