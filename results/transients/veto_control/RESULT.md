# No-injection control at vetoed SMACS injection positions: result (2026-10-16)

Plan: `PLAN.md`, committed before the measurement. Code: `scripts/arc_veto_control.py` (unit-tested in
`tests/test_arc_veto_control.py`) and `arc_diff_morph.py --known`, default parameters. All numbers are `derived`
and come from `compare.json` (per-position rows and strata) and `morph_noinj.json`.

## Inputs
- **No-injection maps:** the stored command of `results/transients/arc_chain/runs/smacs/summary.json` re-run with
  `--maps` (SMACS ERO o001 → VENUS o057, F150W + F444W, 24″ box). **Deviation:** `hits.ecsv` is not byte-identical
  to the stored file; all 73 rows have the same positions, bands, signs and flags, and `max_abs_sigma` differs by
  ≤ 5 × 10⁻¹² (checked in the session, not stored). The maps file was deleted after use.
- **Positions:** every injected star of the 3 stored `--inject-grid 80` runs (seeds 0–2,
  `results/transients/arc_chain/runs/smacs_inj*`) that is a hit and not `near_bright_core`: 67 positions
  (`targets.csv`; verdict = the injected run's morphology verdict). Off the mask: 11 rosette, 8 dipole, 5 extended,
  19 pass. In the mask: 22 pass, 2 edge.
- **Measure:** a pre-existing negative lobe = a lobe peak ≤ −3 (|S| ≥ 3, ≥ 4 px, within 0.8″; the morphology
  code's lobe definition) in either band of the no-injection map ("or" over F150W and F444W).

## Outcome against the frozen rule
The plan names no stratum, so the whole box is the primary reading; the out-of-mask stratum is reported beside it.

| Stratum | Vetoed (rosette + dipole) with a negative lobe | Passed with one | Rule outcome |
|---|---|---|---|
| whole box, both mags | 19 / 19 (1.00) | 23 / 41 (0.561) | **H1** |
| whole box, 27.5 AB | 13 / 13 (1.00) | 0.64 of 25 | H1 |
| whole box, 28.0 AB | 6 / 6 (1.00) | 0.438 of 16 | H1 |
| out of mask, both mags | 19 / 19 (1.00) | 0.789 of 19 | mixed (excess 21 pp < 30 pp) |
| in the mask | no vetoed hit | 0.364 of 22 | undetermined |

- **H2 (injection-induced lobes) is not supported:** every vetoed position already has a negative lobe without
  the injection (H2 needed ≤ 30 %).
- **H1 holds over the whole box but not as a full explanation off the mask.** Off the mask, 79 % of the passed
  positions also have a pre-existing negative lobe, so the lobe's presence is necessary for every veto observed
  but does not alone decide it.

## Post hoc: lobe strength (not in the plan)
`noinj_neg_ratio_post_hoc` = strongest pre-existing negative lobe / the injected hit's largest |S_hit| (lists in
`compare.json`, strata `*/out_mask`). Off the mask, 16 of 19 vetoed positions have a ratio ≥ 1.0, against 9 of 19
passed ones (counted from the stored lists; medians ≈ 1.41 and ≈ 0.99, not stored). The distributions overlap, so
strength separates the two groups only partly. In the mask, no ratio exceeds 0.533.

## Contact sheet (`sheet_noinj_offmask.png`, F150W no-injection S, looked at)
By-eye reading, a judgement, not a rule output. Of the 19 vetoed off-mask positions:
- most sit on a compact ± residual at the centre (a compact source the injection grid placed the star on; the grid
  uses pixels above the epoch-1 90th percentile, and the bright-core flag uses the 99.9th);
- several sit on diffraction-spike stripes or inside the residual field of a bright neighbour.

Passed positions are mostly quiet; some carry vertical stripe patterns or a weaker central ± residual.

## Reading (narrowest statement)
The off-arc efficiency loss in `review_list/RESULT.md` (R2) is not made by the injection procedure: the negative
lobes that can veto an injected star are already in the no-injection difference image at every vetoed position.
Whether the morphology rules then veto depends also on lobe strength and geometry, which this control does not
separate. Nothing here changes a rule or threshold, and nothing here is a candidate.

## Limits
- One box, one pair, two bands; 19 vetoed and 41 passed positions. The 5 `extended` verdicts are not in the test.
- The lobe that exists without the injection is not shown to be the same lobe that vetoed the injected star.
- The cause of the pre-existing residuals (compact-core PSF mismatch, spikes, stripes) is by eye only.

## Next
A morphology-rule change for crowded off-arc positions (for example judging the lobe against the no-injection
map) is a rule change: the owner's call (owner, 2026-10-16: no veto tuning without a decision).
