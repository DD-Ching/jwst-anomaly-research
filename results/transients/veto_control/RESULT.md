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
- **Measure (as planned):** a pre-existing negative lobe = a lobe peak ≤ −3 (|S| ≥ 3, ≥ 4 px, touching the 0.8″
  disk; the morphology code's lobe definition) in either band of the no-injection map ("or" over F150W and F444W).
  This includes the region under the position itself, which the plan did not exclude (see the post hoc split).
- **Deviation:** the plan names per-position `extent_arcsec`; it is in `morph_noinj.json`, not in `compare.json`,
  and is not used.

## Outcome against the frozen rule
The plan names no stratum. Its question is about off-mask vetoes, and every vetoed hit is off the mask, so the
**out-of-mask stratum is the relevant reading**. The whole-box row is diluted by 22 in-mask passed positions with
no vetoes beside them; it is reported, not used.

| Stratum | Vetoed (rosette + dipole) with a negative lobe | Passed with one | Rule outcome |
|---|---|---|---|
| **out of mask, both mags** | 19 / 19 (1.00) | 15 / 19 (0.789) | **mixed** (excess 21 pp < 30 pp) |
| out of mask, 27.5 AB | 13 / 13 | 9 / 12 (0.75) | mixed |
| out of mask, 28.0 AB | 6 / 6 | 6 / 7 (0.857) | mixed |
| whole box, both mags (diluted) | 19 / 19 | 23 / 41 (0.561) | H1 |
| in the mask | no vetoed hit | 8 / 22 (0.364) | undetermined |

- **H2 as defined (≤ 30 % of vetoed positions with pre-existing negative structure) is not met:** 19 / 19 vetoed
  positions carry negative structure without the injection.
- **H1 as defined is not met off the mask:** passed positions carry such structure almost as often.

## Post hoc split: under the star or around it (not in the plan)
`noinj_centre_negative_post_hoc` = the no-injection S at the position is ≤ −3 in some band (`compare.json`,
stratum `all/out_mask`):

| Off the mask | Centre negative | Lobe present, centre not negative | No lobe |
|---|---|---|---|
| vetoed (19) | 11 | 8 | 0 |
| passed (19) | 8 | 7 | 4 |

For 11 vetoed positions the counted "lobe" may be the residual under the star, which a positive star partly fills,
so the measure does not show that a lobe around the star vetoed it. `noinj_neg_ratio_post_hoc` (also stored) uses
global region peaks that may lie outside the disk and may mix bands; it is not interpreted.

## Contact sheet (`sheet_noinj_offmask.png`, F150W no-injection S, looked at)
By-eye reading, a judgement, not a rule output. Of the 19 vetoed off-mask positions:
- most sit on a compact ± residual at the centre (a compact source the injection grid placed the star on; the grid
  uses pixels above the epoch-1 90th percentile, and the bright-core flag uses the 99.9th);
- several sit on diffraction-spike stripes or inside the residual field of a bright neighbour.

Passed positions are mostly quiet; some carry vertical stripe patterns or a weaker central ± residual.

## Reading (narrowest statement)
Every vetoed off-mask injection sits on negative difference-image structure that exists without the injection,
but so do most passed ones, and this control does not separate the two groups (frozen outcome: mixed). It does not
show what makes the morphology code veto one and pass the other. Nothing here changes a rule or threshold, and
nothing here is a candidate.

## Limits
- One box, one pair, two bands; 19 vetoed and 41 passed positions. The 5 `extended` verdicts are not in the test.
- The lobe that exists without the injection is not shown to be the same lobe that vetoed the injected star.
- The cause of the pre-existing residuals (compact-core PSF mismatch, spikes, stripes) is by eye only.

## Next
To separate the groups, measure per position the negative lobes outside the region under the star (location-aware,
with local peaks; `arc_diff_morph` `Regions.local_peak`) on both the injected and the no-injection maps; this needs
the injected maps re-made (≈ 1.5 min per seed). A morphology-rule change for crowded off-arc positions (for example judging the lobe against the no-injection
map) is a rule change: the owner's call (owner, 2026-10-16: no veto tuning without a decision).
