# W1 shear stacked at `radial` centres: result (2026-10-17)

Plan: `PLAN.md` (committed before any run). Code: `scripts/w1_stack_radial.py`. Labels: shapes `observed`;
S values, stacks and nulls `derived`; thresholds ASSUMPTIONs. Every number below is in one of the three
`summary.json` files named next to it.

## Run as planned (`summary.json`)
- 30 `radial` centres in four fields; **19 have a valid S** (≥ 5 sources in the aperture): Abell S1063 5,
  MACS0417 14; the MACS1423 centre and both WHL0137 centres have too few shear sources.
- Stacked mean S = **0.72**, S_× = −0.02; rotation null (2000 draws) median 0.006, 95th percentile 0.45;
  **one-sided p = 0.0035** → frozen outcome "vet the centres".
- Sign check passed as written (mean S with the cluster shear kept is negative in all four fields), but it is weak
  (−0.007 to −0.15): off-centre 10″ apertures see the cluster's shear as nearly uniform, so they barely respond to
  it. **Post hoc:** this check does not test the sign; the sign convention is instead covered by a synthetic test
  (`tests/test_w1_stack_radial.py`: injected radial shear gives S > 3, tangential S < −3; written after the plan,
  **deviation**) and by the D-053 injections that use the same estimator.

## Vetting (both steps **post hoc**: PLAN names "vet the centres" but not these tests)
1. **Circularity: the arcs that define the centres are also shear sources** (`vet_no_radial_arcs/summary.json`).
   The `radial` screen's selected arcs (`exotic_screens.radial_candidates`, 85 used rows over the four fields) point
   at the centres by construction. Without them: mean S = **0.51**, S_× = 0.004, p = **0.0265** (null 95th
   percentile 0.45) → rule outcome "exploratory note".
2. **Neighbours within 2″ of those arcs also dropped** (radius chosen post hoc; `vet_no_radial_arcs_2as/summary.json`,
   225 rows in all): mean S = 0.48, S_× = 0.10, p = 0.071 (null 95th percentile 0.53). The neighbours are not
   fragments of the arcs: their mean cos 2(PA − PA_arc) is −0.11 ± 0.14 (S1063), −0.23 ± 0.10 (MACS0417), −0.27 ± 0.10
   (MACS1423). The stacked mean moves only 0.51 → 0.48; the p change comes mostly from the null's 95th percentile
   (0.45 → 0.53). This step is a diagnostic, not a test.

## Conclusion
Most of the planned excess (mean S 0.72 → 0.51) is the circularity in step 1. What remains (p 0.027, 19 centres,
two fields, 14 of them in MACS0417 with overlapping apertures) is an exploratory note under the frozen rule, not a
candidate: the `radial` centres get no independent support from the shear of other sources at the p < 0.01 level.
Not looked at: a sheet of the highest-S MACS0417 centres (only `results/radial_more/sheet_radial_top.png`, which
shows each field's top centre, was looked at; MACS0417's top centre has S −0.16).

## Ordinary explanations not tested here
- Cluster-model error: `radial` arcs are selected as *anti* to the model's stretch; where the model's shear is
  locally wrong, the residual shear of other sources after model removal can correlate with them. Test: the same
  stack with a second published model per field.
- Spatially coherent systematics (PSF residuals, faint streaks): the rotation null removes them; S_× near 0
  (−0.02 to 0.10) gives no sign of them, but no streak veto was applied (ASSUMPTION in PLAN).

## Reach (rough)
In the D-053 injections the median best-position S of a 2 × 10¹² M☉ W1 lens is 2.05–2.44 in these two fields,
against 0.94–1.06 at 2 × 10¹⁰ M☉ (`results/w1_streak_veto/inject_summary_{abells1063,macs0417}.json`). If every
valid centre held such a lens, the stack would sit well above the null 95th percentile (0.45); this comparison
uses a peak statistic and is not a limit. A limit needs injections at the centres through this stack.
