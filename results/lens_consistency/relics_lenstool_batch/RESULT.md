# RELICS Lenstool v1 2017-12 map batch: amplitude test (H_K; 2026-10-25)

Frozen plan `PLAN.md` (pushed in cf7ec7d8 before any deflection map was read). Numbers come from `summary.json`,
written by `scripts/relics_lenstool_batch.py`. Maps are observed (published deflection fields at D_LS/D_S = 1,
sha256-checked against `jwst_anomaly.lensmodel`); R, classes and gates are derived. No imaging was read.

## Frozen result
**Gate V fails, so H_K is undetermined** (the statistic is not trusted; plan "If V fails").

| Cluster | Role | R (median \|α_LT\|/\|α_ref\|) | 16–84 % | Class | f_obs | 1/R | Gate V |
|---|---|---|---|---|---|---|---|
| MACS J0553 | validation | 1.543 | 1.48–1.78 | offset | 0.665 | 0.648 | pass |
| Abell 3192 | validation | 1.247 | 0.93–1.48 | intermediate | 0.692 | 0.802 | **class fails** (tolerance passes) |
| RXC J2211 | validation | 1.036 | 0.97–1.14 | near 1 | 1.017 | 0.965 | pass |
| RX J0152 | validation | 0.889 | 0.78–1.00 | near 1 | 1.045 | 1.124 | pass |

Abell 3192 misses the offset class by 0.003 in R (cut 1.25, ASSUMPTION); the rule is applied as frozen.

Test clusters, numbers only (not read, gate V failed): RXC J0600 R 1.016 (near 1), MACS J0417 1.129 (near 1, against
the JWST-era CANUCS model), WHL0137 1.037 (near 1), SPT0615 0.605 (intermediate).

## Post hoc (not in the plan)
- **Why the statistic is weak.** The region |α_ref| ≥ 0.5 max covers 84–99 % of each map overlap (`region_fraction` in
  `summary.json`; RXC J0600 lowest, 0.84), not the strong-lensing annulus the plan assumed. So R is close to a
  map-wide amplitude ratio over the 120–240″ overlaps (`overlap_arcsec2`), most of which lies outside the
  multiple-image region (interpretation; the share of weakly constrained area was not measured). A3192's wide 16–84 % range (0.93–1.48)
  shows the ratio is not one number across that region.
- **What the numbers would say if the statistic were trusted** (it is not): both batch maps near 1 would reject H_K.
  This is a reading of a failed statistic and is not used.
- **SPT0615** Lenstool / GLAFIC 0.605: the Lenstool map has the smaller amplitude, the opposite of the offset maps.
  Not tested against images; no flag.

## Failed approach
- A region defined as a fraction of the reference model's peak |α| does not isolate the multiple-image region of a
  cluster deflection map. A region tied to the critical curves (e.g. |μ_ref(z = 2)| above a threshold, from maps whose
  convention gate P has checked) or to the catalogued image positions is needed.

## Deviation
- PLAN.md lists WHL0137's second model as GLAFIC v2; the registered RELICS product is GLAFIC v1
  (`WHL0137_RELICS_GLAFIC`). Used as registered.

## Reading
- The batch association of addendum bi stays untested; the addendum-bd per-cluster amplitude check stays the only
  rule before a RELICS Lenstool v1 map is used.
