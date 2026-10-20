# `arcs`: shape-noise null calibrated on weak-shear rows

Plan: `PLAN.md` (frozen before the weak-row orientation classes were read; the strong-shear excess of #406 was
known). Script: `scripts/arcs_null_calibration.py`; every number below is in `summary.json`. No new data: the
weak-shear (`strong_shear` False, g_min < 0.2) background rows of the same deciding-run `arcs.ecsv` tables as
`results/arcs/pooled-tally/`, with the same cuts (e ≥ 0.5, S/N, size, radius).

## Frozen outcome: **"calibration explains the pooled excess"** (derived), on a point estimate only
The rule fires because pooled B reaches P 0.250 at the isotropic-fraction point estimate f̂ = 0.06. That estimate's
95 % interval is 0–0.58 and includes f = 0 (the uncalibrated null, P 0.0285). At σ̂_ε (0.27), pooled B stays at
P 0.0438. **Post hoc reading:** the label overstates the evidence. The weak rows agree with the σ_ε 0.26 null and
cannot tell it apart from the calibrations that remove the excess. In substance this is "weak rows lack power". The
plan's rule should have asked for an interval that excludes the uncalibrated null (failed approach, below).

## Weak-shear rows (observed counts vs derived expectations)
| Pool | Rows | Anti | Expected (σ_ε 0.26) | Two-sided P | σ̂_ε (95 %) | f̂ at σ_ε 0.26 (95 %) |
|---|---|---|---|---|---|---|
| A (six fields) | 75 | 9 | 11.06 | 0.67 | 0.24 (0.18–0.46) | 0.00 (0.00–0.36) |
| B (A without MACS1423) | 61 | 9 | 9.14 | 1.0 | 0.27 (0.19–0.60, upper at the grid edge) | 0.06 (0.00–0.58) |

S1063 (not in the pool): 14 weak rows, 2 anti vs 1.68 expected. Per field: `weak_by_field`.

## Strong-shear rows under each calibration (P(≥ observed), Poisson)
| Setting | Pooled A: 10 anti vs | Pooled B: 6 anti vs | S1063: 4 anti vs (under B's setting) |
|---|---|---|---|
| σ_ε 0.26 (as run) | 3.22, P 0.0019 | 2.27, P 0.0285 | 0.52, P 0.0021 |
| σ̂_ε | 2.50, P 0.0003 (σ 0.24) | 2.53, P 0.0438 (σ 0.27) | 0.58, P 0.0029 |
| σ_ε 95 % upper | 9.21, P 0.44 (σ 0.46) | 7.53, P 0.76 (σ 0.60) | 1.55, P 0.072 |
| f̂ | 3.22, P 0.0019 (f 0) | 4.22, P 0.250 (f 0.06) | 0.97, P 0.017 |
| f 95 % upper | 21.14, P 0.998 (f 0.36) | 21.06, P 1.0 (f 0.58) | 4.86, P 0.71 |

S1063 under A's settings: `propagated.A.*.S1063` (P 0.0009–0.40).

## Reading (narrowest statement)
- The weak-shear rows show no excess of `anti` over the σ_ε 0.26 null (9 vs 11.06 and 9 vs 9.14). Both simple
  calibration errors, a wider σ_ε and an isotropic contaminant fraction, are allowed by these rows over ranges that
  run from "no change" to "removes the pooled and the S1063 excess". The weak rows therefore neither confirm the
  null at strong shear nor explain the #406 excess.
- At the best-fit σ_ε, and at f = 0, the pooled and S1063 excesses remain (P < 0.05). At f̂ = 0.06 the pooled
  excess is gone and S1063 stays at P 0.017. Neither point estimate is preferred by the data over the null.
- Leverage is low because P(anti) at g ≈ 0.1 changes only from 0.17 to 0.22 between σ_ε 0.26 and 0.40
  (`anti_probability`, not stored), and only 61–75 rows are available.

## Failed approach
- The frozen rule took a point estimate as decisive without asking its interval to exclude the uncalibrated null.
  Future calibration plans should label "explains" only when the 95 % interval excludes the as-run value.

## Limitations
- Gaussian σ_ε per component (ASSUMPTION); no PSF, flexion or blend model. The isotropic model sets P(anti) = 1/3
  for contaminants, which ignores member alignment with the cluster.
- Weak-shear rows sit at larger radii than strong-shear rows: contamination, depth and model error may differ
  between the two.
- Small counts; the σ_ε interval of pool B reaches the grid edge (0.60).
