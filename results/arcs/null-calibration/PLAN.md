# `arcs`: shape-noise null calibrated on weak-shear rows (plan, frozen before the weak-row orientations are read)

Handoff of #406 (`results/arcs/pooled-tally/RESULT.md`): the σ_ε 0.26 shape-noise null under-predicts `anti` in the
strong-shear background rows of six independent fields by ×2.6–3.1; the cause is untested.

## What was seen before freezing
- The pooled strong-shear counts and their expectations (#406); S1063's 4 anti.
- For the weak-shear rows below: only their numbers and `g_min` percentiles (75 rows in the six fields, g_min median
  0.11–0.15; 14 rows in S1063), and `anti_probability` at g 0–0.3 for σ_ε 0.20–0.40. **Not** their orientation
  classes.

## Hypotheses (labels per methodology)
- H0 (assumption under test): the null with Gaussian σ_ε = 0.26 per component is calibrated for these catalogues.
- Hσ (ordinary): the intrinsic or measured ellipticity scatter is wider (σ_ε ≳ 0.27–0.32 per #406), e.g. PSF,
  noise or non-Gaussian tails.
- Hiso (ordinary): a fraction f of "background" rows carry orientations unrelated to the model shear (photo-z
  contamination, blends, clump pairs, model error); for them P(anti) = 1/3 (uniform axis angle).
- Exotic readings are not tested here.

## Sample and data
- W: rows of each field's deciding-run `arcs.ecsv` (paths as in `results/arcs/pooled-tally/summary.json`) with
  `strong_shear` False and `population` "background" (same cuts, so the same e ≥ 0.5, S/N, size and radius
  selection as the strong-shear rows). No new data.
- Calibration pool: the six #406 fields (variant A) and A without MACS1423 (variant B). S1063 is not in the pool
  (it is the target); its weak rows are reported separately.

## Measurements (all derived)
1. K_W (observed anti) vs λ_W(0.26) = Σ anti_probability(round(g_min, 2), σ_ε, 0.5, 60°); two-sided Poisson P.
2. σ̂_ε: maximum of the Bernoulli likelihood Π p_i^a_i (1 − p_i)^(1 − a_i) over σ_ε ∈ [0.15, 0.60] (grid 0.01), with
   the 95 % profile interval (Δ(−2 ln L) ≤ 3.84). ASSUMPTION: Gaussian per component.
3. f̂ at σ_ε 0.26: p_i' = (1 − f) p_i + f/3; MLE over f ∈ [0, 1] (grid 0.01) and 95 % interval as above.
4. Propagation to the strong-shear rows: pooled anti (A, B) and S1063, P(≥ observed) under σ̂_ε, its 95 % upper
   bound, f̂ and its 95 % upper bound (Poisson λ = Σ p_i', as in #406).

## Frozen outcome rule (α = 0.05, ASSUMPTION)
- **"calibration explains the pooled excess"**: pooled B has P ≥ 0.05 at σ̂_ε or at f̂.
- **"null calibrated; excess stands"**: otherwise, if pooled B stays P < 0.05 at the 95 % upper bounds of both σ_ε
  and f.
- **"weak rows lack power"**: otherwise.
S1063 is reported under each calibration; it gets no separate label.

## What would change the plan
- Explains → the #406 excess is a null-calibration effect; S1063 is re-read under the calibrated null, and the
  `arcs` default σ_ε (or an f term) is revisited in a DECISIONS entry.
- Excess stands → weak rows rule out the two simplest calibration errors; the strong-shear excess needs a
  g-dependent cause (model error near critical curves, clump pairs), and S1063 keeps its standing.
- Lacks power → no further `arcs` null calibration from these tables; it needs deeper catalogues.
