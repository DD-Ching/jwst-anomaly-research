# RELICS Lenstool v1 amplitude: map ratio over the high-magnification region (H_K re-test; frozen 2026-10-25)

Frozen before any deflection map of this run was read. Follows D-034 addendum bj
(`../relics_lenstool_batch/RESULT.md`, "Failed approach": a peak-|α| region covers 84–99 % of the overlap).
Hypothesis, ordinary alternatives, clusters, model pairs and "what changes the plan" are unchanged from
`../relics_lenstool_batch/PLAN.md`; only the region S (step 2) changes. The previous run's test-cluster numbers
were seen (they were computed with the failed region); this plan does not depend on them.

## Statistic (frozen)
1. Same 1″ sky grid over the overlap of the two map footprints; bilinear |α| in each map's own WCS.
2. **Region S (new):** |μ_ref(z_s = 2)| ≥ 10 (ASSUMPTION: the neighbourhood of the critical curves, where
   multiple images form). μ_ref is computed from the reference model's own deflection maps: α scaled by
   D_LS/D_S(z_l, 2) (flat ΛCDM with each model's registered H0, Ωm), Jacobian by central differences on the
   map's native pixel grid, μ = 1 / det(I − ∂α/∂θ); then sampled (nearest pixel) onto the sky grid. Cells within
   3″ of the overlap edge are excluded (edge differences).
3. R = median over S of |α_LT| / |α_ref| (both at D_LS/D_S = 1). Also reported: ratio of means, 16–84 %, |S| in
   arcsec², S as a fraction of the overlap, and R for |μ_ref| ≥ 5 and ≥ 20 (sensitivity only, not used by a gate).
4. Classes unchanged (ASSUMPTION): **offset** R ≥ 1.25; **near 1** 0.85 ≤ R ≤ 1.15; else **intermediate**.
5. A cluster with |S| < 100 arcsec² is **unmeasured** (no class).

## Gates and outcomes (unchanged rules)
- **Gate V:** each of MACS J0553 (f 0.665), Abell 3192 (0.692), RXC J2211 (1.017), RX J0152 (1.045) falls in its
  known class (offset when f ≤ 0.8, else near 1), and |1/R − f| ≤ 0.12. Also required: S covers ≤ 50 % of the
  overlap in every validation cluster (the region must be a sub-region, the bj failure mode). If V fails, H_K is
  **undetermined** and the test clusters are reported, not read.
- **If V passes:** H_K **supported** when RXC J0600 and MACS J0417 are both offset and WHL0137 and SPT0615 both
  near 1; **rejected** when any batch map is near 1 or any other-maker map is offset; otherwise undetermined.
  MACS J0417's reference is JWST-era (CANUCS), read with that label.
- **What changes the plan:** supported → "check 2017-12 batch maps by default" enters the addendum-bd check; rejected
  → the per-cluster check stays the only rule and the batch line closes; undetermined → the batch line is paused
  (no third region variant).

## Tool check (added 2026-10-25 before any map was read)
- **Gate T (per reference model):** κ derived as ½ tr(∂α/∂θ) from the same differences must match the model's
  published κ map (D_LS/D_S = 1, same native grid, edges excluded): median |κ_derived − κ_pub| ≤ 0.05 over pixels
  with κ_pub ≥ 0.3 (ASSUMPTION). A reference failing T leaves its cluster **unmeasured** (an axis or unit convention
  the Jacobian does not handle); a validation cluster unmeasured fails gate V.
