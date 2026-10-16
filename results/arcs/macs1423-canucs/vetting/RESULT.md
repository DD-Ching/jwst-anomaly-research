# MACS1423 anti rows, vetting step 1: ICL-subtracted shapes and measurement noise (2026-10-16)

Plan: `PLAN.md` (committed before both runs). Code: #273's `arcs_bgsub_shapes.py` and `arcs_vet_noise.py` at
`43fbbbd`, unchanged, run next to this branch's `scripts/`. Labels: pixels `observed`; injected sources
`simulated`; predictions `model_prediction`; classes and P values `derived`. Thresholds are ASSUMPTIONs (the S1063
values).

## Results
- **ICL-subtracted shapes** (`../macs1423-canucs-maps/bgsub_shapes.json`): known-case check 40 / 43 pipeline-aligned
  rows stay aligned (median axis change 4.3°). Of the 4 screen anti rows, 2197, 2322 and 4234 stay anti (offsets
  70.2°, 77.6°, 71.0° from the z = 2 prediction); **1818 becomes `mixed` (54.0°)**. A fifth row, 1758, is anti
  against the z = 2 prediction only (see the deviation below).
- **Measurement-noise injections** (`../macs1423-canucs-maps/vetting/noise_injections.json`, 60 per row, 55 rows):
  gate passes (bright injections aligned in 100 % for every row). P_anti: 1818 0.017, 2197 0, 2322 0, 4234 0 →
  each **not noise-explicable**; 1758 also 0. Sample (frozen rule): 4 rows anti under the subtracted shapes
  (1758, 2197, 2322, 4234) vs 0.017 expected from measurement noise, P(≥ 4) 0 to the json's precision →
  **excess remains**.
- **Deviation / post hoc (`posthoc.json`):** the #273 scripts classify against the z_s = 2 prediction only. 1758
  (z_phot 6.38, g_min 10.0) is aligned at its own photo-z (5°) and anti only at z = 2 (85°), so the frozen sample
  count includes one row that the screen calls aligned and leaves out 1818. Post hoc, counting only the screen anti
  rows that stay anti after subtraction (3), against the screen's shape-noise expectation (0.949):
  P(≥ 3) = 0.071, above the 0.05 replication threshold. The two numbers answer different questions (measurement
  noise vs intrinsic shape noise); neither is a combined null.

## Reading
- Measurement noise of a truly aligned Gaussian source of the row's flux, size and ellipticity does not produce
  the anti axes of 2197, 2322 and 4234. ICL subtraction moves 1818 out of `anti`. After that, the screen-level
  replication is marginal (post hoc P 0.071 vs shape noise). Nothing else is tested yet.
- **Untested ordinary explanations, in order of the S1063 chain:** deblending fragments and clump pairs (2197 looks
  like one knot of a clumpy chain whose overall axis roughly follows the prediction; by eye, not measured);
  intrinsic shape beyond σ_ε 0.26 for bright rows (4234 is m 24.2 at the g_min threshold 0.20, where the model's
  shear is weakest); photo-z failure (2197 and 2322 both have z_phot 3.42 and lie 6″ apart); lens-model error (a
  preliminary model with 5 spectroscopic systems); second band.
- The noise injections are Gaussian profiles with the row's catalogue ellipticity; they do not test clumpy
  morphology.
- Not a candidate: unvetted rows stay on the draft PR (trial rule T2).
