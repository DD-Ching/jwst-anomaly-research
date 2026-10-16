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

# Step 2: clump-pair census of the 55 strong-shear background rows (2026-10-16)

Plan: `PLAN_census.md` (committed before the run, `1c58b28`). Code: `scripts/arcs_macs1423_census.py` calling
#273's `arcs_vet_census.py` unchanged at `43fbbbd`. **Deviation (no rule change):** the first launch stopped on an
import because #273's deblend module also loads `arcs_bgsub_shapes.py`; that file (also at `43fbbbd`) was added to
the plan's file list and the script's docstring, and the census was run once after that. Outputs: `census.json`, `census.ecsv` (derived;
injections simulated; P values model_prediction). Exploratory (the census idea is post hoc on S1063).

## Results
- **Deblend and gate:** cutout pixel scale 0.0312″. 45 / 55 rows pass the own-flux gate (≥ 50 %, so the
  image-level check of the plan passes); 13 of the 45 are `two-knot` (f̂ 0.289, 95 % upper 0.443); by screen class
  they are 9 aligned, 2 mixed, 2 anti (`census.ecsv`; the json key `two_knot_among_screen_aligned` = 11 counts every
  non-anti row, as in #273). All 55: 39 `single`, 13 `two-knot`, 3 `undecided`.
- **Screen anti rows:** 2197 `two-knot` (ΔBIC 151, separation 0.28″, knot S/N 8.4 / 6.7, gate passed); 2322
  `two-knot` (ΔBIC 150, 0.16″, S/N 10.5 / 17.5, gate passed); 4234 `single` (ΔBIC −2460, gate passed); 1818
  `single` but gate failed (pair injections recovered 20 %), so it gets the mixture probability.
- **Sample test (frozen outcome):** expected anti 2.608 under the clump-aware null (0.946 under shape noise only).
  **Primary k = 4: P(≥ 4) = 0.259 → `anti excess explained by clump orientation (exploratory)`.** Secondary k = 3
  (screen anti still anti after ICL subtraction): P 0.494. **Deviation:** the "still anti" part uses `class_bgsub`,
  which is classified at z_s = 2 only (as `posthoc.json`); it was not recomputed at each row's photo-z. Sensitivity: mixture with f̂ 0.301 / 0.533, with the 95 %
  upper f 0.526 / 0.739 (k = 4 / 3). Shape-noise-only check: 0.0126 / 0.0648.
- **Looked at (`census_anti.png`, F200W, 2″ × 2″, 1″-median-subtracted; `scripts/arcs_macs1423_census_sheet.py`):**
  2197 two compact knots at the end of a faint chain; 2322 two touching compact knots; 4234 one bright elongated
  galaxy with a faint tail; 1818 barely visible next to the BCG halo.

## Reading
- Under the S1063 clump-aware null, the MACS1423 anti count is not an excess. The two z_phot 3.42 rows are both
  two-knot sources whose catalogue axis carries little shear information; 4234 (single, m 24.2, at the g_min
  threshold) stays anti with P_sn 0.069 on its own; 1818 is not testable by the deblend here.
- With MACS0417 (P 0.195), two independent fields with JWST-era models now give no replication of the S1063 anti
  excess under their pre-registered or clump-aware nulls; MACS0647 gave no result. The S1063 excess itself has
  P 0.031 under the same clump-aware null (#273) and stays a single-field, exploratory result.
- Not tested here: the PA conventions of the deblend's `pa_sep` and the catalogue `pa_obs` were not cross-checked,
  so this step does not say whether the knot axes of 2197 and 2322 follow the prediction.
- Not a candidate.
