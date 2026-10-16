# Clump-pair census of the 24 strong-shear background rows (Abell S1063 `arcs`), step 4: plan (2026-10-16)

Written before any measurement in this step. Parent: `RESULT.md` "Step 3 / Next". Script:
`scripts/arcs_vet_census.py` (reuses `arcs_vet_deblend.py`: `load`, `fit_models`, `deblend_class`, `gate`; and
`lens_consistency.anti_probability`, `arcs_vet_clumps.clump_anti_fraction`, `arcs_vet_noise.excess_p`). Data: GLIMPSE
`jw03293-o001_t001` F200W level-3 i2d, 8″ S3 byte-range cutouts, 1″ median filter subtracted. Labels: pixels
observed, fits and classes derived, injections simulations, anti probabilities model_prediction.

**This sample test is post hoc** (the clump idea came from looking at 1316 and 1557 after the screen), so its result
is labelled **exploratory**. No held-out field is claimed.

## Question
Does a null that lets clump pairs (or blends of two compact sources) point anywhere explain the S1063 anti count,
which the catalogue shape-noise null (σ_ε 0.26, `../arcs.json`) puts at P(≥ 4) = 0.0021?

## Ordinary explanation tested
Faint sources whose F200W axis is the separation of two compact knots (clumps of one galaxy, or two blended
sources). Their intrinsic "ellipticity" is maximal, so shear aligns them much less than the shape-noise null
assumes (`clumps_1557.json`: P_anti 0.07–0.14 at g 0.4–1.9).

## Rules (frozen)
- **Per row (all 24 rows of `../bgsub_shapes.ecsv`):** the step-3 deblend unchanged (F200W, 0.4″ fit radius, M1 one
  elliptical Gaussian vs M2 two round Gaussians, `deblend_class`: `two-knot` if ΔBIC ≥ 10, both knot S/N ≥ 3,
  separation ≥ 0.10″; `single` if ΔBIC ≤ −10; else `undecided`). Fit centre: the brightest 1″-filtered pixel
  within 3 px of the catalogue position (as step 3).
- **Gate per row (own fitted fluxes, as PLAN_1557.md):** 20 single injections at the fitted M1 shape and 20 pair
  injections at the fitted M2 knot fluxes and σ, separation = the observed M2 separation clipped to 0.10–0.30″
  (rounded to 0.01″); ≥ 80 % correct for each kind, else the row is `gate-failed`.
- **Per-row anti probability:**
  - `two-knot` with the gate passed: P_cp(g_min) = `clump_anti_fraction(g_min)` (isotropic point pairs, κ < 1).
  - `single` or `undecided` with the gate passed: P_sn(g_min) = `anti_probability(g_min, 0.26, 0.5, 60)` (the
    screen's own null).
  - `gate-failed`: the mixture f̂ P_cp + (1 − f̂) P_sn, with f̂ = (two-knot rows) / (gate-passed rows).
- **Primary test:** k = 4 (the screen's anti count, comparable with P = 0.0021), Poisson-binomial P(≥ 4).
  **Secondary:** k = 3 (bgsub anti: 1316, 1557, 2184) and k = 2 (2184 removed as noise-explicable at core flux,
  step 2). Sensitivity: the mixture model for every row with f̂ and with the upper end of its 95 % Clopper–Pearson
  interval.
- **Outcome:**
  - primary P ≥ 0.05 → `anti excess explained by clump orientation (exploratory)`: 1316 and 1557 are reclassified
    as ordinary-likely clump pairs (not candidates), and the S1063 `arcs` anti line closes pending a reopen input
    (z_spec showing either is a foreground or member source would not reopen it; a z_spec plus a single-component
    morphology at higher resolution would).
  - primary P < 0.05 → `excess remains under the clump-aware null`; next is a z_spec for 1316 and 1557 (MUSE).
- **Also reported, not used in the outcome:** the two-knot fraction among aligned rows and among anti rows, and the
  aligned fraction among two-knot rows against the clump model's expectation (Σ 1 − P_cp − P_mixed is not modelled;
  only P_anti is).

## What changes the plan
If most gates fail (fewer than 12 gate-passed rows) the census is `undetermined` and the outcome is not computed.
Large rows (semimajor > 0.4″ fit radius) are kept: the gate decides whether the deblend is valid for them.
