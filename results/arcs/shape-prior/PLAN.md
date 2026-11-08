# Plan: `arcs` anti null with an empirical blank-field shape prior (exploratory)

Written 2026-11-08 before any blank-field shape was read for this purpose. Follows
`../anti-bands/RESULT.md` "Next" (does intrinsic shape noise at the observed e explain an anti axis?).

## Hypothesis and ordinary explanations
- **Question:** the `arcs` anti null (`lens_consistency.anti_probability`) draws intrinsic complex ellipticity
  (a − b)/(a + b) from a Gaussian, σ_ε 0.26 per component (ASSUMPTION). The pooled anti excess (variant B 6 vs
  2.273, A 10 vs 3.222, S1063 4 vs 0.52; `../pooled-tally/`) and the 3 open single-component rows (WHL0137 2194, RXC
  J0600 581, MACS0308 2157; e 0.531–0.795) could be a Gaussian tail that is too thin: real galaxies (edge-on discs)
  populate high |ε| more often.
- **Ordinary explanation tested:** non-Gaussian intrinsic shape distribution. Not tested: PSF, flexion, blends
  (addendum 9/11), membership (addendum 10), model error (addendum 8).

## Method (frozen)
1. Prior: CEERS blank field `jw01345-o001_t021_nircam_clear-f200w_cat.ecsv` (jwst 3.0.0, sha256 pinned in
   `data/manifests/reference_sample.ecsv`), read with `lens_consistency.load_shapes`, cuts identical to the `arcs`
   runs except the ellipticity cut: `is_extended`, semimajor ≥ 2 px, S/N ≥ 10. |ε| = e / (2 − e) with e = 1 − b/a;
   phases uniform (ASSUMPTION: blank-field shear and PSF anisotropy negligible; same pipeline, filter and PSF).
2. `P_emp(anti | g)`: as `anti_probability` (same lensing map, e_obs ≥ 0.5, angle ≥ 60°, same g_min rounding),
   with |ε| resampled from the prior (|ε| < 0.95 kept, as the Gaussian code does).
3. **Gate (known case):** the weak-shear background rows of pool B (`../null-calibration/`: 61 rows, 9 anti) must
   give a two-sided Poisson P ≥ 0.05 under P_emp, and P_emp at g = 0 must equal 1/3 to 0.01 (isotropy check).
   Gate fails → outcome `tool fails; undetermined`.
4. **Outcome rule** (α 0.05, ASSUMPTION): pooled B strong-shear rows, λ_emp = Σ P_emp(g_min). P(≥ 6 | λ_emp) ≥ α →
   `pooled excess explained by the empirical shape prior`; < α → `excess persists under the empirical prior`.
   A and S1063 are reported, not used for the label.
5. Descriptive (not a label): for each of the 3 open rows, P_emp(anti | g_min) and the prior's fraction with
   |ε| ≥ that row's observed |ε|.

## What would change the plan
Explained → the anti line closes as shape noise (the Gaussian σ_ε 0.26 null is retired for `arcs`). Persists → the
rows remain open; the next test is the deciding models' local shear error.
