# W1 shear screen: aperture-coverage cut on grid centres (plan, fixed before any run; 2026-10-17)

**Why.** TASKS.md "JWST focus" 1 (W1) lists "drop edge apertures with < 50 % coverage" as one way to lower the
B-mode floor. Injection recovery requires S > the field's largest |S_×| (D-053), so the floor sets the limits. In the
committed five-field runs (`results/w1_shear_s1063/rerun_4fields/`, `inject_summary_abells1063.json`) four of the
five S_max centres lie 75–107″ from the model centre, near the catalogue edge (`s_max_xy`), which suggests that edge
apertures with few sources set part of both tails.

**Hypothesis (method, not physics).** Apertures that are partly outside the catalogue footprint have few, one-sided
sources, so their S and S_× have heavier-than-rotation tails. Removing them lowers b_max and raises recovery more
than it loses screened area.

**Ordinary explanations of the B-mode floor that this does not address:** PSF anisotropy, blends, streaks.

**Method.**
- Coverage of a grid centre (`exotic_screens.aperture_coverage`, derived): the fraction of 1″ lattice points with
  r_min ≤ r ≤ R (R = 10″ aperture, r_min = 1″) around the centre that lie in the footprint, i.e. have any
  catalogue row (all rows, before cuts) within 4″ (`FOOTPRINT_RADIUS`, the definition the limits already use).
- `--min-coverage c` (new option, default 0 = unchanged behaviour): centres with coverage < c are invalid
  (NaN) in S, S_×, in every rotation-null draw and in the injection trials. c = 0.5 (ASSUMPTION, the TASKS figure).
- Everything else unchanged: same catalogues, photo-z, Gaia stars, seeds (screen 1, injections 50, 200 per mass),
  `inject_shear.py` defaults, footprint area for the limits.
- Runs: diagnostic first (coverage at each field's S_max and |S_×| max centre, c = 0); then the five fields at
  c = 0.5: real screen + injections + combination.
- Check: c = 0 must reproduce the committed recoveries exactly (all five fields), else stop and report.

**Frozen rules (computable; `scripts/w1_edge_cut.py decide`).**
1. Real fields at c = 0.5: the D-053 rule per field (p_rot < 0.05 and S_max > b_max). A pass goes to
   `/vet-candidate`; nothing is interpreted before that.
2. Adopt c = 0.5 as the shear default only if the combined 95 % upper limit (`upper_limit_deg2`) is lower at
   **each** of 2 × 10¹², 8 × 10¹² and 2 × 10¹³ M☉ than at c = 0. Otherwise c = 0 stays the default and the cut is
   recorded as a failed approach.
3. Either way, the per-field b_max change is reported.

**What changes the plan.** A real pass at c = 0.5 (vet before anything else). A c = 0 mismatch (stop).
