# W1 shear limits as a volume density over lens redshift: plan (frozen 2026-10-17, before any run)

**Question.** The twelve-field W1 shear null (D-053, streak veto D-097; `results/w1_streak_veto/RESULT.md`) gives
surface-density limits for lenses *at each cluster's redshift*. The only literature limit to compare with,
Takahashi & Asada (2013, arXiv:1303.1301: n < 10⁻⁴ h³ Mpc⁻³ for |M| > 10¹² M☉), is a comoving volume density. The
earlier comparison (docs/exotic_limits.md) spread their limit over 0 < z < 1 by hand. Here the efficiency is
measured as a function of the W1 lens redshift, so the null becomes a volume density limit directly.

**Hypothesis being limited.** W1 (negative-mass point lens, docs/exotic_lensing.md) lenses with |M| = 2 × 10¹²,
8 × 10¹², 2 × 10¹³ M☉ distributed uniformly in comoving volume (`hypothesis`). Ordinary explanations do not enter:
the screens are null; this is a limit, not a candidate search.

**What would change the plan.** If the limit is ≤ 10⁻⁴ h³ Mpc⁻³ at some mass, the JWST cluster fields add a new
constraint at that mass (owner judgement before any outside use). If it is weaker by ≥ 10×, the W1 shear line is
closed as a constraint and further W1 field additions are not worth their cost (T4).

**Run (frozen).** `scripts/inject_shear.py` unchanged except the new `--lens-z Z` option: the W1 lens sits at Z
instead of the cluster redshift; θ_E per row uses Z; only screen rows also behind Z (`lensable_mask`, z160 > Z + 0.1,
or no photo-z, which counts as z_s = 2) are painted. Everything else as `results/w1_streak_veto` (seed 50, 200
lenses per field and mass, streak masks regenerated with `w1_streak_mask.py`, unchanged code, D-053 screen
defaults). Grid: Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5}; masses 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉; the twelve headline fields of
`results/w1_streak_veto/combined_12fields_veto.json`. The seeds per mass repeat at every Z (common random numbers).

**Gate Z (input reproduction, computable).** Abell 2744 run with `--lens-z 0.3072` (its cluster redshift) must give
14 / 72 / 146 recoveries of 200 at the three masses, as stored (`results/w1_streak_veto/inject_summary_abell2744.json`).
Fail → stop, no grid run.

**Volume limit (computable; `scripts/w1_volume.py`).** For each mass,
V_eff = Σ_f A_f ∫ ε_f(z) (dV_c / dz dΩ) dz, with A_f the screened area (sr), ε_f(z) linear between grid points and
0 outside [0.1, 1.5] (conservative), dV_c/dz dΩ from Planck18. Zero detections: n₉₅ = −ln 0.05 / V_eff, quoted in
Mpc⁻³ and h³ Mpc⁻³ (h = 0.6766). Compare with 10⁻⁴ h³ Mpc⁻³ as a ratio (`derived`).

**ASSUMPTIONs (stated in advance).**
- Single plane: the W1 lens is added at the observed image positions on top of the cluster; the cluster's lensing
  of a foreground W1 lens's images (Z < z_cl) and the cluster's lensing of the W1 lens position (Z > z_cl) are
  neglected. The screen itself is unchanged (it measures the same rows with the same cluster shear removed).
- Rows without photo-z are lensable at z_s = 2 for every Z ≤ 1.5 (as in the cluster-z runs).
- The comparison with Takahashi & Asada treats their |M| > 10¹² M☉ bound as applying at each of our masses.
