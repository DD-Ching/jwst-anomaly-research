# W1 lenses injected at the `radial` centres through the stack: plan (2026-10-17, frozen before any run)

**Question.** How massive must W1 lenses be, and at what fraction of the `radial` centres, for the W1 shear stack
(`results/w1_stack_radial/RESULT.md`, vet step 1 sample) to exceed what was observed? This turns the exploratory
stack (mean S 0.51, p 0.027) into a limit. Hypothesis tested: some `radial` centres hold W1 lenses. Ordinary
explanation of the observed excess: chance plus cluster-model error (RESULT.md addendum). A limit that excludes
f = 1 at a mass the radial screen also reaches would say the centres are not all W1 lenses of that mass.

**Sample (as the quoted vet step 1).** Abell S1063 and MACS0417 (the fields with valid centres), primary models,
radial arcs dropped from the shear sample and from the lensable background, shear-screen defaults
(`inject_shear.SHEAR_DEFAULTS`). Centres: every `radial` peak with a valid S on the real data.

**Injection.** For each mass M in `inject_radial.MASSES` (2e10, 2e11, 2e12, 8e12, 2e13 M☉) and each valid centre c,
one lens painted at c with `inject_radial.paint_lens` and the image shapes of `inject_shear.ShearInjector` (same
chain as D-050/D-053: blends, detection floor, spike and correction cuts), and S evaluated at c only:
S_inj(M, c). ASSUMPTION (conservative): a lens changes only its own centre's S; overlapping apertures of
neighbouring centres keep their real S.

**Stack for a fraction f** ∈ {0.1, 0.25, 0.5, 1.0}: n = max(1, round(f N)) centres drawn at random without replacement
(1000 draws, seed 316; one draw at f = 1) take S_inj, the others keep S_obs; stacked mean over all N centres.

**Rules (computable; `w1_stack_inject.excluded`, `w1_stack_inject.efficiency`; thresholds ASSUMPTIONs).**
- (M, f) is **excluded** when ≥ 95 % of the draws give a stacked mean > the observed stacked mean (per (M, f) cell,
  pooled stack of both fields).
- **Efficiency** of (M, f): fraction of draws whose stacked mean exceeds the real null's 95th percentile (rotation
  null, 2000 draws, as in the stack).
- **Gate Z (tool):** an injection with θ_E = 0 must reproduce S_obs at every valid centre to 1e-9; otherwise no
  limit.
- **Gate R (response):** at M = 2e13, f = 1, the stacked mean must exceed the null 95th percentile; otherwise the
  stack does not respond to W1 lenses at these centres and the result is "no reach", not a limit.

**What changes the plan.** Gate failure → record and stop. Excluded cells at f = 1 only for M ≥ 8e12 → the stack adds
little to the radial screen's own limits. Exclusion down to 2e11 → the stack is a useful independent W1 constraint at
the centres.
