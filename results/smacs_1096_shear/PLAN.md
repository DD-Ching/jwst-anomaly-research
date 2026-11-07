# SMACS 0723 `f200w_1096`: model shear and magnification at the arc (plan, frozen 2026-11-07 before any evaluation)

Open question from `docs/candidates/jw02736-o001_t001_nircam_f200w_1096.md`: the source's major axis lies 19° from
tangential at r = 46.8″, but the Mahler+2022 ICLv2 shear and magnification there were never evaluated.

- **Hypothesis H-1096:** the orientation of 1096 is inconsistent with the ICLv2 prediction at its redshift.
- **Ordinary explanations:** weak cluster shear plus intrinsic shape (random orientation when g ≪ σ_ε).
- **Inputs:** ICLv2 (`lensmodel.SMACS0723_MAHLER22_ICLV2`, pinned); position and PA from the jwst 3.0.0 F200W
  catalog counterpart (`results/smacs_jwst3`), NED z = 1.3618 (observed). Also z_s = 1, 2, 4 for robustness.
- **Measure (model_prediction):** κ, γ, reduced shear g, signed μ, predicted stretch PA; offset Δ to the observed PA.
- **Frozen rule (ASSUMPTIONS):** with the D-024 shape-noise model (σ_ε = 0.26 per component, Gaussian, e ≥ 0.5),
  compute P(offset ≥ Δ). `inconsistent` if P < 0.01; `consistent` otherwise. If g < 0.05 the orientation test has
  no power: `no_power`.
- **What changes the plan:** `inconsistent` → vet the shape in other bands and the MCMC samples before any claim.
