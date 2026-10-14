# Observational parity test (TASKS Supporting 1, "an observational parity test")

Plan (2026-10-14, cloud routine):
- **Hypothesis:** for a catalogued multiple-image pair (A, B), the lens model's local map
  T = A_B⁻¹ A_A (A = lensing Jacobian ∂β/∂θ) predicts whether B is a mirror image of A (det T < 0) or not.
  A pair whose resolved morphology prefers the other handedness is a parity residual.
- **Ordinary explanations:** near-critical images (|μ| large, parity model-dependent), unresolved or symmetric
  sources (no handedness information), PSF, blends with members, wrong image identification.
- **Measurement:** fit an affine map B ≈ A ∘ M to PSF-matched cutouts in both handedness families
  (det M > 0 and det M < 0); Δχ² between the best fits is the observable. Report only pairs where a
  mirrored synthetic of A is classified correctly (per-pair validation).
- **Validation first:** synthetic pairs (a resolved galaxy mirrored / not mirrored, at the pair's S/N), then
  catalogued pairs far from critical curves, where the model parity is reliable.
- **Would change the plan:** a reliable pair (validated, |μ| < 10 for both images) preferring the
  non-predicted handedness at Δχ² above its own synthetic-null spread.

Acceptance test (2026-10-14, cloud routine, claimed): subtract member + ICL light around each stamp, fit only
inside the source's own segmentation footprint, and require ≥ 95 % `consistent` or `untestable` on catalogued
non-near-critical Abell 2744 pairs before any field run counts. ASSUMPTION: 95 % threshold.
