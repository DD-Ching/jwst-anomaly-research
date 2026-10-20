# W3 negative-mass microlens inside a macro-lensed arc (Chang-Refsdal type): plan

Frozen before any computation (2026-10-20, cloud routine). Everything here is `simulated` / `model_prediction`.

## Question
docs/exotic_lensing.md (W3 row and notes) predicts that a compact source behind an isolated negative-mass point
lens vanishes inside an umbra (β < 2 θ_E, zero lensed flux). Inside a cluster arc the microlens sits on a macro
model with convergence κ and shear γ (**assumption** in the doc: isolated-lens formulas). Does the umbra survive
macro magnification and shear, and how do depth and width change?

## Hypothesis and ordinary explanations
- Hypothesis H: with macro (κ, γ) the negative-mass lens still produces a region of the source plane with no image
  near the lens (a demagnified "hole" in the macro image), so a W3 event in an arc is a dimming of that macro image.
- Competing expectation: shear could fill the umbra (rays from the shear-stretched image plane reach β ≈ 0), so the
  "vanishing" becomes a partial dimming. Ordinary look-alikes for any dimming (positive-mass microlensing near a
  saddle, SN fading, persistence) are unchanged and are not tested here.

## Method
Inverse ray shooting with numpy (β = Aθ − s θ/|θ|², A = diag(1 − κ − γ, 1 − κ + γ), θ in units of θ_E(|M|),
s = +1 positive mass, −1 negative mass). Reuse: `exotic_sim.total_magnification` (isolated lens) for validation.
No new dependency (lenstronomy's point mass is parametrised by θ_E², not a sign; a 30-line numpy shooter is
simpler than a custom profile).

## Validation gates (all must pass before any macro result is read)
- G1: κ = γ = 0, both signs: shooter μ vs `total_magnification` on a radial profile, |Δμ| / μ < 3 % where μ < 20.
- G2: positive mass, κ = 0, γ = 0.3: caustic from det J = 0 (analytic critical curve mapped to the source plane)
  agrees with the ray-shooting μ ridge to within two source pixels.
- G3: mean μ over a region ≫ θ_E equals the macro 1 / |(1 − κ)² − γ²| within 2 %, for both signs.

## Measurements (cases: macro μ_macro ≈ 10 and 100, minimum (κ < 1 − γ) and saddle (κ > 1 − γ), both signs)
- area of the source-plane region with μ / μ_macro < 0.1 (in θ_E²) and its shape;
- light curve along straight tracks (u₀ = 0, 0.5, 1) for a point source and ρ = 0.01, 0.1;
- min μ / μ_macro along each track.

## What would change the plan
- If the hole area at μ_macro = 100 is < 10 % of the isolated value (π·2² θ_E²... i.e. umbra radius 2), W3
  amplitudes in arcs must be revised in docs/exotic_lensing.md and the dimming screens' injections.
- If it is ≥ 50 %, the isolated-lens injections stay a fair approximation (recorded as such).
