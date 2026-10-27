# SN Eos leading pair: does any public archival epoch test it? Plan (frozen 2026-10-27 before the epoch table is evaluated)

Open item from `results/transients/sneos/RESULT.md`: "whether any other archival epoch covers that window was not
checked". The VENUS AstroLensPy v2.1 model predicts two SN Eos images that **lead** 101.1 by 1,321 d (pair a,
μ +14.48) and 1,238 d (pair b, μ −11.92) (`results/transients/sneos/summary.json`, model_prediction).

## Hypothesis and ordinary explanations
- H: the leading pair was bright at some public epoch; a missing point source there would be a lens-model residual.
- Ordinary: the model's Δt near a cluster-member critical curve is uncertain (no ±σ propagated); the SN light curve is
  unknown, so a non-detection where the predicted phase is outside 101.1's observed phases tests nothing.
- What would change the plan: an epoch passing rule T below gets forced photometry next cycle; if none passes, the
  leading-image line closes (outcome U, untestable on public data).

## Inputs
- One MAST `Observations.query_criteria` per predicted position (radius 0.0003°, dataproduct image, intent science,
  HST and JWST), plus 101.1 for its own epochs. Metadata only (`observed`); stored as `epochs.json`.

## Rule T (per epoch × filter × image; all three must hold)
1. **Band:** filter pivot wavelength > 1216 Å × (1 + z), z 5.133 (red of Lyα; the IGM absorbs bluer light).
2. **Phase:** the 101.1-equivalent date t′ = t + |Δt| lies inside the window W where 101.1 is known bright:
   W = [2024-03-11, 2025-09-01] (HST PID 17504 detection to VENUS 6882 o038; `sneos_hst`, `sneos_ratio`).
   ASSUMPTION: Δt uncertainty allowance ±20 % of |Δt|; the epoch passes if any t′ in that range is inside W.
3. **Reference:** 101.1 has a measurement in the same filter (so a predicted flux exists without a colour model).
- Sensitivity variant (reported, not deciding): W extended to 2026-09-28 (12299 o001).
- Outcome: P if ≥ 1 row passes, U if none.
