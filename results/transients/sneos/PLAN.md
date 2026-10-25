# SN Eos (MACS J1931 system 101) predicted images and delays: plan (frozen before any map value is read)

Input: `results/transients/epoch_inventory/RESULT.md` "Next". Model: VENUS AstroLensPy v2.1 (Zenodo
10.5281/zenodo.18646525), already unit-checked and gate-K-tested (`results/lens_consistency/macs1931_gatek/`):
map values in pixels (× 0.04″), α normalised at z_s 5.079. This stage makes predictions only; forced photometry in
6882 o038 / 9493 o001 / 12299 o001 is the next stage with its own plan.

## Hypothesis and ordinary explanations
- H: the pinned model predicts SN Eos images beyond 101.1 / 101.2 whose delays put them inside, or close to, the
  three epochs (2025-09-01, 2025-10-08, 2026-09-28). Then a forced-photometry test can confirm or falsify the model
  on a transient, which static-arc tests cannot.
- Ordinary outcomes: the extra images lead or trail by years (only an "absent where predicted bright" test of the
  host is possible), or the model's delay map is not reproducible from the public files (then no delay test).
- What changes the plan: any extra image with |Δt| ≤ 400 d relative to 101.1 and |μ| ≥ 2 (ASSUMPTION) is a
  photometry target next cycle; if none, the SN Eos line drops to "watch for new epochs".

## Rules (ASSUMPTIONs are labelled)
1. **Gate P (potential reading).** `Lens_pot` must be the potential of `Alpha_x/y`: fit one scalar s so that
   s·∇ψ (centred finite differences, per pixel) matches α (raw map units) over the pixels within 40″ of 101.1 with
   |α| ≥ 1 map unit. P passes when the median |s·∇ψ − α| / |α| ≤ 0.05 (ASSUMPTION). s is reported, and the delay
   uses ψ in arcsec² = s·ψ_map·(0.04″)² at D_LS/D_S(5.079), rescaled by (D_LS/D_S)(5.133)/(D_LS/D_S)(5.079).
2. **Source.** Mean of the back-traced 101.1 and 101.2 at z 5.133 (catalogue z; the gate-K convention).
3. **Images.** `lensmodel.find_images` on a 0.25″ grid over `lens_consistency.grid_half_width`. Catalogued
   images are matched to predicted ones one-to-one by nearest distance ≤ 0.5″ (ASSUMPTION; 101.1–101.2 are 0.9″
   apart). Gate S: both catalogued images matched; otherwise no delays are reported as predictions.
4. **Delays.** Δt_i = (1+z_l)/c · D_l D_s / D_ls · [|θ_i − β|²/2 − ψ(θ_i)] − same for the image matched to 101.1,
   flat ΛCDM H0 70, Ωm 0.3 (the gatek cosmology); ψ bilinearly interpolated. The formula is validated in a unit
   test on a singular isothermal sphere (analytic Δt) before use.
5. **Gate M (reported magnifications).** |μ| from `find_images` at the matched 101.1 / 101.2 vs `Mu_zs5.13` at the
   same pixels; ratio in [0.67, 1.5] (ASSUMPTION) for both, else μ is reported but not used to choose targets.
6. Reported, not gated: predicted count vs the five SN images Allingham et al. 2026 (arXiv:2602.14074) state; the
   ±1σ maps are not used in this stage (no uncertainty on Δt yet: a stated limit).

Data: map zips streamed to the session scratchpad, md5-checked against the record, deleted after the run.
