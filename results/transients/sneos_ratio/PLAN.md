# SN Eos 101.2 / 101.1 flux ratio vs the lens model: plan (frozen before any pixel is read)

Input: `results/transients/sneos/RESULT.md` (merged #540). In the pinned VENUS AstroLensPy v2.1 model 101.2 trails
101.1 by 0.32 d, so the SN flux ratio of the two images is the magnification ratio whatever the light curve:
R_model = |μ_101.2 / μ_101.1| = 28.700 / 24.698 = 1.162 (`sneos/summary.json`; the published `Mu_zs5.13` map gives
28.698 / 24.745 = 1.160). This is a lens-model flux-ratio test on a transient, which the static arcs cannot give.

## Hypothesis and ordinary explanations
- H0: R_obs agrees with R_model within the stated systematic.
- Ordinary reasons for a disagreement, tested or stated before any reading: host-arc light (removed by differencing
  epochs), PSF/registration errors (injection gate), microlensing by intracluster stars (expected at μ ≈ 25; no
  test here, it is a stated explanation), the macro-model's μ error near the critical curve (no ±1σ yet).
- A disagreement would be a lens-model residual to report, never an exotic signal; the strongest wording is
  "R_obs differs from this model's prediction".

## Data
`_i2d` SCI / ERR sections by S3 byte range (`lens_consistency.image_stamper`) from
`jw06882-o038_t041`, `jw09493-o001_t001`, `jw12299-o001_t001` NIRCam `clear-f277w / f356w / f444w` (9 mosaics; no
whole file stored). Positions: catalogued 101.1 / 101.2 (`results/lens_consistency/macs1931_inputs/m1931_images.ecsv`).

## Rules (ASSUMPTIONs labelled; per band, then combined)
1. **Difference.** For each band and each later epoch L ∈ {9493 o001, 12299 o001}: reproject the L stamp onto the
   6882 o038 stamp grid (`reproject_interp`), D = E_6882 − E_L (the SN is brightest in the earliest epoch; if not,
   the sign is taken from the data and reported). Registration check on the 6882 grid: cross-correlation peak offset
   of the two stamps (6″ × 6″, SN pixels within 0.4″ masked) ≤ 0.03″ (ASSUMPTION), else that pair is skipped.
2. **Photometry.** Circular apertures r = 0.20″ (ASSUMPTION; images are 0.9″ apart) at the catalogued positions
   on D; local background from the 0.35–0.6″ annulus excluding the other image's aperture; σ from ERR in quadrature.
   R_band,L = F_101.2 / F_101.1, its σ by first-order propagation.
3. **Gate D (detection).** Both images detected in D at S/N ≥ 5 (ASSUMPTION) in that band and pair; others are
   reported, not used.
4. **Gate I (injection).** Two PSF copies (an isolated star in the same 6882 band stamp, or the WebbPSF-free
   empirical fallback: the brighter SN image itself is not used) with ratio 1.16 and with ratio 0.70 injected at
   3 positions on the host-free side of the pair, same 0.9″ separation and orientation, at the measured SN
   flux; R recovered within 10 % (ASSUMPTION) in ≥ 5 / 6 trials per band, else that band is not used.
5. **Combination.** Inverse-variance mean of ln R over the bands and pairs that pass D and I; σ_ln from the scatter
   or the propagated error, whichever is larger.
6. **Decision.** "Differs" when |ln R_obs − ln 1.162| > 3 σ_ln AND > 0.30 (ASSUMPTION: a floor for microlensing and
   macro-μ error at μ ≈ 25; no published floor was read). "Agrees" when both bounds hold the other way
   (|Δ| ≤ max(3σ_ln, 0.30)). "Undetermined" when no band passes D and I.
7. Reported, not gated: the absolute SN fluxes per epoch; any point source within 3″ that appears in D and is not
   101.1 / 101.2 (contact sheet, then `/vet-candidate` if real).

Literature check before interpretation: whether Allingham et al. 2026 or the SN Eos discovery paper already report
this ratio (to be read and cited by arXiv number only after the measurement is frozen).
