# W1 shear screen: does PSF anisotropy measured from stars set Abell 2744's B-mode floor? (plan, 2026-10-17)

Claimed by cloud routine session 019HJb2RoNpPSjebceJ6ZqNp (#298). TASKS "JWST focus" 1 (W1); follows
`results/w1_floor/RESULT.md` "Next". Written before any screen was re-run (frozen at this commit).

## Hypothesis and ordinary alternatives
- **H (instrumental):** the catalogue shapes carry an anisotropic PSF that the screen's isotropic Gaussian
  deconvolution (`intrinsic_ellipticity`, σ from the narrowest minor axes) leaves in; a PSF pattern varying on the
  ~10″ aperture scale makes the 4.21 B-mode floor at (88, −42).
- **Alternatives:** shape noise above the rotation null's tails (D-053), light gradients of the bright galaxies
  around the peak (`results/w1_floor`), or real structure.
- **What would change the plan:** if the correction lowers the floor (rule R2), run injections with it; if a
  spatial map is infeasible (R1), the next step is a model PSF (STPSF per detector and dither), not catalogue stars.

## Inputs (as committed for W1)
F200W GLIMPSE catalogue `jw02561-o001_t003` (sha256 `9953c869…`), DJA v7.2 zout (`43662686…`), Bergamini+23
model. Stars: Gaia DR3 (VizieR I/355) G < 21 within 5′ of 3.586, −30.40 (`gaia_stars.py`).

## Star selection (computable; `w1_psf_floor.select_stars`)
A catalogue row is a PSF star when **all** hold: a Gaia DR3 source with G < 21 lies within 0.3″; catalogue
S/N > 50; catalogue mag ≥ 19.5 (brighter stars' isophotal moments are spike-dominated; ASSUMPTION).
PSF distortion χ★ = (a² − b²)/(a² + b²) e^{2i PA} per star (observed moments).

## Rules (computable, functions in `scripts/w1_psf_floor.py`, unit-tested before the run)
- **R1 map feasibility:** a star-based spatial PSF map resolves the aperture scale only if the median distance
  from a used source to its nearest PSF star is ≤ 10″ (the aperture radius) and ≥ 10 PSF stars exist (ASSUMPTIONs).
- **R2 adoption for injections:** the anisotropic correction (constant χ_PSF = the median of χ★ by component, and,
  if R1 passes, an inverse-distance map of the 5 nearest stars) is carried to injections only if the field's
  B-mode maximum b_max = max |S_×| falls to ≤ 3.79 (10 % below 4.21; ASSUMPTION). Otherwise not adopted.
- Correction: the PSF's second moments keep the screen's trace 2σ² and take χ_PSF's anisotropy:
  q_nn = σ²(1 + χ1), q_ee = σ²(1 − χ1), q_en = σ² χ2; subtracted from each row's observed moments.

## Outputs
`results/w1_psf/psf_abell2744.json` (stars, χ★, R1 inputs, b_max and S_× at (88, −42) with and without the
correction, rotation-null summary), `RESULT.md`. All `derived`; thresholds ASSUMPTIONs.
