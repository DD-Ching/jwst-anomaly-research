# Plan: DJA v7.4 colours vs CEERS DR1.0 (D-010/D-022 validation target)

Written before any measurement (2026-11-06, cloud routine).

**Question.** Do the DJA v7.4 0.5″-diameter aperture colours used for the CEERS control (`dja05`, D-022) agree
with the PSF-matched CEERS DR1.0 catalog (Cox et al. arXiv:2510.08743) within the photometric errors?

**Hypothesis (H0).** For matched, isolated, unflagged sources the colour difference
Δc = c_DJA − c_CEERS has |median| ≤ 0.05 mag and a robust (NMAD) normalised scatter ≤ 1.5 per colour
(ASSUMPTION thresholds).

**Ordinary reasons to expect a failure.**
- DJA apertures are matched but not PSF-homogenised (D-010): in a fixed 0.5″ aperture the long-wavelength bands lose
  more light, so colours involving F277W–F444W are biased blue→red for compact sources. Predicted size: the
  aperture-flux fraction difference between F115W and F444W PSFs (~0.1–0.3 mag; model_prediction, not measured here).
- CEERS DR1.0 applies residual aperture corrections (2.5–24 %, magnitude dependent, the same factor in all bands)
  and Galactic extinction corrections (small in EGS, E(B−V) ≈ 0.01).
- Different detection images and deblending; mismatches in crowded regions.

**Method.** Stream both catalogs (DJA `ceers-full-grizli-v7.4-fix_phot.fits`, 250.5 MB, D-022 reason; CEERS
`ceers_cat_v1.0.fits.gz`, 204 MB compressed, the only public PSF-matched CEERS catalog), keep only the needed
columns in scratch, delete the files after use. Match within 0.2″ one-to-one; keep CEERS STAR_FLAG = 0,
BAD_REGION_FLAG = 0, S/N ≥ 10 in both bands of each colour. Compare five adjacent colours
(F115W−F150W, F150W−F200W, F200W−F277W, F277W−F356W, F356W−F444W) in the 0.5″ aperture of both catalogs,
split by DJA `flux_radius` (compact vs extended). Also restrict to the t021 footprint used by the control.

**What changes the plan.**
- H0 holds: D-022 colours are validated; close the M2 item.
- A size-dependent offset in the long-wavelength colours: the ranking colours carry a PSF term; record it as a
  known systematic (D-010 "Revisit if" holds) and propose PSF-matched colours or DJA's aperture corrections
  as the fix, since SMACS and CEERS share the same bias only if their size distributions match.
