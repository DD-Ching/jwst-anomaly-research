# SN Eos HST F110W 2024 − 2012 difference-image ratio: plan (frozen before any pixel is read)

Input: `results/transients/sneos_ratio/LITERATURE.md` (#543). Coulter et al. 2026 Table S1 photometry gives a 101.2 /
101.1 HST F110W 2024 ratio of 1.805 (σ_ln 0.039, stat), while our JWST difference images give 1.382 (σ_ln 0.083, #541)
and the published models give R_model 1.066–1.174 (pinned VENUS map 1.162). This run measures the HST ratio from our
own difference images.

## Questions and ordinary explanations
- Q1 (lens model): does R_HST agree with R_model 1.162 under #541's frozen decision rule (`sneos_ratio.decide`)?
- Q2 (published photometry): is R_HST consistent (within 3 σ_ln, combined in quadrature with the comparison value's σ)
  with the Table S1 value 1.805, with our JWST value 1.382, both, or neither (`sneos_hst.compare`)?
- Ordinary reasons for a ratio far from the model: microlensing by intracluster stars (chromatic and time-variable;
  no test here, stated), host residuals from a 12-year template gap (differencing; empirical aperture noise), PSF /
  registration errors (registration check, injection gate), macro-μ error. Nothing here is an exotic test; the
  strongest wording is "R_HST differs from this model's prediction".
- What would change the plan: R_HST ≈ 1.8 (our pixels reproduce Table S1) makes a band- or time-dependent ratio the
  next question (microlensing first); R_HST ≈ 1.1–1.4 points to the Table S1 photometry, not the SN.

## Data (metadata from a MAST query this session; not stored)
- Science: PID 17504 WFC3/IR F110W pipeline `_drz.fits`, visits `ifaa53010, ifaa52010, ifaa50010, ifaa49010,
  ifaa51010, ifaa54010` (MJD 60380–60387).
- Template: PID 12456 (CLASH, 2012) WFC3/IR F110W `ibt2a2020_drz.fits` (1,509 s, the deeper of the two); the second
  template `ibt2a9030` is used only for a reported consistency check on the first science visit.
- SCI sections by anonymous S3 byte range (`s3://stpubdata/hst/public/<4>/<id>/<id>_drz.fits`); no file stored.
- Positions: catalogued 101.1 / 101.2 (`results/lens_consistency/macs1931_inputs/m1931_images.ecsv`).

## Rules (ASSUMPTIONs labelled; per visit, then combined)
1. **Units.** Both images are scaled to flux density with their own PHOTFLAM (e-/s × PHOTFLAM); stamps 6″ × 6″ on the
   science grid.
2. **Difference and registration.** Template resampled onto each science stamp grid through both WCSs (bilinear,
   `sneos_ratio.resample`); cross-correlation offset with the SN pixels (≤ 0.4″) masked (`sneos_ratio.register_offset`).
   If the offset is ≤ 0.3″ (ASSUMPTION) the template is shifted by it (`scipy.ndimage.shift`, order 1) and the offset
   re-measured; the visit is registered when the residual is ≤ 0.03″ (ASSUMPTION, as #541). D = science − template.
3. **Photometry.** Circular apertures r = 0.25″ (ASSUMPTION; ≈ 2 drizzled pixels, the JWST 0.20″ is under 2 WFC3/IR
   pixels) at the catalogued positions, local background from the 0.40–0.60″ annulus excluding the other image's
   aperture (`sneos_ratio.aperture_pair` with a constant error map). Aperture σ = the 1.4826 × MAD of aperture sums
   at 40 random positions in D at least 0.6″ from both SN images and inside the stamp (ASSUMPTION: empirical, includes
   host residuals), in quadrature with SN Poisson noise √(max(F, 0) / (EXPTIME · PHOTFLAM)) × PHOTFLAM.
4. **Gate D.** Both images at S/N ≥ 5 (ASSUMPTION) in D.
5. **Gate I.** Two copies of an empirical PSF (nearest Gaia DR3 star with significant parallax or proper motion,
   isolated, within the science mosaic; implementation of `sneos_ratio.gaia_psf_star` on the science WCS; if none is
   inside the mosaic, the gate is reported as not run and the visit is not used) with ratios 1.16 and 0.70 at 3
   positions 1.8″ from the pair midpoint, at the measured 101.1 flux; ≥ 5 / 6 recovered within 10 % (ASSUMPTION).
6. **Combination.** `sneos_ratio.combine` over visits passing registration, D and I.
7. **Decisions.** Q1 by `sneos_ratio.decide` (floor 0.30 in ln R, as #541). Q2 by `sneos_hst.compare`.
8. Reported, not gated: per-visit fluxes; the second-template check; any point source in D within 3″ other than the
   SN images (contact sheet).
