# g09725 short-wavelength difference, 2561 o002 → o006 (frozen plan, 2026-10-26)

Written and pushed before any SW pixel was read. Open item of `docs/candidates/a2744_2561_flank_g09725.md`
("SW-band difference"). g09725 (RA 3.61170, Dec −30.45765) brightened by 1.0–1.5 mag in F277W / F356W / F444W
(6–9σ) between 2023-07-31 and 2024-07-30 (`../hits.ecsv`, group 36). Off-core: no lensing question; this is
vetting of an ordinary-transient reading.

## Question, hypotheses, what changes the plan
- **Question:** did the source also change in F115W / F150W / F200W, and what is the spectral shape (fν per band)
  of the difference flux Δf over 1.15–4.4 µm?
- Ordinary readings (both `hypothesis`): a supernova in a faint host (host photo-z 0.76 / 0.92 or 4.35), AGN
  variability. Instrumental readings already tested in the LW bands (persistence, dithers, spikes, edges).
- **Changes the plan:** SW brightening at ≥ 5σ → the change spans 1.15–4.4 µm and the Δf SED is recorded for the
  owner's typing (no type is assigned here); SW not detected → a red-only change, recorded as a limit on
  Δf_SW / Δf_F356W. Neither outcome alone separates SN from AGN; the verdict stays `inconclusive` unless an SW
  band shows an instrumental cause (dipole, edge, NaN).

## Data (observed)
`_i2d` level-3 mosaics `jw02561-o002_t001` and `jw02561-o006_t007`, `nircam_clear-<band>`, bands F115W, F150W,
F200W (test) and F277W, F356W, F444W (control), from `s3://stpubdata/jwst/public/jw02561/L3/t/`. Only 6″ × 6″
(o002) and 8″ × 8″ (o006) sections are read by S3 byte range; nothing is stored.

## Rules (frozen; functions in `scripts/g09725_sw.py`, tested in `tests/test_g09725_sw.py`)
1. **Stamp and registration.** o006 is resampled bilinearly onto the o002 section's pixel grid. The residual
   offset (cross-correlation peak, pixels within 0.4″ of g09725 masked) must be ≤ 0.03″ (ASSUMPTION, as in
   `sneos_ratio`); otherwise the band is `unregistered`.
2. **Photometry.** D = o006 − o002 (MJy sr⁻¹). Aperture radius 0.20″ (ASSUMPTION) at the catalogue position of
   g09725, background = median of D in a 0.35–0.60″ annulus; Δf in µJy = aperture sum × PIXAR_SR × 10¹².
3. **Noise.** σ = max(σ_prop, σ_emp): σ_prop from the ERR arrays (quadrature, both epochs); σ_emp = 1.4826 ×
   MAD of the same aperture measurement at 16 positions on a 1.6″ ring around g09725 (a robust empirical scale).
4. **Per band:** `brighter` if Δf / σ ≥ 5, `fainter` if ≤ −5, else `none` with a 5σ limit 5σ. Applied per band.
5. **Control gate C:** ≥ 2 of the 3 LW bands registered and `brighter`. If C fails, the outcome is U (the stamp
   method does not reproduce the known LW change).
6. **Outcome** (C passed): `SW_detected` if ≥ 1 SW band is `brighter`; `SW_none` if all three SW bands are
   registered and `none`; `SW_fainter` if any SW band is `fainter` and none `brighter`; otherwise U.
7. The Δf of all six bands (µJy, with σ) is written to `summary.json`; ratios to F356W are derived from it.

## Contact sheet
o002, o006 (resampled) and D per band, 3″ half-width; looked at before any conclusion is written.
