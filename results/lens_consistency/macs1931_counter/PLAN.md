# MACS J1931 counter-image route: proposal (not frozen)

Status: **stage V frozen 2026-10-24** (section "Stage V: frozen rules" below, committed before any run; claude/macs1931-stagev). Stage 3 stays a proposal. Nothing here has been run; no predicted position, map value or pixel was read for it.
Follows `../macs1931_gatek/RESULT.md` (gate K passes; spec-z systems 1, 5, 6, 7, 100, 101 have more predicted than
catalogued images) and the footprint step (`RESULT.md` here).

## Why a new reference is needed
The extra-image test (`../bullet_extra/PLAN.md` route) needs an arc-aware reference catalogue and passes gate C2 per
field first. MACS J1931 has none public: no DJA field, no CANUCS field (D-034 addendum aa). Pipeline isophotal
segments were rejected as a reference in cluster cores (D-034 addendum 2026-10-22: MACS0417 gate C2 10 / 27).
That addendum names the two alternatives: a deblended arc catalogue, or an arc-aperture sum on member-subtracted
images.

## Proposed route (to be frozen in its own commit before any run)
1. **Own deblended catalogue.** On an S3 byte-range cutout of the GO 6882 F200W level-3 mosaic covering the
   catalogued images plus a margin, subtract a smooth background that carries the cluster light: photutils
   `Background2D` with a box small enough to follow the bCG / member halos (ASSUMPTION, to be fixed: of order 1–2″).
   Then detect sources with `detect_sources` and split them with `deblend_sources`; measure Kron totals and shapes
   with `SourceCatalog`. Map them into the `m0417_arcref` catalogue form (`ra`, `dec`, `a`, `b`, `pa`, `flux`,
   `snr`, `ext`, `ell`, `scale`). photutils is the same library the JWST pipeline's source catalogue step uses; the
   change from the rejected pipeline segments is the background, which keeps cluster light out of the arc segments.
   Adding photutils as an explicit dependency needs `/reuse-check` and a D-NNN entry first.
2. **Validation on a known case (stage V), before MACS J1931.** Run the same builder on the MACS0417 GO pipeline
   F200W mosaic and pass it to the unchanged frozen MACS0417 route. Gate C2 there is 37 / 40 with CANUCS DR1
   (passes) and 10 / 27 with pipeline segments (fails). The builder is accepted for MACS J1931 only if gate C2
   (shared-pair rule) passes on MACS0417 with it; a failure ends the route without a MACS J1931 reading.
3. **MACS J1931 stage**, only after V passes: gate K-WLS repeated (`../macs1931_gatek`, spec-z systems only),
   frame tables with the WLS source, local noise on the same F200W mosaic (constants of `../bullet_extra/PLAN.md`
   rule 5), gate C2 (shared pairs), flags = predicted targets classed `absent` while C2 passes. Photo-z and mixed-z
   systems (8, 10, 11, 13, 15, 16, 17) stay out until their z is settled (`../macs1931_gatek/RESULT.md` reading).

## Hypothesis, ordinary explanations, what would change the plan
- Question (science, stage 3): do the VENUS maps predict a magnified, unmatched counter-image of a spec-z system
  where the F200W mosaic shows no source at the predicted flux?
- Ordinary explanations for any flag: model μ / position error near critical curves, low surface brightness or
  deblending, cluster-member light left by the background filter, a reference flux that includes unrelated light,
  wrong system z. Each goes to `/vet-candidate` before any wording beyond "flag". Exotic readings are not considered.
- Stage V fails → try the other alternative (arc-aperture sums on `member_sub` residuals) as a new plan, or wait
  for a public VENUS catalogue (the VENUS HLSP on archive.stsci.edu could not be checked: TLS fails through the
  cloud proxy; a local session can check it).

## Stage V: frozen rules (2026-10-24, before any run)
Script `scripts/m1931_deblend.py` (builder + stage V driver); tests `tests/test_m1931_deblend.py`. All numbers are
ASSUMPTIONs fixed here, not tuned on the result.
1. **Image.** MACS0417 GO 1208 F200W level-3 mosaic `jw01208-o006_t003_nircam_clear-f200w_i2d.fits` (the mosaic
   whose pipeline catalogue failed C2), SCI and ERR extensions read as one S3 byte-range section each (not stored).
   Section: the bounding box of every control and predicted position of `../macs0417_frame` plus 15″ on each side,
   clipped to the mosaic. Pixels that are non-finite or have ERR ≤ 0 are masked. Units MJy/sr → nJy / pixel (PIXAR_SR).
2. **Background** (the change from the pipeline segments): photutils `Background2D`, box 64 px (2.0″ at 0.0312″),
   filter 3 × 3, `SExtractorBackground`, sigma clip 3σ, coverage mask from rule 1. Subtracted from SCI.
3. **Detection and deblending:** convolve with a normalised 2D Gaussian kernel, FWHM 2.0 px (size 5); `detect_sources`
   at 3.0 × `background_rms` per pixel with npixels 25; `deblend_sources` nlevels 32, contrast 0.001, mode
   exponential.
4. **Measurement:** `SourceCatalog` on the background-subtracted SCI with the ERR map and the convolved image.
   Mapping to the `m0417_arcref` catalogue form: `ra`/`dec` = `sky_centroid`; `a`/`b` = `semimajor_sigma` /
   `semiminor_sigma` × the section's pixel scale; `pa` = the major-axis direction (`orientation`, CCW from +x) carried
   through the WCS to East of North, mod 180°; `flux` = `kron_flux` (Jy); `snr` = `kron_flux / kron_fluxerr`;
   `ell` = `ellipticity`; `ext` = `semimajor_sigma` > 1.5 px (F200W PSF σ ≈ 0.9 px); rows with non-finite or
   non-positive Kron flux or error are dropped. `depth` = `m0417_canucsref.depth_mag` (median AB of rows with S/N 4–6).
5. **Route:** `m0417_arcref.run` unchanged (frame tables, offset, rules 1–6 as frozen there), output
   `../macs0417_ownref/`. Gate C2 with `m0417_arcref.gate_c2_shared` (D-034 addendum i) on its control table.
6. **Stage V verdict:** passes when the shared-pair gate C2 passes (all ≥ 8 testable and ≥ 0.8; arc-referenced ≥ 3 and
   ≥ 0.8). The MACS0417 predicted-target classes are written but **not read** (that field's flags belong to
   `../macs0417_canucsref`). Fails → the builder is not used for MACS J1931; next is the arc-aperture alternative as
   a new plan. A contact sheet of the controls with the deblended segments is looked at before the result is written.
