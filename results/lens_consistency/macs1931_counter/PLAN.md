# MACS J1931 counter-image route: plan (stage V frozen and run; stage 3 frozen)

Status: **stage V frozen and run 2026-10-24**; stage 3 frozen and run 2026-10-24, stopped at frame check F (18 / 51 < 20; `../macs1931_stage3/RESULT.md`); V1 not valid (depth undefined), V2 amendment passes gate C2 at the
threshold (28 / 35, arc 14 / 14; `../macs0417_ownref/RESULT.md`). Stage 3 frozen 2026-10-24 ("Stage 3: frozen rules" below), before any MACS J1931 predicted position, map
value or pixel was read for it (stage V read MACS0417 pixels and frame tables only).
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
   Required addition before freezing (code review of stage V): a footprint rule, i.e. a target whose position is off
   the section or on a masked pixel is classed `off_image` and is neither a control failure nor a flag (stage V
   counted the off-footprint control 15.1 as `absent`).

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

## Stage V2: amendment (frozen 2026-10-24 after the V1 run, before any V2 number)
V1 (rules 1–6 above) ran and its depth is undefined: the own catalogue has no row at Kron S/N 4–6, so
`depth_mag` is NaN and `m0417_arcref.classify` can return neither `absent` nor `undetectable` (every
non-matching control becomes `ambiguous`). V1's gate C2 (28 / 35) is therefore **not a valid gate reading**
(`../macs0417_ownref/RESULT.md`). Plan defect, not a builder result. V2 changes only these, rules 1–4 unchanged:
- **V2-a depth:** `depth` = AB magnitude of 5 × the median `kron_fluxerr` of the kept rows (a 5σ Kron-total limit of
  this catalogue; ASSUMPTION). Recorded with the minimum Kron S/N and the number of rows at S/N 4–6.
- **V2-b chain = the stage 3 chain:** `m0417_arcref.run` with the MACS0417 frame offset (as every MACS0417 run) and the
  containing-source rule **on**, then `m1423_canucsref.apply_local` on the background-subtracted section (nJy per
  pixel, masked pixels NaN) with source positions = the kept rows of the own catalogue, constants and seed unchanged
  (0.3″ apertures, 200 draws at 1.5–4″, 0.5″ exclusion, seed 20261022, 5σ / 10σ).
- **V2-c verdict:** `gate_c2_shared` on `class_local`; passes as in rule 6. Predicted MACS0417 classes are written, not
  read. Output `../macs0417_ownref_v2/`. Fails → the builder is not used for MACS J1931 (rule 6's fallback).

## Stage 3: frozen rules (2026-10-24, before any stage 3 run)
Frozen before any MACS J1931 counter-image prediction, own-catalogue row or F200W pixel was read. Script
`scripts/m1931_stage3.py`; tests `tests/test_m1931_stage3.py` (rule functions tested before the run). Numbers are
ASSUMPTIONs fixed here. Output `results/lens_consistency/macs1931_stage3/`.
1. **Model and images.** `macs1931_gatek.load_models` (md5-verified VENUS v2.1 maps, reading `P_pixels`, z_s 5.079
   normalisation) and `macs1931_gatek.load_images`. Systems: the 9 spec-z systems only (`z_spec` true); z = the
   catalogue z (ASSUMPTION: z fixed). Photo-z and mixed-z systems are out.
2. **Gate K-WLS repeated:** `a2744_uncover.imageplane` on the spec-z back-trace, 0.25″ grid, `m0417_images.gate_k`.
   An outcome different from `../macs1931_gatek/summary.json` "gate_K_numbers_wls" pass / fail → stop.
3. **Frame tables:** `a2744_uncover.wls_backtrace` (WLS source per system), then
   `lens_consistency.predict_counter_images` with its defaults, shapes = the JWST pipeline
   `jw06882-o038_t041_nircam_clear-f200w_cat.ecsv` (`lens_consistency.load_shapes`), depth
   `lens_consistency.depth_mag(shapes)`. Controls = the rows of class `observed`. No pipeline forced photometry
   (the bullet route reports it only; omitted here).
4. **Reference catalogue:** `m1931_deblend.build` (stage V rules 1–4 unchanged, 64 px box = 2.0″ at 0.0312″) on the
   SCI / ERR section of `jw06882-o038_t041_nircam_clear-f200w_i2d.fits` bounding every control and predicted position
   plus 15″, read by S3 byte range (not stored). Depth `m1931_deblend.kron_depth` (V2-a).
5. **Frame check F:** `bullet_extra.frame_check_f` of all 51 used images against the own catalogue (0.5″ match,
   |median dx| and |median dy| ≤ 0.1″), with the minimum number matched set to **20** (ASSUMPTION; 51 images here
   vs 135 in the Bullet). Fails → stop, no stage 3 reading. Frame offset (0, 0).
6. **Chain (= V2-b):** `m0417_arcref.run` with offset (0, 0) and the containing-source rule **on**, then
   `m1423_canucsref.apply_local` on the background-subtracted section (masked pixels NaN), source positions = the
   kept rows of the own catalogue, constants and seed unchanged.
7. **Footprint rule (new, required by the stage V review):** a control or target whose position falls outside the
   section, or whose nearest section pixel is masked (NaN), gets `class_local` = `off_image`. `off_image` is not in
   `m0417_arcref.TESTABLE`, so it counts neither for nor against gate C2, and it is never a flag. Function
   `m1931_stage3.footprint_class`.
8. **Gate C2:** `m0417_arcref.gate_c2_shared` on `class_local` after rule 7 (≥ 8 testable and ≥ 0.8; arc-referenced
   ≥ 3 and ≥ 0.8).
9. **Flag:** a predicted target (classes other than `observed` and the skipped classes) with `class_local` =
   `absent` while F, K-WLS and C2 pass. Contact sheet of controls and flags before RESULT.md; every flag goes to
   `/vet-candidate` in a later cycle, no interpretation here. Stage V's zero margin is stated with any flag.
10. Reported, not rules: the run without the containing-source rule, plain `gate_c2`, classes without rule 7.

Expected: F passes (the image list comes from the same VENUS mosaics); C2 passes or fails narrowly (stage V margin 0);
0–3 flags with ordinary readings near critical curves. Budget: one cycle; an unfinished run is rerun unchanged.
