# Abell 2744: UNCOVER-reference extra-image test, stage 2b (plan, frozen)

Frozen 2026-10-23 before any stage-2b table exists and before any DR3 SUPER flux, F200W pixel or prediction of this
stage was read. Follows `PLAN_1b.md` (gate K-WLS passed) and `PLAN_2a.md` (`theta_J2000` = H_a). Same route as
`../abell370_buffalo_extra/PLAN.md` (stages 2–3), with the deviations listed below. Script:
`scripts/a2744_uncoverref.py --stage 2b`.

**Disclosed before freezing.** The plan author has read the image-class counts of the earlier D-041 run on these maps
(`../abell2744-uncover2/images_predicted_uncover2.ecsv`: WLS source, all 181 images, F277W with a 1″ halo filter;
137 observed, 17 candidate, 5 missing, 10 demagnified, 6 no_flux_ref, 6 outside, 1 faint), not its forced classes in
detail. Footprint check (positions only, no flux): GO 2561 `o001_t003` F200W covers all 181 catalogued images and all
182 of those earlier predicted positions with finite non-zero pixels; `o002_t001` and `o006_t007` cover none,
`o003_t006` 20 / 181. So one tile suffices; the multi-tile extension named in the stage-2a handoff is not needed.

## Hypothesis, ordinary explanations, what would change the plan
- Question 1 (method): does the UNCOVER DR3 SUPER reference with the pipeline-mosaic noise pass gate C2 on Abell 2744?
- Question 2 (field): are the magnified, unmatched images the UNCOVER v2.0 maps predict (WLS source) present in F200W
  at the predicted flux?
- Ordinary explanations of an `absent`: model μ or position error near critical curves (UNCOVER's own rms 0.84″),
  low surface brightness, deblending, cluster-member or ICL light (the pipeline mosaic is not bCG-subtracted), a
  reference image whose catalogue flux includes unrelated light. Exotic readings are not considered.
- C2 fails → the DR3 SUPER reference with this noise rule is not validated; no flag is read. C2 passes and a flag →
  contact sheet, then `/vet-candidate`; no interpretation in this stage.

## Inputs
- model_prediction: `abell2744-uncover2-maps` (UNCOVER v2.0 deflection maps; sha256 in `summary.json`).
- observed: the 172 spectroscopic images of `uncover2_images.ecsv` (the set gate K-WLS validated); z from the
  stage-1 `imageplane_abell2744_abell2744-uncover2-maps.ecsv` `z_used`. The 9 photometric-z images are left out.
- observed: JWST pipeline level-3 GO 2561 `jw02561-o001_t003_nircam_clear-f200w_i2d.fits` (stamps and the SCI
  extension by S3 range reads) and its `_cat.ecsv` (shapes); not stored in git.
- observed: UNCOVER DR3 `UNCOVER_DR3_LW_SUPER_catalog.fits` (sha256 pinned in stage 2a), streamed, deleted after use.

## Frozen rules
1. **Frame tables:** `m0417_extra.run_field` unchanged except that every traced row's source is its family's WLS
   source (`a2744_uncover.wls_backtrace`, the estimator gate K-WLS validated), passed through a new `bt_hook`
   argument (default unchanged). Its pipeline forced test and gate C are reported, not used.
2. **Reference catalogue** (`a2744_uncoverref.to_cat_super`, unit-tested before the run): `id` → label; `ra`, `dec`;
   `a_image`, `b_image` × 0.04″ (PLAN_2a rule 4); PA east of north = deg(`theta_J2000`) + 90 mod 180 (H_a, north-up
   grid measured in 2a); flux = `f_f200w` × 10 nJy; S/N = `f_f200w` / `e_f200w`; ext = not `flag_star`; rows with
   non-finite `f_f200w` dropped. Depth by `m0417_canucsref.depth_mag`.
3. **Stage 3:** `m0417_arcref.run` with this catalogue, frame offset (0, 0) and the containing-source rule on (as
   Abell 370's frozen run), rules otherwise unchanged.
4. **Local noise** (`m1423_canucsref` rule 2, constants unchanged: 0.3″ apertures, 200 draws at 1.5–4″, 0.5″
   source exclusion around DR3 SUPER positions, seed 20261022, 5 σ / 10 σ) on the pipeline F200W SCI converted to
   nJy per pixel (× `PIXAR_SR` × 1e15; exact-zero pixels, i.e. no coverage, set to NaN so apertures touching them are
   redrawn). Deviation: the pixel scale is the mosaic's own (from its WCS; about 0.031″,
   not 0.04″), passed as a new `pix` argument (default unchanged), so all radii stay in arcsec. The mosaic is not
   bCG-subtracted: near members σ_loc is larger, which can only move a test to `undetectable_local` or `ambiguous`
   (conservative against false flags; costs completeness near bright members).
5. **Gate C2** on the classes after rule 4 with the shared-pair counting rule (`m0417_arcref.gate_c2_shared`, D-034
   addendum i), thresholds unchanged (≥ 8 testable, ≥ 80 % recovered / confused overall and on ≥ 3 arc-referenced).
6. **Flag:** a target (not `observed` / skipped classes) with `class_local` = `absent` while C2 passes.
7. Reported, not rules: C2 and flags without the local-noise rule, without the containing-source rule, and with plain
   `gate_c2`.

## Data volume
The SCI extension (12020 × 8940 float32, ≈ 430 MB) is read once into memory by a range stream for the random
apertures; ~170 positions × 8″ cutouts with 64 kB read-ahead would cost several GB. Nothing is written to disk
beyond the 9 MB pipeline catalogue and the DR3 catalogue (deleted after the run).

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged.
