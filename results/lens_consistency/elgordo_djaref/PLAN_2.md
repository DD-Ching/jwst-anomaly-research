# El Gordo reference-catalogue extra-image test: stage 2 (plan, frozen)

Frozen 2026-10-23 before any DJA flux row, F200W pixel or stage-2 prediction was read. Follows `PLAN_1.md` (gate K
passed; rule 2 failed as written). Same route as Abell 2744 stage 2b
(`../abell2744_uncoverref/PLAN_2b.md`), with the deviations listed below. Script:
`scripts/elgordo_djaref.py --stage 2`.

**Disclosed before freezing.** Read so far: the DJA table header (column names, `TUNIT`, `APER_k` / `ASEC_k`,
`MINKRON`, `TOTCFILT`, `ZP` = 28.9 while `TUNIT` of the flux columns says `uJy`: the two disagree, so rule 2a
decides the unit from the rows). Footprint check (positions only, array bounds of the SCI header): GO 1176
`jw01176-o241_t012_nircam_clear-f200w_i2d.fits` (10076 × 4579) contains all 56 images. Stage 1's
`extra_predicted` counts per family (15 families) were read; no predicted position or flux.
The stage-1 handoff proposed `flux_aper_k × tot_corr`. The grizli 1.14.3 source (`pipeline/photoz.py`,
`prep.get_kron_tot_corr`; SOURCES.md) shows that `tot_corr` only corrects a Kron flux to total for a point source;
the total is `{band}_flux_aper_k × (flux_auto / flux_aper_k) × {band}_tot_corr`, with the ratio from the detection
image. Rule 2 uses that definition (deviation from the handoff wording, made before any row was read).

## Hypothesis, ordinary explanations, what would change the plan
- Question 1 (method): does the DJA v7.0 reference with grizli total fluxes and the pipeline-mosaic local noise pass
  gate C2 on El Gordo?
- Question 2 (field): are the magnified, unmatched images the Caminha+2023 model predicts present in F200W at the
  predicted flux?
- Ordinary explanations of an `absent`: model μ or position error (stage-1 rms 0.754″), low surface brightness,
  deblending, cluster-member or ICL light (pipeline mosaic not bCG-subtracted), a reference flux that includes
  unrelated light, the point-source assumption of `tot_corr` for extended images. Exotic readings are not considered.
- Rule 2a or 2b fails → stop before stage 3; nothing read. C2 fails → the reference is not validated; no flag is
  read. C2 passes and a flag → contact sheet, then `/vet-candidate`; no interpretation in this stage.

## Inputs
- model_prediction: `elgordo-caminha23` (sha256 in `lensmodel.ELGORDO_CAMINHA23`), frame offset of D-034 rule 1.
- observed: the 56 spectroscopic images of `caminha23_images_jwstframe.ecsv`; z from the stage-1
  `imageplane_elgordo_elgordo-caminha23.ecsv` `z_used`.
- observed: GO 1176 level-3 `jw01176-o241_t012_nircam_clear-f200w_i2d.fits` (stamps and the SCI extension by S3
  range reads) and its `_cat.ecsv` (shapes); not stored in git.
- observed: DJA v7.0 `elgordo-grizli-v7.0-fix_phot.fits` (17,170,560 bytes), read into memory, sha256 recorded,
  not stored.

## Frozen rules
1. **Frame tables:** `m0417_extra.run_field` unchanged (mean back-traced source, the estimator gate K validated);
   its pipeline forced test and gate C are reported, not used.
2. **Reference catalogue** (`elgordo_djaref.to_cat_dja`, unit-tested before the run), aperture k = 1 (0.5″):
   - 2a, unit (`flux_unit_jy`): u = median over rows with finite positive `flux_auto` and finite `mag_auto` of
     `mag_auto + 2.5 log10(flux_auto)`; u within 0.05 of 23.9 → µJy (1e-6 Jy); within 0.05 of 28.9 → 10 nJy;
     otherwise stop.
   - 2b, grid (`grid_rotation_deg`): least-squares affine map of (`x_image`, `y_image`) → (ΔRA cos δ, ΔDec) on all
     rows (ξ east-positive). Rotation = deg(atan2(∂ξ/∂y, ∂η/∂y)), the angle of +y from north toward east. Passes
     when |rotation| ≤ 1° and the determinant is negative (east to the left, standard parity); otherwise stop.
     Pixel scale = √|det| of that map.
   - flux = `f200w_flux_aper_1 × flux_auto / flux_aper_1 × f200w_tot_corr` × unit; S/N =
     `f200w_flux_aper_1 / f200w_fluxerr_aper_1`; rows dropped when the flux is not finite or ≤ 0, `flux_aper_1` ≤
     0, or `f200w_mask_aper_1` > 0.2 × π (APER_1 / 2)² (grizli's own mask rule).
   - `a_image`, `b_image` × pixel scale; PA east of north = `cr.sky_pa`(deg(`theta_image`)) (SEP convention,
     radians ccw from +x, as measured for SEP on Abell 2744 in stage 2a); ext = true for every row (ASSUMPTION: the
     table has no star flag; deviation from Abell 2744's `flag_star`). Depth by `m0417_canucsref.depth_mag`.
3. **Stage 3:** `m0417_arcref.run` with this catalogue, frame offset (0, 0) and the containing-source rule on, rules
   otherwise unchanged.
4. **Local noise** (`m1423_canucsref` rule 2, constants unchanged) on the pipeline F200W SCI converted to nJy per
   pixel (`a2744_uncoverref.to_njy_per_pixel`), source exclusion around all DJA rows, pixel scale from the mosaic
   WCS (as Abell 2744 stage 2b).
5. **Gate C2** with the shared-pair rule (`m0417_arcref.gate_c2_shared`), thresholds unchanged.
6. **Flag:** a target with `class_local` = `absent` while C2 passes.
7. Reported, not rules: C2 and flags without the local-noise rule, without the containing-source rule, and plain
   `gate_c2`.

## Data volume
SCI extension 10076 × 4579 float32 ≈ 185 MB, read once into memory by a range stream; DJA table 17 MB in memory;
pipeline `_cat.ecsv` cached in the session and deleted after use.

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged.
