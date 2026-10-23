# SMACS 0723 reference-catalogue extra-image test: stage 2 (plan, frozen)

Frozen 2026-10-23 before any DJA flux row, F200W pixel, frame table or stage-2 prediction was read. Follows
`PLAN_1.md` (gate K passed, rule 2 passed). Same rules as El Gordo stage 2b (`../elgordo_djaref/PLAN_2.md` with
`PLAN_2b.md`'s rule 2a), with the deviations listed below. Script: `scripts/smacs_djaref.py --stage 2`.

**Disclosed before freezing.** Read so far: stage 1's outputs (gate K, the per-family `extra_predicted` counts, the
image-plane `z_used` per family), the DJA table header, and a footprint check (positions only, from the SCI header
WCS): ERO GO 2736 `jw02736-o001_t001_nircam_clear-f200w_i2d.fits` (10276 × 4762) contains all 60 images. No
predicted position or flux, no row and no pixel.

## Hypothesis, ordinary explanations, what would change the plan
- Question 1 (method): does the DJA v7.4 reference with grizli total fluxes and the pipeline-mosaic local noise pass
  gate C2 on SMACS 0723?
- Question 2 (field): are the magnified, unmatched images the Mahler+2022 ICLv2 model predicts present in F200W at
  the predicted flux?
- Ordinary explanations of an `absent`: model μ or position error (stage-1 rms 0.321″ on the fixed-z set; the free
  systems' z are fitted, not measured), low surface brightness, deblending, cluster-member or ICL light (pipeline
  mosaic not bCG-subtracted), a reference flux that includes unrelated light, the point-source assumption of
  `tot_corr` for extended images. Exotic readings are not considered.
- Rule 2a or 2b fails → stop; nothing read. C2 fails → the reference is not validated; no flag is read. C2 passes
  and a flag → contact sheet, then `/vet-candidate`; no interpretation in this stage.

## Inputs
- model_prediction: `smacs0723-iclv2` (sha256 in `lensmodel.SMACS0723_MAHLER22_ICLV2`), no frame offset.
- observed: the 60 images of `mahler22_images.ecsv`; z from stage 1's `imageplane_smacs0723_smacs0723-iclv2.ecsv`
  `z_used` (fixed z for systems 1, 2, 3, 5, 7; `fitted_z` for the rest).
- observed: GO 2736 level-3 `jw02736-o001_t001_nircam_clear-f200w_i2d.fits` (stamps and the SCI extension by S3
  range reads) and its `_cat.ecsv` (shapes); not stored in git.
- observed: DJA v7.4 `smacs0723-grizli-v7.4-fix_phot.fits` (38,626,560 bytes), read into memory, sha256 recorded,
  not stored.

## Frozen rules
Rules 1–7 of `../elgordo_djaref/PLAN_2.md`, with rule 2a as in `../elgordo_djaref/PLAN_2b.md`
(`flux_unit_jy(t, with_tot_corr=True)`), unchanged in every constant and threshold, applied to the inputs above.

**Deviations (made before any read):**
- D1, all 60 images. Unlike Abell 2744 stage 2b (photometric-z images left out), the free-z systems take part
  with their `fitted_z`, as in MACS0417 (`../macs0417_extra/PLAN.md`); PLAN_1.md's stage-2 inputs were defined
  this way.
- D2, z at the grid edge. A family whose `z_used` is ≥ 10.0 (the `fitted_z` grid's upper edge; system 8, its free z
  unconstrained, D-024) is not a flag target: its targets are classified and reported (`flags_z_edge`) but not
  flagged. Function `split_z_edge`, unit-tested.
- D3, labels. Each flag carries whether its family's z is fixed (rule 0) or fitted; both are flags.

## Data volume
SCI extension 10276 × 4762 float32 ≈ 196 MB, read once into memory by a range stream; DJA table ≈ 38.6 MB in
memory; pipeline `_cat.ecsv` (3.4 MB) cached in the session and deleted after use.

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged.
