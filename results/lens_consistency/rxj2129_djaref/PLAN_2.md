# RXJ2129 reference-catalogue extra-image test: stage 2 (plan, frozen)

Frozen 2026-10-23 before any DJA flux row, F200W pixel, frame table or stage-2 prediction was read. Follows
`PLAN_1.md` (gate K passed, rule 2 passed). Same rules as El Gordo stage 2b (`../elgordo_djaref/PLAN_2.md` with
`PLAN_2b.md`'s rule 2a), with the deviation listed below. Script: `scripts/rxj2129_djaref.py --stage 2`.

**Disclosed before freezing.** Read so far: stage 1's outputs (gate K, per-family rms, the per-family
`extra_predicted` counts), the DJA table header, and a footprint check (positions only, array bounds of the SCI
header WCS; not stored): GO 2767 `jw02767-o002_t001_nircam_clear-f200w_i2d.fits` (4539 × 4520) contains all 31
images; it is the only level-3 tile of program 2767 on S3. No predicted position or flux, no row and no pixel.

## Hypothesis, ordinary explanations, what would change the plan
- Question 1 (method): does the DJA v7.0 reference with grizli total fluxes and the pipeline-mosaic local noise pass
  gate C2 on RXJ2129?
- Question 2 (field): are the magnified, unmatched images the Dhanasingham+2026 Lenstool model predicts present in
  F200W at the predicted flux?
- Ordinary explanations of an `absent`: model μ or position error (stage-1 rms 0.449″; 0.78″ / 0.66″ for the
  S3-1 / S3-2 knots), low surface brightness, deblending, cluster-member or ICL light (pipeline mosaic not
  bCG-subtracted), a reference flux that includes unrelated light, the point-source assumption of `tot_corr` for
  extended images. Exotic readings are not considered.
- Rule 2a or 2b fails → stop; nothing read. C2 fails → the reference is not validated; no flag is read. C2 passes
  and a flag → contact sheet, then `/vet-candidate`; no interpretation in this stage.

## Inputs
- model_prediction: `rxj2129-lenstool-maps` (sha256 in `lensmodel.RXJ2129_LENSTOOL_MAPS`), no frame offset.
- observed: the 31 images of `dhanasingham26_images.ecsv`, every family at its spectroscopic z (stage-1 `z_used`).
- observed: GO 2767 level-3 `jw02767-o002_t001_nircam_clear-f200w_i2d.fits` (stamps and the SCI extension by S3
  range reads) and its `_cat.ecsv` (shapes); not stored in git.
- observed: DJA v7.0 `rxj2129-grizli-v7.0-fix_phot.fits` (9,152,640 bytes), read into memory, sha256 recorded,
  not stored.

## Frozen rules
Rules 1–7 of `../elgordo_djaref/PLAN_2.md`, with rule 2a as in `../elgordo_djaref/PLAN_2b.md`
(`flux_unit_jy(t, with_tot_corr=True)`), unchanged in every constant and threshold, applied to the inputs above.
SMACS stage 2's D2 (z at the grid edge) does not apply: no z here is fitted.

**Deviation D1, knot group (made before any read; `merge_knot_targets`, unit-tested).** Families 3, S3-1 and S3-2
are one source galaxy (S3-1 / S3-2 are knots in system 3's arcs, all at z 1.5194). For flags only:
- Two predicted targets of *different* families of this group are joined when their search circles overlap
  (separation < sum of their `search_arcsec`); joins are transitive. Each joined set is one target.
- A joined set's class is the first member class in the order recovered / confused, then ambiguous, then
  absent, then the first member's class (so a detected knot or arc anywhere in the region keeps the region from
  being flagged; this is the reverse of the shared-pair order, which is conservative for controls).
- Flags (rule 6) = targets outside any joined set with `class_local` = `absent`, plus joined sets whose class is
  `absent`. Every per-row class is reported unchanged (`target_classes`); joined sets are listed (`knot_sets`).
- Controls are counted by rule 5 unchanged (the shared-pair rule already merges knot controls of different
  families that share one source pair).

## Data volume
SCI extension 4539 × 4520 float32 ≈ 82 MB, read once into memory by a range stream; DJA table ≈ 9 MB in memory;
pipeline `_cat.ecsv` cached in the session and deleted after use; Lenstool members (≈ 300 MB) re-staged by range
reads from the Zenodo zip into the session cache and deleted before the session ends.

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged.
