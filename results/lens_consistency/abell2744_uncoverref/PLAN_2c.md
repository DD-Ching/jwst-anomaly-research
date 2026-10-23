# Abell 2744: UNCOVER-reference extra-image test, stage 2c (plan, frozen)

Frozen 2026-10-23 before any F277W / F356W flux of the DR3 SUPER catalogue or any F277W / F356W pixel of GO 2561 was
read (only the catalogue's column names and the S3 file list were read). Follows `PLAN_2b.md` (gate C2 passed in F200W,
no flag; 4.2 at 9.7 σ and 78.3 at 8.7 σ `ambiguous`, just under the 10 σ absent threshold). Script:
`scripts/a2744_uncoverref.py --stage 2c`.

## Hypothesis, ordinary explanations, what would change the plan
- Question: are the stage-2b `ambiguous` targets 4.2 and 78.3 (and every other target) present in the long-wavelength
  bands at the flux the UNCOVER v2.0 maps predict from their reference images, measured with the stage-2b rules?
- Ordinary explanations of an `absent`: as `PLAN_2b.md` (model μ or position error, low surface brightness,
  deblending, BCG / ICL light, unrelated light in the reference flux), plus colour: an image whose reference flux in
  this band is blended with a redder neighbour. Exotic readings are not considered.
- A frozen flag → contact sheet, then `/vet-candidate`; no interpretation in this stage. Present in a band → the
  stage-2b ambiguity is resolved as present for that target. Gate C2 fails in a band → that band is not read.

## Inputs
- model_prediction / derived: the committed stage-2b frame tables `frame/predicted_abell2744.ecsv` and
  `frame/control_abell2744.ecsv` (WLS source, F200W pipeline frame), not recomputed.
- observed: UNCOVER DR3 `UNCOVER_DR3_LW_SUPER_catalog.fits` (sha256 pinned in stage 2a), streamed, deleted after use.
- observed: GO 2561 `jw02561-o001_t003_nircam_clear-{f277w,f356w}_i2d.fits` SCI extensions (≈ 827 MB files; the
  SCI extension alone is read once per band by a range stream into memory, not stored).

## Frozen rules
Every rule is applied **per band**, for F277W and F356W independently, with the stage-2b code paths:
1. **Reference catalogue:** `to_cat_super(t, band)`: as `PLAN_2b.md` rule 2 with `f_<band>` / `e_<band>` in place of
   `f_f200w` / `e_f200w`; shapes, PA and star flag unchanged (they are band-independent in the SUPER catalogue).
2. **Stage 3:** `m0417_arcref.run` on the committed frame tables, frame offset (0, 0), containing-source rule on.
3. **Local noise:** `m1423_canucsref.apply_local` (constants unchanged, seed 20261022 restarted per band) on that
   band's pipeline SCI converted by `to_njy_per_pixel`, with the mosaic's own pixel scale from its WCS.
4. **Gate C2 (per band):** `m0417_arcref.gate_c2_shared`, thresholds unchanged.
5. **Flag (per band):** a target with `class_local` = `absent` in a band whose gate C2 passes.
6. **Target outcome** (`a2744_uncoverref.outcome_2c`, unit-tested before the run), per target over the two bands:
   `flag` if flagged (rule 5) in at least one band; else `present` if `class_local` is `recovered` or `confused` in
   at least one band whose C2 passes; else `undetermined`. Reported for every target; 4.2 and 78.3 are the question.
7. Reported, not rules: per-band C2 without the local-noise rule; the per-band class of 3.2.

Two bands double the chance of a false `absent` against stage 2b's one band (ASSUMPTION accepted: a flag only starts
vetting, never a claim).

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged.
