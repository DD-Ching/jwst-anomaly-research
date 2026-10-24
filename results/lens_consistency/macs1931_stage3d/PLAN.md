# MACS J1931 stage 3d: continuum-restricted frame check F, validated on MACS0417 first (frozen before any run)

**Why.** Stages 3, 3b and 3c stopped at frame check F (18, 18, 17 of 51 within 0.5″; rule ≥ 20). Their medians were
≤ 0.05″, and the 3c contact sheet (by eye, not measured) showed about 20 of 51 positions with no visible F200W emission.
F is a check of the catalogue *frame* (registration), so stage 3d counts only images that carry measurable F200W
light at their position. Completeness for the extra-image test stays the job of gate C2 (unchanged).

**Not blind.** The threshold below was chosen after the 3c post hoc reading (about 31 of 51 visible by eye, 17
matched). It is therefore not an independent test of that reading; the MACS0417 validation (rule 3d-2) and gate C2
remain the quality gates. Written before any stage 3d code runs on either field.

## Rules (functions in `scripts/m1931_stage3d.py`, unit-tested before the run)
- **3d-0 builder:** stage V2 builder unchanged (`m1931_deblend.build`, background box 64 px), i.e. the box of stage 3b
  and of the MACS0417 V2 gate that passed (28 / 35).
- **3d-1 continuum criterion (per image):** on the builder's background-subtracted section, the sum in a 0.3″-radius
  circular aperture at the image position (`m1423_canucsref.aperture_sum`) divided by the local noise
  (`m1423_canucsref.local_noise`: σ-clipped std of 200 random 0.3″ apertures 1.5–4″ away, ≥ 0.5″ from every
  catalogue row, seed 20261022) is ≥ 5 (ASSUMPTION; the existing local-noise rule's MIN_SIGMA). An image off the
  section or with a NaN aperture is `off_section` and not counted. The criterion uses no catalogue match.
- **F_c (per field, combined with "and"):** among the n_c continuum-detected images, n_m are matched within 0.5″ by
  the stage 3 matcher (`bullet_extra.frame_check_f` on that subset). F_c passes iff n_c ≥ 10 **and**
  n_m ≥ 0.5 × n_c **and** |median dx| ≤ 0.1″ **and** |median dy| ≤ 0.1″ **and** the stage 3c chance control on the
  same subset passes (5″ shifts along 8 position angles, n_m ≥ 2 × mean shifted). ASSUMPTIONs: 5σ, 0.3″, 10, 0.5.
- **3d-2 validation (MACS0417, first):** F_c on the MACS0417 gold images (`m0417_extra.load_images`, moved by the
  measured `macs0417_frame` offset as in every MACS0417 own-reference stage) with the 3d-0 catalogue on the stage V2
  section. **If F_c fails on MACS0417, stage 3d stops** and nothing runs on MACS J1931.
- **3d-3 (MACS J1931):** F_c on the 51 used images (3b section, all images inside) with the 3d-0 catalogue. If F_c
  fails: no reading, stop; the field's own-reference route is then closed.
- **3d-4 (only if 3d-3 passes):** stage 3 rules 6–10 unchanged (V2-b chain with and without the containing-source
  rule, footprint rule, local-noise rule, gate C2 with shared pairs, flags = `absent` predicted rows when C2 passes);
  needs the VENUS maps (`--maps`). Flags are unvetted until `/vet-candidate`.

## Run
`python scripts/m1931_stage3d.py validate`, then `python scripts/m1931_stage3d.py f` (refuses unless 3d-2 passed),
then, only if F_c passed, `python scripts/m1931_stage3d.py run --maps DIR` (wiring `use_fc`, committed before its run:
stage 3 `run` with the 3b section and the 64 px box, rule 5 replaced by the stored F_c after a row-count check).
