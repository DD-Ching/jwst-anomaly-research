# MACS J0308 5293 o005 → VENUS 6882 o015: D-027 two-epoch grid (result, 2026-10-27)

Plan `PLAN.md` (pushed before any pixel was read). Tool `scripts/transient_grid.py`, unchanged. Numbers:
`summary.json`, `hits.ecsv` (27 rows, 27 groups), post hoc `posthoc.json` (`scripts/m0308_grid_posthoc.py`, tested
in `tests/test_m0308_grid_posthoc.py`). Imaging (observed): GO 5293 o005 (2025-02-04) and VENUS 6882 o015
(2025-08-23) level-3 `_i2d` in F300M / F410M. Four mosaics (bytes in `posthoc.json` `mosaic_bytes`) and four
`_cat.ecsv` streamed from MAST into the gitignored data root and deleted after the run. Every number below is
derived; thresholds are the D-027 ASSUMPTIONs.

## Frozen result
- **Outcome U: gate G1 fails.** Injected −0.75 mag brightenings were recovered 0.518 (500; 0.759 where the
  fraction set the flux), against the frozen ≥ 0.85. Per band at −0.75 mag: F300M 0.654, F410M 0.58; with two
  bands, both must pass. At −0.44 mag: 0.038.
- **G2 passes:** noise scales F300M 1.57, F410M 1.50.
- 2,215 positions from the four pipeline catalogues, 2,035 measured in ≥ 1 band, 2,030 in both. 27 hits in 27
  groups, 5 of them mixed-sign.
- **Known case** (`posthoc.json` `known_case`): no grid row lies within 0.5″ of the listed SN Ares position. The
  nearest row (`g02256`, from the epoch-2 F410M catalogue) is 0.678″ away and brightens at +3.5σ (F300M) and +4.4σ
  (F410M) from an epoch-1 flux consistent with zero. That is below the 5σ hit threshold. The rows nearest
  "Ares_Image2" (1.366″) and "Ares_Image3" (0.399″) do not change (|σ| ≤ 1.3). SN Ares is not recovered as a hit,
  so the run has no positive control.

## Groups (contact sheet e1 / e2 / e2 − e1 per band, 3″ stamps on the epoch-1 pixel grid, looked at; not stored)
Reading from the contact sheet. The post hoc `classify` rules (`posthoc.json` `classification`, every matching
label kept) agree; groups 4 and 24 get two labels each. The epochs differ in roll by about 180° (`pa_aper` 75.12° vs
255.0°), so star spikes and wings flip between the epochs.

| Groups | Reading |
|---|---|
| 7, 23, 24, 25, 26 | PSF wings and spikes of a bright star 1.41–2.32″ away; dipole residuals at the star |
| 0, 2, 5, 6, 8–22 (and 4, 24) | Covered by one exposure in epoch 1 in both bands (`CON` 1, relative weight 0.16–0.35). Epoch-1 stamps look noise-dominated at the display stretch |
| 1, 3, 4 | A compact source in epoch 1 only, at the **same offset from a bright star** (8.32–8.35″, position angle 120.0–120.2°, three different stars of 6939.5–11892.8 µJy in epoch-1 F410M, `bright_star`). Read as an optical artefact of the star (post hoc, below); group 4 is also single-exposure |

No group is left without an ordinary explanation. Because G1 fails, this is not a sensitivity statement.

## Post hoc (not in the plan; `posthoc.json`)
- **Groups 1, 3, 4 (fixed-offset star image).**
  - Three different bright stars (positions in `bright_star`) give the same sky-frame offset vector at the
    epoch-1 roll. An astrophysical source has no reason to sit at a fixed offset from three unrelated stars.
  - In the epoch-1 `_cal` exposures (`detector_level`), groups 1 and 3 are missing in dither 1 (aperture sums
    −0.04–0.25) and present in dithers 2 and 3 (1.20–2.07), in both the F410M and the later F300M sequence. A
    variable sky source could not turn off, on, off and on again in step with the dither pattern. Group 4 falls on
    a detector only in dither 3, so this test does not apply to it; it rests on the offset vector alone.
  - Steady control rows near each group are seen in every dither (`detector_level_controls`, 0.57–2.58). So the
    dither-1 gap is not a WCS offset of that exposure.
  - The same positions are detected in the epoch-1 SW mosaics F115W / F150W / F200W (`sw_epoch1`, aperture S/N
    29–125 against the stamp's pixel rms).
  - In epoch 2, no source is visible at the same offset rotated with the roll (position angle 300°; stamps
    looked at, not stored).
  - The optical path is not identified (ghost or scattered light), so the reading is "instrumental, fixed offset
    from a bright star". It is not a named artefact type. The stored vector (8.33″, 120.1°) is a sky-frame offset
    measured at `PA_APER` 75.12°; an instrument-frame artefact rotates with the roll, so on another pair it must
    be rotated by the change in `PA_APER` before use.
- Single-exposure groups (labelled `single_exposure` only): 18 of the 19 have epoch-1 relative weight 0.32–0.35, the edge strip covered by one
  dither (`coverage`). Group 6 has F300M weight 0.16.
- G1 failure, not diagnosed. Only by-band rates are stored (above). Whether single-exposure regions or the
  shallower epoch 1 (1,095 s vs 2,061 s per band, `PLAN.md`) drive it is not counted.

## Reading and limits
- Undetermined. The automatic two-band rule recovered about half of the −0.75 mag injections, and the known SN
  was not a hit. This pair cannot bound transients in the MACS J0308 core with the frozen thresholds.
- The 27 hits contain no unexplained group. That is a statement about these hits only, not a null.
- The fixed-offset star image is a new artefact family for the grids: three fading groups at 8.3″ from bright
  stars (two of them confirmed by the dither test). A later grid should test this offset, rotated by its epoch's
  `PA_APER` − 75.12°, before reading a fading compact source near a bright star.
- Only LW medium bands; SW not gridded; catalogued positions only.
- `grid.ecsv` is not kept (as on the earlier grids); the command in `PLAN.md` regenerates it in about 10 s once
  the mosaics are local. The products are jwst 3.0.0 builds (header `CAL_VER`; not stored).

## Deviations
- None in the frozen run. Notes: the post hoc checks above, the `classify` thresholds (3″ star distance, the
  8.33″ / 120.1° offset ± 0.2″ / 1°, one exposure; none of the 27 is off the footprint) and the control rule are post hoc ASSUMPTIONs fixed after
  the contact sheet was seen.
