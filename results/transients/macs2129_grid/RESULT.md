# MACS J2129 4903 o004 → VENUS 6882 o040: D-027 two-epoch grid (result, 2026-10-27)

Plan `PLAN.md` (pushed before any pixel was read). Tool `scripts/transient_grid.py`, unchanged. Numbers:
`summary.json`, `hits.ecsv` (16 rows, 13 groups), post hoc `posthoc.json` (`scripts/m2129_grid_posthoc.py`, tested
in `tests/test_m2129_grid_posthoc.py`). Imaging (observed): GO 4903 o004 (2024-10-01) and VENUS 6882 o040
(2025-10-28) level-3 `_i2d` in F277W / F444W (bytes in `posthoc.json` `mosaic_bytes`) and four `_cat.ecsv`,
streamed from MAST into the gitignored data root and deleted after the run. Every number below is derived;
thresholds are the D-027 ASSUMPTIONs.

## Frozen result
- **Outcome U: gate G1 fails.** Injected −0.75 mag brightenings were recovered 0.78 (500; 0.928 of the 349 where the
  fraction set the flux), against the frozen ≥ 0.85. Per band at −0.75 mag: F277W 0.878, F444W 0.824; both bands
  must pass. At −0.44 mag: 0.208.
- **G2 passes:** noise scales F277W 1.47, F444W 1.25.
- 2,524 positions from the four pipeline catalogues, 847 measured in ≥ 1 band, 843 in both. 16 hits in 13 groups,
  0 mixed-sign. All 16 are fades except three epoch-2-catalogue rows that brighten and share groups 7, 9, 10 with a
  fading row.
- No known transient in this core is pinned, so there is no positive control beyond the injections.

## Groups (contact sheet e1 / e2 / e2 − e1 per band, 60-pixel stamps (about 3.8″) on the epoch-1 pixel grid, looked at; not stored)
The post hoc `classify` rules (`posthoc.json` `classification`) agree with the reading. The two epochs differ in
roll by 4.0° (`pa_aper` 77.16° vs 73.16°), so diffraction spikes move slightly between them.

| Groups | Reading |
|---|---|
| 0, 1, 7, 8, 9, 10, 11 | On diffraction spikes / wings of two bright stars 2.54–2.84″ away (`star_wings`); the difference stamps show spike residuals |
| 2, 3, 4, 5, 6, 12 | A compact epoch-1-only spot (1–2 pixels at the display stretch; not stored) at a position covered by one exposure in epoch 1 in both bands (`CON` 1, relative weight 0.331–0.343; `single_exposure`). Absent in the four-exposure epoch 2 |

No group is left without an ordinary explanation. Because G1 fails, this is not a sensitivity statement.

## Post hoc (not in the plan; `posthoc.json`)
- The six single-exposure spots appear in both F277W and F444W. A cosmic ray hits one exposure of one filter, so it
  does not explain a spot in both. A detector-fixed defect (warm or hot pixel) seen through the same single dither
  position in both filters would; not tested at detector level here (no `_cal` read), so the reading stays
  "single-exposure coverage, origin not identified", as in #576.
- The #576 fixed-offset star image, rotated to this epoch-1 roll (8.33″ at 122.14°, `ghost_pa_rotated_deg`),
  matches no group.
- G1 failure, not diagnosed. Of the grid rows with a finite epoch-1 flux, 124 (F277W) and 121 (F444W) lie on
  single-exposure epoch-1 coverage, 29 / 31 on two exposures and 806 / 804 on three (`grid_e1_n_exp`). Whether
  these rows carry the missed injections is not counted (injection positions are not stored by the tool).

## Reading and limits
- Undetermined. The automatic two-band rule recovered 0.78 of the −0.75 mag injections, short of the frozen 0.85.
  This pair cannot bound transients in the MACS J2129 core with the frozen thresholds.
- The 16 hits contain no unexplained group. That is a statement about these hits only, not a null.
- With #576, two consecutive VENUS pairs against a shallower GO epoch 1 with single-exposure strips fail G1
  (0.518 in `../macs0308_grid/summary.json`, 0.78 here). A further grid of this kind should first count the recovery on rows with ≥ 2 epoch-1 exposures,
  which needs the tool to store injection positions (a tool change, so a new plan).
- Only LW wide bands; SW not gridded; catalogued positions only. Products are jwst 3.0.0 builds (header `CAL_VER`;
  not stored). `grid.ecsv` is not kept; the command in `PLAN.md` regenerates it in about 6 s once the mosaics are
  local (not stored).

## Deviations
- None in the frozen run. Notes: the post hoc checks above and the `classify` thresholds (3″ star distance, the
  rotated 8.33″ / 122.14° offset ± 0.2″ / 1°, one exposure) are post hoc ASSUMPTIONs carried over from #576.
