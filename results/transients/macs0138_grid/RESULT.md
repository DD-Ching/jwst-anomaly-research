# MACS J0138 2345 o005 → 6549 o003: D-027 two-epoch grid with SN Encore as positive control (result, 2026-10-27)

Plan `PLAN.md` (pushed before any pixel was read). Tool `scripts/transient_grid.py` as on main, with `--coverage`.
Numbers: `summary.json`, `hits.ecsv` (7 rows, 4 groups), `injections.ecsv`, post hoc `posthoc.json`
(`scripts/m0138_grid_posthoc.py`, which reuses the tested `a1689_grid_posthoc` functions; test
`tests/test_m0138_grid_posthoc.py`). Imaging (observed): GO 2345 o005 (2023-11-17, SN Encore discovery epoch) and
DD 6549 o003 (2024-01-08) level-3 `_i2d` in F150W / F444W (bytes in `posthoc.json` `mosaic_bytes`; the epoch-1 F150W
mosaic is 1.29 GB, streamed from S3 for the field-wide grid and deleted after the run) and four `_cat.ecsv`. Every
number below is derived; thresholds are the D-027 ASSUMPTIONs.

## Frozen result
- **Outcome U: gate G1 fails.** Injected −0.75 mag brightenings were recovered 0.50 (500; 0.825 of the 263 where
  the fraction set the flux), against the frozen ≥ 0.85. Per band: F150W 0.602, F444W 0.718. At −0.44 mag: 0.062.
- **G2 passes:** noise scales F150W 1.52, F444W 1.31.
- 2,346 positions from the four catalogues, 871 measured in ≥ 1 band, 819 in both. 7 hits in 4 groups (3 hits
  mixed-sign).
- **Positive control (SN Encore): not recovered.** No hit group lies on the MRG-M0138 arcs. By eye (frozen plan:
  identification by eye against arXiv:2509.12301 Figure 1), images 1a / 1b sit on the southern arc. There, e2 − e1
  shows a dipole along the whole arc in both bands (bright on one side, dark on the other), and that dipole is much
  stronger than any point-like fade. No catalogue position on the arc passes in either band. The sheet was looked
  at and is not stored. Reading (by eye, not quantified): a relative sub-pixel registration offset between the
  epochs leaves a residual on bright arcs that this aperture-at-catalogue-position tool cannot separate from an SN
  ≈ 0.1″ from the host light.

## Groups (contact sheets e1 / e2 / e2 − e1 per band, 3″ stamps on the epoch-1 grid, looked at; not stored)
The epochs differ in roll by 30.1° (the difference of `PA_APER` 65.52° and 35.41° in `posthoc.json` `pa_aper`; the difference is not stored).

| Group | Reading |
|---|---|
| 0, 2 | Both sides of one bright point source with diffraction spikes about 1″ away (by eye, not stored); the difference stamps show spike residuals at the new roll. `classify` says `unexplained` because this source is below the post hoc F444W 99.5th-percentile star cut (nearest source above it 31.9″ / 33.3″). Epoch-1 coverage 2 exposures (`WHT` 0.51) |
| 1 | Spike residuals of a bright star 1.35″ away (`star_wings`); three hits, fade and brightening on opposite sides |
| 3 | A bright compact knot on the radial arc north-east of the BCG, not on the MRG-M0138 arcs: the difference is a dipole across the knot, and the hit sits on its fading side, 0.43″ from the knot centre (`star_wings` label from the flux cut; it is the knot, not a star). Fade 0.61 / 0.44 mag in F150W / F444W at 5.6σ / 5.5σ; epoch-2 F150W on 2 exposures |

All four groups have an instrumental reading (spike residuals at a new roll, or a registration dipole on a bright
compact source). Because G1 fails, this is not a sensitivity statement.

## Post hoc (not in the plan; `posthoc.json` `injection`)
- Recovery of the frac-1.0 injections by their expected calibrated σ (smaller band), as in #579:

  | expected σ | < 3 | 3–4 | 4–5 | 5–6 | 6–8 | ≥ 8 |
  |---|---|---|---|---|---|---|
  | n | 3 | 67 | 194 | 81 | 153 | 2 |
  | recovered | 0.333 | 0.149 | 0.222 | 0.556 | 0.987 | 0.000 |

- Where the expected σ is ≥ 5 in both bands (236 / 500), recovery is 0.831. Here even the 5–6σ bin is only 0.556.
  So unlike #579 (0.939, `results/transients/a1689_grid/posthoc.json`), restricting G1 to expected σ ≥ 5 would also fail the 0.85 gate on this pair.
- Floor-set injections (237 / 500) reach 4.98σ (F150W) and 4.92σ (F444W) at the median errors and are recovered
  0.140 (≥ 2 epoch-1 exposures, `split_frac_1.0`). Nominal injections on ≥ 2 exposures: 0.830 (253).
- Epoch-1 coverage is not the cause: 489 / 500 injections sit on ≥ 2 epoch-1 exposures, recovery 0.497.

## Reading and limits
- Undetermined. This pair cannot bound transients in the MACS J0138 core with the frozen thresholds, and the
  catalogue-position grid misses the one known transient in the field.
- Failed approach: a pair with a known lensed SN does not give the grid a working positive control when the SN sits
  on a bright arc. Catalogued arc positions get a registration dipole, and the SN is not at a catalogue centroid.
  A positive control for this tool needs a transient that is separate from bright extended light (e.g. SN Eos in
  #571), or a PSF-matched, registration-refined difference image (SFFT was trialled earlier: `results/transients/sfft_trial`).
- The SN Encore image positions are by eye from a published figure (ASSUMPTION); no coordinate is stored. The 3″ star
  reading for groups 0 / 2 is by eye.
- Catalogued positions only; F150W / F444W only; `grid.ecsv` is not kept (the `PLAN.md` command regenerates it in
  about 5 s once the mosaics are local).

## Deviations
- None in the run. Everything under "Post hoc", the group readings and the positive-control reading are post hoc or by
  eye, as the plan allowed for the control identification.
