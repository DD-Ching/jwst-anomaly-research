# Abell 1689 1840 o019 → 5782 o001: D-027 two-epoch grid (result, 2026-10-27)

Plan `PLAN.md` (pushed before any pixel was read). Tool `scripts/transient_grid.py` as on main, with `--coverage`.
Numbers: `summary.json`, `hits.ecsv` (15 rows, 12 groups), `injections.ecsv`, post hoc `posthoc.json`
(`scripts/a1689_grid_posthoc.py`, tested in `tests/test_a1689_grid_posthoc.py`). Imaging (observed): GO 1840 o019
(2023-02-25) and GO 5782 o001 (2025-07-12) level-3 `_i2d` in F277W / F356W (bytes in `posthoc.json`
`mosaic_bytes`) and four `_cat.ecsv`, streamed from S3 into the gitignored data root and deleted after the run.
Every number below is derived; thresholds are the D-027 ASSUMPTIONs.

## Frozen result
- **Outcome U: gate G1 fails.** Injected −0.75 mag brightenings were recovered 0.42 (500; 0.678 of the 307 where
  the fraction set the flux), against the frozen ≥ 0.85. Per band at −0.75 mag: F277W 0.536, F356W 0.422; both
  bands must pass. At −0.44 mag: 0.094.
- **G2 passes:** noise scales F277W 1.41, F356W 1.38.
- 3,460 positions from the four pipeline catalogues, 1,387 measured in ≥ 1 band, 1,384 in both. 15 hits in 12
  groups, 0 mixed-sign: 9 fades from epoch-1 catalogues, 6 brightenings from epoch-2 catalogues.
- Epoch-1 coverage is not the cause here: 495 / 500 injections sit on ≥ 2 epoch-1 exposures, and their recovery is
  0.420 (`summary.json` `by_e1_n_exp_min`).
- No known transient in this core is pinned, so there is no positive control beyond the injections.

## Groups (contact sheet e1 / e2 / e2 − e1 per band, 40-pixel stamps on the epoch-1 pixel grid, looked at; not stored)
The epochs differ in roll by 171.6° (the difference of `PA_APER` 287.6° and 116.0° in `posthoc.json` `pa_aper`; the difference is not stored), so diffraction spikes
land in different places.

| Groups | Reading |
|---|---|
| 0, 2, 3, 4, 5, 6, 7, 9, 10, 11 | On diffraction spikes / wings of bright stars 1.28–2.34″ away (`star_wings`, post hoc cut 3″); the difference stamps show spike residuals at the new roll |
| 8 | On a diffraction spike of a bright star 3.11″ away, just outside the 3″ cut (`classify` gives `unexplained`; by eye the spike crosses the position in both bands and both epochs at different angles) |
| 1 | At the edge of the epoch-1 mosaic: one epoch-1 exposure (`CON` 1, relative weight 0.259), and half the epoch-1 stamp has no data, beside the halo of a bright source (`single_exposure`) |

No group is left without an ordinary explanation. The #576 fixed-offset star image, rotated to this roll, matches no
group. Because G1 fails, this is not a sensitivity statement.

## Post hoc: why G1 fails here (not in the plan; `posthoc.json` `injection`)
- For each frac-1.0 injection we computed the calibrated significance it would reach with no noise draw: the
  injected flux over the error model of `transient_forced.compare` (both `ERR` terms, the 0.05 mag systematic floor,
  epoch 1 on the epoch-2 zero point), divided by the band's noise scale. Each injection keeps the smaller of the two
  bands' values.
- Recovery follows that expected significance (`frac_1.0_z_min_bins_all`):

  | expected σ (smaller band) | < 3 | 3–4 | 4–5 | 5–6 | 6–8 | ≥ 8 |
  |---|---|---|---|---|---|---|
  | n | 207 | 47 | 32 | 34 | 167 | 13 |
  | recovered | 0.010 | 0.000 | 0.219 | 0.794 | 0.988 | 0.692 |

- Where the expected σ is ≥ 5 in both bands (214 / 500), recovery is 0.939. The misses are injections the 5σ
  threshold cannot reach without a favourable noise draw; this run gives no sign of a tool failure.
- Floor-set injections (193 / 500) reach 3.36σ (F277W) and 2.04σ (F356W) at the median errors and are recovered
  0.010. Of the 307 nominal injections, 93 have expected σ < 5 and 34 lie at 5–6σ
  (`frac_1.0_z_min_bins_nominal`). For a source of equal flux in both epochs the 0.05 mag systematic floor caps a
  −0.75 mag brightening at about 10.9σ before the noise scale, about 7.7–7.9σ after it (an estimate from the error
  model, not stored), and the epoch-1 errors are the larger ones (median 0.036 / 0.043 against epoch-2 0.020 / 0.013
  in F277W / F356W).
- The ≥ 8σ bin (13 injections, 0.692) is small; those are bright positions, and what removes them is not checked.

## Reading and limits
- Undetermined. This pair cannot bound transients in the Abell 1689 core with the frozen thresholds.
- With #576 / #577 / #578, three pairs fail G1, and on this one the failure tracks the expected significance of each
  injection: G1 at "frac 1.0" counts faint catalogue positions where a −0.75 mag brightening is not 5σ. A G1 defined
  on injections with expected σ ≥ 5 (or a floor set in calibrated σ, #578) would be a D-027 amendment; this PR does
  not make it, the owner decides.
- The 15 hits contain no unexplained group. That is a statement about these hits only, not a null.
- The bin edges and the "smaller band" summary are post hoc. The 3″ star cut and the classify labels are post hoc
  ASSUMPTIONs carried over from #576 / #577. Only LW wide bands; catalogued positions only. `grid.ecsv` is not kept;
  the command in `PLAN.md` regenerates it in about 7 s once the mosaics are local.

## Deviations
- Pair: the previous handoff named 5782 o001 → 6882 o048; the plan chose 1840 o019 → 5782 o001 and says why. No
  deviation from the frozen plan in the run itself. Everything under "Post hoc" is post hoc.
