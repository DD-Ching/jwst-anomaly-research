# SMACS 0723 arc difference imaging — result

Box: 24″ × 24″ centred between SN Atalanta's images A and B, (110.824442, −73.452577); it holds the Mahler+22
system 5 arc (PLAN.md, written before the run). `arc_diff.py --kernel-px 4 --kernel-stamps --inject 100
--min-bands 2`, `OMP_NUM_THREADS=1`, `--procs 4` (full commands in each `summary.json`). Every number is `derived`
unless labelled; thresholds are ASSUMPTIONs. The frozen "on an image" rule is the tested function
`classify` in `scripts/smacs_arcdiff_images.py` (known control ≤ 0.15″; else morph `pass` and ≤ 1.0″ from a
Mahler+22 ICLv2 image → vet). Sheets were looked at before writing this.

## Deviations from PLAN.md
- **Post hoc:** `arc_diff_morph.py --all-hits` (`morph_all.json`, `image_distances_all.json`) in addition to the
  default run, which judges only the hits not flagged `near_bright_core`.
- **Post hoc:** a 12″ box re-centred on h061 (`h061_box/`), because h061's `edge` veto came from this run's box
  boundary.
- **Added:** a second, independent pair (GO-4043 o001 → DD 12774 o002), not in the PLAN.

## Pair 1: ERO o001 (2022-06-07) → VENUS o057 (2026-06-05), F150W + F444W (`o001_o057/`)
- **Positive control passed:** Atalanta A = h021 (S +58.2 F150W, +49.5 F444W), B = h024 (+52.8, +54.3); both
  `pass` in `arc_diff_morph.py`; on the sheet both are clean, isolated point brightenings in both bands.
- Injection recovery per band (`summary.json`): F150W 0.82 / 0.60 / 0.37 at 27.5 / 28.0 / 28.5 AB; F444W
  0.76 / 0.62 / 0.30 (0.68 / 0.40 / 0.17 and 0.61 / 0.44 / 0.18 away from bright cores).
- **73 hits** (27 not `near_bright_core`). Default morph on those 27: 4 pass, 9 dipole, 8 extended, 3 rosette,
  2 edge, 1 sign_flip. All-hits morph: 5 pass.
- **Frozen rule** (`image_distances.json`, `image_distances_all.json`): 2 `known` (Atalanta), **0 `vet`**. The
  other passes are `off_image`: h039 (3.71″ from 18.1), h063 (1.37″ from 6.1), h019 (4.91″ from 17.2; all-hits
  run only). On the sheet h039 is a one-band dipole on a compact galaxy (F150W), h063 a one-band dipole on a
  compact source beside a diffraction spike (F444W), h019 a mixed-sign residual on a bright compact core.
- **h061** (0.15″ from image 12.1, `edge`, 0.25–0.44″ from the box edge): in the re-centred box its position has
  S −16.2 (F150W) and +5.1 (F444W), opposite signs, so it is not a hit there. The nearest hit there, h032
  (0.16″ away, 0.27″ from 12.1), is a `dipole`; on the sheet a positive/negative pair across a compact knot in
  F150W (dipole ratio 1.38) with a weak one in F444W (0.91). Classified ordinary by the frozen rule (not a pass).
- **Outcome:** 0 candidates on the catalogued images in this box for this pair, beyond SN Atalanta. Not a limit
  below the injection depth above.

## Pair 2: GO-4043 o001 (2024-05-11) → DD 12774 o002 (2026-07-06), F115W + F444W (`o4043_o12774/`)
- **Positive control failed under the 2-band rule:** Atalanta S = +2.5 / +2.4 (F115W, A / B) and +15.4 / +19.0
  (F444W); not hits. Forced photometry of B (`vet/lightcurves.ecsv`): F115W 13.7 ± 16.7 nJy in DD 12774, F444W
  87.4 ± 7.2 nJy. Injection recovery is low: F115W 0.26 / 0.11, F444W 0.34 / 0.10 at 27.5 / 28.0 AB. Per the PLAN this
  pair is not trusted for a null; its hits are recorded but give no limit.
- 38 hits; all-hits morph: 1 pass, 22 extended, 11 edge, 3 dipole, 1 rosette. The frozen rule gives 0 `vet` (`image_distances_all.json`);
  the two hits ≤ 1″ from an image (h018, h020 near 12.1) are diffraction-spike stripes on the sheet.
- **h000** (the one pass; 1.81″ from image 6.1; S 469): a point source with a full PSF only in the 2024 GO-4043
  mosaics. Forced photometry in all 19 covering mosaics (`vet/lightcurves.ecsv`): F090W 2579 nJy, F115W 2688,
  F444W 258 on MJD 60441.7–60441.9, versus 7–49 nJy in F090W–F444W in 2022, 2026-06 and 2026-07. It holds the same sky position in three
  filters taken over ≈ 0.2 d. It lies 0.081″ from AT 2024jfz (`h000_tns.json`; TNS position as recorded on #237,
  not re-queried). **Known transient, not a candidate.** Not checked: whether it is multiply imaged (counter-images
  in the lens model).

## Regenerate
Commands are in `o001_o057/summary.json`, `o4043_o12774/summary.json` and `h061_box/summary.json`; then
`arc_diff_morph.py --hits … --maps … [--all-hits] --out …`, `smacs_arcdiff_images.py MORPH --known
known_controls.csv --out …`, `arc_diff_sheet.py MAPS MORPH OUT [RA,DEC …]` and
`hit_lightcurve.py --pos h000 110.83411267801681 -73.45173415397286 --pos Atalanta_B 110.823747 -73.452975
--query-only h000 --out o4043_o12774/vet`. The maps (`--maps`) are not stored.
