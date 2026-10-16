# SMACS 0723 arc difference imaging — result

Box: 24″ × 24″ centred between SN Atalanta's images A and B, (110.824442, −73.452577); it holds the Mahler+22
system 5 arc (PLAN.md, written before the run). `arc_diff.py --kernel-px 4 --kernel-stamps --inject 100
--min-bands 2`, `OMP_NUM_THREADS=1`, `--procs 4` (full commands in each `summary.json`). Every number is `derived`
unless labelled; thresholds are ASSUMPTIONs. The frozen "on an image" rule is the tested function
`classify` (rows also carry `posthoc_2arcsec`, and each file `controls_recovered`) in `scripts/smacs_arcdiff_images.py` (known control ≤ 0.15″; else morph `pass` and ≤ 1.0″ from a
Mahler+22 ICLv2 image → vet). Sheets were looked at before writing this.

## Deviations from PLAN.md
- **Post hoc:** `arc_diff_morph.py --all-hits` (`morph_all.json`, `image_distances_all.json`) in addition to the
  default run, which judges only the hits not flagged `near_bright_core`.
- **Post hoc:** a 12″ box re-centred on h061 (`h061_box/`), because h061's `edge` veto came from this run's box
  boundary.
- **Added:** a second, independent pair (GO-4043 o001 → DD 12774 o002), not in the PLAN.
- **The frozen 1.0″ radius fails the known case** (found in review): Atalanta A and B lie 1.69″ and 1.59″ from
  their catalogued images 5.1 and 5.2 (`image_distances_all.json`), so a transient placed like Atalanta would be
  `off_image`. **Post hoc:** every morph pass within 2.0″ of an image (the smallest round radius holding both
  Atalanta images; ASSUMPTION) was given multi-epoch forced photometry (`hit_lightcurve.py`, every covering
  public mosaic, conservative random-aperture errors): h063 (pair 1) and h000 (pair 2), plus h061.

## Pair 1: ERO o001 (2022-06-07) → VENUS o057 (2026-06-05), F150W + F444W (`o001_o057/`)
- **Positive control passed:** Atalanta A = h021 (S +58.2 F150W, +49.5 F444W), B = h024 (+52.8, +54.3); both
  `pass` in `arc_diff_morph.py`; on the sheet both are clean, isolated point brightenings in both bands.
- Injection recovery per band (`summary.json`): F150W 0.82 / 0.60 / 0.37 at 27.5 / 28.0 / 28.5 AB; F444W
  0.76 / 0.62 / 0.30 (0.68 / 0.40 / 0.17 and 0.61 / 0.44 / 0.18 away from bright cores).
- **73 hits** (27 not `near_bright_core`). Default morph on those 27: 4 pass, 9 dipole, 8 extended, 3 rosette,
  2 edge, 1 sign_flip. All-hits morph: 5 pass.
- **Frozen rule** (`image_distances.json`, `image_distances_all.json`): 2 `known` (Atalanta), **0 `vet`**. The
  other passes are `off_image`: h039 (3.71″ from 18.1), h063 (1.37″ from 6.1), h019 (4.91″ from 17.2; all-hits
  run only). All three stay morph `pass` (the function's verdict). By eye only (not a rule): h039 shows a
  positive/negative pair on a compact galaxy in F150W, h063 a negative residual with a positive side lobe on a
  compact source in F444W, h019 a mixed-sign residual on a bright compact core.
- **h063 under the post hoc 2.0″ radius** (`o001_o057/vet/lightcurves.ecsv`): forced flux F150W 142.9 ± 4.6 nJy
  (2022) → 126.9 ± 9.4 (2026-06) → 129.2 ± 14.4 (2026-07), F444W 59.9 ± 2.6 → 60.8 ± 4.6 → 49.1 ± 7.2; the
  2022 → 2026-06 changes are −16.0 ± 10.4 and +0.9 ± 5.3 nJy (1.5σ and 0.2σ). No flux change ≥ 2σ in forced
  photometry; not vetted further.
- **h061** (0.15″ from image 12.1, `edge`, 0.25–0.44″ from the box edge): in the re-centred box the S values at
  its position change (first box: negative in both bands, max |S| 11.4, `hits.ecsv`; re-centred box: F150W −16.2,
  F444W +5.1), so the difference measurement
  there depends on the box (kernel and registration fit); it is not a hit in the re-centred box. The nearest hit
  there, h032 (0.164″ away, just outside the 0.15″ match radius; 0.27″ from 12.1), is a `dipole` (ratio 1.38
  F150W, 0.91 F444W). Forced photometry at h061 (`o001_o057/vet/lightcurves.ecsv`): 2022 → 2026-06 F150W
  87.2 ± 4.6 → 83.1 ± 9.4 nJy (−0.4σ), F444W 40.4 ± 2.6 → 49.8 ± 4.6 (+1.8σ): no change ≥ 2σ **in this pair**.
  But F444W in 2024 (GO-4043) is 76.0 ± 2.2 nJy, +35.6 ± 3.4 over 2022 (≈ 10σ; F115W 59.0 ± 3.4 in 2024 vs
  41.4 ± 18.6 in 2026-07). By eye a diagonal stripe crosses the aperture in the 2024 F444W stamp
  (`vet/stamps_h061.png`), and that epoch's pair-2 hits h018 (0.51″ away) and h020 (0.91″) are stripes on the
  sheet; the spike geometry is **not computed**, so **h061 stays open (unvetted)**.
- **Outcome:** beyond SN Atalanta, the one morph pass within 2.0″ of a catalogued image (h063) shows no ≥ 2σ
  forced-photometry change; h061 (not a pass) is open because of its 2024 F444W excess (above). No candidate so
  far; h061 is undetermined (follow-up section below). Not a limit below the injection
  depth above.

## Pair 2: GO-4043 o001 (2024-05-11) → DD 12774 o002 (2026-07-06), F115W + F444W (`o4043_o12774/`)
- **Positive control failed under the 2-band rule:** Atalanta S = +2.5 / +2.4 (F115W, A / B) and +15.4 / +19.0
  (F444W); not hits. Forced photometry of B (`vet/lightcurves.ecsv`): F115W 13.7 ± 16.7 nJy in DD 12774, F444W
  87.4 ± 7.2 nJy; whether F115W misses it because the band is shallow (F115W forced error 16.7 nJy vs 13.1 nJy in
  F150W in the same visit) or because the SN is faint in F115W is not separated. Injection recovery is low: F115W 0.26 / 0.11, F444W 0.34 / 0.10 at 27.5 / 28.0 AB. Per the PLAN this
  pair is not trusted for a null; its hits are recorded but give no limit.
- 38 hits; all-hits morph: 1 pass, 22 extended, 11 edge, 3 dipole, 1 rosette. The frozen rule gives 0 `vet` (`image_distances_all.json`);
  the two hits ≤ 1″ from an image (h018, h020 near 12.1) are diffraction-spike stripes on the sheet.
- **h000** (the one pass; 1.81″ from image 6.1; S 469): a point source with a full PSF only in the 2024 GO-4043
  mosaics. Forced photometry in all 19 covering mosaics (`vet/lightcurves.ecsv`): F090W 2579 nJy, F115W 2688,
  F444W 258 on MJD 60441.7–60441.9, versus 7–49 nJy in F090W–F444W in 2022, 2026-06 and 2026-07. It holds the same sky position in three
  filters taken over ≈ 0.2 d. It lies 0.081″ from AT 2024jfz (`h000_tns.json`; TNS position as recorded on #237,
  not re-queried). **Known transient, not a candidate.** Not checked: whether it is multiply imaged (counter-images
  in the lens model). Under the post hoc 2.0″ radius it is the one pass in this pair (1.81″ from 6.1).

## h061 follow-up (PLAN.md addendum 1, 2026-10-16): undetermined, not a candidate
- **Frozen test** (`h061_spike/`, `scripts/stripe_phot.py`, rule = unit-tested `classify_excess`): the strongest
  stripe through h061 in the 2024 GO-4043 F444W mosaic is at PA 65.5° E of N (`summary.json`); after subtracting
  its across-stripe profile in every F444W mosaic, 2024 is 64.3 ± 5.0 nJy versus 43.8 ± 10.7 (2022),
  51.7 ± 11.8 (2026-06) and 33.3 ± 14.0 (2026-07) (`stripe_phot.ecsv`); z = 2.39 (`excess_z`) → **`undetermined`**
  (2 ≤ z < 3).
- **The injection check is weak** (found after the run): recovered / alone = 1.000 by construction, because the
  profile excludes |along| < 0.8″, so a source never enters its own model; it tests only that, not the stripe
  removal (`summary.json` `injection.note`).
- **Post hoc, the sheet** (`h061_spike/sheet.png`, looked at): h061 lies 4.70″ from a bright star
  (`h061_ctrl/near_star.json`) whose diffraction spikes cross the box at a different angle in each visit; in 2024
  one spike crosses h061.
- **Post hoc, calibrated controls** (`h061_ctrl/`, `hit_lightcurve_controls.py pick` + `compare`, 60 compact
  sources from the 2022 F444W box): 2022 → 2024 F444W z_cal 5.04, control tail fraction 0.033
  (`control_z.ecsv`; one of the two controls that extreme, c57, is h061's own knot). Split by distance to the star
  (`near_star.json`, c57 excluded): within 7″, 4 of 17 controls (0.24) have |z| ≥ h061's 14.2 (median z −10.2,
  max |z| 26.0); beyond 7″, 0 of 42 (max |z| 9.1). The near-star controls mostly fade (negative z) while h061
  brightens, so the sign is not matched; only the size is.
- **Post hoc, along the 2024 spike line** (`h061_spike/line_profile.ecsv`, `h061_ctrl/along_line.json`):
  2024 − 2022 raw flux at 0.4–2.0″ from h061 has median +5.0 nJy, robust std 5.4, but 7 of 34 positions (0.21)
  reach |Δ| ≥ h061's 35.6 nJy (max 66.7): the stripe is not uniform along its length at that level.
- **Same visit, other bands:** F090W 2022 → 2024 z_cal 0.19 (`control_z.ecsv`); F115W has no 2022 epoch.
- **Outcome:** h061 is **undetermined** under the frozen rule. Post hoc, its F444W change is matched in size by
  sources near the same star and by positions along the same 2024 spike, so it is not distinguishable from
  bright-star spike systematics with these tests; **not a candidate**. Reopen if a model of the star's 2024 F444W
  PSF (spikes included) leaves a residual at h061, or if image 12.1's counter-images show a change.

## Regenerate
Commands are in `o001_o057/summary.json`, `o4043_o12774/summary.json` and `h061_box/summary.json`; then
`arc_diff_morph.py --hits … --maps … [--all-hits] --out …`, `smacs_arcdiff_images.py MORPH --known
known_controls.csv --out …`, `arc_diff_sheet.py MAPS MORPH OUT [RA,DEC …]` and
`hit_lightcurve.py --pos h000 110.83411267801681 -73.45173415397286 --pos Atalanta_B 110.823747 -73.452975
--query-only h000 --out o4043_o12774/vet` (and `--pos h063 … --pos h061 … --pos Atalanta_A …
--query-only h063 --query-only h061 --out o001_o057/vet`). The maps (`--maps`, written to the session scratchpad,
shown as `<scratchpad>` in the stored commands) are not stored. h061 follow-up: the command in
`h061_spike/summary.json`; `hit_lightcurve_controls.py pick --uri <2022 ERO F444W i2d> --ra 110.82204720838848
--dec -73.44918269902084 --n 60` (`h061_ctrl/pick.txt`), `hit_lightcurve.py --pos h061 … <pick.txt positions>
--query-only h061 --out h061_ctrl`, `compare` (pairs in `control_z.ecsv`), then `scripts/h061_near_star.py
h061_ctrl`. Only the stamps of h061 and the four near-star controls with the largest |z| are kept.
