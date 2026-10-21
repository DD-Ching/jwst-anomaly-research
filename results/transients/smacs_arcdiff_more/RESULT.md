# SMACS 0723 arc difference imaging, remaining arc boxes (ERO o001 → VENUS o057) — result

Seven 24″ boxes cover all 45 Mahler+22 ICLv2 images outside the system 5 box of `../smacs_arcdiff/` (PLAN.md,
written before the run; addenda 1–3 written before each follow-up step). Pair jw02736-o001 (2022-06-07) →
jw06882-o057 (2026-06-05), F150W + F444W, same `arc_diff.py` settings as pair 1 there; commands in each
`bN/summary.json`. Numbers are `derived`; thresholds are ASSUMPTIONs. Sheets were looked at before writing this.

## Validation
- **b2 known case passes:** SN Atalanta A = h021 (S +58.6 F150W / +47.6 F444W), B = h023 (+53.1 / +53.4), both
  morph `pass`, clean point brightenings on the sheet. The control is off-centre in a new box, so it also checks
  the box placement.
- **Injection gate passes in all 7 boxes** (≥ 0.50 at 27.5 AB in both bands). Recovery at 27.5 / 28.0 AB:

| box | images | hits | F150W | F444W | rule hits (pass, ≤ 2″ of an image) |
|---|---|---|---|---|---|
| b1 | 17 (4, 5.3, 8.1, 9, 10.3, 11, 12.3, 15.2/3, 17.2, 18.3, 19.3, 21.2/3) | 21 | 0.84 / 0.52 | 0.92 / 0.76 | h019 |
| b2 | 9 (1.1, 2.1, 3.1, 6.1/2, 12.1, 13.1, 15.1, 17.3) | 73 | 0.81 / 0.60 | 0.67 / 0.52 | h032 (+ Atalanta) |
| b3 | 7 (1.2/3, 2.2/3, 3.3, 6.3, 17.1) | 9 | 0.81 / 0.47 | 0.91 / 0.70 | — |
| b4 | 7 (4.3, 7, 25) | 9 | 0.95 / 0.75 | 0.95 / 0.65 | — |
| b5 | 3 (8.2, 26) | 8 | 0.94 / 0.52 | 0.97 / 0.71 | — |
| b6 | 1 (13.3) | 17 | 0.65 / 0.19 | 0.95 / 0.50 | — |
| b7 | 1 (21.1) | 18 | 0.98 / 0.77 | 0.98 / 0.71 | — |

All other morph passes (b1 5, b5 3, b6 3, b7 4) lie 2.6–22″ from any catalogued image and are recorded only; on
the b1 sheet the off-image passes are negative residuals on compact galaxy cores.

## Rule hits
- **b1 h019** (1.51″ from 17.1, arc_diff sign −). Forced photometry (`b1/vet/`): +73.6 ± 16.1 nJy F150W, +33.4 ± 5.1
  F444W on a 1.71 / 0.99 μJy source → control test (Addendum 1, 25 controls of 0.5–5 μJy, `b1/h019_ctrl/`):
  z_cal 2.08 (F150W), 1.09 (F444W). **Frozen outcome: undetermined** (between 2 and 3). Post hoc: its F150W flux
  ratio 1.043 equals the control median 1.042 (robust σ 0.019); F444W 1.034 vs 0.993 (σ 0.074). The change is the
  common ≈ 4 % flux-scale difference between the two mosaics. The arc_diff difference, which fits the flux scale,
  has the opposite sign. Not a candidate.
- **b2 h032** (0.09″ from image 1.1, arc_diff sign −). Forced photometry (`b2/vet/`): F150W ratio 1.053, F444W 1.153 →
  control test (34 controls, `b2/h032_ctrl/`): z_cal 0.44 / **5.65** → open. Step R (Addendum 3): the catalogue
  registration between the mosaics is (0.009, 0.020)″ F150W and (0.002, 0.022)″ F444W at h032 (`b2/h032_reg/`,
  M0 held-out rms ≈ 0.009″). Corrected photometry: F150W ratio 1.020 (z_cal 0.02), F444W 0.943 (z_cal −2.23).
  **Frozen outcome: undetermined.** A 0.02″ correction moved the F444W ratio by 0.21, so the aperture flux on the
  arc's edge depends on registration about ten times more steeply than on the compact controls. Step C
  (`b2/sys1_counter/`, corrected): catalogue 1.1 itself gives F444W 1.063, counter-images 1.2 0.984 and 1.3 0.931,
  i.e. ±7 % scatter on unchanged arcs. By eye (sheet) the F444W S map at h032 is an elongated +/− pair across the
  arc's width, the pattern of a sub-pixel shift. Post hoc reading: registration residual on a steep arc gradient.
  Not a candidate.

## Outcome
No candidate in the remaining SMACS 0723 arc boxes for this pair. Both rule hits are undetermined as frozen;
post hoc each has an ordinary explanation (flux-scale offset; registration on an arc edge). Point-source depth
on the arcs: injection recovery ≥ 0.50 at 27.5 AB in every box (≥ 0.65 F150W, ≥ 0.67 F444W), 28.0 AB in 5 of 7
for both bands. Not a limit below the injection depth.

## Limitations
- Aperture forced photometry on extended arcs is dominated by registration (h032), so it cannot confirm or reject
  a faint change on an arc ridge; the arc_diff S maps and morphology are the better test there.
- Control z from compact sources does not carry the gradient sensitivity of arc positions.
- One epoch pair; transients shorter than the sampling or outside 2022-06 → 2026-06 are not covered.

## Regenerate
`bN/summary.json` holds each `arc_diff.py` command (maps in the session scratchpad, not stored); then
`arc_diff_morph.py --hits bN/hits.ecsv --maps <maps> --all-hits --out bN/morph_all.json`,
`smacs_arcdiff_images.py bN/morph_all.json --known ../smacs_arcdiff/known_controls.csv --out
bN/image_distances_all.json`, `arc_diff_sheet.py <maps> bN/morph_all.json bN/sheet.jpg` (b2's sheet rescaled to
75 % to stay under 1 MB). Follow-ups: `hit_lightcurve.py` / `hit_lightcurve_controls.py pick|compare` /
`reg_offsets.py` with the parameters in PLAN.md addenda 1–3; only the hits' stamps are kept.
