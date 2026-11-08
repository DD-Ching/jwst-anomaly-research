# `arcs`: RXC J0600 581 under the RELICS GLAFIC v3 realisations (exploratory, 2026-11-08)

Plan: `PLAN.md` (frozen and committed before any realisation value was read). Script
`scripts/arcs_shear_error_581g.py`; resampling and gate functions unit-tested before the run
(`tests/test_arcs_shear_error_581g.py`); rule functions are addendum 13's, imported unchanged. Rule numbers are in
`summary.json`; the post hoc κ numbers are in `posthoc_kappa.json`. Maps streamed into memory from the RELICS HLSP
(DOI 10.17909/T9SP45; sha256 per file in both JSON files); only the 3 κ gate maps touch disk, briefly. Observed axis
and photo-z from the stored GLAFIC run (`../rxc0600-relics/rxc0600-relics-glafic/arcs.ecsv`); predicted axes and κ
`model_prediction`; classes and labels `derived`; thresholds ASSUMPTIONs (`summary.json` `params`).

## Frozen outcome: **undetermined** (route gate V2 fails)
| Gate / quantity | Value | Threshold | Result |
|---|---|---|---|
| V1 best fit pixel × 0.1″ vs arcsec maps (full map) | max 3.8 × 10⁻⁶″ | ≤ 1 × 10⁻³″ | passes |
| V2 best fit resampled to the 0.4″ realisation grid, 12 stored strong-shear background rows | classes equal; max Δoffset_min 16.29° | equal, ≤ 2° | **fails** |
| V3 median-over-realisations α vs best fit at the 12 rows | median rel. diff 0.0462 | ≤ 0.10 | passes |
| Map gate, κ, map000 / 049 / 098 | median rel. diff 4.3 / 2.8 / 4.1 × 10⁻⁵ | ≤ 0.05 | passes |
| Control retention (11 GLAFIC `aligned` rows) | median 1.0 (range 0.909–1.0) | ≥ 0.70 | passes |
| 581 `anti` share f | 1.0 (99 / 99; f_aligned 0.0) | robust ≥ 0.68 | not used: any failed gate → `undetermined` |

Secondary (not rules): offset_min 16 / 50 / 84 % 77.868 / 78.8 / 80.1°; g_min 0.241 / 0.264 / 0.296. Resampled best
fit: 581 `anti`, offset_min 79.9°, g_min 0.236. Stored GLAFIC run (0.1″ maps): 63.6°, g_min 0.394.

## Why V2 fails (post hoc diagnostics, not rules; every number in this section is post hoc)
- Per row (`summary.json` `validation.V2_rows_posthoc`, added after the first run; rules unchanged, second run gives
  identical realisation sha256 and outcome): 11 of 12 rows reproduce within 0.01° in offset_min and 0.003 in g_min.
  Only 581 moves: 63.62° → 79.91°, g_min 0.394 → 0.236.
- `posthoc_kappa.json`: the GLAFIC best-fit κ has a local peak of 5.24 at 0.2″ from 581 (κ 1.18 at the row; 2″-box
  median 0.80). For the other 11 rows the box maximum lies on the box edge (0.75–1.23, at 2.0–2.8″), so none has an interior peak. The
  Lenstool v1 best fit has no local peak there either (2.48 at the box corner, 2.83″, rising toward the cluster). The
  0.4″ realisation grid cannot resolve a galaxy-scale halo at 0.2″, so the realisations cannot test this row.
- Reading (hypothesis, not tested): GLAFIC places a galaxy-scale halo at 581, the usual way cluster members enter
  such models. 581 is a bright (m 21.2, `../rxc0600-relics/RESULT.md`) near edge-on disc with a bulge (contact sheet `../anti-bands/anti_bands.png`,
  looked at again). Its photo-z (0.70) puts it behind the lens, but the photo-z membership test failed its gate
  (addendum 10). If 581 is a member or sits on one, its axis is its own shape and not a lensing prediction.
- The local halo makes 581 *less* anti in GLAFIC (63.6° with it at 0.1″ vs 79.9° smoothed). Lenstool, without it, puts
  581 at 74.6–76.3° (16–84 %, `../anti-shear-error/RESULT.md`, addendum 13). The `anti` class therefore does not come from that halo.

## Reading (narrowest statement)
- The GLAFIC realisations cannot test 581 at 0.4″: the frozen label is `undetermined`, and the 99 / 99 `anti` share
  is not a robustness result.
- New ordinary route for 581: the GLAFIC member list (or the RELICS GLAFIC model's input catalogue), if public,
  decides whether the model treats 581 as a member. Spectroscopy still decides membership directly (none public,
  addendum 15).
- Not a lensing anomaly claim and not a candidate. The open set stays 581 (Lenstool `robust`, GLAFIC `undetermined`)
  and 2157 (`robust`, membership `undetermined`).
