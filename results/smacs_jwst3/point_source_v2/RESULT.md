# SMACS point-source screen v2 (locus 20–24 mag): not validated (V1a and V1b fail)

Plan: `PLAN.md` (window, gates and decision rule pushed before the run; `decide_v2` and `held_out_by_group`
unit-tested before the run). Run `20261107T141219Z-86cb61bb` of `configs/smacs_jwst3_gaia_only.yaml` at main
`ae36074` (2938 ranked, count not stored; `f200w_1816` rank 14 and `f200w_384` rank 19, as #765 / #766). Measurement:
`scripts/smacs_point_source.py <run> --phot <DJA v7.4 phot> --out . --plan v2` → `summary.json` (phot
sha256-checked). Contact sheet of the run looked at.

## Numbers (derived; ASSUMPTION thresholds from the plan)
- Locus: 101 Gaia-split stars in ref_mag [20, 24]; `r_psf` 2.1293 px, 1.4826 MAD 0.1988 px, cut 2.7257 px
  (v1: 2.9605 px, `../point_source/summary.json`). Locus 5th / 50th / 95th percentiles 1.934 / 2.129 / 2.956 px.
- **V1a fails:** parity held-out recovery 0.6909 / 0.9348 (gate ≥ 0.90 both ways).
- **V1b fails:** spatial held-out recovery (split at the locus median `ra` 110.770329°) 0.9608 (east cut tested on
  west) / 0.76 (west cut tested on east).
- **V2 passes:** 0 / 19 matched ICLv2 images at or below the cut (60 listed).
- **V3 descriptive** (5 matched stars with ref_mag in (24, 27], fewer than 10): 0.8 point-like; no gate.
- Verdict (`decide_v2`): "not validated (V1a, V1b fail)"; `targets_labelled` false. By the plan, `f200w_1816`
  (2.276 px, ref_mag 25.897) and `f200w_384` (1.906 px, 26.305) stay for normal vetting.
- Descriptive: point-like galaxy-stratum rows by ref_mag 0 / 13 (18–22), 2 / 90 (22–24), 59 / 213 (24–26),
  44 / 256 (26–27), 61 / 438 (27–28), 68 / 1024 (28–31) (`galaxy_point_like_by_ref_mag`).

## Pre-run expectation that was wrong
The plan expected V1a to pass because only 5 of v1's outliers lay in [20, 24]. Dropping the brighter stars
narrowed the MAD (v1 0.2681 in `../point_source/summary.json` → 0.1988 px), so the cut fell to 2.7257 px and 11 locus stars lie above it
(`post_hoc_locus_outliers`).

## Post hoc (deviation: added after V1a and V1b failed)
- Locus radius by ref_mag (`post_hoc_locus_by_ref_mag`): median 2.649 px with 8 / 19 above the cut at 20–20.7;
  2.162 px, 0 / 22 at 20.7–22; 2.032 px, 3 / 60 at 22–24. The broadened core thus extends to about 20.7 mag in
  this sample (a saturation or non-linearity reading is a hypothesis; no pipeline flag was tested).
- Outliers by V1b half (`post_hoc_outliers_east_west`): 9 east, 2 west. Whether that reflects the halves' magnitude
  mix or a PSF change across the mosaic was not tested.

## Reading
- A constant-width locus (median + 3 × MAD) on DJA `flux_radius` fails held-out star recovery in both v1 (19–24)
  and v2 (20–24). Moving the window again (e.g. to 20.7) would be a third post-hoc choice on the same 101 stars,
  so it is not done here.
- Independently of the gates, only 5 matched Gaia-split stars lie in (24, 27], all at ref_mag ≤ 24.657, while the
  targets are at 25.9 and 26.3; this field cannot validate a size label at their magnitude (V3 needs ≥ 10). The size route does not
  decide them; they stay for normal vetting by another route.
- The catalogued images stay well above any tried cut (V2 0 / 19 in both versions): size does not remove known
  lensed images, but that alone does not make a screen.

## Limitations
- One field, one catalog version; `flux_radius` from DJA's detection image, not per band.
- The faintest matched Gaia-split star in (24, 27] has ref_mag 24.657 (`v3_faint_stars`); faint-end behaviour is
  untested.
- Data: DJA phot streamed to scratch and deleted after the run (Cloud disk decision).
