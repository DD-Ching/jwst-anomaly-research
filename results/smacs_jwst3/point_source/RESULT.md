# SMACS point-source screen on DJA `flux_radius`: not validated (V1 fails); lensed images are kept (V2 passes)

Plan: `PLAN.md` (rule, validation gates and targets pushed before the run). Run `20261107T131039Z-cd7a1929` of
`configs/smacs_jwst3_gaia_only.yaml` at main `0d7ef22` (ranks reproduce #765: `f200w_1816` rank 14, `f200w_384`
rank 19). Measurement: `scripts/smacs_point_source.py <run> --phot <DJA v7.4 phot> --out .` → `summary.json` (phot
sha256-checked). Rule functions (`locus_cut`, `locus_mask`, `held_out_recovery`, `decide`) were unit-tested before
the run. Contact sheet of the run looked at.

## Numbers (derived; ASSUMPTION thresholds from the plan)
- Locus: 113 of 133 matched Gaia-split stars in ref_mag [19, 24]; `r_psf` 2.1563 px, 1.4826 MAD 0.2681 px, cut
  2.9605 px. Locus 5th / 50th / 95th percentiles 1.939 / 2.156 / 5.563 px.
- **V1 fails:** held-out recovery 0.7333 (cut from even `id`s, tested on odd) and 0.9057 (the reverse); the gate
  needs ≥ 0.90 both ways.
- **V2 passes:** 19 of the 60 listed ICLv2 images match DJA at 0.2″; 0 / 19 are at or below the cut (image radii 5th /
  50th / 95th percentiles 3.724 / 6.636 / 15.447 px).
- Verdict (`decide`): "not validated (V1 fails)". By the plan, both targets stay for normal vetting. Their values,
  for the record: `f200w_1816` 2.276 px, `f200w_384` 1.906 px (both below the cut; not a label, since the screen is
  not validated).
- Descriptive: point-like fraction of galaxy-stratum rows by ref_mag: 0 / 13 (18–22), 2 / 90 (22–24), 62 / 213
  (24–26), 48 / 256 (26–27), 82 / 438 (27–28), 107 / 1024 (28–31); 8 of the top 20 are below the cut (`top20`).

## Post hoc (deviation: added after V1 failed)
- The 14 locus stars above the cut (`post_hoc_locus_outliers`): 9 have ref_mag < 20 (radii 3.5–7.1 px), 3 lie at
  20.4–20.7 and 2 at 22.8 / 23.5. A reading (hypothesis, not tested): F200W saturation or a broadened core starts
  fainter than the plan's 19 mag, so the plan's locus window was too bright. How the outliers divide between the
  two parity halves was not stored, so the asymmetric recovery (0.733 vs 0.906) is not explained here.

## Reading
- The stellar locus core is narrow (MAD 0.2681 px around 2.1563 px), and all 19 matched catalogued images lie above
  the cut. The frozen star-recovery gate fails, however, so no screen is adopted on this run.
- The contact sheet shows `f200w_1816` (#14) and `f200w_384` (#19) as compact, isolated sources (visual, not stored).
  Neither is a lens-model residual candidate on current evidence, but that is not a screen result.

## Limitations
- V2 checks only catalogued images (19 matched, mostly extended arcs/knots); it does not bound faint compact lensed
  galaxies or high-z galaxies, which can also be PSF-sized at 26–28 mag.
- `flux_radius` comes from DJA's detection image, not per band; one field, one catalog version.
- Data: DJA phot streamed to scratch and deleted after the run (Cloud disk decision).
