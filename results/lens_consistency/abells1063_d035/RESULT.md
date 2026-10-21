# Abell S1063 D-035 revisit: CATS v4.1 map with the Chan+2020 redshifts (result)

Plan: `PLAN.md` (frozen before the run). Script: `scripts/s1063_d035.py`. Numbers: `summary.json`, per image
`imageplane.ecsv` (model_prediction). Post hoc: `zcompare_posthoc.json` (`scripts/s1063_d035_zcompare.py`).

- **Route:** the CATS v4.1 files load from a cloud session once the trust store holds Mozilla's roots (certifi has
  emSign Root CA - G1; the cloud bundle does not). TLS verification stays on; the files match their sha256 pins.
- **Frame:** 49 of 117 Chan images have a CATS `arcs.txt` image within 1″, median separation 0.054″ (no frame offset).
- **R1 (whole list, ≤ 0.72″):** fails on both subsets. S (12 z_spec systems, 36 images): rms 18.5″, all solved.
  A (117 images): rms 11.7″, 116 solved.
- **R2 (per system, ≤ 0.5″, S only):** no system passes. Five S systems (1, 2, 10, 11, 52) get a single predicted
  image from the mean back-traced source (`n_predicted` 1; rms 14–37″). The other seven have rms 0.80–1.83″
  (`per_system`).
- **Frozen outcome:** **map only**; the S1063 image list stays gated off and `images` is not run on it.
- **Post hoc (deviation from the plan, diagnostic only):** for the 33 matched z_spec images, the CATS `arcs.txt`
  redshift equals Chan's MUSE z_spec within 0.005 (`zcompare_posthoc.json`; all 72 CATS images carry z > 0). So the
  missing-redshift reason D-035 gives for S1063 does not explain these systems: the same z and positions give rms
  0.8–1.8″ or a single predicted image. Not yet separated: the mean-source estimator (Lenstool optimizes the
  source; a mean of traced positions can fall outside the caustic for near-critical systems), the map resolution
  (D-035's κ difference is largest for S1063, 1.8e-2), and the map normalization.

## Next
- Check the route on a cluster where it should work (MACS1149 passed D-035 at 0.67″ / 0.63″) with the systems that
  fail here replaced by a source fitted in the image plane; only a route that reproduces ≤ 1.5 × quoted on S1063's
  own z_spec systems reopens D-035. Until then, S1063 counter-image tests stay map-only.
