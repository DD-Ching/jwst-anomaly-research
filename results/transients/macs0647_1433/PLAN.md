# MACS0647 1433 o010 → o020: D-027 two-epoch grid, F200W only (frozen plan, 2026-10-26)

Written and pushed before any pixel was read. Pair from the frozen epoch-inventory rule
(`results/transients/epoch_inventory/RESULT.md`): program 1433 o010 (2022-09-23, MJD 59845.6) → o020 (2023-01-08,
MJD 59952.7), 107 d, same pointing centre (101.9583, +70.2486), sharing only F200W. Both epochs cover the
strong-lensing core of MACS J0647.7+7015, so this is the first two-epoch test on its multiply imaged systems
(Meena et al. 2023 list, `results/lens_consistency/macs0647_inputs/meena23_images.ecsv`, 72 images).

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in F200W between the epochs.
- **H_lens (what would matter for the mission):** a change on a catalogued multiple image (within 0.5″ of a
  Meena+2023 position) without the same change on its counter-images (expected for an intrinsic variable only when
  the delay is long or the variability fast; read with the model delay, not as an anomaly by itself).
- Ordinary, tested in this order for every hit group: diffraction spikes / star wings (roll change between the
  visits), persistence, detector-fixed artefacts (one epoch only), proper motion (stars), mosaic edges / chip
  gaps, galaxy-edge or bright-core residuals, AGN or supernova in a host (an ordinary transient is a result, not an
  anomaly).
- What changes the plan: a hit group surviving the instrumental tests goes to `/vet-candidate` and the PR gets
  `needs-human`; with 0 survivors only the recovery-weighted sensitivity below is reported.

## Run 1: catalogue grid (tool unchanged)
`scripts/transient_grid.py --band F200W <o010 i2d> <o020 i2d> --catalog <o010 cat> --catalog <o020 cat>
--epoch2-catalogs 1 --min-bands 1 --inject 500 --inject-frac 1.0 0.5 --seed 0`, other options at their defaults
(`--min-sigma 5 --min-dmag 0.3`, recentre in origin, mixed signs). One shared band, so `--min-bands 1` as in the
single-band Abell 2744 grid (D-027 addendum 2026-10-15j). Thresholds are the D-027 ASSUMPTIONs.

## Run 2: difference-image pass on the core (D-027 addendum 2026-10-15j requires it before an SN-on-host null)
`scripts/arc_diff.py --band F200W <o010 i2d> <o020 i2d> --ra 101.96636 --dec 70.24716 --half-arcsec 55
--min-bands 1 --kernel-px 4 --inject 300`, other options at their defaults, centred on the median position of the
72 Meena+2023 images. The box holds the multiple-image region (images span 101″ × 83″).

## Gates
- **G1 (run 1):** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else run 1 is undetermined (U).
- **G2 (run 1):** noise scale ≤ 2.0, else U.
- **G3 (run 2):** recovery of the 27.5 AB injections ≥ 0.70, else run 2 is U (the Sunrise box gave 0.83–0.86).

## Outcomes
- **N:** gates pass and every hit group (run 1) and hit (run 2) has an ordinary explanation (contact sheet: e1, e2,
  difference, looked at before any conclusion).
- **C:** ≥ 1 group without an instrumental explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails (the other run is still reported).

## Data
Two `_i2d.fits` (1.29 / 1.30 GB, 2.6 GB in all) streamed from MAST into the gitignored data root and deleted after
the run (cloud-disk decision): the catalogue grid touches every row, and F200W is the only shared band. Reason
logged here, as the cloud-disk decision asks. Catalogues: the pipeline `_cat.ecsv` (2.1 / 2.8 MB). Run 2 reads its
box from the same local mosaics.

## Amendment (2026-10-26, before any pixel was read)
The first version set `--half-arcsec 30` and S3 reads. Computing the box centre from the Meena+2023 table showed the
images span 101″ × 83″, so a 60″ box would miss part of the multiple-image region: the half-width is now 55″, and
run 2 reads the mosaics already streamed for run 1. No other change.
