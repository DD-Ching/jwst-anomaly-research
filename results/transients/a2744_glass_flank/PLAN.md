# Abell 2744 GLASS flanking pair 1324 o001 → o007: D-027 two-epoch grid (frozen plan, 2026-10-26)

Written and pushed before any pixel was read. Pair chosen by the frozen #547 rule (largest unused overlap,
`results/transients/a2744_epoch_footprint/RESULT.md`): GLASS 1324 o001 (2022-06-28/29) → o007 (2023-07-07), 373 d,
outside the 1′ strong-lensing core (low magnification; a general two-epoch rate, not a lens-model test).

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in ≥ 2 of F277W, F356W, F444W
  between the epochs (supernova, AGN, or anything unexplained).
- Ordinary, tested in this order for every hit group: diffraction spikes / star wings (roll change), persistence,
  detector-fixed artefacts, proper motion (stars), mosaic edges / chip gaps, galaxy-edge background, AGN or
  supernova in a host (ordinary astrophysical transients are a result, not an anomaly).
- What changes the plan: any hit surviving the instrumental tests goes to `/vet-candidate` and the PR gets
  `needs-human`; a field with 0 survivors gives only the recovery-weighted sensitivity below.

## Run (tool unchanged)
`scripts/transient_grid.py` on the three shared LW bands (SW bands not used: 1.7–1.8 GB per mosaic per epoch,
download budget), pipeline `_cat.ecsv` of both epochs and all three bands, defaults (`--min-sigma 5 --min-dmag 0.3
--min-bands 2`, recentre in origin, mixed signs allowed), `--inject 500 --inject-frac 1.0 0.5` (seed 0).
Thresholds are the D-027 ASSUMPTIONs used on every earlier grid.

## Gates
- **G1 injection:** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else outcome undetermined (U).
- **G2 noise:** each band's noise scale ≤ 2.0 (earlier grids 1.18–1.44), else U.

## Outcomes
- **N:** G1, G2 pass and every hit group has an ordinary explanation (contact sheets looked at: e1, e2, difference
  per band).
- **C:** ≥ 1 group without an instrumental explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails.

## Data
Six `_i2d.fits` (F277W/F356W/F444W × 2 epochs, 289–345 MB each, 1.9 GB in all) streamed from MAST into the
gitignored data root and deleted after the run (cloud-disk decision); the whole mosaic is needed because the grid
touches every row (tool docstring). Catalogues: the pipeline `_cat.ecsv` files.
