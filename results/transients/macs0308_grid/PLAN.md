# MACS J0308 5293 o005 → VENUS 6882 o015: D-027 two-epoch grid (frozen plan, 2026-10-27)

Written and pushed before any pixel was read. Pair from the extended epoch inventory
(`results/transients/epoch_inventory_ext/RESULT.md`): GO 5293 o005 (2025-02-04) → VENUS 6882 o015 (2025-08-23),
200.2 d, both footprints covering the full 1′ core disc. First two-epoch run on MACS J0308.9+2645.

## Why this pair (information per hour)
The largest unused core overlap in the inventory, and it carries a known transient: VENUS reported the strongly
lensed core-collapse SN "SN Ares" from this cluster's imaging (discovered 2025-08-25, JWST programme 9478 target
list, `www.stsci.edu/jwst-program-info/download/jwst/pdf/9478/`, read 2026-10-27; not a refereed source). The
first image ("Ares", RA 47.2225746°, Dec 26.76629°) is ≈ 40″ from the cluster centre, inside both footprints.
So the run is a positive-controlled transient search through a lensing core. It is not a lens-model test:
no MACS J0308 lens model is pinned here.

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in both of F300M and F410M between
  the epochs (supernova, AGN, caustic-crossing star, or anything unexplained).
- Ordinary, tested in this order for every hit group: diffraction spikes / star wings (roll change), persistence,
  detector-fixed artefacts (single-exposure coverage), proper motion (stars), mosaic edges / chip gaps,
  galaxy-edge background, AGN or supernova in a host (ordinary astrophysical transients are a result, not an
  anomaly).
- What changes the plan: any hit group surviving the instrumental tests and not SN Ares goes to `/vet-candidate`
  and the PR gets `needs-human`; with 0 survivors only the recovery-weighted sensitivity below is reported.

## Bands
The only LW filters in both epochs are the medium bands F300M and F410M (5293 has no F277W / F356W / F444W).
SW mosaics (F115W / F150W / F200W, shared) are not used: about 1.7 GB per mosaic per epoch (as in the Abell 2744
flanking grid). With two bands, `--min-bands 2` means both must pass: a stricter hit rule than the three-band
grids, stated as a limit. 5293 has 1095 s per band, VENUS 2061 s, so epoch 1 sets the depth.

## Run (tool unchanged)
`scripts/transient_grid.py --band F300M <e1> <e2> --band F410M <e1> <e2>`, the four pipeline `_cat.ecsv` (epoch 1
F300M, F410M; epoch 2 F300M, F410M; `--epoch2-catalogs 2 3`), defaults (`--min-sigma 5 --min-dmag 0.3
--min-bands 2`, recentre in origin, mixed signs allowed), `--inject 500 --inject-frac 1.0 0.5` (seed 0).
Thresholds are the D-027 ASSUMPTIONs used on every earlier grid.

## Gates
- **G1 injection:** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else outcome undetermined (U).
- **G2 noise:** each band's noise scale ≤ 2.0, else U.

## Known case (reported; a positive-control check, not a gate)
- SN Ares image "Ares": the row nearest the listed position (≤ 0.5″) and its per-band significance and Δm.
  Expected: brighter in epoch 2 (epoch 2 is 2 d before the reported discovery date, so probably the discovery data; epoch 1 is 6.5 months earlier),
  so a hit if it was bright enough in both LW medium bands and catalogued. If it is not a hit, the reason is
  reported (no row within 0.5″, or below 5σ / 0.3 mag in a band).
- "Ares_Image2" (47.2342638°, 26.76035°) and "Ares_Image3" (47.2386733°, 26.76029°) are listed targets of the same
  programme without stated status; their nearest rows are reported too. Nothing is concluded from them.

## Outcomes
- **N:** G1, G2 pass and every hit group has an ordinary explanation or is SN Ares (contact sheets looked at: e1,
  e2, difference per band).
- **C:** ≥ 1 other group without an instrumental explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails.

## Data
Four `_i2d.fits` (F300M / F410M × 2 epochs; MAST sizes 324.4 MB × 2, 306.5 MB × 2, 1.26 GB in all) and four
`_cat.ecsv`, streamed from MAST into the gitignored data root and deleted after the run (cloud-disk decision;
the grid touches every row of the mosaic). Each single file is > 200 MB: the reason is the field-wide grid,
which reads the whole mosaic (as #571 / #573).
