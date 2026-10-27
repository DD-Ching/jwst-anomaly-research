# Abell 1689 1840 o019 → 5782 o001: D-027 two-epoch grid (frozen plan, 2026-10-27)

Written and pushed before any pixel was read. Pair from the extended epoch inventory
(`results/transients/epoch_inventory_ext/RESULT.md`): GO 1840 o019 (MJD 59990.16, 2023-02-25) → GO 5782 o001
(MJD 60860.01, 2025-07-12), 869.8 d, core overlap 2.80 arcmin² in F277W / F356W inside the 1′ disc. First two-epoch
run on Abell 1689.

## Why this pair (information per hour)
The handoff named 5782 o001 → 6882 o048, but that pair shares only F150W / F356W (one SW mosaic; 2.0 / 1.9 arcmin² of
core). This pair shares two wide LW bands (F277W, F356W) over a larger core (2.80 arcmin²), and neither epoch is
VENUS, so it is not another GO-vs-VENUS depth mismatch like #576 / #577. No known transient in this core is pinned
here, so the run has no positive control beyond the injections. It is a transient search through a lensing core, not
a lens-model test: no Abell 1689 lens model is pinned here.

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in both F277W and F356W between the
  epochs (supernova, AGN, caustic-crossing star, or anything unexplained).
- Ordinary, tested in this order for every hit group: diffraction spikes / star wings (roll change), persistence,
  detector-fixed artefacts (single-exposure coverage, `--coverage` counts), the fixed-offset bright-star image of
  #576 (8.33″ at 120.1° measured at `PA_APER` 75.12°, rotated by this epoch's `PA_APER` − 75.12°), proper motion
  (stars), mosaic edges / chip gaps, galaxy-edge background, AGN or supernova in a host (ordinary transients are a
  result, not an anomaly).
- What changes the plan: any hit group surviving the instrumental tests goes to `/vet-candidate` and the PR gets
  `needs-human`; with 0 survivors and G1 / G2 passing, only the recovery-weighted sensitivity is reported.

## Bands
F277W and F356W, both epochs. `--min-bands 2`: both bands must pass.

## Run (tool as on main, with `--coverage`)
`scripts/transient_grid.py --band F277W <e1> <e2> --band F356W <e1> <e2>`, the four pipeline `_cat.ecsv` (epoch 1
F277W, F356W; epoch 2 F277W, F356W; `--epoch2-catalogs 2 3`), defaults (`--min-sigma 5 --min-dmag 0.3
--min-bands 2`, recentre in origin, mixed signs allowed), `--inject 500 --inject-frac 1.0 0.5` (seed 0),
`--coverage`. Thresholds are the D-027 ASSUMPTIONs used on every earlier grid.

## Gates (D-027, unchanged)
- **G1 injection:** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else outcome undetermined (U).
- **G2 noise:** each band's noise scale ≤ 2.0, else U.
- Reported alongside, not a gate (#578 found that G1 mixes two quantities): recovery split by epoch-1 exposure count
  (`summary.json`) and by nominal vs floor-set injected flux (from `injections.ecsv`).

## Outcomes
- **N:** G1, G2 pass and every hit group has an ordinary explanation (contact sheets looked at: e1, e2, difference
  per band).
- **C:** ≥ 1 group without an instrumental or ordinary explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails.

## Data
Four `_i2d.fits` and four `_cat.ecsv` streamed from MAST into the gitignored data root and deleted after the run
(cloud-disk decision). The field-wide grid reads every row of the mosaic (as #571 / #573 / #576 / #577); mosaic bytes
are recorded in the result.
