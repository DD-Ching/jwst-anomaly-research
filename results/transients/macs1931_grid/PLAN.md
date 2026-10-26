# MACS J1931 VENUS 6882 o038 → 12299 o001: D-027 two-epoch grid (frozen plan, 2026-10-26)

Written and pushed before any pixel was read. Pair from the epoch inventory
(`results/transients/epoch_inventory/RESULT.md`): VENUS 6882 o038 (2025-09-01) → GO 12299 o001 (2026-09-28,
public the same day), 392.6 d, the longest-baseline unused MACS J1931 pair. 12299 targets SN Eos image 101.1,
≈ 8″ from the VENUS centre, so the pair covers the strong-lensing core (a lens-field two-epoch search, not only a
flanking-field rate). The 9493 o001 epoch (2025-10-08, 773 s per band) is not used here.

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in ≥ 2 of F277W, F356W, F444W
  between the epochs (supernova, AGN, caustic-crossing star, or anything unexplained).
- Ordinary, tested in this order for every hit group: diffraction spikes / star wings (roll change), persistence,
  detector-fixed artefacts, proper motion (stars), mosaic edges / chip gaps, galaxy-edge background, AGN or
  supernova in a host (ordinary astrophysical transients are a result, not an anomaly).
- What changes the plan: any hit surviving the instrumental tests goes to `/vet-candidate` and the PR gets
  `needs-human`; a field with 0 survivors gives only the recovery-weighted sensitivity below.

## Run (tool unchanged)
`scripts/transient_grid.py` on the three shared LW bands, pipeline `_cat.ecsv` of both epochs and all three bands,
defaults (`--min-sigma 5 --min-dmag 0.3 --min-bands 2`, recentre in origin, mixed signs allowed),
`--inject 500 --inject-frac 1.0 0.5` (seed 0). Thresholds are the D-027 ASSUMPTIONs used on every earlier grid.
12299's F356W / F444W exposures are short (601 s / 515 s vs 1245–1675 s in VENUS), so epoch 2 sets the depth there.

## Gates
- **G1 injection:** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else outcome undetermined (U).
- **G2 noise:** each band's noise scale ≤ 2.0 (earlier grids 1.18–1.44), else U.

## Known case (reported, not a gate)
SN Eos 101.1 / 101.2 (positions in `results/transients/sneos/summary.json`). Its 2025 S/N was ≈ 4 in at least one
band (`results/transients/sneos_ratio/RESULT.md`), so a 5σ grid need not recover it; the row nearest each image
(≤ 0.5″) is reported with its per-band significance. A hit there is the known SN, not a candidate.

## Outcomes
- **N:** G1, G2 pass and every hit group has an ordinary explanation (contact sheets looked at: e1, e2, difference
  per band).
- **C:** ≥ 1 group without an instrumental explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails.

## Data
Six `_i2d.fits` (F277W/F356W/F444W × 2 epochs) streamed from MAST into the gitignored data root and deleted after
the run (cloud-disk decision); the whole mosaic is needed because the grid touches every row (tool docstring).
Sizes are logged in the result. Catalogues: the pipeline `_cat.ecsv` files.
