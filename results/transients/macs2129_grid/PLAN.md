# MACS J2129 4903 o004 → VENUS 6882 o040: D-027 two-epoch grid (frozen plan, 2026-10-27)

Written and pushed before any pixel was read. Pair from the extended epoch inventory
(`results/transients/epoch_inventory_ext/RESULT.md`): GO 4903 o004 (MJD 60586.74, 2024-10-01) → VENUS 6882 o040
(MJD 60978.81, 2025-10-28), 392.1 d, core overlap 2.22–2.27 arcmin² per shared filter inside the 1′ disc. First
two-epoch run on MACS J2129.4-0741 (not RX J2129.7+0005, which the `arcs` / W1 runs used).

## Why this pair (information per hour)
The longest single-year baseline in the inventory with shared wide LW bands (F277W, F444W), so it avoids the
two-medium-band insensitivity that made the MACS J0308 grid fail G1 (#576). No known transient in this core is
pinned here, so the run has no positive control beyond the injections. It is a transient search through a
lensing core, not a lens-model test: no MACS J2129 lens model is pinned here.

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in both F277W and F444W between the
  epochs (supernova, AGN, caustic-crossing star, or anything unexplained).
- Ordinary, tested in this order for every hit group: diffraction spikes / star wings (roll change), persistence,
  detector-fixed artefacts (single-exposure coverage), the fixed-offset bright-star image of #576 (8.33″ at
  120.1° measured at `PA_APER` 75.12°, rotated by this epoch's `PA_APER` − 75.12°), proper motion (stars), mosaic
  edges / chip gaps, galaxy-edge background, AGN or supernova in a host (ordinary transients are a result, not an
  anomaly).
- What changes the plan: any hit group surviving the instrumental tests goes to `/vet-candidate` and the PR gets
  `needs-human`; with 0 survivors only the recovery-weighted sensitivity is reported.

## Bands
F277W and F444W, both epochs (`_i2d` 131.9 MB × 2 in 4903, 306.5 MB × 2 in 6882; 877 MB in all). SW F115W / F200W
are shared too but not used in the frozen run (two more mosaics per epoch; only for a post hoc check if a group
needs it). `--min-bands 2`: both bands must pass.

## Run (tool unchanged)
`scripts/transient_grid.py --band F277W <e1> <e2> --band F444W <e1> <e2>`, the four pipeline `_cat.ecsv` (epoch 1
F277W, F444W; epoch 2 F277W, F444W; `--epoch2-catalogs 2 3`), defaults (`--min-sigma 5 --min-dmag 0.3
--min-bands 2`, recentre in origin, mixed signs allowed), `--inject 500 --inject-frac 1.0 0.5` (seed 0).
Thresholds are the D-027 ASSUMPTIONs used on every earlier grid.

## Gates
- **G1 injection:** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else outcome undetermined (U).
- **G2 noise:** each band's noise scale ≤ 2.0, else U.

## Outcomes
- **N:** G1, G2 pass and every hit group has an ordinary explanation (contact sheets looked at: e1, e2, difference
  per band).
- **C:** ≥ 1 group without an instrumental or ordinary explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails.

## Data
Four `_i2d.fits` and four `_cat.ecsv` streamed from MAST into the gitignored data root and deleted after the run
(cloud-disk decision). The two 6882 mosaics are > 200 MB each: the field-wide grid reads every row of the mosaic
(as #571 / #573 / #576).
