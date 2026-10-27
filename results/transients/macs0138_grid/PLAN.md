# MACS J0138 2345 o005 → 6549 o003: D-027 two-epoch grid with SN Encore as positive control (frozen plan, 2026-10-27)

Written and pushed before any pixel was read. Pair from the extended epoch inventory
(`results/transients/epoch_inventory_ext/RESULT.md`, the only `core` pair on MACS J0138): GO 2345 o005 (MJD 60265.21,
2023-11-17, the SN Encore discovery epoch) → DD 6549 o003 (MJD 60317.22, 2024-01-08), 52.0 d, core overlap
1.004 / 0.983 arcmin² in F150W / F444W inside the 1′ disc.

## Why this pair (information per hour)
#576, #577 and #579 failed G1 and had no positive control. This pair contains a known transient: SN Encore, a
z = 1.95 lensed SN Ia in MRG-M0138, discovered in F150W of exactly this epoch 1 and followed by 6549 (Pierel et al.,
arXiv:2509.12301, §2; image 1a was between peak and the second IR maximum at discovery). Over 52 d (≈ 18 rest-frame
days) the SN images should fade. So this run tests the grid on a real transient on top of a bright lensed host, which
injections into catalogue positions do not. It is a transient search through a lensing core, not a lens-model test.

## Hypothesis and ordinary explanations
- **H_var:** a catalogued position changed by ≥ 0.3 mag at ≥ 5σ (calibrated) in both F150W and F444W.
- **Positive control (pre-declared):** a fade hit group on an MRG-M0138 image carrying SN Encore 1a or 1b. The SN sits
  ≈ 0.1″ from the host light (arXiv:2509.12301 §2.2), so the catalogue position is the host image. No SN coordinates
  are stored in this repository and none were found in the arXiv source; the host images are identified by eye
  on the contact sheet against that paper's colour figure (ASSUMPTION; identification is by eye, not by a frozen
  coordinate match). Recovered / not recovered is reported either way; a miss in F444W (rest ≈ 1.5 µm, near the
  second IR maximum) is possible and is not by itself a tool failure.
- Ordinary explanations, in this order, for every other group: diffraction spikes / star wings (roll change),
  persistence, single-exposure coverage (`--coverage`), the #576 fixed-offset bright-star image (8.33″ at 120.1° at
  `PA_APER` 75.12°, rotated by this epoch's `PA_APER` − 75.12°), proper motion, mosaic edges / chip gaps, galaxy-edge
  background, AGN or supernova in a host (ordinary transients are a result, not an anomaly), SN Requiem images
  (faded by 2016; image 2d predicted ≈ 2026–2027, after both epochs).
- What changes the plan: any group surviving these goes to `/vet-candidate` and the PR gets `needs-human`.

## Bands
F150W and F444W, both epochs (the only shared filters). `--min-bands 2`.

## Run (tool as on main, with `--coverage`)
`scripts/transient_grid.py --band F150W <e1> <e2> --band F444W <e1> <e2>`, the four pipeline `_cat.ecsv` (epoch 2
catalogues as `--epoch2-catalogs 2 3`), defaults (`--min-sigma 5 --min-dmag 0.3 --min-bands 2`), `--inject 500
--inject-frac 1.0 0.5` (seed 0), `--coverage`. Thresholds are the D-027 ASSUMPTIONs of every earlier grid.

## Gates (D-027, unchanged)
- **G1 injection:** recovery at −0.75 mag (frac 1.0) ≥ 0.85, else outcome U.
- **G2 noise:** each band's noise scale ≤ 2.0, else U.
- Reported alongside, not gates: recovery by epoch-1 exposure count, and by expected calibrated σ (as #579's post hoc;
  still post hoc here, because the G1 amendment is the owner's call).

## Outcomes
- **N:** G1, G2 pass and every group has an ordinary explanation (contact sheets looked at).
- **C:** ≥ 1 group without an ordinary explanation → `/vet-candidate`, `needs-human`.
- **U:** a gate fails. The positive-control result is reported in every outcome.

## Data
Four `_i2d.fits` and four `_cat.ecsv` streamed from MAST S3 into the gitignored data root and deleted after the run
(cloud-disk decision); mosaic bytes recorded in the result.
