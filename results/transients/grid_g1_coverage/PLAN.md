# Why did the MACS J0308 / J2129 grids fail G1? Injection recovery by epoch-1 exposure count (frozen plan, 2026-10-27)

Written and pushed before the tool change was run on any pixel. Inputs: the two frozen grids
`../macs0308_grid/` (#576, G1 0.518) and `../macs2129_grid/` (#577, G1 0.78).

## Hypothesis and ordinary explanations
- **H_cov:** the −0.75 mag injections are missed mainly at positions with one epoch-1 exposure (`CON` bit count 1),
  where the `ERR` array understates the true pixel noise less than the noise scale assumes, or where single-exposure
  artefacts inflate the band's robust noise scale.
- Alternatives: (a) recovery is low at every coverage (the noise scale, set by the whole grid, is too large for the
  ordinary rows); (b) misses follow the injected flux (faint hosts near the floor), not the coverage.
- What changes the plan: if recovery on rows with ≥ 2 epoch-1 exposures in every band is ≥ 0.85 at −0.75 mag in
  both grids, later grids report G1 on that subset and restrict the search to it (an amendment, new D-027
  addendum). If not, H_cov is rejected and the GO-vs-VENUS depth mismatch stays undiagnosed; no rule change.

## Tool change (`scripts/transient_grid.py`)
- `injections.ecsv`: one row per injection and fraction (pos_id, ra, dec, frac, injected flux and per-band pass sign,
  recovered, nominal), so a post hoc split needs no re-run.
- `--coverage`: the epoch-1 and epoch-2 exposure count (bits set in the mosaic `CON` extension) at every grid
  position, stored as `<band>_e1_n_exp` / `<band>_e2_n_exp`, and the injection recovery split by the minimum
  epoch-1 count over bands (1, ≥ 2) in `summary.json`.
- Same seed, thresholds and positions as the frozen runs, so the overall recovery must reproduce 0.518 / 0.78
  exactly (a check of the change). Hits are unchanged by construction.

## Run
Both frozen commands from the two `PLAN.md` files with `--coverage` added. Mosaics streamed from MAST into the
gitignored data root and deleted after each run (877 MB for J2129, the J0308 set as in #576).

## Outcome
- **S (supported):** both reproduce, ≥ 2-exposure recovery ≥ 0.85 in both grids.
- **R (rejected):** both reproduce, ≥ 2-exposure recovery < 0.85 in either.
- **U:** a reproduction check fails (tool change wrong) or too few injections (< 100) in a subset.
