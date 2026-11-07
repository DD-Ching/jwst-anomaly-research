# Plan: what does the DJA eazy `z_phot = -1` (no-fit) flag mark among SMACS ranked sources? (TASKS Next M2)

Written 2026-11-07 before any fetch or measurement. Follows `results/smacs_jwst3/zphot_stratify/RESULT.md` (post hoc:
`z_phot = -1` rows are 4 / 20 of the top 20 vs 156 / 2938 ranked).

## Hypotheses
- H-star: the no-fit flag marks uncatalogued stars / compact non-galaxy sources missed by the Gaia-only star split.
- H-cov (ordinary, catalog effect): eazy did not fit because too few usable filters (footprint edge, masked or
  saturated photometry); the flag then marks a data-quality condition, not a source class.
- H-mix: both, in proportions the measurement states.

## Measurement (derived from observed DJA / MAST catalogs)
- Fetch the DJA v7.4 zout member (SOURCES.md sha256) and phot catalog; match ranked sources at 0.2″ (as before).
- For no-fit rows vs fitted rows of the ranked stratum: number of usable filters (`nusefilt` in the zout, if present),
  compactness (pipeline/DJA size, e.g. flux_radius in F200W) and the DJA phot `flag`/S/N.
- Re-run the same comparison for the 4 top-20 no-fit rows individually.

## Decision rule (ASSUMPTION thresholds)
- If ≥ 80 % of no-fit rows have fewer usable filters than the 5th percentile of fitted rows → H-cov: the flag is a
  coverage flag; no star rule; record it.
- Else if no-fit rows have median F200W flux_radius below the 10th percentile of fitted rows and ≥ 2 of the 4 top-20
  no-fit rows are compact → H-star supported: propose a compactness screen item (not adopted here).
- Otherwise H-mix / undetermined; record numbers.

## What changes the plan
H-star would add a compactness check to the D-019–D-021 screen before SMACS top-k vetting.
