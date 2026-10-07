# Methodology

## Epistemic labels

Every table, figure and written claim is tagged with one of these labels (`schema.Provenance`):

| Label | Meaning | Example |
|---|---|---|
| observed | Archive products as delivered, including pipeline catalogs | `_cat.ecsv` fluxes |
| derived | Deterministic computation from observed data | colors, ellipticity |
| simulated | Synthetic or injected data | injection-recovery sources |
| model_prediction | Output of a fitted or learned model | anomaly scores |
| assumption | Taken as given without verification | "pipeline WCS accurate to 0.05″" |
| hypothesis | Proposed interpretation awaiting tests | "candidate is a lensed arc" |

## Anomaly scores are not discoveries

A high anomaly score means only "unusual relative to this sample under this feature set and model".
The expected order of explanations for a top-ranked source is:

1. Processing or instrument artifacts: image edges, low exposure/weight, diffraction spikes,
   saturation, persistence, snowballs, ghosts or wisps, bad deblending.
2. Catalog effects: mismatched cross-band associations, NaN magnitudes from negative fluxes,
   blended neighbours.
3. Known astrophysical populations that are merely rare in the sample: stars and brown dwarfs,
   high-redshift dropouts, dusty or line-dominated galaxies, AGN, mergers, lensed arcs in cluster fields.
4. Only then, genuinely unexplained objects. Even these are reported as "unexplained under tests X, Y, Z",
   never as exotic physics.

Each candidate's vetting record lists which tests were run, their outcomes and their evidence links.

## Baselines

Every new model is compared against the simple baseline (robust z-scores, Isolation Forest, LOF on
catalog features) using injection-recovery on simulated outliers and the stability of rankings across
seeds. A model that doesn't beat the baseline doesn't replace it.
