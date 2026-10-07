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
   Check SIMBAD/NED/Gaia (`crossmatch.py`) and published models, e.g. SMACS 0723 lens models
   (docs/landscape.md), before calling anything unexplained.
4. Only then, genuinely unexplained objects. Even these are reported as "unexplained under tests X, Y, Z",
   never as exotic physics.

Each candidate's vetting record lists which tests were run, their outcomes and their evidence links.

## Baselines

Every new model is compared against the simple baseline (robust z-scores, Isolation Forest, LOF on
catalog features) using injection-recovery on simulated outliers and the stability of rankings across
seeds. A model that doesn't beat the baseline doesn't replace it.

## Tracked metric: top-k composition

Every run reports, per ranked stratum, how its top k (the stored candidates) splits up. The numbers are also
saved under `samples[].topk` in `run_record.json`.
- **Known objects:** a SIMBAD/NED/Gaia match within the cross-match radius.
- **Lens-related:** matches with lens object types.
- **Catalogued stars.**
- **Cutout-flagged:** any image quality flag.
- **`spikes`:** D-018.
- **Screened out:** D-019; these are counted before selection.

Use it to judge a change to the pipeline. A drop in flagged or star fractions means less contamination. Known
and lens-related fractions show recovery of real, catalogued populations. The uncatalogued remainder is the
discovery set, and its members still need `/vet-candidate`. Baseline: run `20261007T183007Z-4342c5c2` (SMACS NIRCam).
- Galaxy top 20: 35% known, 10% lens-related, 0% stars, 15% cutout-flagged, 2 screened.
- Star top 10: 100% spikes.
