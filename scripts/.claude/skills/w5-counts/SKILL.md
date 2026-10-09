---
name: w5-counts
description: Recipe for the W5 count-deficit search (a negative-mass lens empties a disc of radius ~θ_E in the background galaxy counts) on survey count maps - server-side HEALPix counts, predicted profile, matched filter, cross-region null, vetting, injections and the sky-density limit - with its failed-approach rules. Use when changing src/jwst_anomaly/countmap.py or scripts/w5_counts.py, or adding a survey or depth to W5.
---

# W5 count deficits

Rules: root CLAUDE.md, `scripts/CLAUDE.md`, `src/jwst_anomaly/CLAUDE.md`. Decisions: D-047 (prediction), D-054,
D-063 (this search). Results: docs/exotic_limits.md "W5 count deficits".

## Commands
`python scripts/w5_counts.py fetch | numcounts | predict | screen | vet | sheet | inject | limit`
(tracked outputs in `results/w5_counts/`; count chunks cached under `$JWST_ANOMALY_DATA`).

## Method
- Counts, not catalogues: `jwst_anomaly.countmap` (`signatures.CountMapSurvey`) aggregates galaxies per `nest4096`
  pixel on the Data Lab server, with a masked fraction and depth per pixel; astropy-healpix for geometry.
- Prediction: `exotic_sim.count_ratio` with the number-count slope measured from the same sample
  (`model_prediction`); ~10 % of galaxies inside θ_E are missing, so small θ_E is untestable (0.09 galaxies per lens at
  10″ in DR10).
- Matched filter on several θ_E; significance from a cross-region null built from *vetted* peaks of the other region
  (galaxy clustering makes Poisson errors 1.1–5.8× too small); detection needs N_false < 0.01 (ASSUMPTION).
- Vetting, cheapest first: mask / bright star / large galaxy (HyperLEDA D25) / depth and tile edges / core reach /
  cosmic variance via the null. Look at the contact sheet.
- Limit: injections through screen + vetting at each θ_E; n₉₅ = 3 / (area × efficiency) per deg²; convert to mass
  only with a stated geometry.

## Failed approaches (rules)
- Data Lab ADQL rejects sub-selects, CASE, SIGN and GROUP BY on expressions: group by `nest4096` only.
- Chunked RA/Dec fetches split edge pixels: sum them, or every chunk border looks like a deficit; drop region-border
  pixels from the screen.
- An unbounded bright-end count slope blows the profile up at x → 0: fix it (0.6).
- A null from the other region's raw peaks inherits its artefacts: vet the null peaks first.
- Veto radii growing with θ_E removed 90 % of random positions at 32′: veto only mimics that can empty ≥ 10 % of the
  core.
- Sky over-subtraction around large galaxies (NGC 1398) makes deficits: the large-galaxy veto is required.
- A straight-edged drop in *all* sources (not only galaxies) is a depth/tile edge: the depth-edge test; re-run
  injections after adding any test.
- HyperLEDA returns sexagesimal unless `_RAJ2000` / `_DEJ2000` are requested.
- Data Lab served ~1 chunk/min after a 30-min 502 outage: probe the service and resume from cached chunks.
- Euclid Q1 counts are clustering-limited at θ_E ≤ 4′ (D-065: ×1.0–1.4 S/N over DR10); do not port the count screen
  for small θ_E, use shapes. IRSA TAP: spatial cuts only with `CONTAINS(POINT, CIRCLE)`; RA/Dec ranges time out.
