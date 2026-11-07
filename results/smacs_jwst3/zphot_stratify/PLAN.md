# Plan: DJA photo-z as a stratifier of the SMACS galaxy ranking (TASKS Next (M2), D-022 addendum 2)

Written 2026-11-07 before any fetch or re-score. Run: `configs/smacs_jwst3_gaia_only.yaml` (jwst 3.0.0 catalogs,
Gaia-only star separation; CDS unreachable from the cloud session). Photo-z: DJA v7.4 SMACS eazy `zout` (SOURCES.md),
matched at 0.2″ as in `scripts/science_control_features.py`.

## Question
The ranked SMACS galaxy stratum mixes z ≈ 0.39 cluster members and background galaxies, whose feature distributions
differ strongly (D-022 addendum 2). A single robust-z / isolation-forest / LOF score over a mixed population can rank
minority-population members as "anomalous" because they belong to the minority, not because they are odd.
Does scoring each z_phot stratum separately change the top 20, and in which direction?

## Hypotheses and ordinary explanations
- H0 (ordinary): members are ~10–20 % of the ranked stratum and their colour/size offsets push bright members into
  the overall top 20 (population effect, not anomaly).
- H1: the overall top 20 is already dominated by background sources whose oddness survives within-stratum scoring.
- Catalog effects (photo-z failures, unmatched sources) put rows in the "other" stratum; reported, not interpreted.

## Measurement (derived)
- Strata by z_phot (ASSUMPTION ranges, as D-022 addendum 2): member 0.25–0.55, background ≥ 0.8, other (incl.
  unmatched or NaN).
- Baseline: re-score the run's galaxy stratum (`rank.score_anomalies`, random_state 0) and check it reproduces the
  pipeline ranks exactly.
- Stratified: score the member and background strata separately with the same procedure.
- Report: stratum counts; member fraction of the overall top 20 vs of the stratum; overlap of the overall top 20
  with the union of the per-stratum top 20s (scaled to the stratum's share); rank of each Mahler+2022 ICLv2
  catalogued multiple image (0.3″ match, ASSUMPTION) in the overall and background rankings (descriptive only).

## Decision rule (ASSUMPTION thresholds)
- Member fraction of the overall top 20 ≥ 2× the member fraction of the ranked stratum, and ≥ 4 of 20: population
  effect confirmed (H0) → z_phot stratification becomes the rule for SMACS top-k statements (D-022 addendum) and a
  TASKS item wires it as a pipeline option.
- Otherwise: no rule change; record the numbers.

## What changes the plan
If background ranks of catalogued multiple images rise sharply (median rank improves by > 2×), the background
stratum is the better place to look for lens-model residual candidates; this would move the next cycle to wiring
the option before any new top-k vetting.

## Data
Pipeline catalogs (small, S3) and the SMACS DJA photometry/zout (59 MB tarball, one member used) to the scratch data
root; deleted after the run (Cloud disk decision).
