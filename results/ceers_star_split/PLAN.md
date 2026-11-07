# Plan: CEERS t021 with stars split (TASKS "Then" 8, D-022)

Written 2026-11-07 10:20Z before any fetch. Re-runs `configs/reference_sample.yaml` unchanged; since #743 the Gaia
chain falls back from CDS XMatch to the ESA TAP `field` backend (ESA TAP answered 200 at 10:12Z), and SIMBAD to the
CfA mirror. CDS still resets connections (10:12Z).

## Question
#750 ranked CEERS as one mixed stratum because classification failed. Does splitting CEERS stars (D-012) change the
science-vs-control verdicts of `results/science_control_features/RESULT.md`?

## Expectation and ordinary explanations
- CEERS (high latitude, small field) holds few Gaia / SIMBAD stars; the stellar locus needs 10 reference stars and
  may not apply (configured, D-022 comment). Expected: < 2 % of the 3156 ranked CEERS rows move to a star stratum.
- Stars are compact and blue-ish in NIR colours, so their removal moves bright-bin morphology (CI, sharpness) most.

## Measurement
`scripts/science_control_features.py` on the new run (same frozen thresholds as #750: binned KS D ≤ 0.1 comparable,
> 0.2 differ strongly; ASSUMPTION), plus the CEERS star count by basis and a feature-by-feature D change.

## What changes the plan
- Any feature verdict changes category, or any binned D moves by > 0.05 (ASSUMPTION): record a D-022 addendum and
  re-state the control-comparison rule.
- CEERS classification fails again: record the failing backend and stop.
- Otherwise: Then 8 closed for the star split; record the bound.

## Data
SMACS and CEERS pipeline catalogs and the DJA catalogs (CEERS 250.5 MB, stated reason in D-022) streamed to the
scratch data root and deleted after the run (Cloud disk decision).
