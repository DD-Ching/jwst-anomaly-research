# Plan: forced-photometry vetting of review-list entries SMACS h039 and Sunrise o052 h005

Frozen 2026-10-19 before any measurement (cloud routine). Every number will be `derived`; thresholds are ASSUMPTIONs.

## Why
`../RESULT.md` leaves two review entries open: h039 (SMACS ERO o001 → VENUS o057, ± pair across a compact core)
and h005 (Sunrise o010 → VENUS o052, negative on a galaxy core, near the map edge). Both were called artifacts
by eye only. SMACS h063 was closed with forced photometry (`../../smacs_arcdiff/RESULT.md`); this applies the
same test to the two open entries.

## Hypotheses
- H0 (ordinary, expected): the difference-image residual is a subtraction artifact (PSF/kernel mismatch or
  registration on a compact core); the source's aperture flux is unchanged between the pair's epochs.
- H1: the source really changed (an AGN, SN or other variable; ordinary unless later vetting says otherwise).
- Nothing here can test an exotic-lensing hypothesis: both entries are 3.7″ / 10.5″ from the nearest catalogued
  image and outside the arc mask.

## Method
`scripts/hit_lightcurve.py` (unchanged: 0.15″ aperture, 0.35–0.6″ annulus, 300 random-aperture errors) at the
review positions in `../smacs_o001_o057.json` and `../sunrise_o010_o052.json`, every covering public NIRCam
mosaic. The test uses only the pair's own mosaics: SMACS jw02736 o001 (2022-06) and jw06882 o057 (2026-06);
Sunrise jw02282 o010 (2022-07-30) and jw06882 o052 (2025-07-08); bands F150W and F444W.

## Frozen rule (per entry; Δ = flux_B − flux_A, σ = √(err_A² + err_B²))
- **no ≥ 2σ change**: |Δ|/σ < 2 in both bands. The entry is closed as a difference-imaging residual.
- **change**: |Δ|/σ ≥ 2 in both bands with the same sign. The entry stays open and goes to `/vet-candidate`.
- **inconclusive**: anything else (one band only, or opposite signs).
- **untested**: a pair mosaic is missing in either band at the position.
Other epochs are reported post hoc only. The stamps are looked at before any conclusion is written.

## What would change the plan
A `change` reopens the entry as an unvetted transient (ordinary first). An `untested` means the review stays open.
