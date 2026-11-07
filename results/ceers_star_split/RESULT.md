# CEERS t021 with stars split: 6 bright stars leave the ranked stratum; no verdict changes, background match unchanged

Plan: `PLAN.md` (pushed before the run). Run `20261107T101020Z-c9a64ce3` of `configs/reference_sample.yaml`, unchanged,
at main `29b31d7`. Comparison: `scripts/science_control_features.py <run> --zout <SMACS DJA v7.4 zout>` (zout sha256
as in SOURCES.md) → `summary.json`, `feature_cdfs.png` (looked at). Reference: `../science_control_features/` (#750).

## Tool check (observed)
- CDS XMatch refused both `vizier:I/355/gaiadr3` and `simbad`; CDS SIMBAD TAP reset (10:18–10:26Z). The default
  chains fell back to the ESA Gaia TAP `field` backend and the CfA SIMBAD mirror; CEERS classification took 124 s,
  most of it CDS retries.

## CEERS strata (derived)
- 7582 merged sources, 3156 pass the quality gate (as in #750). 6 are stars (4 SIMBAD/Gaia astrometry, 2 Gaia
  position only); all 6 pass the gate, so the galaxy stratum is 3150. The stellar locus (D-015) did not apply: too few
  reference stars (6 < `min_ref_stars` 10). One source was returned to "other" by the extended veto (D-025).
- All 6 stars have `ref_mag` 17.6–21.1, so only the ≤ 24 bin changes (CEERS 55 → 49 rows).

## Frozen rule on the new stratum (ASSUMPTION thresholds as #750)
- **No feature changes verdict** (overall or excluding the bright bin).
- Largest binned D change: `ref_log_ci_70_50` ≤ 24, 0.363 → 0.308 (Δ 0.055 > the plan's 0.05 trigger); next
  `ref_log_ci_50_30` 0.372 → 0.338, `ref_log_nn_dist` 0.244 → 0.217, area 0.298 → 0.321. Every D outside the ≤ 24
  bin moves by ≤ 0.010.
- **Background match (post hoc z_phot split, D-022 addendum 2) unchanged:** SMACS background (z_phot ≥ 0.8) vs CEERS
  moves by ≤ 0.010 in every 24–28 bin. The member-vs-CEERS bright-bin D values move by ≤ 0.05 (e.g. area 0.310 →
  0.357); members still differ strongly in colour and morphology.

## Post hoc (looked at the plot; not a rule)
The CEERS ≤ 24 CI curves keep a low-CI step: 11 / 49 ranked rows have `ref_log_ci_70_50` < 0.2 (0.135–0.163) with
F150W−F200W −0.39 to +0.03, compact and blue as stars are (`ref_mag` 20.1–23.5). Hypothesis: uncatalogued stars
fainter than the Gaia limit, which the locus cannot catch without 10 reference stars. Not tested; it bears only on
the bright bin, which #750 already called weak (KS D ≈ 0.23 reachable by chance at n ≈ 50).

## What this changes
- TASKS Then 8 "re-run CEERS with the star split": done. D-022 addendum 3 restates the control rule unchanged: SMACS
  z_phot ≥ 0.8 in `ref_mag` bins 24–28; the ≤ 24 bin is not used for control statements (small n, possible
  unsplit CEERS stars).

## Limitations
- Star split by catalogues only; CEERS stellar locus not applicable (above).
- CDS unavailable: SIMBAD came from the CfA mirror (D-006 addendum; it may lag CDS).
- Data: pipeline catalogs, DJA catalogs (CEERS 250.5 MB, D-022) and the SMACS zout streamed to scratch, deleted after
  the run (Cloud disk decision).
