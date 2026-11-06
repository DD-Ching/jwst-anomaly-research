# SMACS 0723 jwst 3.0.0 galaxy-stratum re-run (Gaia-only star separation): six records stay in the top 31

Plan: `../PLAN_galaxy_stratum.md` (pushed before the run). Run `20261106T142054Z-68e24ff4`,
`configs/smacs_jwst3_gaia_only.yaml` (reference sample, SMACS NIRCam only, `classify.services: [gaia]`,
`crossmatch.services: [ned, gaia]`). Code: branch `claude/smacs-jwst3-galaxy-stratum` with the Gaia `field` backend.
Tables: `candidates.csv`, `summary.json` (`scripts/smacs_jwst3_compare.py --out`), `screened.ecsv` (D-019–D-021),
`rescore.txt` (`scripts/smacs_jwst3_rescore.py`).

## Tool check first (observed)
- CDS (SIMBAD TAP, XMatch) still drops connections from the cloud session (relay log `ws_closed_mid_exchange`).
  The ESA Gaia TAP answered in 1.7 s.
- The Gaia `field` backend (one cone around all targets, local association) returned 7 rows in 13 s for the five
  test targets. It gives the same best match, `is_star`, `n_gaia` and match type as the per-target `tap` backend on
  every target, both live and from the recorded fixture `tests/data/xmatch_raw_gaia_field.ecsv`.
- On SMACS (6879 sources), classification took 48 s. The previous run was stopped after 10 min of per-target cones.

## Strata (derived)
- 3077 sources pass the quality gate. Of them, 139 are stars and 2938 galaxies are ranked. Star basis among
  quality-passing stars: 36 SIMBAD/Gaia astrometry, 97 stellar locus (D-015), 6 Gaia position only.
- Spike and weight screening ran (galaxy stratum only): 6 of the top 28 were screened out (5 spikes, D-019/D-020;
  1 low weight, D-021).

## Records (derived; match radius 0.15″, ASSUMPTION)
| Record | 2.0.1 rank | 3.0.0 galaxy rank | Δ | Paired re-score: without → with the 6 Gaia-position stars |
|---|---|---|---|---|
| f200w_2804 | 6 | 6 | 0 | 5 → 5 |
| f200w_1032 | 8 | 8 | 0 | 7 → 8 |
| f200w_438 | 11 | 11 | 0 | 11 → 13 |
| f356w_1243 | 12 | 15 | +3 | 15 → 17 |
| f277w_829 | 3 | 7 | +4 | 6 → 7 |
| f200w_2915 | 21 | 27 | +6 | 27 → 31 |

The three "catalog effect" records behave as in `../RESULT.md`: 1096 and 2559 rank 551 and 687 (sharpness NaN), and
2925 has no counterpart.

- **Hypothesis (±3) holds for 4 of 6; it fails for 829 (+4) and 2915 (+6).** The "changes the plan" condition (a
  record outside the top 50) is not met: every record stays in the top 31.
- **Missing SIMBAD override:** with the 6 Gaia-position-only stars returned to the galaxy stratum, record ranks
  move by +0 to +4. Four of the six would rank 9, 10, 11 and 23 (f200w_3209, 3999, 3469, 3794), so the galaxy top 25
  is uncertain by at most these four sources. A star rule with SIMBAD may keep them as stars (catalogued stars) or
  return them (nearer catalogued galaxies). The plan's threshold was "> 5 such stars among the top 20 neighbours"
  (ASSUMPTION). At most 4 can enter the top 25, so the result counts as determined, with this bound.
- **Reading:** the jwst 3.0.0 catalogs and the changed star rule move the six records by 0 to +6 ranks. None
  drops out, so no recorded verdict changes. The rank differences are not stored per feature, and old ranks come from
  runs of different code versions (`../RESULT.md`).

## New entries (contact sheet seen; not candidates, not vetted)
15 of the galaxy top 25 shown are not among the 2.0.1 records: f200w_2063, 318, 3446, 3177, 1816, 3774, 1552,
1788, 384, 429, 3650, 1110, 2297, 3400, 2695. The 2.0.1 runs' full top-20 lists were not kept, so "new to the top
20" cannot be decided. By eye:
- f200w_3177 (#10): no source at the target. It lies in the spike and halo of a bright star, so it is likely an
  artefact the D-019 screen did not catch.
- f200w_3446 (#9) and 1552 (#17): faint, low-contrast targets.
- The rest are compact sources or galaxies near brighter neighbours.

None is called a candidate before `/vet-candidate`.

## Limitations and failed approaches
- Without SIMBAD, the "nearer SIMBAD galaxy wins" rule and SIMBAD star types are absent. The paired re-score
  above bounds the effect.
- Re-scoring from the saved `features.ecsv` does not reproduce the pipeline's ranks: 2822 of 2938 ranks differ,
  by at most 239, while the record ranks move by ±1. So only paired re-scores, with the same procedure on both sides,
  are compared. The cause (ECSV round trip or input order feeding the isolation forest and LOF) is not isolated.
- NED and Gaia cross-matches of the top 20 ran; SIMBAD lens flags are missing, so "lens-related" comes from NED only.
