# SMACS 0723 (program 2736) on MAST's jwst 3.0.0 catalogs: six of nine candidate records reproduce

Plan: `PLAN.md` (pushed before any fetch). Run `20261106T121217Z-c191b5fb`, `configs/reference_sample.yaml`, main at
`5fbb95d`. Tables: `candidates.csv` (per record), `summary.json` (per sample). Script: `scripts/smacs_jwst3_compare.py`.

## Data (observed)
- MAST now serves all ten 2736 level-3 catalogs (6 NIRCam, 4 MIRI) built with jwst 3.0.0 / photutils 3.0.0; the
  manifest records the new sizes and sha256. NIRCam rows per band rise 22–30 % over the 2.0.1 builds (both row sets in SOURCES.md; ratio
  not stored).
- Merged sources (derived, `summary.json`): SMACS NIRCam 6879 (5254 on 2026-10-07 from 2.0.1, CHANGELOG
  2026-10-07), MIRI 908 (530, same entry).
  **Control:** CEERS t021, whose catalogs are byte-identical to 2026-10-07, still merges to 7582 sources with 3156
  passing the gate (both as recorded then), so code changes since then do not explain the SMACS changes.

## Candidate records (derived; match radius 0.15″, ASSUMPTION)
| Record | Recorded verdict (short) | Old rank | New rank | Counterpart |
|---|---|---|---|---|
| f200w_2804 | inconclusive: needs NIRSpec | 6 | 5 | 0.002″ |
| f277w_829 | spectroscopically ruled-out dropout | 3 | 4 | 0.002″ |
| f200w_1032 | lensed arc | 8 | 9 | 0.004″ |
| f356w_1243 | red high-z candidate | 12 | 11 | 0.001″ |
| f200w_438 | Balmer-break galaxy (hypothesis) | 11 | 12 | 0.001″ |
| f200w_2915 | quiescent galaxy | 21 | 20 | 0.001″ |
| f200w_1096 | catalog effect (sharpness, EE colours) | 16 | 520 | 0.054″; `ref_mag` and sharpness NaN |
| f200w_2559 | catalog effect (sharpness, EE colours) | 8 | 728 | 0.003″; sharpness NaN |
| f200w_2925 | catalog effect (sharpness, EE colours) | 3 | — | none (nearest 0.22″, 7 sources within 1″) |

- The six records with an astrophysical or open verdict reproduce at 0.001–0.004″ with rank changes of ±1 and
  scores within 0.007 of the recorded ones (difference not stored; both in `candidates.csv`). A pipeline-version (processing) origin is therefore not supported for
  them at the catalog level; their verdicts stand unchanged.
- The three "catalog effect" records no longer rank: for 1096 and 2559 the ranking feature (`ref_sharpness`) is now
  NaN by the D-013 point-like gate (`daofind_max_ci` 1.8), a **code** change made after their vetting, not a jwst
  3.0.0 effect; 2925 has no counterpart within 0.15″ and 7 sources within 1″, consistent with a
  changed deblending (images not checked). This is consistent with
  their recorded verdicts.

## Limitations
- **Deviation:** CDS (SIMBAD TAP and XMatch) reset every connection from this cloud session (`curl` exit 35, proxy
  log `ws_closed_mid_exchange`), so `classify_sources` failed by design and each sample was ranked as one mixed
  stratum, without D-012 star separation or D-019 spike screening. New ranks are therefore not the galaxy-stratum
  ranks the records quote; they are reported only as a stability check. The run's mixed top 20 (contact sheet seen:
  5 spike-flagged panels, one target on a bright halo) is not a candidate list and was not vetted.
- Old ranks come from runs of different code versions (each record names its run and commit).
- Not done: a full galaxy-stratum re-run and a per-source feature comparison against the 2.0.1 catalogs (no longer
  served by MAST; earlier copies were not kept).
