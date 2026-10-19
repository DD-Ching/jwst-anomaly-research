# W3 dimming: D-039 persistence test on the compact single-epoch vanish flags (partial run)

Run 2026-10-19 16:10–16:42 UTC under the frozen `PLAN.md`. Screen code, thresholds and the D-039 test unchanged.
Catalogues from S3 (`fetch_cats.py`, sha256 checked against the committed manifests); Gaia mask from the ESA TAP
(`gaia_esa.py`, as in #376). Every number below is in `summary.json`, the per-field `verdicts.ecsv` /
`exposures_e*.ecsv`, or the committed manifests.

## Deviations (post hoc)
- **NEXUS ran on the 17 targets only** (`run.py ... targets`): with 151 monitored controls its first epoch had not
  finished after 13 min (each `_cal` open streams several MB; not stored). No NEXUS controls.
- **MACS0416 and Abell 2744 did not finish** in the session: MACS0416 epochs 0–1 and Abell 2744 epoch 0 are done
  (160 and 1,200 `_cal` files per epoch); their verdicts were rebuilt from the stored per-exposure tables with
  `partial.py` (same `classify`). None of their 20 targets was detected in a finished epoch, so **their targets are
  not tested yet**.

## Targets (derived; NEXUS, all 17 in epoch 2 = `jw05105-o014`, 2,000 `_cal` files, 0 failed reads)
| verdict | n | uids |
|---|---|---|
| `on_sky` | 12 | d10423 d13073 d15524 d18183 d20291 d22021 d23045 d24470 d29073 d29638 d30717 d46213 |
| `inconclusive` (one clean detection, F200W only) | 4 | d11436 d12360 d12528 d30428 |
| `undetected` (no exposure at S/N ≥ 5) | 1 | d30295 |
| `persistence` | 0 | |

- Mosaic F200W S/N of the targets: 10.3–32.0. Most `on_sky` targets are detected in several filters of the same
  visit (F090W–F444W; `nexus/exposures_e2.ecsv`), so they are sky sources in that epoch, not afterimages.
- **Open ordinary explanation (not tested here):** every NEXUS target is detected only in epoch 2 (`o014`), whose
  F200W catalogue is 65 MB against 12–27 MB for the other NEXUS epochs (`data/manifests/dimming_nexus.ecsv`). A depth
  or footprint difference of that epoch, not a vanish, is the first thing `forced` must test.

## Controls (derived; partial: 61 monitored (source, epoch) pairs, MACS0416 51, Abell 2744 10)
- `on_sky` 61 / 61 (S/N 10–20: 5 / 5; 20–50: 56 / 56); `persistence` 0 / 61. No pair below S/N 10 or above 50.
- Frozen rule 1 evaluates to "usable" on these pairs, but the set is partial (no NEXUS; 3 of 7 MACS0416 / Abell 2744
  F200W epochs), so **no swap decision is taken from it**.

## Failed approach
- `pkill -f persistence_run/run.py` killed the calling shell as well (scripts/CLAUDE.md warns about this); stop
  workers by PID.

## Next (in order)
1. `dimming_screen.py forced` + exposure-level check on the 12 NEXUS `on_sky` targets, with the `o014` depth and
   footprint test first; survivors to `/vet-candidate`.
2. Finish MACS0416 epochs 2–3 and Abell 2744 epochs 1–2 (`run.py`; about 8 min per 1,200-file epoch with
   controls, not stored), which also tests their 20 targets; NEXUS controls on a subset of monitored sources.
3. Only then the veto swap in `efficiency_table` with the per-bin `on_sky` efficiency, and `combine`.
