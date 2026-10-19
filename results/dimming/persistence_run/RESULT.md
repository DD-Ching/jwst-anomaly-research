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

## Continuation 2026-10-19 18:08–19:05 UTC (all 37 targets tested; Abell 2744 epoch 2 controls missing)
Same frozen `PLAN.md`, catalogues re-fetched from S3 (sha256-checked), Gaia mask from the ESA TAP. `run.py` gained
`epochs=` to run only the unfinished epochs; `partial.py` rebuilds `verdicts.ecsv` from all stored per-exposure
tables (same `classify`; reads `.ecsv.gz`, used where a table is over 1 MB) and records `complete`.
- **MACS0416 complete** (epoch 2 = `jw01208-o004`, 320 `_cal`; epoch 3 = `jw01176-o213`, 160 `_cal`; 0 failed
  reads). All **12 targets `on_sky`**, every one detected only in epoch 2 (CANUCS, the field's only non-1176 F200W
  epoch; F200W catalogue 5.5 MB against 3.4–3.8 MB for the three 1176 epochs). Mosaic S/N 11.5–22.2; 9–45
  detected exposures each, 0–6 suspect.
- **Abell 2744 epoch 1** (`jw02561-o002`, 1,080 `_cal`, 0 failed reads): all **8 targets `on_sky`**, every one only
  in this epoch. Mosaic S/N 17.6–33.3; 21–108 detected exposures each (d07198: 19 suspect, 81 clean). Epoch 2
  (360 `_cal`, controls only: no target is detected there) did not finish in the session.
- **All 37 targets tested:** 32 `on_sky` (MACS0416 12, Abell 2744 8, NEXUS 12), 4 `inconclusive`, 1 `undetected`,
  0 `persistence`. In each field every target sits in one epoch only (MACS0416 e2, Abell 2744 e1, NEXUS e2). The
  persistence test therefore removes none of them; the NEXUS 12 were all explained by depth / footprint in `forced`
  (#378), so the same test comes first for the other 20. No flag is interpreted.
- **Controls:** 194 monitored pairs (MACS0416 121, Abell 2744 73): `on_sky` 194 / 194 (S/N 10–20: 19 / 19;
  20–50: 175 / 175), `persistence` 0 (`summary.json`; rule 1 "usable"). Every pair has mosaic S/N 11–33, so the
  efficiency below S/N 11 is **unmeasured**; NEXUS and Abell 2744 e2 controls are still missing.

## Next (replaces the list above)
1. `forced` with the epoch depth / footprint test first (as `results/dimming/nexus_o014/`) on the 12 MACS0416 and
   8 Abell 2744 `on_sky` targets; survivors to `/vet-candidate`.
2. Then the swap: given 0 `persistence` among 37 targets and 194 controls, the D-039 veto removes nothing that
   `forced` does not; decide whether the swap is still worth building once (1) is in.

## Control completion 2026-10-19 23:11–23:33 UTC (frozen addendum of 23:12; NEXUS partial by the fixed order)
Same code; catalogues from S3 (sha256-checked), Gaia mask from the ESA TAP. Done: Abell 2744 e2 (360 `_cal`) and
NEXUS e0 (800), e1, e3, e4 (120 each), in the frozen order 0, 1, 3, 4, …; NEXUS e5–e7 and the e2 controls were not
reached (stopped at the time budget, not by outcome). 0 failed reads. NEXUS `counts.json` lists e2 in
`epochs_done` because its `exposures_e2.ecsv` is the earlier targets-only run; `controls_missing_epochs` names
the epochs still without controls (2, 5, 6, 7). Per-field figures from `summarise.py` on each
field alone; pooled in `summary.json`.

| field | monitored pairs | S/N 10–20 on_sky | S/N 20–50 on_sky | persistence | rule 1 |
|---|---|---|---|---|---|
| MACS0416 | 121 | 8 / 8 | 113 / 113 | 0 | usable |
| Abell 2744 | 127 | 19 / 19 | 107 / 108 (1 `undetected`) | 0 | usable |
| NEXUS (e0, e1, e3, e4) | 190 | 14 / 16 | 129 / 174 (43 `inconclusive`, 2 `undetected`) | 0 | **not usable** (0.74 < 0.8) |
| pooled | 438 | 41 / 43 | 349 / 395 | 0 (Wilson 95 % upper 0.9 %) | usable as written |

- **Frozen reading:** the false-afterimage rate is 0 / 438 in every field. MACS0416 and Abell 2744 pass rule 1
  (247 / 248 `on_sky`, Wilson 0.989–0.998). NEXUS fails it, so its swap needs its own efficiency (addendum rule).
- **Post hoc (not frozen): the NEXUS loss tracks exposures per source, not mosaic S/N.** NEXUS pairs with ≤ 4
  NIRCam exposures covering the source (the 120-file wide epochs e1, e3, e4) are `on_sky` 49 / 96 (Wilson
  0.46–0.56); pairs with > 4 exposures 94 / 94 (0.989–1). Every `inconclusive` pair has ≤ 1 clean detection: D-039
  needs two, and a visit's first exposure has no prior exposure to clear it. MACS0416 / Abell 2744 pairs all have
  ≥ 8 exposures. A mosaic-S/N efficiency would therefore mis-state the swap in shallow-tiling epochs.

## Next (replaces the list above)
1. Freeze an efficiency model by exposures per (source, epoch) (e.g. ≤ 4 / > 4, from this table), and check it on
   NEXUS e5–e7 (not yet run; held out from this reading) before using it.
2. `forced` + epoch depth / footprint test on the 12 MACS0416 and 8 Abell 2744 `on_sky` targets (handoff of #379,
   still open): with the swap they are unvetoed flags, so the zero-event limit needs them explained first.
3. Then fold the swap into `inject` with the per-(source, epoch) efficiency and re-run `combine`.
