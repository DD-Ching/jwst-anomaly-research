# W3 dimming: k = 1 NaN-tolerant aperture on unseen NEXUS e2 controls: not tested (R3); > 4 class holds (2026-10-20)

Plan frozen first (`PLAN.md`, committed before any data read). Runs: `persistence_run/run.py nexus <gaia>
results/dimming/nan_k1/nexus/<dir> epochs=<k> max_nan=<K>` for e0 (K = 1, 800 `_cal` files) and e2 (`jw05105-o014`,
2,000 files, K = 1 and K = 0 in parallel), all in about 20 min, 0 failed reads, `_cal` data streamed from S3 (0 MB
stored). Scored by `score.py` (`score.json`). Every number is derived; thresholds are ASSUMPTIONs.

## Gates
- **G0 (code identity): pass.** e0 against the committed k = 0 table: 1,006 / 1,006 finite fluxes identical, no
  row missing, 66 pairs without a NaN exposure keep their verdict; k = 1 makes 36 more exposures finite. e2, k = 1
  against this run's k = 0: 1,453 / 1,453 identical, 36 pairs unchanged, 253 exposures made finite.
- **Reproduction:** the e2 k = 0 target verdicts equal the committed targets-only run (`persistence_run`) for all 17.

## Unseen check (NEXUS e2 monitored pairs, mosaic S/N ≥ 10; all 151 are ≥ 10)

| class | k = 0 `on_sky` / n (Wilson 95 %) | k = 1 `on_sky` / n (Wilson 95 %) |
|---|---|---|
| `n_exposures` ≤ 4 | 8 / 8 (0.676–1.000) | 8 / 8 (0.676–1.000) |
| > 4 | 143 / 143 (0.974–1.000) | 143 / 143 (0.974–1.000) |
| `persistence` (all monitored) | 0 / 151 (0–0.025) | 0 / 151 (0–0.025) |

- **Decision: R3 (not tested).** e2 has only 8 ≤ 4 pairs (< 10, the frozen minimum), and they are all `on_sky`
  without the fill. The swap keeps ε(≤ 4) = 0.412; k = 1 is neither adopted nor rejected.
- **> 4 class: holds.** First NEXUS held-out test of ε(> 4) = 0.961 (handoff 3 of `expcount_eff`): 143 / 143,
  interval above 0.961. NEXUS monitored `persistence` over all controls run so far: 0 / 430.

## Post hoc (not frozen)
- e0 at k = 1: ≤ 4 9 / 10 `on_sky` (unchanged from k = 0), > 4 94 / 94.
- The ≤ 4 deficit is a property of the wide 120-file NEXUS epochs (e1, e3–e7: two-exposure pairs with
  DO_NOT_USE aperture pixels); the deep epochs (e0, e2) lose almost nothing at ≤ 4, with or without the fill.
- Target nexus-d30295 (12 covering exposures, mosaic S/N 32): `undetected` at k = 0 becomes `inconclusive` (one
  clean detection) at k = 1. Not a limit input here (targets keep k = 0); recorded for its vetting
  (`results/dimming/nexus_o014/RESULT.md`).

## Failed approach
- Using NEXUS e0 / e2 as the unseen ≤ 4 sample: they hold 10 and 8 such pairs, too few to test a k = 1 efficiency,
  and both are already near 1 at k = 0. No unseen ≤ 4 sample exists in the three fields (MACS0416 and Abell 2744
  have none).

## Next
1. A k = 1 adoption needs new wide-epoch data: a future NEXUS wide epoch (or JADES) as the held-out sample, frozen
   with the 0.584 value of `PLAN.md`.
2. Split the held-out `inconclusive` pairs by cause (unchanged).
