# W3 dimming: D-039 `on_sky` efficiency by exposure count, held-out check on NEXUS e5–e7 (2026-10-20)

Plan frozen first (`PLAN.md`, commit before any e5–e7 verdict). Run: `persistence_run/run.py nexus ... epochs=K`
for K = 5, 6, 7 (120 `_cal` files each, about 2 min each, 0 failed reads), tables merged into
`persistence_run/nexus/` with `partial.py`; scored by `score.py` (`score.json`). Provenance: derived.

| split | class | on_sky / n | Wilson 95 % | swap efficiency (frozen) | rule 1 |
|---|---|---|---|---|---|
| train (e0, e1, e3, e4) | `n_exposures` ≤ 4 | 49 / 96 | 0.412–0.608 | 0.412 | — |
| train | > 4 | 94 / 94 | 0.961–1.000 | 0.961 | — |
| held out (e5, e6, e7) | ≤ 4 | 33 / 89 | 0.278–0.475 | 0.412 | **holds** (interval contains 0.412) |
| held out | > 4 | 0 / 0 | — | 0.961 | not tested (no pairs) |

- **Frozen reading:** the ≤ 4 class holds as written; the > 4 class is untested on NEXUS held-out data (e5–e7 are
  all wide 120-file epochs). Monitored `persistence` 0 / 89 (Wilson upper 4.1 %; rule 3 not triggered).
  Combined over all NEXUS controls run so far: 0 / 279 `persistence`.
- **Post hoc (not frozen):**
  - The held-out point estimate (0.37) sits below the swap value 0.412; the pooled ≤ 4 rate is 82 / 185
    (Wilson 0.374–0.515). The rule passes on width, not on agreement, so the swap result should be quoted with the
    pooled lower bound 0.374 as a sensitivity case.
  - The ≤ 4 class is almost all two-exposure pairs (train 39 / 83, held out 29 / 81 `on_sky` at n = 2; n = 4:
    9 / 10 and 4 / 5). S/N is 10–33 in both splits, so S/N cannot be separated from exposure count here.
  - 11 held-out `undetected` pairs (train 2) are monitored sources at mosaic S/N 13–33 not found in their own one
    or two covering exposures; not examined (stamps not looked at). They lower the ≤ 4 rate and could be a
    position (detector-gap / edge) effect rather than an exposure-count one.
- **Not done:** NEXUS e2 (`o014`, 2,000 files) controls, the last epoch in the 23:12 addendum's order; started and
  stopped at the time budget (`nexus/counts.json` lists it as missing).
- **Failed approach:** running e2 after e5–e7 in the same 40-minute run (2,000 files; the 120-file epochs took
  2 min, e2 did not finish in 11 min).

## Next
1. Fold the swap into `inject` with the per-(source, epoch) efficiency by exposure count (0.412 / 0.961, pooled
   0.374 as sensitivity), counting covering exposures from the `_cal` footprints; re-run `combine`.
2. Look at the 11 held-out `undetected` stamps before (1) if the swap gain depends on NEXUS ≤ 4 pairs.
3. NEXUS e2 controls (2,000 files, a run of its own) would test the > 4 class on NEXUS.
