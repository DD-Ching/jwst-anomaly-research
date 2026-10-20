# W3 dimming: D-039 `on_sky` efficiency by exposures per (source, epoch)

Frozen 2026-10-20 01:20 UTC, before any NEXUS e5–e7 control verdict exists (handoff 1 of
`results/dimming/persistence_run/RESULT.md` "Next"). Code (`persistence_run/run.py`, `persistence_check`) unchanged.

- **Hypothesis.** The D-039 test's `on_sky` efficiency on known on-sky sources depends on the number of NIRCam
  `_cal` exposures covering the source in that epoch (`n_exposures` in `verdicts.ecsv`), not on mosaic S/N: a pair
  needs two clean detections, and a visit's first exposure cannot be cleared.
- **Ordinary alternatives:** the e0–e4 dependence is a field-position or epoch effect that happens to correlate
  with exposure count (e.g. only the 120-file wide epochs e1/e3/e4 are shallow), or an S/N effect at fixed count.
- **Model (training data: NEXUS e0, e1, e3, e4 monitored controls; MACS0416 / Abell 2744 all have ≥ 8):**

  | class | train on_sky | point | Wilson 95 % | efficiency used by the swap (ASSUMPTION) |
  |---|---|---|---|---|
  | `n_exposures` ≤ 4 | 49 / 96 | 0.510 | 0.412–0.608 | 0.412 |
  | `n_exposures` > 4 | 94 / 94 | 1.000 | 0.961–1.000 | 0.961 |

- **Held-out check:** `run.py nexus <gaia> results/dimming/persistence_run/nexus epochs=5,6,7` (monitored + targets,
  frozen order 5, 6, 7 until the time budget ends; a partial set is not chosen by outcome), merged with `partial.py`.
  Only monitored pairs with mosaic S/N ≥ 10 are scored (same population as the training table).
- **Decision rules (ASSUMPTION).**
  1. Per class with ≥ 10 held-out pairs: the model holds if the held-out Wilson 95 % interval contains the
     swap efficiency or lies above it (the swap does not over-state efficiency). A class whose held-out upper
     bound is below the swap efficiency fails, and the swap then uses the pooled (train + held-out) Wilson lower
     bound for that class.
  2. Classes with < 10 held-out pairs: not tested; recorded as such.
  3. Monitored `persistence` verdicts in e5–e7: > 2 % (Wilson lower bound above 2 %) makes the test unusable for
     the NEXUS swap (as in the 23:12 addendum).
- **What would change the plan:** rule 1 failing in the ≤ 4 class (exposure count is not the driver; look at
  position / epoch), or any `persistence` among monitored sources.
