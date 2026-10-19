# W3 dimming: D-039 persistence test on the 37 compact single-epoch vanish flags, and its `on_sky` efficiency

Frozen 2026-10-19 16:13 UTC, before any verdict was read (handoff of #376,
`results/dimming/persistence_scope/RESULT.md` Next 1–2). Screen code, thresholds and the D-039 test are unchanged.

- **Hypotheses.** (H0) The 37 compact single-epoch vanish flags (MACS0416 12, Abell 2744 8, NEXUS 17) are
  detector afterimages, cosmic-ray residues or faint sources at the catalogue limit, i.e. not on the sky in a way a
  W3 full vanish needs. (H1) Some are real sources present in one epoch and absent in the others.
- **Ordinary explanations still open after an `on_sky` verdict:** SN / AGN / variable star (TNS, colour, host),
  catalogue-limit flicker (mosaic S/N near the detection threshold in the other epochs), edges. Those go through
  `forced` and `/vet-candidate`; this run does not interpret any flag.
- **Measurement.** `run.py FIELD`: rebuild each field's screen (catalogues from S3 by sha256, Gaia mask from the
  ESA TAP as in #376), select (a) the 37 targets (`flagged & single_epoch & ~other catalogue vetoes & point_like
  & vanish`) and (b) every monitored source; for each F200W epoch, run `persistence_check.run` on all NIRCam `_cal`
  exposures of that observation (all visits, all filters, so the 3 h lookback sees the earlier exposures) for the
  targets detected there and the monitored sources detected there.
- **Efficiency (derived).** For monitored sources (known on sky: detected in ≥ 2 epochs), the fraction of
  (source, epoch) pairs called `on_sky`, in mosaic F200W S/N bins [0, 10), [10, 20), [20, 50), [50, ∞).
- **Decision rules (ASSUMPTION; computed by `summarise.py`).**
  1. The test is usable for the swap if the monitored `on_sky` fraction is ≥ 0.8 in every S/N bin with ≥ 10 pairs
     at S/N ≥ 20, and the `persistence` verdict occurs in ≤ 2 % of monitored pairs (false-afterimage rate).
     Otherwise record the efficiency curve; the swap then uses it per bin (no all-or-nothing).
  2. A target goes to `forced` + `/vet-candidate` iff its verdict is `on_sky`. `persistence` targets are
     explained (D-039). `undetected` / `inconclusive` targets are not cleared and stay vetoed.
- **What would change the plan:** a monitored `persistence` rate > 2 % (the test would remove real sources and
  needs recalibration before any swap), or ≥ 1 `on_sky` target (vetting work before the swap).
