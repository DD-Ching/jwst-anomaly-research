# W3 dimming: why NEXUS monitored pairs come out `undetected` (frozen 2026-10-20 05:11 UTC)

Handoff 1 of the 2026-10-20 swap_inject entry (`results/dimming/expcount_eff/RESULT.md` "Next" 2).

**Seen before this plan (observed, from the committed `persistence_run/nexus/exposures_e*.ecsv`):** all 13
monitored `undetected` NEXUS pairs (11 held out e5–e7, 2 train e4) have NaN flux in *every* covering exposure,
i.e. no per-exposure measurement, not a low S/N. (The `undetected` flag d30295 in e2 has 4 finite exposures, all
S/N < 5, and 8 NaN ones.) NaN-flux exposures are 70 / 1076 of the
monitored exposures in e0 and 11–26 per epoch (about a quarter to a third) in e1 and e3–e7. Nothing below has
been computed yet.

**Hypothesis.** (H_nan) The ≤ 4-exposure `on_sky` deficit (D-039 swap value 0.412) is mostly a measurement
coverage effect: `aperture_flux` returns NaN when any aperture pixel is no-data (array edge, DQ DO_NOT_USE pixels
set to NaN in `_cal`), and a pair needs two finite, clean detections. (H_cnt) It is an exposure-count / S/N effect
among measurable exposures. Ordinary explanations of a NaN aperture: array edge, bad-pixel / reference-pixel
masks, cosmic-ray (jump) flags set to DO_NOT_USE, saturation, a read failure.

**Measurement.**
1. `stamps.py` (range reads, 31 × 31 px SCI + DQ around each position, no full download) for every NaN-flux
   exposure of the monitored pairs in e4–e7 and of d30295; `reason()` classifies each one (function + unit
   test, written before the run): `edge` (centre closer than r_px to the border, `measure_pixel`'s rule),
   `aperture_no_data` (interior, ≥ 1 NaN aperture pixel; the OR of the DQ bits of those pixels is recorded),
   `annulus` (aperture finite, < 10 finite annulus pixels), `not_reproduced` (finite now).
2. `score.py` (tables only, no download), per class (≤ 4 / > 4 exposures) and split (train e0, e1, e3, e4 /
   held out e5–e7): `on_sky` / n over all pairs (reproduces `expcount_eff`), and over the *measurable* pairs
   (≥ 2 exposures with finite flux), with Wilson 95 % intervals.

**Decision rules (ASSUMPTION; computed by `score.py`'s `decide()`).**
- R1 (H_nan supported), applied to held-out ≤ 4 only: the Wilson 95 % lower bound of the measurable-pair
  `on_sky` fraction is ≥ 0.80 (the D-039 rule-1 level). Then the swap's ≤ 4 efficiency is best modelled as
  P(≥ 2 measurable) × P(`on_sky` | measurable), and the next PR changes `persistence_eff_by_exposures` to use it
  (positions of injected copies get their own measurability from the same `_cal` footprints).
- R2 (H_cnt): measurable lower bound < 0.80 in held-out ≤ 4 → the deficit is not only coverage; keep the
  exposure-count model and record the measurable rate.
- The `undetected` stamps change nothing in the limits by themselves (the swap already counts them as misses);
  they are explained when every exposure has an `edge` or `aperture_no_data` reason.

**What would change the plan:** a `not_reproduced` or read-failure reason (the measurement itself is unstable;
re-run before any modelling), or a real visible source in an `aperture_no_data` stamp with only isolated flagged
pixels (then a NaN-tolerant aperture is the fix, not a new efficiency model).
