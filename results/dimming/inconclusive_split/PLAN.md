# Plan (frozen before any number is read): split the held-out NEXUS `inconclusive` pairs by cause (2026-10-20)

Handoff item 2 of `undetected_stamps`, `nan_aperture`, `nan_k1`, `o015_k1`.

## Question
`undetected_stamps` left 12 measurable held-out (e5–e7) ≤ 4-exposure monitored pairs `inconclusive` (≥ 2 finite
fluxes, mosaic S/N ≥ 10) and the held-out ≤ 4 measurable rate at 33 / 45 (Wilson lower 0.590 < 0.80, R2). Why does
D-039 give fewer than two clean detections for these visible sources?

## Ordinary explanations (to be counted, one cause per pair, first match in this order)
- **C1 `n_prior`:** a covering exposure fails only because it is the first of its visit (no prior to compare).
- **C2 low S/N:** at least one finite-flux exposure is below the per-exposure S/N 5 detection rule.
- **C3 other:** anything else (listed per pair, not explained here).

The cause is read from the committed `persistence_run/nexus` tables by `split.py`, using the per-exposure columns
the D-039 test wrote. Provenance: derived. Thresholds are the frozen D-039 ones (ASSUMPTIONs).

## What would change the plan
- If C1 explains ≥ 6 of 12 (ASSUMPTION), the ≤ 4 deficit is mainly a test-design effect of visit structure; the
  next step is a D-039 rule question for the owner/next run (not changed here).
- If C2 explains ≥ 6 of 12, the deficit is depth: ε(≤ 4) is honest as it stands and the line stops for ≤ 4.
- Otherwise list the C3 pairs with stamps as the next step.
No limit changes in this run either way.
