# W3 dimming: k = 1 NaN-tolerant aperture, frozen check on unseen NEXUS controls (frozen 2026-10-20 07:20 UTC)

Handoff 1 of `results/dimming/nan_aperture/RESULT.md`. Written before any k = 1 run and before any e2 control read.

**Seen before this plan (derived, committed):** k = 1 passed G1 as written (1 masked pixel: median flux error 3.7 %,
S/N call 42 / 42). Post hoc k = 1 (approximated from the k = 2 tables) on NEXUS e1, e3–e7, ≤ 4-exposure monitored
pairs, mosaic S/N ≥ 10: 115 / 175 `on_sky` (Wilson 0.584–0.723), against k = 0 73 / 175. NEXUS e0 at k = 0: ≤ 4
class 9 / 10 `on_sky` (k = 1 can move at most one pair: no power), > 4 class 94 / 94. MACS0416 and Abell 2744 have
no ≤ 4 monitored pairs. **The only NEXUS controls never measured at any k are epoch 2 (`jw05105-o014`, 2,000
`_cal` files)**, the epoch holding all 17 NEXUS targets.

**Hypothesis.** (H_k1) The k = 1 fill raises the ≤ 4-exposure `on_sky` efficiency to about the seen 0.58–0.72 on
data it was not tuned on, without inventing illuminators. Ordinary alternatives: the seen gain is specific to the
wide 120-file epochs (e2 is a deep 2,000-file epoch with other dithers and detector positions), or e2 has few ≤ 4
pairs, so nothing is tested.

**Proposed value.** ε_k1(≤ 4) = 0.584, the Wilson lower bound of the seen pooled k = 1 rate (ASSUMPTION, same
convention as `expcount_eff`: lower bound of the train interval).

**Measurement.**
1. G0 (code identity on an epoch never run at k ≥ 1): `run.py nexus <gaia> results/dimming/nan_k1/nexus/e0
   epochs=0 max_nan=1`; every exposure with finite committed k = 0 flux keeps the identical flux, no row is missing,
   and pairs without a NaN k = 0 exposure keep their verdict.
2. Unseen check: `run.py nexus <gaia> results/dimming/nan_k1/nexus/e2_k1 epochs=2 max_nan=1` and, in parallel,
   the same with `max_nan=0` into `e2_k0` (the first k = 0 e2 controls). Run budget: both stopped at 07:52 UTC if
   unfinished; an unfinished run writes nothing and is reported as not done (no partial chosen by outcome).
3. `score.py`: monitored pairs with mosaic S/N ≥ 10, `on_sky` / n by class (≤ 4 / > 4), Wilson 95 %, k = 0 and
   k = 1 side by side; `persistence` counts.

**Decision rules (ASSUMPTION).**
- R1 (adopt k = 1 for the next PR): G0 passes; e2 k = 1 has n(≤ 4) ≥ 10; its ≤ 4 Wilson interval contains 0.584
  or lies above it, and its lower bound is > 0.412; e2 monitored `persistence` Wilson lower bound ≤ 2 %. Then the
  next PR re-runs exact k = 1 on every NEXUS epoch and the 37 targets in all three fields, sets ε(≤ 4) = 0.584
  in the swap and re-runs `inject` / `combine`.
- R2 (not adopted): any condition of R1 fails. The swap keeps 0.412.
- R3 (not tested): e2 k = 1 unfinished, or n(≤ 4) < 10. The swap keeps 0.412; the > 4 class is reported.
- The > 4 class on e2 (first NEXUS held-out test of 0.961, handoff 3 of `expcount_eff`): holds if the k = 0
  Wilson interval contains 0.961 or lies above it. Reported; nothing changes a limit in this PR.

**What would change the plan:** G0 failing; a filled exposure making a monitored pair `persistence`; e2 ≤ 4
k = 1 rate below 0.412 (the fill helps only the epochs it was seen on).
