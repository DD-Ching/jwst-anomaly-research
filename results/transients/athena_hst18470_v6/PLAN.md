# SN Athena reappearance in HST 18470 visit 6: plan (written before measuring, 2026-11-24)

## Question
Has the next image of SN Athena (MACS J0417 h002, D-034 / #323) appeared in HST 18470 visit 6 (`ifu906010`,
WFC3/IR, MJD 61368.23 = 2026-11-24 05:30 UTC; public per MAST CAOM TAP `dbo.ObsPointing`, 2026-11-24 19:10Z,
31 rows)? It is the first public visit since visit 5 (`results/transients/athena_hst18470_v5/`).

## Predictions, ordinary explanations, method and rule
Unchanged from `results/transients/athena_hst18470_v5/PLAN.md` (frozen): same two positions, same photometry
(`scripts/athena_hst18470.py`), visit 6 compared with visit 1. **Rule (ASSUMPTION):** `appeared` if
Δ = F₆ − F₁ ≥ 5 σ_d at either position **and** a point source is visible at that position in the visit-6 stamp
(looked at); `not seen` if Δ < 3 σ_d at both; otherwise `undetermined` (function `decide` / `combine`).
Visit 6 must be F105W like visits 1–5; if it is another filter, the outcome is `undetermined` (filter mismatch)
and no difference is interpreted.

## What would change the plan
`appeared` → the discovery team's window holds and our RELICS glafic window (2027-09-14 to 2027-11-22, #309) is
too late (D-034 addendum; an ordinary model error, not an anomaly); notify the owner. `not seen` → the window
stays open to early Mar 2027; repeat at the next public visit.
