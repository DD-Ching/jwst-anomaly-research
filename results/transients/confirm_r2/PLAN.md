# Confirmation round 2: plan (fixed before any new comparison)

Owner instruction 2026-10-18 (after #327). Budget 08:29Z → **12:29Z** (session clock), ≥ 30 min kept for review
and records.

## A. SMACS box 1 h005 control comparison (rules of `../batch_venus4/confirm/PLAN.md`, unchanged)
- **Stall diagnosis (done before this plan, 08:35Z):** `hit_lightcurve.covering_products` issues one MAST
  `query_criteria` per position; with 28 controls in 6 position groups, a portal request in the second group waited
  > 150 s for an HTTP status line (stack dump: `http.client._read_status` inside astroquery `discovery_portal`).
  Not data reads or computation. Fix: new `--group-arcsec` option; with `--group-arcsec 90 --query-only h005` the
  controls share h005's mosaics and one MAST query (all controls are within 30″ of h005).
- **Controls:** the 28 positions already picked (`../batch_venus4/confirm/smacs_h005_pos30.txt`, seed 1), minus any
  within 0.3″ of h005 (the entry itself), then de-duplicated: in listed order, a control within 0.5″ of an earlier
  kept control is dropped. ≥ 20 remaining required, else not decided.
- Decisive pair F150W 2736 o001 → VENUS o057; secondary F444W same pair and F444W VENUS → 12774. Rule: tail < 0.05 ⇒
  exceeds matched-control scatter; ≥ 0.05 ⇒ consistent with it (not proof of an instrumental cause).
- Cap: 30 min wall; if it stalls again, the reproducible case is saved and A is closed as not decided.

## B. One discriminating test: chosen in `B_choice.md` before any new data are read.
