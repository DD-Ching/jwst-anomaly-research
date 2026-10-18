# Confirmation round 2: result (plan `PLAN.md`, test choice `B_choice.md`, both committed before the new data)

Every number is `derived` from the files named (`simulated` for the simulation). Session-clock times.

## A. B4-SJ0723.3-7327-box1-h005: consistent with matched-control scatter
- **Stall cause:** one MAST portal `query_criteria` per position; a request in the second position group waited
  > 150 s for an HTTP status line (stack dump in the session; not stored). Not data reads or computation. Fix:
  `hit_lightcurve.py --group-arcsec` (tested), so controls around a hit share its single query; the run then took
  ~2 min (322 measurements).
- Controls: 28 picked → 22 after removing the entry itself (c14) and 5 duplicates within 0.5″
  (`smacs_h005_controls.txt`).
- Decisive F150W 2736 → VENUS: ratio 1.051, control tail **0.591** (22 controls) → consistent with the matched-control
  scatter. Secondary: F444W same pair tail 0.227; F444W VENUS → 12774 tail 0.318 (`smacs_h005/cz_*.ecsv`). Not proof
  of an instrumental cause.

## B. Is the B1-PLCKG165-box2-h002 fade achromatic? Not decided; the data cannot separate the hypotheses
- H1 intrinsic transient (colour evolves) vs H2 magnification of a persistent compact source by a point mass, an
  Ellis wormhole (n = 2, ε > 0) or a negative-mass lens (n = 1, ε < 0), achromatic in all three cases.
- **Frozen rule: not decided.** Too few usable controls (S/N ≥ 5 at both epochs) in F090W: 10 / 15 / 15 with the 16″
  pick (`g165_h002/calibration.json`); deviation (recorded before any target ratio): 30″ pick, 14 / 31 / 32
  (`g165_h002_30/calibration.json`). F090W is the band that limits.
- **Identifiability (`simulation.json`, 2,000 draws per class, real epochs/bands/errors):** equal-ratio false
  rejection 0.013 / 0.017; power against SN templates at p < 0.01: hsiao (Ia) 0.52 / 0.42, nugent-sn1bc 0.59 / 0.51,
  nugent-sn2p 0.26 / 0.21, nugent-sn2l 0.21 / 0.18 (16″ / 30″ calibration). Below the plan's 0.5 for type II, so even
  a passing control gate would mostly give an uninformative "not distinguished".
- Post hoc (gate bypassed; not a decision; `posthoc_ungated.json`): calibrated ratios F090W 0.62 / 0.66, F150W 0.87,
  F200W 0.98 / 0.95; χ² p = 0.033 / 0.30. The blue band fades faster (the direction of a supernova), not significantly.
- **Conclusion:** existing public data do not distinguish an intrinsic transient from compact-deflector magnification
  for h002, and the test says nothing about which deflector. Missing: a second bright-epoch with F090W calibrated to a
  few per cent (or a spectrum); 4446 o003 (2023-05-10) gives no finite flux at h002 in our photometry (cause not checked), and the later
  epochs (2025, 2026) are non-detections. No such data are public; this path is **paused**.
- Recorded, not decisive: H2 needs a persistent source; the 2025 F150W non-detection (−2.0 ± 5.3 nJy) implies a
  magnification ≳ 15 at the first epoch under H2.

## Not done (by the plan's conditions)
- No new search boxes: no calibrated search can measure the H1/H2 difference (it needs multi-band repeat photometry,
  not more two-epoch area).

## Cost
- A: ~2 min photometry after the fix. B: control picks + photometry ~3 min; simulation 25 s; total wall for this
  round's computation < 15 min (session clock 08:33 → 08:43).
