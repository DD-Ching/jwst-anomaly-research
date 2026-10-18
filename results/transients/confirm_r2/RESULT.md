# Confirmation round 2: result (plan `PLAN.md`, test choice `B_choice.md`, both committed before the new data)

Session-clock times in the plans disagree with git timestamps by a few minutes; the commit order is the record.

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

## B. Is the B1-PLCKG165-box2-h002 fade achromatic? Not decided; low power even if it had been
- H1 intrinsic transient (colour evolves) vs H2 magnification of a persistent compact source by a point mass, an
  Ellis wormhole (n = 2, ε > 0) or a negative-mass lens (n = 1, ε < 0), achromatic in all three cases for a
  uniform-colour source. The test cannot exclude a chromatic caustic crossing of a source with a colour gradient and
  cannot separate the three deflectors.
- **Frozen rule: not decided** (too few usable controls). 16″ pick: 10 / 15 / 15 in F090W / F150W / F200W, all
  below 20 (`g165_h002/calibration.json`). Deviation (30″ pick, recorded before any target ratio; commit order
  d299909 plan < a40052c photometry < 72fdec7 deviation): 14 / 31 / 32, F090W below 20 (`g165_h002_30/`). h002 is
  near the edge of the 4446 coverage: 4446 o002 is finite for 16 of 24 positions in the 16″ set (33 / 63 in the 30″
  set) and o003 for 13 / 24 (31 / 63); o003 has no finite flux at h002.
- **Identifiability (`simulation_plan.json`, `simulation_excess.json`; 2,000 draws per class at the real epochs,
  bands and errors; phases from each template's B-band peak).** "plan": the calibration scatter as the plan states
  it (this double-counts the controls' own noise, so it is conservative); "excess" (post hoc, after review): scatter
  minus the controls' median ratio error in quadrature. Power at p < 0.01, 16″ / 30″ calibration:
  - plan: hsiao (Ia) 0.53 / 0.43, nugent-sn1bc 0.35 / 0.29, nugent-sn2p 0.21 / 0.17, nugent-sn2l 0.17 / 0.15;
  - excess: hsiao 0.62 / 0.54, nugent-sn1bc 0.41 / 0.36, nugent-sn2p 0.29 / 0.22, nugent-sn2l 0.21 / 0.17;
  - equal-ratio false rejection 0.006–0.011.
  The plan's "power < 0.5" condition was not written per class; read per class (interpretation), every template
  except hsiao is below 0.5 in both models and both calibrations, and hsiao is below 0.5 in the 30″ plan case. So
  even a decided "not distinguished" would have been weak for most supernova types.
  (An earlier version of these numbers, with multiplicative calibration noise and phases from template zero, is in
  git history at 364f546; superseded.)
- Post hoc (gate bypassed; not a decision; `posthoc_ungated.json`): calibrated ratios F090W 0.62 / 0.66, F150W 0.87,
  F200W 0.98 / 0.95; χ² p = 0.033 / 0.30. The blue band fades faster (the direction of a supernova), not significantly.
- **Conclusion:** the frozen test did not decide, and with these two epochs it has low power against most supernova
  templates. The existing JWST epochs covering h002 (#296 photometry) and the 2016 HST limit do not give a usable
  third bright epoch; other archives were not searched. Missing: a bright-epoch repeat with F090W calibrated to a few
  per cent, or a spectrum. This path is **paused**.
- Recorded, not decisive (not stored; 240 / (3 × 5.3)): H2 needs a persistent source; the 2025 F150W non-detection
  (−2.0 ± 5.3 nJy) implies a magnification ≳ 15 at the first epoch under H2.
- Photometry note: with `--group-arcsec 90` the forced-photometry error is one robust value over the group box
  (here ~60″), not local; h002's F150W error is 5.4 nJy (16″ set) vs 6.1 nJy (30″ set).

## Not done (by the plan's conditions)
- No new search boxes: no calibrated search can measure the H1/H2 difference (it needs multi-band repeat photometry,
  not more two-epoch area).

## Cost
- Not stored (session observations): A ~2 min photometry after the fix; B picks + photometry ~3 min, each simulation
  < 1 min; all computation of this round < 20 min wall.
