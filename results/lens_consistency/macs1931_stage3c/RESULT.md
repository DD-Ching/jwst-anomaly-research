# MACS J1931 stage 3c: wider background box; MACS0417 gate passes, frame check F fails again (D-034 addendum af)

Plan: `PLAN.md`, frozen in e1f21a8 before any build with the new box. Code: `scripts/m1931_stage3c.py`
(`m1931_deblend.build(..., bkg_box)`, `m1931_stage3.run(..., bkg_box, extra_f)`).

| Rule | Outcome (stored in) |
|---|---|
| 3c-0 builder | background box 256 px (0.0312″/px, ≈ 8.0″); everything else as stage V2 |
| **3c-1 gate (MACS0417)** | **passes**: gate C2 shared pairs 24 / 27 = 0.89, arc-referenced 10 / 10 (`../macs0417_ownref_3c/summary.json`; V2 with the 2″ box: 28 / 35, 14 / 14) |
| 3c-2 rules 1–4 (MACS J1931) | gate K-WLS repeats (rms 0.630″, 26 / 26); 3b section (51 / 51 images inside); 826 segments → 1,164 after deblending; depth 26.66 AB (`summary.json`) |
| **3c-2 frame check F** | **fails**: 17 / 51 within 0.5″ (rule ≥ 20); medians dx 0.033″, dy −0.021″; chance control passes (shifted counts 4, 5, 7, 6, 2, 1, 3, 1; mean 3.6; 17 ≥ 7.25) |
| 3c-3 | not run: F failed, so no gate C2 and no flag reading |

Reading (frozen rule): **no stage 3c reading for MACS J1931**; not a null. The wider box did not raise the match count
(stage 3b: 18 / 51 with the 2″ box).

## Post hoc (not rules; `posthoc.json`)
- Own catalogue within 0.5 / 1.0 / 1.5″: 17 / 19 / 26 of 51 (3b: 18 / 20 / 28); spec-z 11 / 12 / 15 of 26. The JWST
  pipeline catalogue reaches 18 / 24 / 33. 15 / 51 positions lie inside an own deblended segment (3b: 14).
- Looked at (contact sheet of the 51 positions, 6″ boxes, segment outlines and catalogue rows; 1.1 MB, not stored):
  - About 20 positions show no visible F200W emission within ~0.5″: 5.1, 5.2, 6.1–6.3, 7.1, 7.2, 11.1, 11.2, 12.3,
    16.1–16.3, 17.1, 100.1–100.4, 101.1, 101.2 (a reading by eye, not measured); 13 of these 20 are spec-z images (systems 5, 6, 7, 100, 101).
  - Thin faint arcs (4.1, 9.1, 12.2, 13.x) are visible; 10.1 is now segmented.
  - With the 8″ box, cluster and bright-galaxy light stays in the image and forms large connected segments
    (around 1.1, 1.2, 6.3, 12.1, 14.x, 100.x, 101.x).
- So, by eye, a single-band F200W detection cannot reach 20 / 51 here unless most visible positions match. Both box
  sizes give 17–18. The limit is the images' F200W emission and the arcs' thinness, not the background box alone
  (not measured).

## Limitations
- Only the background box was varied (one value). The detection threshold, minimum area and band were not.
- F's minimum (20) and the chance-control ASSUMPTIONs were not changed after the run.
- The gate on MACS0417 has fewer testable pairs with the wider box (27 vs 35).
