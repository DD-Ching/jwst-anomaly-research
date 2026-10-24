# MACS J1931 stage 3d: continuum-restricted F passes; gate C2 passes; 0 of 11 extra-image predictions testable (D-034 addendum ag)

Plan: `PLAN.md`, frozen in 1314a97; rule functions and tests in bce45c1, both before any run. The rule 3d-4 wiring
(`use_fc`) was committed after the 3d-2 / 3d-3 runs and before its own run (labelled in PLAN.md). Code:
`scripts/m1931_stage3d.py`. Builder: stage V2 / 3b (64 px background box).

| Rule | Outcome (stored in) |
|---|---|
| **3d-2 validation (MACS0417)** | **passes**: 46 / 56 images continuum-detected (7 faint, 3 off the section); 41 / 46 matched within 0.5″; medians dx 0.022″, dy −0.006″; chance mean 3.75 (`../macs0417_ownref_3d/summary.json`, `images_fc.ecsv`) |
| **3d-3 F_c (MACS J1931)** | **passes**: 18 / 51 continuum-detected (33 faint); 17 / 18 matched; medians 0.040″ / 0.017″; chance mean 1.125 (`summary.json`, `images_fc.ecsv`) |
| 3d-4 rule 2 | gate K-WLS repeats (rms 0.630″); catalogue 1,153 rows, the same count as 3d-3; the stage 3 F alone would fail (18 / 51) (`run_summary.json`) |
| **3d-4 gate C2 (shared pairs)** | **passes**: 11 / 11, arc-referenced 7 / 7 (`run_summary.json`) |
| 3d-4 targets | 11 predicted extra images: 10 `no_reference`, 1 `undetectable_local` (system 1, predicted 1.11 σ); **0 testable, 0 flags** (`predicted_macs1931.ecsv`) |

Reading (frozen rules): no flag. This is **not a null** for MACS J1931: no predicted extra image could be tested.
- The 10 `no_reference` targets belong to systems 5, 6, 7, 100 and 101.
- None of their 13 observed images has F200W continuum: aperture S/N −1.2 to 1.4 (`images_fc.ecsv`). So no flux can
  be predicted from a sibling.

## Looked at
Contact sheet of the 11 targets (4″ boxes, background-subtracted F200W; 30 kB, not stored): no source at any target
position. The nearest objects are compact sources 0.5–1.5″ away and a star with diffraction spikes. One system-101
panel lies on a noisy section edge.

## Limitations
- F_c's threshold was chosen after the 3c post hoc reading (PLAN.md "Not blind"). The MACS0417 validation and gate C2
  are the quality gates.
- F200W only. Systems 5, 6, 7, 100 and 101 have spec-z but no continuum in this band. This stage cannot test their
  extra images in any band.
- Gate C2 has 11 testable controls (7 arc-referenced). That is fewer than on MACS0417 (35 in V2).

## Next
The F200W own-reference route for MACS J1931 is finished: no target was testable. A test of the five continuum-free
spec-z systems needs a line-emission reference (e.g. a public MUSE cube) and its own frozen plan. Otherwise move to
another field with public JWST-era maps.
