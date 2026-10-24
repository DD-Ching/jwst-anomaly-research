# MACS0417 stage V: own deblended F200W reference catalogue against the frozen gate C2 (result)

Plan: `../macs1931_counter/PLAN.md`, "Stage V: frozen rules" (ff990f1) and "Stage V2: amendment" (12f91c2), both
committed before the run they govern. Script `scripts/m1931_deblend.py` (rules unit-tested in
`tests/test_m1931_deblend.py` before each run). Numbers: `summary.json` and the ECSV tables here (V1) and in
`../macs0417_ownref_v2/` (V2).

## Frozen outcomes
| Step | Result | Rule |
|---|---|---|
| Catalogue (both runs) | 1,232 detected segments → 1,577 after deblending, all kept; section y 1200–4948, x 5277–9354; 9.1 % masked | V rules 1–4 |
| V1 depth | **undefined**: 0 rows at Kron S/N 4–6 (minimum S/N 8.78) | V rule 4 |
| V1 gate C2 (shared pairs) | 28 / 35 (0.80), arc-referenced 14 / 14: **not a valid reading** (with a NaN depth `classify` cannot return `absent` or `undetectable`; 7 controls became `ambiguous`) | V rule 6 |
| V2 depth | 27.78 AB (5 × median Kron error) | V2-a |
| **V2 gate C2 (shared pairs, local classes)** | **28 / 35 = 0.80, arc-referenced 14 / 14: passes** (at the threshold; no merged groups) | V2-c |
| V2 control classes | 25 recovered, 3 confused, 7 absent, 10 inconsistent_reference, 9 no_reference | — |

- **Reading (frozen rule):** the builder is accepted for the MACS J1931 stage. CANUCS DR1 gave 37 / 40 on the same
  controls (`../macs0417_canucsref`), JWST pipeline segments 10 / 27 (`../macs0417_arcref`).
- **Margin is zero:** 28 / 35 is exactly 0.80, so one more non-ok testable control would fail the gate. Any MACS J1931
  flag from this builder carries that caveat.
- MACS0417 predicted-target classes are written, not read (V2-c; the field's flags belong to `../macs0417_canucsref`).

## The 7 absent V2 controls (`../macs0417_ownref_v2/control_macs0417.ecsv`)
- 1a.1, 1b.1, 1c.1, 1m.1: knots of the bright lensed spiral; best flux ratio 0.238–0.305, just below the 1 / 3 limit
  (predicted 18.8–19.1 AB). The same family-1 pattern as the pipeline route (`../macs0417_arcref`: shared segment
  holds four families). Contact sheet: the spiral is split into several deblended segments; which piece holds
  which knot was not measured.
- 8.1: best ratio 0.029; its reference 8.2 is `confused` (ratio 34.7), on a compact galaxy (as in `../macs0417_arcref`).
- 15.1, 15.2: no source in the radius, 25.3σ and 12.0σ local. Contact sheet: 15.1 is off the F200W footprint; route
  rule 6 has no footprint check (known defect, `../macs0417_arcref/RESULT.md`), so it counts as a failure here. That
  direction can only lower the fraction.

## Looked at
Contact sheet of the 54 controls (6″ boxes) on the section with the deblended segment outlines (not stored): arcs 3.x,
12.x and 6.x are single elongated segments and recovered; family 1 is split into several segments; 15.1, 16.2 and 16.3
lie off the image.

## Limitations and deviations
- **Deviation:** V2 is an amendment after V1 (plan defect: the route's depth proxy needs S/N 4–6 rows, which a 3σ /
  25-pixel detection does not produce). V2 was frozen before any V2 number; its depth and chain (containing rule on,
  local noise) are the planned MACS J1931 stage 3 chain.
- One field, one band, one model (CANUCS v1 maps via `../macs0417_frame`); no injections. The builder's parameters
  (2″ background box, 3σ, 25 px, contrast 0.001) were fixed before the run and not varied.
- The first V1 run stopped on a NaN in the summary JSON (the undefined depth); the stored V1 files are from the
  re-run with non-finite values written as null (same classes).
- After the runs (code review, no rule change): photutils pinned `<4`, an empty-catalogue guard, the row filter
  written once (`keep_rows`), the cache under `$JWST_ANOMALY_DATA`, a tighter PA test and a row-filter test. The
  footprint gap is carried into the stage 3 plan as a required rule.
