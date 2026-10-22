# MACS1423: extra-image test with CANUCS DR1 reference fluxes and a local-noise rule (result)

Plan: `PLAN.md` (frozen and committed before the run, 6c114b7; rule-2 function and tests committed before the run).
Numbers: `summary.json`, `control_macs1423.ecsv`, `predicted_macs1423.ecsv`; post hoc `posthoc.json`. Scripts:
`scripts/m1423_canucsref.py` (runs `m0417_arcref.run` on MACS1423, then rule 2; tests
`tests/test_m1423_canucsref.py`), `scripts/m1423_canucsref_posthoc.py`.

## Frozen outcome: gate C2 fails on the arc-referenced controls; no flag is read
| Step | Result | Rule |
|---|---|---|
| Gate C2, all controls | 9 / 11 testable recovered or confused (0.82) | n ≥ 8 and ≥ 0.8: passes |
| Gate C2, arc-referenced controls | 2 / 4 (0.50) | n ≥ 3 and ≥ 0.8: **fails** |
| Control classes (14 units) | 9 recovered, 2 absent (3.2, 4.3), 3 inconsistent_reference (4.1, 4.2, 4.4) | — |
| Predicted extra images (2) | system 1 `missing` (|μ| 17.3) `absent`; system 1 `candidate` (|μ| 8.8) `recovered` | flag needs C2 |

- Rule 2 (local noise) changed no class: every testable prediction is 36–2,266 σ_loc (`pred_over_sigma`), so the
  outcome equals the run without it (`summary.json` "without_local_noise_rule").
- Depth (S/N 4–6 median) 27.92 AB; pixel scale 0.0400″; 12,118 rows with a finite F200W flux.
- **Answer to plan question 1:** the CANUCS DR1 reference does not pass gate C2 on MACS1423 as frozen, so MACS0417's
  pass (37 / 40, 10 / 10) is not shown to transfer. Per the plan, other CANUCS fields need a new gate before this
  reference is used there.
- The 2 `demagnified` central predictions stay excluded (rule 8), as planned.

## Post hoc (not frozen rules; `posthoc.json`; contact sheet of the four items looked at, not stored)
- **3.2 `absent`:** its position lies inside CANUCS 4101312 (ellipticity 0.84, the long arc on the sheet), whose
  centroid is 1.02″ away, just outside the 1.0″ search radius; that source holds 1.31× the predicted flux. Same
  failure mode as the MACS0417 pipeline-segment run (own centroid outside the radius).
- **4.3 `absent`:** inside CANUCS 4100458 (0.07″), flux 0.324× predicted, just below the 1/3 bound. Its reference is
  4.5, the central image (|μ| 1.26) beside the BCG (σ_loc 37 nJy there vs 4–9 elsewhere); a reference holding BCG
  residual light, or a model μ error for the central image, would both raise F_pred. Not separated here.
- **System 1 `missing` (|μ| 17.3):** no source's ellipse contains the position and the 0.3″ aperture holds 0.5 nJy
  against a predicted 1,291 nJy (σ_loc 8.1); the sheet shows empty sky, with elongated light ≳ 2″ away. This
  prediction is a model copy of system 1 whose catalogued images 1.3–1.5 the CANUCS v1 joint source does not
  reproduce (`../macs0417_extra/summary.json` "unpredicted"; image-plane rms 1.59″ in `../macs0417_images`), so
  model position error is the ordinary reading; the earlier pipeline forced test `recovered` it with flux found
  2.2″ from the position (`../macs0417_extra/predicted_macs1423.ecsv` best_dx/dy). Not a candidate; with C2 failing
  it is not read.

## Limitations
- One model (CANUCS v1 maps); no injections; no limit. Frame offset assumed 0. PA convention as in
  `../macs0417_canucsref` (not re-tested here). Only 4 arc-referenced controls, so the arc gate is coarse.

## Next
- The arc failure (3.2) is the known centroid-outside-radius mode; a rule that also accepts the source whose
  ellipse contains the position would address it, but must be frozen and gated on a field other than MACS0417 and
  MACS1423 (both now seen). MACS1423 extra-image work is otherwise done: its 2 testable predictions are 1 recovered
  and 1 model copy of an unreproduced system.
