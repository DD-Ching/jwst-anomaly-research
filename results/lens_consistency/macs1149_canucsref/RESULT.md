# MACS1149: CANUCS DR1 reference extra-image test with a containing-source rule (result)

Plan: `PLAN.md` (frozen and committed before the run, with the rule function and its tests). Numbers: `summary.json`
(all stages), `control_macs1149.ecsv`, `predicted_macs1149.ecsv`, stage-1 `imageplane_*.ecsv` / `readings_*.ecsv`,
stage-2 `frame/`, comparison `without_contain/`; post hoc `posthoc_contain.json` and `posthoc_contain/`. Scripts:
`scripts/m1149_canucsref.py` (tests `tests/test_m1149_canucsref.py`), `scripts/m1149_canucsref_posthoc.py`.

## Frozen outcome: gate K and gate C2 pass; no flag
| Step | Result | Rule |
|---|---|---|
| Gate K (23 gold images, CANUCS v1 maps) | rms 0.990″, 22 / 23 within 3.0″ (3.2 not reproduced) | ≤ 1.0″ and ≥ 0.8: passes |
| Gate C2, all controls | 15 / 15 testable recovered or confused | n ≥ 8 and ≥ 0.8: passes |
| Gate C2, arc-referenced controls | 5 / 5 | n ≥ 3 and ≥ 0.8: passes |
| Control classes (21 units) | 14 recovered, 1 confused, 2 inconsistent_reference (13.1, 13.3), 4 no_reference (system 9) | — |
| Targets (6) | 3 recovered, 1 confused, 2 inconsistent_reference; 0 absent | flag needs C2 |

- **Answer to plan question 1:** the CANUCS DR1 reference passes gate C2 on a third field (MACS1149), with the
  containing-source rule as frozen. By the plan, reference + rule may be used on Abell 370 and MACS0416.
- **The containing-source rule was not exercised here:** the comparison without it gives the same gate and the same
  class counts (`summary.json` "comparison_without_contain"). So MACS1149 validates the CANUCS DR1 reference on a
  third field, not the new rule itself.
- The local-noise rule changed no class (`without_local_noise_rule` equals the frozen counts).
- Gate K passes by 0.010″ (rms 0.990″ against 1.0″): a narrow pass.
- Frame check (reported, not a gate): 19 gold images with a CANUCS counterpart sit a median (−0.019″, −0.020″) from
  its centroid, so the assumed zero offset holds at the 0.02″ level.
- Depth (CANUCS, S/N 4–6 median) 27.77 AB; 12,680 rows with a finite F200W flux. Stage 2 (pipeline forced test,
  reported only): gate C 7 / 8; targets 1 recovered, 1 inconsistent_reference, 4 no_reference; no flag.

## Targets (sheet of all 6 targets and 6 controls looked at; not stored)
- `missing` system 3 (|μ| 17.8): `inconsistent_reference`, so not tested. The sheet shows empty sky in the 1″
  circle. Its family has the unreproduced 3.2 (gate K), so this prediction is a model copy of an image the model
  does not reproduce. That is the same pattern as MACS1423 system 1 (post hoc reading, not tested).
- `candidate` system 4 (|μ| 6.9): `inconsistent_reference`, not tested.
- `candidate` system 5 (|μ| 1.5): `confused` on a BCG-subtracted cluster member (max F / F_pred 1,154), so it is
  not informative.
- `candidate` systems 1 (|μ| 2.0), 5 (|μ| 5.2) and 110 (|μ| 3.6): `recovered`.
- System 9 controls are `no_reference`: the sheet stamp at 9.1 is blank (no mosaic data there).

## Post hoc: the containing-source rule on MACS0417 and MACS1423 (`posthoc_contain.json`; not gates)
- MACS1423 (with the local-noise rule): C2 all 10 / 11, arc 3 / 4, **still fails**. 3.2 becomes `recovered` (its own arc
  contains the position). 4.3, the BCG-adjacent central-image reference, stays `absent`.
- MACS0417 (its frozen run, no local-noise rule): C2 38 / 40 (was 37 / 40), arc 10 / 10, passes. One of the two
  family-10 flags (64.394629, −11.906917, `missing`) becomes `confused` (one containing source); the
  `candidate` flag that #455 closed as ordinary remains.
- So the rule fixes the arc-centroid failure mode it was designed for (MACS1423 3.2). It does not fix MACS1423's
  BCG-adjacent reference. A position inside a bright extended source becomes `confused`, which can hide an `absent`
  there. That is a cost of the rule, stated rather than measured.

## Deviations
- Gate K used a 0.25″ deflection grid (the `m0417_images` default, as on MACS0417 / MACS1423); PLAN.md does not state
  it (deviation). Post hoc, gate K at the 0.1″ grid of stage 2 is in `posthoc_contain.json`
  "macs1149_gate_K_grid_step": rms 0.990″, 22 / 23, passes (unchanged).

## Limitations
- With rule 4, a control that has a rule-1 counterpart always has that source in its comparison, because its own
  ellipse contains its position by definition. Gate C2 then tests the flux scale (reference × μ) and no longer tests
  detection at the position. Position agreement rests on gate K and the frame check. The `posthoc_absent_controls`
  key "own_segment_centroid_outside_radius" can no longer occur under rule 4 (from /code-review).
- One model per field (CANUCS v1 maps); no injections; no limit. 23 gold images only (the CANUCS model used 30
  gold systems from Schuldt et al. 2024, not transcribed here). MACS1149's controls never needed rule 4, so its
  validation rests on MACS1423 post hoc only. Of the three CANUCS fields, C2 passes on two.

## Next
- Abell 370 and MACS0416 with this reference + rule (their CANUCS DR1 `model/` directories are not yet checked for
  maps; a gold image list is needed per field). Expect few testable targets; most extra predictions near the core are `inconsistent_reference`
  or `confused`.
