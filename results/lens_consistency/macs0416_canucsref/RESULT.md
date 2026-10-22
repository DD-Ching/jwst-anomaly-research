# MACS0416: CANUCS DR1 reference extra-image test with the containing-source rule (result, flag not yet vetted)

Plan: `PLAN.md` (frozen and pushed before the run, ec5834f). Numbers: `summary.json`; tables `control_macs0416.ecsv`,
`predicted_macs0416.ecsv`, stage-1 `imageplane_*.ecsv` / `readings_*.ecsv`, stage-2 `frame/`, comparison
`without_contain/`. Script: `scripts/m0416_canucsref.py` (tests `tests/test_m0416_canucsref.py`).

## Frozen outcome: gate K and gate C2 pass; one flag (D115b), not yet vetted
| Step | Result | Rule |
|---|---|---|
| Gate K (303 model images, `macs0416-canucs`, 0.25″ grid) | rms 0.512″, 303 / 303 within 3.0″ | ≤ 1.0″ and ≥ 0.8: passes |
| Gate C2, all controls | 160 / 172 testable recovered or confused (0.930) | passes |
| Gate C2, arc-referenced controls | 58 / 59 (0.983) | passes |
| Targets | 15 recovered, 8 confused, 1 absent, 6 inconsistent_reference, 8 no_reference, 6 undetectable, 5 undetectable_local | — |
| Flag (rule 8) | `D115b` @ (64.048551, −24.064825), |μ| 26.6, m_pred 24.45, max F / F_pred 0.049 over 3 sources | — |

- **Plan question 1:** the CANUCS DR1 reference with both rules passes gate C2 on a fourth field, with 172 testable
  controls (MACS1149 had 15).
- The flag is the same without the local-noise rule; without the containing-source rule there is a second flag,
  `D1f` @ (64.040935, −24.071016), that the rule turns into a non-`absent` class (`comparison_without_contain`).
- Frame check (reported): 283 images with a CANUCS counterpart sit a median (−0.011″, +0.014″) from its centroid.
- Stage 2 (pipeline forced test, reported only): gate C 188 / 191; no flag.

## First look at the flag (cutout looked at, not stored; not a vetting)
- D115b is the second knot of the D115 source (z 2.2999; D115a has three catalogued images, D115b two). The model
  predicts three extra D115b images; by the positions in `predicted_macs0416.ecsv`, the flagged one lies about 0.6″ from a `recovered` D115a extra prediction
  (|μ| 20) and about 1.6″ from a `recovered` D115b extra prediction (|μ| 22), on BCG-subtraction residuals near a
  cluster galaxy. |μ| 20–30 there means a galaxy-scale critical curve, where the model's positions are least
  reliable. That is an ordinary reading to test, not a conclusion.
- Not yet done: `/vet-candidate` (cheapest first: model position error near the critical curve, e.g. the D115a /
  D115b relative geometry; the CANUCS posterior or an independent model; deblending and BCG-model residuals). Until
  then this is an unvetted screen hit, not a candidate.

## Vetting, 2026-10-22 (`docs/candidates/macs0416-canucsref-D115b.md`)
Verdict `inconclusive`: the extra D115a / D115b images are made by cluster member 8515 (z_spec 0.3907); they vanish
at σ × 0.7–0.8, but the χ² cost of that is unresolved (0.25″ `validate`: 330.8 → 5261 from a D115a match failure; a
0.1″ solve keeps D115a matched). Next: the fine-grid `validate` in the record's "Open questions".
