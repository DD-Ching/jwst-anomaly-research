# W3 ZTF fade-control rate across areas: result

Plan: `PLAN.md` (committed before computing). Scripts: `diag_script.py.txt` (frozen tests → `diag.json`),
`parent_script.py.txt` (post hoc → `parent.json`). Offline; committed tables only. No method change.

## Frozen tests (derived)
- **Areas:** 25 (v16 8, v20 7, v21 2, v22 8), each with n = 1,000 fade-control trials. Pooled rate 9.08 %.
- **Homogeneity:** χ² = 233.8, 24 dof, p = 2.6 × 10⁻³⁶. The per-area rates (3.9–13.9 %) are **not** binomial
  scatter around one rate (σ ≈ 0.9 points per area).
- **Inventory predictors** (Spearman with the rate; pass needs |ρ| ≥ 0.6 and p < 0.0014, ASSUMPTION):
  median exposures −0.14, infobits −0.35 (p 0.084), product ratio −0.22, worst-cell score −0.28, mean necessary
  fraction −0.23, ecliptic latitude −0.07, Galactic latitude −0.03. **None passes.**
- **Outcome C** (as frozen): the area-to-area spread is real, and the pre-fetch inventory metadata does not predict
  it. v16-inventory areas median 7.3 % vs v22 11.3 % (Mann–Whitney p 0.018; descriptive only).

## Post hoc (pattern seen after the frozen tests; not a tested rule)
- Grouping the order-5 partitions by their HEALPix order-4 parent (npix // 4, ≈ 13 deg²): **between parents
  χ² = 207.9 / 11 dof (p 1.7 × 10⁻³⁸); within parents 25.8 / 13 dof (p 0.018)**. About 89 % of the excess spread
  sits between sky regions; siblings agree to within binomial-plus-small scatter. Examples: parent 622 (2488–2491)
  10.6–13.9 %, all fail; 2031 (8124/8126/8127) 5.4–7.4 %, all pass the fade limit; 339 (1356–1358) 3.9–6.6 %.
- The mechanism is undetermined: no stored sampling metric tracks it (above). Seeds differ per area, so a shared
  seed is not the cause. Candidates (hypotheses, untested): the field / CCD-quadrant season pattern that siblings
  share, or the carrier stars' light-curve noise in that region.
- **v22 remaining partitions** (17, `parent.json`): 2345–2347 are siblings of 2344 (8.1 %, the only v22 pass);
  1149, 1852/1854/1855, 4988/4990 and 5052/5053 are siblings of areas at 10.2–13.8 % (all failed);
  2812/2813 and 2876–2879 have no tested sibling.

## What this changes
- Nothing in the frozen calibration, the 10 % limit or any result. The v22 "7 of 8 fail" mostly reflects that
  the v22 queue drew several partitions from two high-rate parents (287, 622 in v20; 287, 1247 in v22).
- **Testable prediction for the next round (to be frozen in its plan, not adopted here):** siblings of a
  fade-passing parent pass the fade control more often than siblings of a failing parent. A v23 queue ordered
  "untested parents and passing-parent siblings first" would test it (2345–2347 predicted ≤ 10 %; 1149 and
  4988/4990 predicted > 10 %) and avoid spending fetch budget on predicted failures. Whether the queue rule
  changes is the next round's frozen-plan decision; the fade limit itself stays an owner decision.
