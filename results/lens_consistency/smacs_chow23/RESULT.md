# SMACS 0723 Chow+2023 maps: input check and fixed-z gate K (2026-10-24)

Plan frozen in 014dc2c (`PLAN.md`) before any map value was read. Script `scripts/smacs_chow23.py` (tests
`tests/test_smacs_chow23.py`); numbers from `summary.json`; tables `imageplane_fixedz.ecsv` (mean source) and
`imageplane_fixedz_wls.ecsv` (WLS source). Provenance: image positions observed (Mahler+22 ICLv2 `arcs.dat`),
predictions model_prediction, verdicts derived.

**Deviation (procedure):** the rule functions and tests were written before the run, but their commit failed on a
line-length lint error and the run went ahead; they were committed after it, changed only by docstring rewrapping
and the model description below. **Plan text correction (post hoc):** the plan calls the model "free-form"; the paper
(arXiv:2310.09790, HTML read 2026-10-24) builds it with GLAFIC (parametric). No rule depends on this.

| Rule | Found | Status |
|---|---|---|
| C coverage | 16 / 16 images inside the 3260² CAR maps | **passes** |
| U unit and x sign | median per-system rms: A 0.767″, A_neg_x 2.394″, P 11.345″, P_neg_x 11.346″ (≤ 1.0″) | **passes**, reading A (arcsec, α_x as given) |
| K fixed-z gate (mean source) | rms 10.18″, 12 / 16 (0.75) within 3″ (needs ≤ 1.0″ and ≥ 0.8) | **fails** |
| K-WLS (reported, not a gate) | rms 1.947″, 15 / 16 within 3″ | — |

- **Per system** (mean source / WLS rms; catalogued → predicted images): 1 0.53″ / 0.41″ (3 → 3); 2 0.70″ / 0.48″
  (3 → 3); 3 0.75″ / 0.67″ (4 → 5); **5 22.98″ / 0.62″ (3 → 1)**; **7 4.82″ / 4.34″ (3 → 1)**. System 5's failure is
  the mean-source estimator (one predicted image; the WLS source reproduces it), as for Bullet and S1063; system 7
  fails with both sources.
- **ICLv2 in the same run** (the harness check): gate K rms 0.321″, 16 / 16, equal to D-034 addendum r; K-WLS 0.293″.
  So the failure is the Chow model on this list, not the harness. Systems with different predicted counts in the
  two models (mean source): 5 and 7 (`predicted_count_differs`).

## Post hoc (not in the plan; `posthoc.json`, `scripts/smacs_chow23_posthoc.py`, same map sha256)
- Free source z per system (grid 0.5–10.0, step 0.01): systems 1, 2, 3, 5 have minima at z 1.44, 1.38, 1.97, 1.47
  (catalogue 1.449, 1.3779, 1.9914, 1.425), rms 0.51–0.76″. **System 7's minimum lies on the grid edge, z 10.0, rms
  1.08″**; at its catalogue z 5.17 the rms is 1.90″ (`summary.json` U, reading A).
- Without system 7: mean source rms 11.06″ (13 images, 0.85 within 3″; system 5 again), WLS 0.566″, 13 / 13.
- The paper lists five spectroscopic systems at 1.38 ≤ z ≤ 2.21 (Golubchik+22, Mahler+23, Noirot+23); a z 5.17 value
  appears in its table as a fitted (geometric) redshift of another candidate system (read by eye, row identity not
  established, not stored). So Mahler+22 system 7 at z 5.17 is a fixed model redshift in ICLv2's `input.par`, not a
  spectroscopic one as rule 0 of D-034 addendum r assumed (that rule is labelled ASSUMPTION there). The paper's own
  rms on its spectroscopic images is 0.492″ (paper, not stored).

## Reading
- The Chow+2023 maps load (unit arcsec, α_x along +x = West, D_ds / D_s = 1, CAR grid) and reproduce the four
  Mahler+22 systems with spectroscopic z at ≤ 0.76″ per system with the WLS source; the frozen gate fails, from
  system 5 (mean-source estimator) and system 7 (no z ≤ 10 brings it below 1.08″ in this model).
- Ordinary explanations of system 7 not tested: different image identification of that system in Chow+23 (their
  numbering differs from Mahler's), the z 5.17 itself being a model value, or the GLAFIC model's mass in that region.
  Nothing here is an anomaly: two published models disagree on one photometric-z system.
- **Route:** stopped as frozen. A next stage could gate on the four spectroscopic Mahler systems with the WLS source
  only as a labelled deviation citing these numbers (addendum k precedent), then compare the two models' predicted
  extra images; the SMACS reference route itself is closed (D-034 addendum s), so that stage would need a validated
  reference first.
