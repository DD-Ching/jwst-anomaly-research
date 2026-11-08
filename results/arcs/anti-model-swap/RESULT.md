# `arcs`: pooled anti rows under a second published model (2026-11-08)

Plan `PLAN.md` (frozen in the first commit, before any second-run class of a primary anti row was read); script
`scripts/arcs_model_swap.py`, rule functions unit-tested before the run (`tests/test_arcs_model_swap.py`). Every
number below is in `summary.json`. Stored tables only (no new data, no network). Provenance: observed axes from the
catalogues, predicted axes `model_prediction`, classes and counts `derived`; thresholds ASSUMPTIONs.

## Frozen outcome: **`mixed`** (rule 3: gate passes, m = 6, k = 4; `persists` needs k ≥ 5, `model-dependent` k ≤ 3)
| Rule | Found | Status |
|---|---|---|
| 1 control: primary aligned rows still `aligned` in the second run | 66 / 76 matched (0.868); per field 0.73–1.00 | **passes** (≥ 0.70) |
| 2 tested anti rows | m = 6 matched of 6 (MACS0417 3, WHL0137 1, RXC J0600 1, MACS0308 1) | — |
| 3 anti in the second run | k = 4 | **`mixed`** |

| Field | Row | Primary g_min | Second-run class |
|---|---|---|---|
| MACS0417 (CANUCS v1 → RELICS Lenstool v1) | 1993 | 3.69 | aligned |
| MACS0417 | 2341 | 0.26 | mixed |
| MACS0417 | 4071 | 1.18 | anti |
| WHL0137 (RELICS Lenstool → GLAFIC) | 2194 | 0.24 | anti |
| RXC J0600 (RELICS GLAFIC → Lenstool) | 581 | 0.39 | anti |
| MACS0308 (RELICS GLAFIC → Zitrin LTM) | 2157 | 0.64 | anti |

MACS1423 (4 anti) and RXC J0018 (0 anti) have no second run: untested.

## Secondary (in the plan; not rules)
- Aligned → anti flip rate of the control: 5 / 76 (0.066). The anti rows keep their class (4 / 6) far more often than
  an aligned row turns anti, so the anti class of these rows is mostly not a one-model artefact (post hoc reading of
  the secondary numbers, small counts).
- Variant B with the second runs swapped in for the four fields: 13 anti vs 2.63 expected (σ_ε 0.26), P 4.1 × 10⁻⁶,
  against 6 vs 2.27 (P 0.029) with the deciding runs. Most of the rise is MACS0417's HST-era RELICS Lenstool run (8
  anti vs 1.42 on 61 rows; its strong-shear set differs from the CANUCS run's 49). The excess is not an artefact of
  which published model decides.

## Reading (narrowest statement)
- H_model (the pooled anti class comes from deciding-model shear errors) is not supported for the four anti rows
  outside MACS0417 (all persist) and is possible for two MACS0417 rows (1993 lies deep inside the critical curve,
  g_min 3.69, where axes are model-sensitive). The frozen label is `mixed`.
- Together with Abell S1063 (4 / 4 anti persist in three HFF models, `../abells1063-cats/second_models.json`), model
  choice does not explain the pooled excess. Still untested ordinary explanations: blends and clump pairs (the
  MACS1423 census explains that field), photo-z contamination of the background cut, and σ_ε at strong shear.
  An anomaly score is not evidence of new physics; nothing here is a candidate.

## Limitations
- Second runs are not independent of the data: they share the catalogue, photo-z and shape measurement; only the
  predicted axis changes. The test addresses model error only.
- m = 6; one row more or less moves the label. Two of the four second runs are HST-era RELICS maps.
- Exploratory line (the pooled excess itself was found post hoc, `../pooled-tally/`).
