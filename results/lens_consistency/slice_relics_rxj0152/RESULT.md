# SLICE JWST-era images vs three RELICS HST-era maps of RX J0152.7−1357 (2026-10-25)

Frozen plan `PLAN.md`, maps pinned, image list, script and rule tests (`tests/test_slice_relics_rxj0152.py`) pushed in
d72854c before any map was read. Numbers come from `summary.json`, `validate/*/validate.json` and (post hoc)
`posthoc.json`; per-image residuals in `imageplane_<model>.ecsv` (provenance `model_prediction`; positions observed,
SLICE Table 8, clump 1, 20 images in 7 families). No imaging was read.

## Frozen result
| | Lenstool v1 | GLAFIC v2 | Zitrin LTM-Gauss v1 |
|---|---|---|---|
| Gate P (κ, μ(z = 2) median ratio) | passes: 1.000 / 1.002 | passes: 1.000 / 1.002 | passes: 1.000 / 1.001 |
| Gate K0 (families 1–6, 18 images) | passes: rms 1.22″, 17 / 18 within 3″ (2.2 at 3.32″) | passes: rms 0.49″, 18 / 18 | passes: rms 0.63″, 18 / 18 |
| Frame check F (median RA·cos δ / Dec) | passes: +0.49″ / −0.47″ | passes: +0.21″ / −0.29″ | passes: +0.14″ / −0.19″ |
| Reading A (family 7, 2 images) | rms 0.22″, 2 / 2 | rms 0.34″, 2 / 2 | rms 0.39″, 2 / 2 |
| Reading B (family 1 at z_spec 3.93, 4 images) | rms 1.51″, 4 / 4 | rms 0.50″, 4 / 4 | rms 1.29″, 4 / 4 |
| Rule 7 amplitude factor, median (range), 7 families | 1.045 (0.805–1.343), `near_1` | 1.028 (0.730–1.390), `near_1` | 1.045 (0.768–1.177), `near_1` |
| Rule 6, more predicted than catalogued | families 4 (2 → 4), 7 (2 → 4) | 1 (4 → 6), 4 (2 → 3) | 1 (4 → 5), 4 (2 → 5), 7 (2 → 3) |

All three models are usable. Joint flags (rule 5): none in reading A, none in reading B. **Q1: `H_none`.**

## Answers to the plan's questions
- **Q1:** the RX J0152 Lenstool v1 map needs no amplitude factor (median 1.045, as GLAFIC 1.028 and Zitrin 1.045). The
  Lenstool v1 offset of SMACS 0723 and MACS J0553 is therefore not a property of every RELICS Lenstool v1 map: three
  of five clusters tested so far need none (RXC J0232, RX J0152; RXC J2211 post hoc on a model that failed F, addendum bc). The derived practical rule of addendum
  bd stays as written (check the amplitude against a second model or z_spec families before reading a Lenstool v1
  map); it is a per-cluster check, not a known global scaling.
- **Q2:** no image misses every usable model by > 3″. The one miss (2.2 in Lenstool v1, 3.32″) is a single-model
  miss inside a family that the other two models fit within 3″. Reading A on family 7 was declared uninformative in
  advance (image separation 1.45″ < 3″ match radius); it is listed, not read.

## Post hoc (not in the plan)
`scripts/slice_relics_rxj0152.py --posthoc` lists, at each model's fitted z, the predicted images that match no
catalogued image (`posthoc.json`):
- **Family 4** (SLICE model z 1.69; pair 4.1 / 4.2 at predicted |μ| 26–27 (Lenstool), 11–13 (GLAFIC), 18 (Zitrin),
  i.e. a pair straddling a critical curve): **all three models predict a further positive-parity image
  16.7–22.2″ from the pair and 21.4–23.2″ from the cluster centre** (SLICE Table 1), at μ +3.81 (Lenstool,
  28.17895, −13.95293), +2.38 (GLAFIC, 28.17805, −13.95268) and +2.39 (Zitrin, 28.17945, −13.95343). The three
  positions lie 2.5–5.6″ apart. Predicted pair-to-this |μ| ratios: 6.9–7.2 (Lenstool), 4.7–5.3 (GLAFIC), 7.5–7.7
  (Zitrin). SLICE Table 8 lists no such image. Lenstool and Zitrin also predict fainter images 4.1″ from the pair
  (μ −2.59; −0.32 and +0.16), not read.
  Ordinary explanations not tested yet: below the imaging depth at that flux ratio, blended with a cluster member,
  outside the SLICE footprint, model position scatter, or a pair that is not a multiple-image system. This is a
  lensing-residual **input** for imaging, not a candidate.
- Family 1 (GLAFIC, Zitrin): the unmatched predictions are central demagnified images (listed in `posthoc.json`;
  not read). Family 7 (Lenstool, Zitrin): unmatched predictions next to the pair (galaxy-scale feature; not read).

## Reading
- Three HST-era RELICS models, built on the six HST families, reproduce SLICE's JWST positions of those families,
  the one spectroscopic family at its z_spec and the one JWST-new family; no flag goes to vetting.
- Practical: on RX J0152, a RELICS Lenstool v1 map gives a usable reading after the amplitude check.
- **Next:** family 4's predicted third image: forced photometry at the three predicted positions in the SLICE
  NIRCam imaging, with the pair's measured flux scaled by the predicted μ ratio as the expected flux (frozen plan
  first; ordinary explanations above in that order).
