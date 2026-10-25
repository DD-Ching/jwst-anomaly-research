# SLICE JWST-era images of MACS J0358.8−2955 vs two RELICS HST-era maps of Abell 3192 (2026-10-25)

Frozen plan `PLAN.md`, map pins, image list, script and rule tests pushed in 56a0b3bb before any map was read.
Numbers come from `summary.json`, `validate/*/validate.json` and (post hoc) `posthoc.json`; per-image residuals in
`imageplane_<model>.ecsv` (provenance `model_prediction`; positions observed, SLICE Table 8, 16 images in 6
families). No imaging was read.

## Frozen result
| | Lenstool v1 | GLAFIC v2 |
|---|---|---|
| Gate P (κ, μ(z = 2) median ratio) | passes: 1.000 / 1.000 | passes: 1.000 / 1.000 |
| Gate K0 (families 1–3, 7 images, fitted z) | passes: rms 0.25″, 7 / 7 within 3″ | passes: rms 0.96″, 7 / 7 |
| Frame check F (median RA·cos δ / Dec) | passes: −0.01″ / −0.05″ | passes: +0.24″ / −0.17″ |
| Reading A (families 4–6, 9 images, fitted z; exploratory) | rms 0.88″, 9 / 9, no miss | rms 1.34″, 9 / 9, no miss |
| Reading B (families 1–3 at z_spec, 7 images) | rms 8.59″, 5 / 7; misses 1.1, 1.2 | rms 1.45″, 7 / 7 |
| Rule 7 amplitude factor, median (range), 6 families | 0.692 (0.674–0.721), `lenstool_like` | 0.999 (0.777–1.030), `near_1` |
| Rule 6, more predicted than catalogued | families 2, 3 (2 → 3), 5 (3 → 5) | all six (1: 3 → 5; 2, 3: 2 → 7; 4: 3 → 5; 5: 3 → 7; 6: 3 → 8) |

Both models are usable. Joint flags (rule 5): none in reading A; none in reading B (1.1 and 1.2 are missed by
Lenstool v1 only; GLAFIC matches them). **Q1: `H_sys`.**

## Answers to the plan's questions
- **Q1:** the Abell 3192 Lenstool v1 map needs a factor 0.674–0.721 in all six families, while GLAFIC v2 does not
  (median 0.999). This is the third RELICS Lenstool v1 map with such an offset (SMACS 0723, MACS J0553, now
  Abell 3192), against two without (RXC J0232, RX J0152; RXC J2211 failed F). The addendum-bd rule (check the
  Lenstool v1 amplitude against a second model or z_spec families before reading it) stays a per-cluster check.
- **Q2:** no image misses every usable model. Lenstool v1's reading-B misses of 1.1 and 1.2 (16.2″ and 15.6″ at
  z_spec 1.261) match at its fitted z 0.76 (0.41″, 0.48″), so they belong to the amplitude offset of Q1, not to a
  position residual. Reading A is exploratory, because the RELICS input lists are not published (PLAN.md).
- Weakness stated in advance and seen: GLAFIC fits families 2 and 3 (z_spec 3.07; close pairs) at z 1.26 / 1.32,
  yet matches them at z_spec too (reading B). These pairs constrain z weakly. Not re-weighted.

## Post hoc (not in the plan)
- **Amplitude as a reference redshift** (`posthoc.json` `amplitude_as_dls_ds_zref`; derived): the Lenstool v1
  median factor 0.692 equals D_LS/D_S at z_ref 2.09 for z_lens 0.425 (flat, H0 70, Ωm 0.3). Such a value would arise if
  the map were normalised at a source z ≈ 2 rather than at D_LS/D_S = 1, as its readme states. This is a
  **hypothesis** about the map product, not tested here. It could be tested on SMACS 0723 and MACS J0553 (with
  the same rule) or against the RELICS team.
- Unmatched predicted images at the fitted z are listed per model in `posthoc.json` (GLAFIC: up to 5 extra per
  family; Lenstool: families 2, 3, 5). GLAFIC's brightest one is family 4: μ 16.1, 4.2″ from the nearest catalogued
  image; Lenstool v1 predicts no such image. It is single-model, at |μ| typical of a merging pair, so it is not read.
  None is a lensing-residual input under the two-model rule.

## Reading
- Two HST-era RELICS models reproduce SLICE's JWST-era positions of the three Hsu et al. 2013 families and of
  the three SLICE-added families. GLAFIC also reproduces the spectroscopic families at z_spec. No flag goes to
  vetting.
- The Lenstool v1 amplitude offset recurs (3 of 5 clusters tested) with one factor across families. A
  normalisation-at-z ≈ 2 hypothesis fits this cluster (post hoc).

## Deviation
- PLAN.md was edited after 56a0b3bb in one caveat sentence, not in any rule: a foreground-structure redshift written
  from memory, with no checked source, was replaced by what the SLICE text and the Hsu et al. 2013 title support.
