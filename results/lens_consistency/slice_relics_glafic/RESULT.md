# SLICE JWST-era images vs RELICS GLAFIC v2 maps: MACS J0553, RXC J2211 (2026-10-25)

Frozen plan `PLAN.md` (pushed in 43fade5 before any GLAFIC map was read); maps pinned, script and rule tests
(`tests/test_slice_relics_glafic.py`) pushed in f7bdcde before the run. Image list, families and rules 1–6 are those
of `../slice_relics_heldout/` (D-034 addendum bc). Numbers come from `summary.json` and `validate/*/validate.json`;
per-image residuals are in `imageplane_<cluster>.ecsv` (provenance `model_prediction`; positions observed, SLICE
Table 8, clump 1). No imaging was read.

## Frozen result
| | MACS J0553 (GLAFIC v2) | RXC J2211 (GLAFIC v2) |
|---|---|---|
| Gate P (κ, μ(z = 2) median ratio) | passes: 1.000 / 1.007 | passes: 1.000 / 1.000 |
| Gate K0 (families not held out, fitted z) | passes: rms 0.88″, 30 / 30 within 3″ | passes: rms 0.40″, 11 / 11 |
| Frame check F (median RA·cos δ / Dec) | passes: −0.05″ / −0.32″ | passes: −0.17″ / −0.01″ |
| Reading A (held-out families, fitted z) | 27 images: rms 1.16″, 26 / 27; **1 flag: 12.2 (3.84″)** | 22 images: rms 2.77″, 16 / 22; **6 flags: 4.1, 4.4, 5.2, 5.3, 8.1, 8.2** |
| Reading B (spectroscopic families at z_spec) | 11 images: rms 1.47″, **11 / 11, no flag** | none (as addendum bc) |
| Rule 7 amplitude factor, median (range) over families | **0.966 (0.831–1.052), "near 1"** | 1.029 (0.996–1.245), "near 1" |
| Rule 9, MACS J0553 11.1 | GLAFIC 0.08″ → **"lenstool_only"** | — |

Frozen labels: 12.2 and the six RXC J2211 rows `lensing_residual_input` (single model each). Inputs, not candidates.

## Answers to the plan's questions
- **Q1:** 11.1 is matched by GLAFIC at 0.08″ (Lenstool v1 3.49″). It is a Lenstool-v1-only miss; closed.
- **Q2:** GLAFIC needs no amplitude factor on MACS J0553 (every family 0.83–1.05; Lenstool v1 0.637–0.708 for the
  same 21 families), and it reproduces all four spectroscopic families at z_spec. So the factor ≈ 0.66 belongs to the
  RELICS Lenstool v1 MACS J0553 product or model (H_L), not to the SLICE list or redshifts (H_S). The seven Lenstool
  reading-B flags of addendum bc are Lenstool-specific.

## Post hoc (not in the plan)
- No image is missed by both models on MACS J0553: 12.2 is 0.49″ in Lenstool v1; GLAFIC fits family 12 at the
  grid's upper edge (z 12.0; SLICE model z 7.25), and also predicts 2 → 3 images for families 12, 13, 14, 17, 19, 20
  (`extra_predicted`).
- RXC J2211: every GLAFIC flag shares its family with another image beyond or near 3″ (4.2 2.68″, 5.1 2.81″,
  8.3 2.64″), so they are family-wide offsets, not isolated misses; families 4 and 5 fit z 5.05 / 8.84 against SLICE
  model z 2.98 / 3.30. Lenstool v1 (frame check F failed) also misses 5.2 (5.00″) but matches 4.1, 4.4, 5.3, 8.1, 8.2
  within 1.8″. Not read further.
- The RELICS lens-model readme v3 states D_LS/D_S = 1 for every method's deflection maps, and both Lenstool v1's and
  GLAFIC's κ / μ(z = 2) products agree with their own deflections. A documented normalisation therefore does not
  explain the Lenstool v1 MACS J0553 factor; its cause (model normalisation through its assumed source z, or a
  product scaling) stays unisolated. With SMACS 0723 (Lenstool v1 ≈ 1.25 × GLAFIC, addendum aw) this is the second
  RELICS Lenstool v1 map set whose amplitude disagrees with GLAFIC in the same direction; RXC J0232 and RXC J2211
  Lenstool v1 need no factor.

## Reading
- With a second HST-era model the MACS J0553 Lenstool flags (11.1 and the seven z_spec misses) do not survive:
  GLAFIC reproduces them. No image of either cluster is missed by two usable models; no flag goes to vetting.
- **Practical rule (derived):** a RELICS Lenstool v1 map gives no reading on its own; check its amplitude against a
  second model (or z_spec families) first.
- **Limit:** both models are HST-era; the JWST-era SLICE Lenstool models (LAM atlas, still redirected to
  "restricted" when checked at about 09:10Z; PLAN.md wrote 09:20Z) remain the reference for RXC J2211 families 4, 5, 8 and MACS J0553 12.
