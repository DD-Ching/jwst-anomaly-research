# SLICE JWST-era images vs RELICS Lenstool v1 maps: RXC J0232, RXC J2211, MACS J0553 (2026-10-25)

Frozen plan `PLAN.md` (pushed in c1e2644 before any map was read); maps pinned, the Table 8 parser, the image list
and the tests (`tests/test_slice_relics_heldout.py`) were pushed in 47c7efd before the run. Rule functions are those
of `macs0647_clash_heldout` / `macs2129_clash_heldout`. Numbers come from `summary.json`, `posthoc.json` and
`validate/*/validate.json`; per-image residuals are in `imageplane_<cluster>.ecsv` (provenance `model_prediction`;
positions observed, SLICE Table 8, clump 1 of each source). No imaging was read.

## Frozen result
| | RXC J0232 | RXC J2211 | MACS J0553 |
|---|---|---|---|
| Images (families) | 26 (7) | 33 (10) | 57 (21) |
| Gate K0, families not held out (fitted z) | passes: rms 0.69″, 26 / 26 within 3″ | passes: rms 0.96″, 11 / 11 | passes: rms 1.32″, 30 / 30 |
| Frame check F (median RA·cos δ / Dec) | passes: −0.02″ / −0.07″ | **fails**: +0.38″ / **−0.52″** | passes: +0.00″ / −0.31″ |
| Usable | yes | no | yes |
| Reading A (held-out families, fitted z) | none held out | no reading (not usable) | families 11–21, 27 images: rms 1.14″, 26 / 27 within 3″; **1 flag: 11.1 (3.49″)** |
| Reading B (spectroscopic families at z_spec) | families 5 and 7 (8 images): rms 0.37″, 8 / 8; no flag | none (family 1 z ambiguous) | families 1, 5, 10, 13 (11 images): rms 18.4″, 4 / 11; **7 flags**: 1.1, 1.3, 5.1, 5.2, 5.3, 10.1, 10.2 |

Frozen labels: 11.1 `lensing_residual_input`; the seven reading-B rows `redshift_residual_input`. These are inputs,
not candidates.

## Post hoc (not in the plan)
- **MACS J0553 reading B fits one global amplitude factor.** The fitted z of every MACS J0553 family is 0.71–1.18
  (`summary.json` `z_fit`), each below that family's SLICE z (SLICE z 1.14–7.25). For the four spectroscopic families, D_LS/D_S at the
  fitted z over D_LS/D_S at z_spec is 0.664, 0.663, 0.677, 0.637 (`posthoc.json`; z_spec 1.14–3.76), i.e. one
  factor ≈ 0.66 for all redshifts. On RXC J0232 the same ratio is 1.005 and 0.982. A single factor across a
  redshift range of 1.14–3.76 is what a map-amplitude or normalisation difference predicts, not a family-specific
  lensing residual, so the seven reading-B flags are read as that ordinary effect. Its cause is not isolated: the
  published κ and μ(z = 2) maps agree with our deflection-derived values (median ratio 1.000 / 1.006,
  `validate/macs0553-relics-lenstool/validate.json`), so the released products are internally consistent and the
  factor sits in the released model itself (its z_lens, cosmology or z reference) or in our D_LS/D_S = 1 reading of it.
- **11.1 is not affected by that factor** (rule 1 fits z per family, which absorbs a global amplitude). It misses by
  3.49″ while 11.2 / 11.3 match at 1.26″ / 0.77″ (`imageplane_macs0553.ecsv`). The threshold is 3″, the HST-era
  model never had family 11, and the SLICE model rms is 0.66″ (paper text); an HST-era model missing one image of
  a new family by 3.5″ is within the ordinary expectation in PLAN.md. Not vetted further here.
- RXC J2211 fails F by 0.016″ on Dec. Reported only: its held-out families give rms 2.21″, 18 / 22 within 3″,
  with 4.2, 5.1, 5.2, 6.2 beyond 3″ (`imageplane_rxc2211.ecsv`); family 5 fits z 9.79. Family 1 fits z 1.04, which
  is closer to the Table 8 value 1.051 than to the text's 1.501 (not used: frozen out of reading B).
- Table 8 clump-1 rows give 33 RXC J2211 images where the text says 32 (not resolved).

## Reading
- Two of three RELICS models are usable references on the SLICE lists; the one frozen held-out flag (MACS J0553
  11.1, 3.49″) is a lensing-residual input for the SLICE Lenstool model, not a candidate.
- **Limit:** each cluster has one HST-era model, so "every usable model misses" is one model. The JWST-era SLICE
  Lenstool models (LAM atlas, closed on 2026-10-25) remain the reference to re-run 11.1 and the MACS J0553 scale
  question against.
