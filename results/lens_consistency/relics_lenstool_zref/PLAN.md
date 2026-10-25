# Plan: is the RELICS Lenstool v1 amplitude offset a z_ref normalisation? (frozen 2026-10-25, before any header was read)

## Question
D-034 addendum bh (post hoc) noted that the Abell 3192 Lenstool v1 median amplitude factor 0.692 equals
D_LS/D_S at z_ref 2.09 (z_lens 0.425; flat, H0 70, Ωm 0.3), i.e. what a deflection map normalised at a source
z ≈ 2 instead of D_LS/D_S = 1 would give. This plan tests that hypothesis on the factors already stored and
on the published map headers. No imaging and no new image lists.

## Inputs (all stored; derived)
Rule 7 median amplitude factors of the RELICS Lenstool v1 maps (z_lens as configured in
`scripts/lens_consistency.py`):

| Cluster | z_lens | Median factor | Source |
|---|---|---|---|
| RXC J0232 | 0.2836 | from `slice_relics_heldout/posthoc.json` `scale` | addendum bc |
| RXC J2211 | 0.397 | same (map fails frame check F; factor still defined) | addendum bc |
| MACS J0553 | 0.412 | same | addendum bc |
| RX J0152 | 0.8269 | 1.045 (`slice_relics_rxj0152` RESULT, rule 7) | addendum be |
| Abell 3192 | 0.425 | 0.692 (`slice_relics_macs0358/posthoc.json`) | addendum bh |

SMACS 0723 (addendum aw) has no per-family factor (post hoc |α| ratio ≈ 1.25 to GLAFIC only); it is listed as
context, not tested.

## Hypotheses and frozen rules
- **H_U (one RELICS-wide z_ref):** every RELICS Lenstool v1 deflection map is normalised at one z_ref ≈ 2.09.
  Prediction per cluster: f_pred = D_LS/D_S(z_lens, 2.09). **Rejected** if any of the five clusters has
  |f_obs − f_pred| > 0.10 (ASSUMPTION: the per-family spreads seen are ≤ ±0.05 around the median).
- **H_C (a per-map normalisation stated in the product):** read the primary and first extension headers of
  each map's `kappa` and `x-arcsec-deflect` FITS (byte-range, header only). **Supported** if a header keyword
  gives a source z or D_LS/D_S whose implied factor matches f_obs within 0.05 on the offset clusters
  (MACS J0553, Abell 3192) and ≈ 1 on the others. **Not supported** if no such keyword exists or it is the
  same for offset and non-offset maps.
- **Ordinary explanation kept open (H_M):** the HST-era models fitted families at assumed (photometric)
  redshifts; a model's overall mass scale is degenerate with those assumed redshifts, so a global factor can
  arise from the model's own inputs. Not testable here (the RELICS Lenstool constraint lists are not
  published with the maps); it is what remains if H_U and H_C both fail.

## What would change the plan
- H_U survives: apply one correction factor to every RELICS Lenstool v1 map before any reading (rule change).
- H_C supported: read that keyword in the map loader for all RELICS maps.
- Both fail: keep the addendum-bd per-cluster amplitude check as the only rule; no correction.
