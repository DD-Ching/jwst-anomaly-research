# RX J0152 family 4: revised forced photometry at the predicted third image (2026-10-25)

Frozen plan `PLAN.md` (pushed in e566d3e with the script and rule tests, before any pixel was read). Numbers:
`result.json` (frozen run) and `posthoc.json` (post hoc, `scripts/rxj0152_fam4_photometry_v2_posthoc.py`).
Imaging (observed): all ten VENUS GO 6882 NIRCam level-3 mosaics, 50″ cutouts read through MAST byte ranges (not
stored). Predicted positions and μ ratios: model_prediction (addendum be). Exploratory: addendum bf already looked
at this region.

**Deviation:** after the run, review fixes changed the code without changing the result (re-run on fresh cutouts gives
a byte-identical `result.json`): 4.1's position is passed to `inject` explicitly, `measure` skips unused work,
`chance_p` returns NaN instead of 0 when no field rate is measured (a 0-field run would otherwise read as
significant; here n_field = 90), and a default temporary cutout directory is deleted after the run.

## Frozen result
- 4.1 isolated (its own F277W segment). F41(F277W) = 33.6 (0.3″ aperture sum); counterpart window [1.50, 18.96].
- **Gate V fails:** 27 / 30 injected single images recovered (≥ 24 needed), but only 2 / 3 at the predicted
  positions: the copy at the Lenstool prediction (and one 1.5″ Lenstool offset, one Zitrin offset) is missed.
- Gate C passes: 7 / 12 catalogued images of families 2, 3, 5, 6 colour-match a sibling (2.1, 2.2, 2.3, 5.1, 5.3 do not).
- **Gate D fails:** the lower window edge 1.50 is below 5 σ_ap(F277W) = 1.67.
- Gate X: 90 field segments pass rule 6; λ = 2.27 expected in the 61.5 arcsec² region; observed 3 (segments 112,
  127, 131), Poisson P(≥ 3) = 0.40.
- **Outcome: O3_undetermined** (gate V failed). No flag. Had V passed, the frozen table would give O1c (P 0.40 > 0.05):
  the three matches are not distinguished from field sources at this tolerance.

## Post hoc (not in the plan)
- **V misses** (`posthoc.json`): in all three missed injections the copy is detected 0.10–0.12″ from the target with
  F277W 4.25–4.31 (expected 4.50–4.70, inside the window) but fails the colour test; the failing band was not stored.
  The fixed 0.3 dex tolerance over up to nine bands rejects some faint true copies.
- **Gate V credit check** (raised by code review): the GLAFIC at-prediction recovery and one GLAFIC offset were
  credited only to segment 131, which passes rule 6 without any injection (`gate_V_recovered_credit`). Counting only
  new counterparts, V is 1 / 3 at the predictions and 25 / 30 overall; the outcome is unchanged.
- Gate C (`gate_C_colour_offsets`): 2.1, 2.3 and 5.3 (F277W 0.84, 0.85, 1.05) fail because fewer than 3 bands are
  detectable against either sibling; 2.2 fails at 0.53 dex (3 bands vs 2.3), 5.1 at 0.86–0.89 dex (also failed in bf).
- No 0.3″ aperture used here contains a non-finite pixel in any band (`aperture_nonfinite_pixels`).
- Region matches vs models: 127 and 131 (F277W 9.86, 10.59) are 2.11″ and 0.97″ from the GLAFIC prediction and
  within a factor 3 of its expected flux; 112 (1.57) is 2.85″ from the Zitrin prediction.
- Contact sheet looked at (F150W / F277W / F444W, not stored): 4.1 / 4.2 form a tangential arc; 127 and 131 are
  compact sources, 112 is faint; no arc-like source at the predicted positions. A μ 2.4–3.8 image of the arc's
  source would be compact, so morphology does not separate the matches from field galaxies.

## Reading and next
Nothing is concluded about the third image. With λ = 2.27 colour-and-flux-matched field sources expected in the
region, aperture colours alone cannot decide between (a) and H_miss even with a working injection gate. A deciding
test needs an independent discriminant: photometric or spectroscopic redshifts of 127 / 131 / 112 against family 4
(SLICE z 1.69, photometric), e.g. from a public JWST photometric catalogue of the field. Without one, this line
stays a lensing-residual input and is paused.
