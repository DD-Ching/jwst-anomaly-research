# RX J0152 family 4: revised forced photometry at the predicted third image (2026-10-25)

Frozen plan `PLAN.md` (pushed in e566d3e with the script and rule tests, before any pixel was read). Numbers:
`result.json` (frozen run) and `posthoc.json` (post hoc, `scripts/rxj0152_fam4_photometry_v2_posthoc.py`).
Imaging (observed): all ten VENUS GO 6882 NIRCam level-3 mosaics, 50″ cutouts read through MAST byte ranges (not
stored). Predicted positions and μ ratios: model_prediction (addendum be). Exploratory: addendum bf already looked
at this region.

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
- Gate C: the failed images are the four faintest in F277W (2.1 0.84, 2.3 0.85, 5.3 1.05, 2.2 1.25) plus 5.1 (1.74;
  also failed in bf). Best sibling offsets over all bands are 0.35–9.4 dex for 2.1, 2.2, 5.1, 5.3; 2.3's are within
  0.16 dex, so it presumably fails the ≥ 3 detectable-band condition (inferred, not stored).
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
