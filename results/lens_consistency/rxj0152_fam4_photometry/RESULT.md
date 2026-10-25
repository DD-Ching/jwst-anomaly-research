# RX J0152 family 4: forced photometry at the predicted third image (2026-10-25)

Frozen plan `PLAN.md` (pushed in ad41f8a before any pixel was read). Numbers: `result.json` (frozen run) and
`posthoc.json` (post hoc). Imaging (observed): VENUS GO 6882 NIRCam level-3 mosaics F115W, F150W, F277W, F444W, 50″
cutouts read through MAST byte ranges (not stored). Predicted positions and μ: model_prediction (addendum be).

**Deviation:** the script and rule tests were committed after the run (a lint failure stopped the pre-run commit);
changes after the run: a shortened docstring and a trailing newline on the JSON writes.

## Frozen result
- Pair 4.1 / 4.2 detected as two F277W segments (203, 206); F277W window for a counterpart [4.51, 60.96] (pixel sums).
- **Gate V fails: 1 / 3** injected copies recovered (GLAFIC position yes; Lenstool, Zitrin no).
- Gate C passes: 11 / 12 catalogued images of families 2, 3, 5, 6 colour-match a sibling (5.1 does not).
- Gate D passes (F277W σ_ap 0.335; faintest accepted flux 4.51 = 13.4 σ_ap).
- Three segments in the region pass rule 5 (119, 127, 131); five region segments in all.
- **Outcome: O3_undetermined** (gate V failed). No flag, no reading of the three segments.

## Post hoc (not in the plan)
- **Why V failed** (`posthoc.json` `gate_V_diagnosis`): the injected stamp copies both pair images (two blobs ≈ 1.9″
  apart, each ≈ half of the expected 13.5). At the Lenstool and Zitrin positions the 2σ segmentation keeps only the
  cores, two segments 0.77–1.23″ from the target with F277W sums 1.56–2.27, below the window floor 4.51. At the GLAFIC
  position the copy merged with existing segment 131. Two design faults: the injection should be one image, not the
  pair, and isophotal segment sums lose most of a faint source's flux, which would also push a real third image
  below the window (biasing toward O2).
- **Chance control** (`chance_control`): 62 non-pair segments in the 2,501 arcsec² cutout pass the flux window and
  colour rule; for the 61.5 arcsec² region that predicts 1.53 by chance, observed 3, Poisson P(≥ 3) = 0.20. At the
  0.3 dex four-band colour tolerance the three segments are not distinguishable from field sources.
- Contact sheet looked at (not stored): 4.1 is a compact source; the three passing segments are compact sources near
  the GLAFIC / Lenstool predictions; no arc-like source at the predicted positions.

## Reading and next
Nothing is concluded about the third image. A revised frozen plan needs: a single-image injection (4.1 scaled by
μ3 / |μ4.1| per model), fixed-aperture forced fluxes instead of isophotal sums, all ten VENUS bands (tighter colour
match), and the chance control above as a frozen gate.
