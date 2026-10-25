# RX J0152 family 4: forced photometry at the predicted third image (frozen plan, 2026-10-25)

Written and pushed before any pixel was read. Script `scripts/rxj0152_fam4_photometry.py`; rule functions tested
in `tests/test_rxj0152_fam4_photometry.py`.

## Question and inputs
- Input (`../slice_relics_rxj0152/posthoc.json`, model_prediction): three HST-era RELICS models (Lenstool v1, GLAFIC
  v2, Zitrin LTM-Gauss v1) all predict a positive-parity image of SLICE family 4 (pair 4.1 / 4.2, observed positions
  from SLICE Table 8) that SLICE does not list, at (28.17895, −13.95293), (28.17805, −13.95268), (28.17945, −13.95343).
- Predicted flux ratio R = (|μ4.1| + |μ4.2|) / μ3 per model: Lenstool 14.05, GLAFIC 10.05, Zitrin 15.10 (derived
  from the image-plane tables of addendum be).
- **Hypothesis H_miss:** the third image is absent at the predicted flux: a lens-model residual input.
- **Ordinary explanations, in order:** (a) the image is there and SLICE did not list it; (b) it is fainter than
  predicted because the pair's μ near the critical curve is model-sensitive (|μ| 11–27 across the models); (c) below
  the imaging depth; (d) blended with a cluster member; (e) 4.1 / 4.2 are not a multiple-image pair.

## Data (observed)
VENUS GO 6882 level-3 NIRCam mosaics `jw06882-o008_t008_nircam_clear-<band>_i2d.fits`, bands F115W, F150W, F277W,
F444W, read as one 50″ cutout per band (MAST download byte ranges, `jwst_anomaly.cutouts.make_cutouts`; not stored)
centred on (28.18100, −13.95480), which contains the pair and all three predictions.

## Measurement (rules)
1. Background: subtract the 3σ-clipped median of each cutout; noise σ_ap = 3σ-clipped std of the sums in 300 random
   0.3″-radius apertures whose centres are off detected segments.
2. Detection: photutils `detect_sources` on the F277W cutout convolved with a 2-pixel-FWHM Gaussian, threshold
   2 × the clipped pixel std, ≥ 10 connected pixels, then `deblend_sources` (defaults). Segment fluxes in all bands
   are sums over the same F277W segment (the four mosaics are resampled to their own grids, so each band's segment
   map is the F277W map reprojected by WCS onto that band's pixel grid, nearest neighbour).
3. Pair flux F_pair(band): sum of the distinct segments containing 4.1 and 4.2 (within 0.3″ of each position;
   none found → the run stops: "pair not detected").
4. Search region: union of 3.0″ circles (ASSUMPTION: the gate-K match radius) around the three predicted positions.
5. Counterpart rule: a segment whose flux-weighted centroid lies in the region, with F277W flux in
   [F_pair / (3 · 15.10), 3 · F_pair / 10.05] (ASSUMPTION: factor 3 around the predicted range) and
   every colour |log10(F_b / F_F277W) − log10(F_pair,b / F_pair,F277W)| ≤ 0.3 (ASSUMPTION) for b in F115W,
   F150W, F444W. Segments with any band flux ≤ 0 fail the colour test.
6. Depth gate D: the faintest accepted flux F_pair / (3 · 15.10) in F277W must be ≥ 5 σ_ap (F277W).
7. Validation gate V (synthetic injection, before reading the result): a copy of the F277W/all-band pair stamp
   (segment pixels of 4.1 / 4.2, background-subtracted) scaled by 1 / 15.10 is injected at each of the three
   predicted positions in turn (one at a time); V passes when rules 2–5 return a counterpart for 3 / 3.
8. Known-case check C (colour rule on catalogued families): for families 2, 3, 5, 6 (three catalogued images each),
   the colour test of rule 5 applied between each image's segment and the family's other images must accept at
   least 6 of the 12 images (an image is accepted if it matches ≥ 1 sibling). C fails → no reading of rule 5.

## Outcomes
- **O1 counterpart:** V and C pass, ≥ 1 segment passes rule 5 → explanation (a); no flag. List them.
- **O2 missing:** V, C and D pass and no segment passes rule 5 → flag `missing_predicted_image` for
  `/vet-candidate`, with (b), (d) and (e) to test there. Not evidence of anything exotic.
- **O3 undetermined:** V, C or D fails.

Post hoc work (contact sheet, other bands, other apertures) is labelled post hoc.
