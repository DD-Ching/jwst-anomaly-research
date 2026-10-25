# RX J0152 family 4, revised forced photometry at the predicted third image (frozen plan, 2026-10-25)

Written and pushed before any pixel of this run was read. Revises addendum bf
(`../rxj0152_fam4_photometry/RESULT.md` "Reading and next"), whose gate V failed for two design faults: the
injection copied the whole pair, and isophotal segment sums lost a faint source's flux. Script
`scripts/rxj0152_fam4_photometry_v2.py`; rule functions tested in `tests/test_rxj0152_fam4_photometry_v2.py`.

**Not held-out:** addendum bf's run looked at the same region (three compact segments near the GLAFIC / Lenstool
predictions, contact sheet viewed). This run is a re-measurement with a fixed design, labelled exploratory.

## Question and inputs
- Input (`../slice_relics_rxj0152/posthoc.json`, model_prediction): Lenstool v1, GLAFIC v2 and Zitrin LTM-Gauss v1
  each predict a positive-parity image of SLICE family 4 that SLICE does not list, at (28.17895, −13.95293),
  (28.17805, −13.95268), (28.17945, −13.95343).
- Predicted single-image ratio R_m = |μ4.1| / μ3 (derived, same file): Lenstool 7.165, GLAFIC 5.325, Zitrin 7.482.
- **Hypothesis H_miss:** the third image is absent at the predicted flux (a lens-model residual input).
- **Ordinary explanations, in order:** (a) the image is there and SLICE did not list it; (b) it is fainter than
  predicted because μ4.1 near the critical curve is model-sensitive (|μ4.1| 12.7–27.3 across the models); (c) below
  the imaging depth; (d) blended with a cluster member; (e) 4.1 / 4.2 are not a multiple-image pair.

## Data (observed)
All ten VENUS GO 6882 level-3 NIRCam mosaics `jw06882-o008_t008_nircam_clear-<band>_i2d.fits` (F090W, F115W,
F150W, F200W, F210M, F277W, F300M, F356W, F410M, F444W; list from a MAST CAOM query, 2026-10-25), one 50″ cutout per
band centred on (28.18100, −13.95480) through MAST byte ranges (`jwst_anomaly.cutouts.make_cutouts`; not stored).

## Measurement (rules)
1. Background: subtract the 3σ-clipped median of each cutout.
2. Detection (unchanged from bf): photutils `detect_sources` on the F277W cutout convolved with a 2-pixel-FWHM
   Gaussian, threshold 2 × the clipped pixel std, ≥ 10 connected pixels, `deblend_sources` defaults. Detection
   only gives candidate positions: the flux-weighted F277W centroid of each segment.
3. **Forced fluxes:** a 0.3″-radius circular aperture (photutils `aperture_photometry`, exact overlap) at each
   position in every band. Aperture noise σ_ap(b): 3σ-clipped std of 300 random aperture sums per band whose
   centres fall off detected segments (seed 20261025).
4. Reference: the forced fluxes F41(b) at SLICE's catalogued 4.1 position. The run stops ("4.1 not isolated") if
   4.1's segment also contains 4.2's position.
5. Search region: union of 3.0″ circles (ASSUMPTION: the gate-K match radius, as in bf) around the three predictions.
   Candidates: segments with centroids in the region, excluding segments within 0.3″ of 4.1 or 4.2.
6. **Counterpart rule:** F277W forced flux in [min_m F41/R_m / 3, 3 · max_m F41/R_m] (ASSUMPTION: factor 3) and the
   colour test: for every band b ≠ F277W where the reference SED scaled to the candidate's F277W flux,
   F41(b) · F(F277W) / F41(F277W), is ≥ 3 σ_ap(b) (ASSUMPTION: "detectable band"), require
   |log10(F(b)/F(F277W)) − log10(F41(b)/F41(F277W))| ≤ 0.3 dex (ASSUMPTION, as in bf); a non-positive flux in such a
   band fails. Fewer than 3 detectable bands fails the test (not informative).
7. **Gate V (single-image injection):** the stamp is 4.1's own segment pixels (background-subtracted; in each band
   the F277W segment mapped by WCS, nearest neighbour), scaled by 1 / R_m and added at (i) each model's predicted
   position (3 injections) and (ii) 9 positions 1.5″ from each prediction at angles 0°, 40°, …, 320° (27 injections),
   one at a time, using that model's R_m. Recovered: rules 2–6 return a counterpart within 1.0″ of the injection
   point. V passes when 3 / 3 of (i) and ≥ 24 / 30 of all (ASSUMPTION: 80 %) are recovered.
8. **Gate C (known cases):** forced fluxes at the catalogued positions of families 2, 3, 5, 6; an image is accepted
   when it passes rule 6's colour test against ≥ 1 sibling (sibling as reference). C passes with ≥ 6 / 12.
9. **Gate D (depth):** the lower window edge min_m F41/R_m / 3 is ≥ 5 σ_ap(F277W).
10. **Chance control X (frozen gate):** n_field = segments outside the region and away from 4.1 / 4.2 (0.3″) that
    pass rule 6; λ = n_field · A_region / (A_valid − A_region), areas from finite F277W pixels;
    P = Poisson P(N ≥ n_obs | λ), n_obs = candidates passing rule 6.

## Outcomes
- **O3 undetermined:** V or C fails.
- **O1 counterpart:** V and C pass, n_obs ≥ 1 and P ≤ 0.05 (ASSUMPTION) → explanation (a); no flag. List them with
  the nearest model and whether the F277W flux lies within a factor 3 of that model's F41/R_m.
- **O1c chance-level:** V and C pass, n_obs ≥ 1, P > 0.05 → a counterpart is not distinguished from field sources;
  undetermined between (a) and H_miss; no flag.
- **O2 missing:** V, C and D pass and n_obs = 0 → flag `missing_predicted_image` for `/vet-candidate`, with (b), (d)
  and (e) to test there. Not evidence of anything exotic.
- V and C pass, n_obs = 0, D fails → O3.

Post hoc work (contact sheet, other apertures, other radii) is labelled post hoc.
