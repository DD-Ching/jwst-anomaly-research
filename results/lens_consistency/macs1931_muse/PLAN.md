# MACS J1931 stage 4: MUSE Lyα line-flux reference for systems 5, 6, 7, 100 and 101 (frozen before any run)

**Why.** Stage 3d (D-034 addendum ag) left 10 of 11 predicted extra images `no_reference`: systems 5, 6, 7, 100 and
101 (spec-z 4.745–5.339) have no F200W continuum. They are line emitters, so a line-flux reference can predict how
bright each extra image should be. Only the cube footprint and the target list were looked at before freezing this
plan (no spectrum or narrow-band image of any position).

**Hypothesis under test.** The VENUS AstroLensPy v2.1 model (K-WLS source, as in stage 3d) predicts images that the
MUSE Lyα data do not show. **Ordinary explanations of an `absent` reading** (to vet before any interpretation):
model error near critical curves (μ of order 10–30 is uncertain by factors), differential magnification of extended
Lyα halos, a predicted position off by more than the search radius, MUSE PSF blending, cube edge noise.
**What changes the plan:** gate L fails → the MUSE route is not validated for this field and no target is read; any
`absent` target → contact sheet, then `/vet-candidate` in a later cycle; all targets `recovered` / `confused` /
untestable → the MACS J1931 extra-image line is closed for this model.

## Inputs
- observed: ESO MUSE combined cube `ADP.2017-03-27T15:26:07.707` (095.A-0525(A), 8,127 s, WFM, 0.2″ spaxels,
  footprint RA 292.94744–292.96750, Dec −26.58529 to −26.56779), read only as SODA `BAND` cutouts (full field,
  ±80 Å around each system's nominal Lyα), ≈ 105 MB each, streamed to the session scratchpad and deleted after use;
  URL, size and sha256 of every cutout go to `summary.json`.
- observed: image list `../macs1931_inputs/m1931_images.ecsv` (positions, spec-z).
- model_prediction: |μ| of observed images from `../macs1931_gatek/imageplane_spec.ecsv` (`magnification`), and
  of targets from `../macs1931_stage3d/predicted_macs1931.ecsv` (`abs_mu`, with its `search_arcsec`).

## Frozen rules (functions in `scripts/m1931_muse.py`, unit-tested before the run)
1. **Cube band:** λ_nom = 1215.67 Å (1 + z). Continuum slices: |λ − λ_nom| > 30 Å inside the cutout.
2. **Line window (per system):** S(λ) = sum over the system's in-cube observed images of the 0.8″-radius aperture
   spectrum minus each spaxel's continuum median. The peak is the max of S within |λ − λ_nom| ≤ 20 Å; the window is
   the peak slice ± 4 slices (11.25 Å). ASSUMPTIONs: 0.8″, 20 Å, ± 4.
3. **Narrow-band image:** NB = 1.25 Å × Σ_window (data − per-spaxel continuum median), 10⁻²⁰ erg s⁻¹ cm⁻².
4. **Aperture flux and local noise:** circular 0.8″ (4-spaxel) aperture sum, NaN if it touches a NaN or the edge.
   σ_loc = σ-clipped std of 200 random apertures 2–6″ away, redrawn within 1.0″ of any listed image (all systems),
   seed 20261024. ASSUMPTIONs.
5. **Blending:** a position within 1.6″ (two aperture radii) of an observed image of **another** system whose window
   overlaps its own is `blended`, untestable (systems 100 and 101 share z = 5.133).
6. **Reference (per system):** the in-cube, unblended observed image with the highest S/N = F / σ_loc; it needs
   S/N ≥ 5, else the system is `no_line_reference`.
7. **Prediction and measurement:** F_pred = F_ref × |μ_t| / |μ_ref|. F_meas = max aperture sum over spaxel
   centres within the search radius (targets: their stage 3d `search_arcsec`; controls: 1.0″, ASSUMPTION).
   Classes, in order: `off_cube` (centre outside the cube or NaN F_meas), `blended`, `undetectable_local`
   (F_pred < 5 σ_loc), `recovered` (1/3 ≤ F_meas / F_pred ≤ 3), `confused` (> 3), `absent` (< 1/3 and
   F_pred ≥ 10 σ_loc), else `ambiguous` (the existing FLUX_RANGE and 5σ / 10σ constants).
8. **Gate L (controls):** every other in-cube observed image of the five systems, predicted from its own system's
   reference. Passes iff ≥ 4 controls are testable (`recovered` / `confused` / `absent` / `ambiguous`) **and** ≥ 75 %
   of them are `recovered` or `confused`. ASSUMPTIONs. Fails → no target reading.
9. **Flag:** a target with class `absent` while gate L passes. Unvetted until `/vet-candidate`.

## Budget and stop
One cycle. Network ≈ 4 cutouts × ≈ 105 MB (systems 100 and 101 share one band), each < 200 MB. If the run does not
finish, the plan stays frozen and the next cycle runs it unchanged.
