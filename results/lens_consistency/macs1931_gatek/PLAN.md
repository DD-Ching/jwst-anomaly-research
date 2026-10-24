# MACS J1931: unit check U and gate K of the VENUS AstroLensPy v2.1 maps (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any value of the deflection maps was read. Follows
`../macs1931_inputs/RESULT.md` (D-034 addendum aa). Script `scripts/macs1931_gatek.py`; tests
`tests/test_macs1931_gatek.py`.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the published best-fit deflection maps, loaded by our map-model code, reproduce the
  spectroscopic multiple images the model was fitted to? A pass makes the maps usable for later model-residual
  stages on this field (no deblended reference catalogue is public, so gate C2 and the extra-image test cannot run
  yet); a fail stops the fitted-model route on this field.
- Ordinary explanations of a fail: our loading (deflection unit, α normalisation redshift, orientation, WCS, frame),
  map resolution, or the model. A catalogued image that is not reproduced is a residual of a fit to these very
  images, so it is **not** a lensing-residual input and no flag is read here.
- What would change the plan: U picks one reading and K passes → the next cycle looks for an image-level test that
  needs no reference catalogue (e.g. the VENUS mosaics on MAST). U picks no reading or K fails → the route stops on
  this field; nothing exotic is read either way.

## Inputs (Zenodo 10.5281/zenodo.18646525, md5 checked against the record, sha256 recorded)
- model_prediction: `Alpha_x.fits.zip`, `Alpha_y.fits.zip` (best model `vA_52`; 3001 × 3001, 0.04″ / px; 67.7 /
  67.0 MB zipped, under the 200 MB rule). Cloud session: kept in the session scratchpad only, deleted after the run.
- model output, used as a check: `sourcepos.txt` (the authors' source positions per system, pixels; the README
  says to use the last two columns, `Bary_Pos_x`, `Bary_Pos_y`); `chires.txt` (their image-plane rms, reported).
- observed: `../macs1931_inputs/m1931_images.ecsv` (51 used images; 26 in the 9 spec-z systems).
- **α normalisation (decided before the maps, from the authors' text):** the record states that every value that
  needs a source redshift is computed for z_s = 5.079, and `Readme_Magnifications_summary.txt` lists
  `k_dlds` 0.89043 / 0.99885 / 1.03257 at z_s 2 / 5 / 9. With flat ΛCDM H0 70, Ωm 0.3 and z_lens 0.3525
  (`Z_CL`), (D_LS/D_S)(z) / (D_LS/D_S)(5.079) gives exactly these three values (derived; recomputed and stored by
  the script as `normalisation_check`). So D_LS/D_S = 1 deflection = map / (D_LS/D_S)(5.079). This check is a
  function (`k_dlds_check`, tolerance 5e-5 per value, ASSUMPTION) and must pass before the maps are read.
- ASSUMPTION: the maps' pixel axes follow the FITS WCS (+i West, +j North), and `sourcepos.txt` pixels are 1-based
  map pixels (a ±1 px ambiguity is 0.04″, far below the U tolerance). Frame offset 0 (the list and the maps share
  the authors' frame; FK5 vs ICRS is far below the tolerances).

## Frozen rules (computable)
1. **Unit check U** (`unit_check`), run separately for reading A (map values in arcsec, README) and reading P (map
   values in pixels, header `BUNIT`; × 0.04″): back-trace each of the 26 spec-z images at its system z, average
   per system, and take the distance of that mean to the authors' `Bary` source of the system (converted to the
   model frame through the map WCS). A reading is **accepted** when the median of the 9 per-system distances is
   ≤ 0.5″ (ASSUMPTION). U **passes** when exactly one reading is accepted; then only that reading is used below.
   If U fails, stop: gate K does not run.
2. **Gate K** (`m0417_images.gate_k`, same code as MACS0417 / MACS1423 / MACS1149 / Bullet), combined over the 26
   spec-z images at their system z with the accepted reading: per-system mean source, image plane solved on a
   0.25″ grid, each image's distance to its nearest predicted image; passes when rms ≤ 1.0″ **and** ≥ 80 % of the
   images lie within 3.0″. The function decides; both conditions must hold.
3. Reported, not gates (fixed now): (a) gate K's numbers with the D-041 WLS family source
   (`a2744_uncover.imageplane`, the Abell 2744 / Bullet precedent); (b) gate K's numbers over the 44 used images
   with a non-zero catalogue z (photo-z / model-optimised systems at the file z; the 7 z = 0 images skipped); (c) the
   per-system rms and predicted image counts above or below the catalogued count; (d) the authors' `chires.txt` rms
   (Bary 0.444″) next to ours; (e) the per-system U distances for both readings.

## Stop rules
- No change of any rule after the first map value is read; any departure is labelled "deviation" or "post hoc".
- No extra-image, flux or flag reading in this stage.
