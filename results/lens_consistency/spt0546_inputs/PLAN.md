# SPT-CL J0546-5345 LTM input check, unit check U and free-z gate K (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any value of the deflection maps was read. Next field of the
fitted-model route (CHANGELOG 2026-10-24 RXCJ2211 handoff: "e.g. SPT-CL J0546 LTM maps, input check first"). Script
`scripts/spt0546_inputs.py`; tests `tests/test_spt0546_inputs.py`.

## Inputs (Zenodo 10.5281/zenodo.16534337, Allingham+2025 LTM model, CC-BY-4.0; md5 checked against the record)
Already seen when this section was written (input listing only):
- The record description: the maps of every quantity that needs a source redshift are computed for z_s = 3.5 (the
  paper's anchor redshift); `alpha_x`, `alpha_y` (and 1σ error maps) are released; 16 systems "including the most likely
  candidates" in `actcl0546_mulfile_v1.mul`; a spec-z catalogue `actcl0546_spec-z.cat` (conf 3 / 4 good / very good,
  6 stars).
- The `alpha_x` FITS header (first 11.5 kB by HTTP range): 2600 × 2600, RA---TAN / DEC--TAN, FK5, |CDELT| =
  1.111e-5° = 0.04″ / px, CRPIX 1073, no BUNIT and no source-redshift keyword.
- The mul file: 59 rows in Lenstool format (z column 0 for every row, i.e. no redshift); ids `N.k` for 10 systems
  (1–10, 45 rows) and `c`-prefixed candidates (14 rows: c2.3, c4.11, c4.12, c5.6, c6.21–23, c7.6–8, systems c11–c15).
- A position match of the mul rows to the spec-z catalogue (≤ 2″, not stored, re-done by the script as check Z): the
  nearest entries are cluster members (z 1.064–1.078) or a z 0 entry; **no lensed image has a spectroscopic
  redshift**. So every system's z is free in what follows; there is no spec-z gate K on this field.

Data handling: the two α maps (27 MB each, under the 200 MB rule) are read into the session scratchpad, md5-checked,
sha256 recorded, and deleted after the run. Nothing but the image list and `summary.json` is stored.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the published LTM deflection maps, loaded by our map-model code, reproduce the model's own
  multiple images? A pass makes the maps usable for later stages on this field (a reference-catalogue test, extra-image
  predictions); a fail stops the route on this field. Nothing exotic is read either way: a catalogued image that is not
  reproduced is a residual of a fit to these very images, not a lensing-residual input, and no flag is read here.
- Ordinary explanations of a fail: our loading (unit, normalisation, orientation, WCS), map resolution or extent, the
  free-z degeneracy, or the model.
- Free z weakens the gates: a per-system z absorbs any overall deflection scale within the D_ds / D_s range of the z
  grid, so U only distinguishes readings that differ by more than that range (A vs P differ by 25×; the grid spans
  D_ds / D_s ≈ 0.23–1.29 of its z_s 3.5 value, derived), and K is weaker than the spec-z gate K of other fields (as addendum w
  showed for MACS0647, free z can hide offsets).

## Frozen rules (computable; each is a function with a unit test written before the run)
Cosmology flat ΛCDM H0 70, Ωm 0.3 (ASSUMPTION: LTM's usual choice; under free z it only relabels z). z_lens 1.07
(record). Normalisation: map value at z_s 3.5, so deflection at D_ds / D_s = 1 is map × unit / (D_ds / D_s)(3.5).
Pixel axes follow the FITS WCS. "Main images" = the 45 non-`c` rows of systems 1–10.
1. **C (coverage)** `coverage`: every main image lies inside the map. Fails → stop.
2. **Z (spec-z)** `specz_systems`: a system is spectroscopic when one of its main images lies ≤ 0.5″ from a conf 3 / 4
   catalogue entry with z > z_lens + 0.1. Reported (expected 0 systems from the listing above); it decides only that no
   spec-z gate runs.
3. **U (unit)** `unit_check`, for reading A (map values in arcsec) and P (map values in pixels, × 0.04″): for each of the
   10 systems, back-trace its main images on the z grid (60 log-spaced z from z_lens + 0.2 to 10, as
   `m0417_images.fitted_z`) and take the smallest rms of the source-plane offsets mapped to the image plane
   (`backtrace_images` `dtheta_arcsec`). A reading is **accepted** when the median of the 10 per-system minima is
   ≤ 0.5″ (ASSUMPTION). U **passes** when exactly one reading is accepted; otherwise stop.
4. **K (free-z gate)** `m0417_images.gate_k` on the 45 main images, each system at its fitted z (`m0417_images.best_z`
   over exact image-plane solutions on the same z grid, grid step 0.25″), accepted reading only: passes when the rms
   ≤ 1.0″ **and** ≥ 80 % of the images lie within 3.0″ of a predicted image. Both must hold; the function decides.
5. **R (reference catalogue)** `dja_check`: the DJA v7 imaging index has a field whose name holds `spt0546`, `spt-0546`,
   `j0546`, `sptcl0546` or `actcl0546`.
Reported, not gates (fixed now): per-system fitted z (and whether it sits on a grid edge), rms, catalogued and predicted
image counts; the U minima of both readings; the `c` candidates' offsets from the nearest image predicted for their
system's source at its fitted z (systems 2, 4–7 only; c11–c15 have no main images and are not used).

## Stop rules
- No change of any rule after the first map value is read; any departure is labelled "deviation" or "post hoc".
- No extra-image, flux or flag reading in this stage.
