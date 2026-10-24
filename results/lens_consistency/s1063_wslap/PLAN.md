# Abell S1063 WSLAP+ JWST 2025: input check and spec-z gate K (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any value of the deflection maps was read. Next field of the
fitted-model route (CHANGELOG 2026-10-24 SPT-CL J0546 handoff: "a field with public JWST-era maps and spectroscopic
systems (Abell S1063 WSLAP+ ...)"). The S1063 image list has been gated off since the CATS v4.1 map failed
(D-035 addenda; "waits on another public S1063 model"). Script `scripts/s1063_wslap.py`; tests
`tests/test_s1063_wslap.py`.

## Inputs
Zenodo 10.5281/zenodo.17726321 ("WSLAP+ lens model for JWST 2025 data", published 2025-11-26, CC-BY-SA-4.0), one zip
`LensModel_AS1063_WSLAPplus_JWST_Case3_2025.zip` (1,016,075,135 B, md5 6a6558557078c1dc2b0c8a2a13d17944). Already seen
when this section was written (listing and headers only):
- Record description: deflection, convergence, shear and potential "derived with WSLAP+ based on BUFFALO arc
  identifications and new redefinitions for some of them from the latest JWST data"; maps "produced for
  Dds(zl,sz)/Ds(zs)=1"; magnification maps for several z. No image list in the record.
- Zip central directory (HTTP range reads): 14 FITS members; `AlphaX_DdsDs1_arcsec.fits` and
  `AlphaY_DdsDs1_arcsec.fits` are 141,122,880 B each (≈ 131 MB deflated). Magnification maps at z 0.73, 1.23, 1.261,
  2, 3, 6, 12 (and 1.261, 0.73 match Chan systems 2 and 10), so the model was fitted to at least some of the systems
  used below: **this is a reproduction check of a fit (method), not an independent test.**
- `AlphaX` header (first block, by range read of the member): 4200 × 4200, BITPIX −64, RA---TAN / DEC--TAN, FK5,
  CDELT ∓1.66667e-5° (0.06″ / px), CRPIX 2100, CRVAL (342.1827257, −44.5313911), no BUNIT, no rotation.

Data handling (cloud disk rule): only the two α members are fetched, each by one HTTP range request over its
compressed bytes, inflated while streaming, CRC-32-checked against the zip directory, written to the session
scratchpad, sha256 recorded, and deleted after the run. Transfer ≈ 262 MB (> 200 MB): stated reason, the deflection
maps are the only released form of the model and nothing smaller carries α; the other 12 members (≈ 750 MB) are not
fetched. Nothing but `summary.json` and the residual table is stored.

Image list: Chan et al. 2020 (ApJ 888, 35), `results/lens_consistency/abells1063_images/chan20_images.ecsv`
(observed positions): the **36 images of the 12 systems with a MUSE z_spec** (systems 1, 2, 5, 6, 7, 10, 11, 12, 13,
14, 19, 52; none flagged `dup_position`), each at its z_spec.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the WSLAP+ maps, loaded by our map-model code, reproduce the spectroscopic multiple images? A
  pass makes S1063 usable for the extra-image stages (gate C2, predictions) with a JWST-era model; a fail stops the
  route on this field. Nothing exotic is read: no extra-image, flux or flag reading in this stage.
- Ordinary explanations of a fail: our loading (unit, normalisation, orientation, WCS), map resolution, the
  mean-source estimator near caustics (the CATS v4.1 failure, D-035 addendum b), different image identifications
  (WSLAP+ used BUFFALO / redefined identifications, not Chan's list), or the model (WSLAP+ is free-form and its
  published rms values are typically ≳ 1″).

## Frozen rules (computable; each is a function with a unit test written before the run)
Cosmology flat ΛCDM H0 70, Ωm 0.3 (ASSUMPTION; only the ratio D_ds / D_s enters). z_lens 0.348
(`lensmodel` HFF table). Normalisation: map value = deflection at D_ds / D_s = 1 (record), so `MapLensModel` takes
the maps as they are (× the reading's unit).
1. **C (coverage)** `coverage`: all 36 images inside the map. Fails → stop.
2. **U (unit)** `unit_check`, readings A (values in arcsec, as the file names state) and P (values in pixels,
   × 0.06″): per system at its z_spec, the rms of `backtrace_images` `dtheta_arcsec`; a reading is accepted when the
   median of the 12 per-system values is ≤ 1.0″ (ASSUMPTION). U passes when exactly one reading is accepted; else stop.
3. **K (spec-z gate)** `m0417_images.gate_k` on the 36 images at z_spec with the accepted reading, exact image-plane
   solutions (`lensmodel.imageplane_residuals`, mean back-traced source, grid step 0.25″): passes when rms ≤ 1.0″
   **and** ≥ 80 % of the images lie within 3.0″ of a predicted image.
4. **K-WLS (declared now, reported, not a gate):** the same with the WLS family source (`a2744_uncover.imageplane`),
   because the mean source failed near caustics on S1063 before. It changes no verdict in this stage; a future stage
   that wants the WLS source must cite this number and freeze its own rule.
5. **R (reference catalogue)** `dja_check`: the DJA v7 imaging index lists a field whose name holds `abells1063`.

Reported, not gates (fixed now): per-system rms, catalogued vs predicted image counts and `shared_match`; the U values
of both readings; the CATS v4.1 numbers of D-035 addendum b for comparison (rms 18.5″ with the mean source).

## Stop rules
- No change of any rule after the first map value is read; any departure is labelled "deviation" or "post hoc".
- No extra-image, flux or flag reading in this stage. A pass hands the field to a separately frozen extra-image plan.
