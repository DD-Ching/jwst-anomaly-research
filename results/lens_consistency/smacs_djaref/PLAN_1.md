# SMACS 0723 reference-catalogue extra-image test: stage 1 (gate K) — frozen plan

Frozen 2026-10-23 before any residual of this run was read. Script: `scripts/smacs_djaref.py`.

**Question.** Does the Mahler+2022 ICLv2 Lenstool model (`smacs0723-iclv2`, D-024) reproduce its own
fixed-redshift image list well enough (gate K) to read extra-image predictions in a later stage with the
reference-catalogue route used on MACS0417 / MACS1149 / MACS0416 / Abell 370 / Abell 2744 / El Gordo
(D-034 addenda b–q)? D-034 (2026-10-08) ran SMACS with single-image references only (0 absent); D-024 back-traced
all 60 `arcs.dat` images at the model's own redshifts (rms 0.319″), which is not the gate-K estimator.

**Hypothesis and ordinary explanations.** A predicted extra image (model_prediction) that is absent from the data
would be a lens-model residual input. Ordinary explanations come first in later stages (model position error as
Abell 2744 78.3c and El Gordo 23c, reference contamination, local noise in BCG / ICL light, photo-z or z errors).

**Inputs.** Model `best.par`, `input.par` and `arcs.dat` pinned by sha256 in `lensmodel.SMACS0723_MAHLER22_ICLV2`;
60 images in 21 systems. No frame offset (the model was built on the ERO NIRCam frame; `lc.MODELS` carries none).

**Rule 0 (spectroscopic set, ASSUMPTION).** A system counts as `z_spec` when `input.par` has no `z_m_limit` line
for it (its redshift is fixed in the fit). `input.par` frees 16 systems (4, 6, 8–13, 15–19, 21, 25, 26), so the
fixed set is systems 1, 2, 3, 5, 7 (16 images). Function `fixed_z_systems`, unit-tested.

**Rule 1 (gate K, `m0417_images.gate_k`, ASSUMPTIONs as frozen there):** each family's source is the mean
back-trace of its catalogued images (`m0417_images.Solver.image_plane`, grid step 0.25″); pass when the
image-plane rms of the fixed-z images is ≤ 1.0″ and ≥ 80 % lie within 3″ of a predicted image. The other
systems get `m0417_images.fitted_z` (recorded, not gated). Own-constraint validation only.

**Rule 2 (stage-2 input check, no flux read):** the DJA v7.4 `smacs0723-grizli-v7.4-fix_phot.fits` reference is
reachable from the session (HTTP range, headers only) and its table header carries every column El Gordo stage 2b
used for the grizli total flux and the containing-source rule (D-034 addenda p–q): `x_image`, `y_image`, `a_image`,
`b_image`, `theta_image`, `flux_auto`, `mag_auto`, `tot_corr`, `flux_aper_1`, `f200w_flux_aper_1`,
`f200w_fluxerr_aper_1`, `f200w_mask_aper_1`, `f200w_tot_corr`, plus the `APER_1` header card. Function
`stage2_inputs`, unit-tested; it passes only when every item is present.

**Outcomes.** Gate K passes and rule 2 passes → stage 2 (frozen plan first): El Gordo stage 2b's rules on the
SMACS reference (`flux_auto × tot_corr` unit check, grid check, frame tables with `m0417_extra.run_field`,
gate C2 with the shared-pair rule, local noise). Gate K fails → the route stops for this model (no prediction
read). Rule 2 fails → stage 2 first freezes a different flux definition. Families with more predicted than
catalogued images are recorded as stage-2 inputs, not read.
