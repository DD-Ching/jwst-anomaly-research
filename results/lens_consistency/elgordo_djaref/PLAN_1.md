# El Gordo reference-catalogue extra-image test: stage 1 (gate K) — frozen plan

Frozen 2026-10-23 before any residual of this run was read. Script: `scripts/elgordo_djaref.py`.

**Question.** Does the Caminha+2023 Lenstool model (`elgordo-caminha23`, D-030; frame offset D-034 rule 1)
reproduce its own spectroscopic image list well enough (gate K) to read extra-image predictions in a later stage
with the reference-catalogue route used on MACS0417 / MACS1149 / MACS0416 / Abell 370 / Abell 2744
(D-034 addenda b–n)? D-034 (2026-10-08) ran El Gordo with single-image references only (0 absent).

**Hypothesis and ordinary explanations.** A predicted extra image (model_prediction) that is absent from the data
would be a lens-model residual input. Ordinary explanations come first in later stages (model position error,
frame offset, reference contamination, local noise in BCG / ICL light, photo-z or z errors).

**Inputs.** Model `best.par` and `arcs.dat` pinned by sha256 in `lensmodel.ELGORDO_CAMINHA23`; 56 images in 23
systems, all with spectroscopic z (docs/fields/elgordo.md). Frame offset (+0.224″, −0.016″) applied as in D-034.

**Rule 1 (gate K, `m0417_images.gate_k`, ASSUMPTIONs as frozen there):** each family's source is the mean
back-trace of its catalogued images (`m0417_images.Solver.image_plane`, grid step 0.25″); pass when the
image-plane rms of the 56 images is ≤ 1.0″ and ≥ 80 % lie within 3″ of a predicted image.
This is own-constraint validation only (the model was fitted to these positions).

**Rule 2 (stage-2 input check, no flux read):** the DJA v7.0 `elgordo-grizli-v7.0-fix_phot.fits` reference is
reachable from the session; its header lists the size, row count and the segment shape columns the containing-
source rule needs (`a_image`, `b_image`, `theta_image`, `x_image`, `y_image`, a total-flux column per band).
Only the FITS headers are read (HTTP range).

**Outcomes.** Gate K passes → stage 2 (frozen plan first): predicted / control tables with
`m0417_extra.run_field`, gate C2 with the shared-pair rule (D-034 addendum i) and the local-noise rule on the DJA
reference. Gate K fails → the route stops for this model (no prediction read), as for Abell 370 CANUCS.
Families with more predicted than catalogued images are recorded as stage-2 inputs, not read.
