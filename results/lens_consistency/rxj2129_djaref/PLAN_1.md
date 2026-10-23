# RXJ2129 reference-catalogue extra-image test: stage 1 (gate K) — frozen plan

Frozen 2026-10-23 (20:20Z) before any residual of this run was read. Script: `scripts/rxj2129_djaref.py`.

**Why this field.** After SMACS 0723 (D-034 addendum s) the handoff named WHL0137. Its inputs add no new test: the
only JWST-era image list (Scofield+2025 Table 3) covers systems 1–2, which D-033 addendum 4 already tested with all
four RELICS models (0 flags), and the one JWST-era model (Scofield M5-G1) has an undocumented deflection unit
(D-033 addendum). RXJ2129 has a JWST-era Lenstool model with public maps (`rxj2129-lenstool-maps`, Dhanasingham+2026,
Zenodo 21460864; pinned, D-033-style map model used for `arcs` on 2026-10-16) **and** the image list it was fitted
to, on the same JWST/Gaia DR3 frame (paper Table "Systems_RA_and_dec", Appendix B; read from the arXiv source of
2604.11882). The Zenodo zip holds only FITS maps (583 members, no image list; listed 2026-10-23).

**Question.** Does the Lenstool model reproduce its own image list well enough (gate K) to read its extra-image
predictions against a reference catalogue in a later stage (D-034 addenda b–s)?

**Hypothesis and ordinary explanations.** A predicted extra image (model_prediction) absent from the data would be
a lens-model residual input. Ordinary explanations come first in later stages (model position error near critical
curves, reference contamination, local noise in BCG / ICL light, z errors).

**Inputs.**
- Images (observed): the paper table transcribed to `dhanasingham26_images.ecsv`: 31 images in 10 families
  (1–7, G, S3-1, S3-2). Ids are rewritten `<family>.<n>` (paper id kept in `paper_id`) so
  `lensmodel.image_family` reads the family. Every family has a spectroscopic z in the table (rows marked "--"
  take their family's z; S3-1 / S3-2 take system 3's z, as the paper states). The SN 2022riv row (a single image)
  is not a family and is left out. The paper's text says "11 systems"; the table has 10 families (recorded, not
  resolved).
- Model: `rxj2129-lenstool-maps` (sha256 pinned in `lensmodel.RXJ2129_LENSTOOL_MAPS`; arcsec at D_LS/D_S = 1,
  z_lens 0.234). No frame offset (same frame by construction). GLAFIC (`rxj2129-glafic-maps`) is not used: its
  maps are normalised at z_s = 3 and `MapLensModel` assumes D_LS/D_S = 1.

**Rule 1 (gate K, `m0417_images.gate_k`, ASSUMPTIONs as frozen there):** each family's source is the mean
back-trace of its images (`m0417_images.Solver.image_plane`, grid step 0.25″); pass when the image-plane rms of
all 31 images is ≤ 1.0″ and ≥ 80 % lie within 3″ of a predicted image. Own-constraint validation only (the
published Lenstool rms is 0.32″).

**Rule 2 (stage-2 input check, no flux read):** the DJA v7.0 `rxj2129-grizli-v7.0-fix_phot.fits` (the only v7
version reachable, 9,152,640 B; v7.1–v7.5 return 403) is read by HTTP range (headers only) and passes when its
table header carries every column SMACS rule 2 named (`smacs_djaref.STAGE2_COLS` and the `APER_1` card).

**Outcomes.** Gate K passes and rule 2 passes → stage 2 under a new frozen plan with El Gordo stage 2b's rules.
Gate K fails → the route stops for this model (no prediction read). Rule 2 fails → stage 2 first freezes another
flux definition (as El Gordo stage 2). Families with more predicted than catalogued images are recorded as
stage-2 inputs, not read.
