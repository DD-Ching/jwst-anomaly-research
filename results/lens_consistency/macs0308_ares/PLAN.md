# MACS J0308.9+2645: SN Ares image order and time delays in the pre-JWST RELICS maps (plan, frozen before the run)

**Input (observed, literature).** SN Ares (Larison et al. 2026, arXiv:2610.33161v1, read 2026-11-08): lensed SN IIb at
z_spec 1.276 ± 0.002, RA 03:08:53.418, Dec +26:45:58.63 (J2000), in image 1.3 of a triply imaged host. Three JWST-era
models (Lenstool, GLAFIC, a third by Allingham et al. in prep.) say Ares is the **first** image to arrive; the other
two arrive after 72–79 yr (paper text; none of these models is public).

**Question.** Do the two public pre-JWST RELICS models of the cluster (GLAFIC v2, Zitrin LTM-Gauss v1; Acebron et al.
2018; already pinned in `lensmodel.MACS0308_RELICS_*`) give the same image order? A model in which another image
arrives **earlier** than Ares predicts that the SN was already visible there in archival imaging (a testable
two-epoch prediction); a model in which Ares leads agrees with the published JWST-era models.

**Method** (reuses `scripts/smacs_at2024jfz.py`): trace the SN position to the source plane at z_s 1.276, re-image
it with `lensmodel.find_images`, Fermat potential φ = |θ − β|²/2 − ψ with ψ differences as line integrals of the
deflection (straight path, checked against an L-shaped path), delay = (φ_k − φ_Ares) × (1+z_l) D_l D_s / (c D_ls).

**Frozen rule** (`decide`, per model, unit-tested before the run):
- gate P: path check ≤ 1 % (as in the AT 2024jfz plan), else `stopped`;
- gate S: an image within 0.3″ of the SN position (ASSUMPTION), else `stopped`;
- `conflicts`: some image arrives ≥ 1 yr (ASSUMPTION) **before** the Ares image;
- `agrees`: ≥ 3 images and every other image arrives ≥ 10 yr (ASSUMPTION) **after** the Ares image;
- otherwise `undetermined` (fewer than 3 images, or a delay between −1 and +10 yr).

**Ordinary explanations of a `conflicts`.** The RELICS models predate JWST and were constrained by HST only
(Acebron et al. 2018); a different image order can follow from a wrong model redshift of the host family, a
missing member halo or the map edge. A `conflicts` is a model disagreement, not an anomaly, until archival imaging
at the predicted earlier epoch is checked.

**Not done here.** No published JWST-era model is public, so no direct comparison of maps. Delays are best-model
values (no MCMC). H0 70, Ωm 0.3 (the RELICS map cosmology).
