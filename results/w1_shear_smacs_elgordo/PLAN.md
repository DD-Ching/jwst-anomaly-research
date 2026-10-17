# W1 shear screen on SMACS 0723 and El Gordo (plan, frozen before the run)

Cloud routine, 2026-10-17 (claim PR #303). D-053 "Revisit if": the DJA photo-z tarballs for these fields
answered 404 on 2026-10-08; on 2026-10-17 both answer HTTP 200 and match the sha256 in SOURCES.md
(SMACS v7.4 `13b5b6f4…`, El Gordo v7.0 `9a4a39ce…`; zout members `82db8910…`, `efd051e4…`). Tarballs deleted
after extracting the zout member.

**Hypothesis (H).** A negative-mass (W1) lens of |M| ≥ 2 × 10¹² M☉ in either field would give a radial
(E-mode positive in the D-053 sign convention) aperture-mass peak.
**Ordinary explanations of a peak.** Diffraction spikes and streaks (radial lines), PSF anisotropy and
catalogue isophotal-moment noise (B-mode floor), blends near bright members, residual cluster shear from model
error (R g removal), edge apertures.
**Rule (D-053 defaults, unchanged; ASSUMPTIONs).** Schirmer Q_TANH, R_ap 10″, 1″ grid, sources 1–10″,
rotation null; a peak counts only if p_rot < 0.05 **and** S > the field's max |S_×|. Grid half-width =
the FIELDS `max_radius` (SMACS 60″, El Gordo 110″). Spike veto: catalogue spike segments plus Gaia DR3
G < 17 within 4′ of the model centre (SMACS (110.826989, −73.454723), El Gordo (15.740693, −49.271992)),
as for Abell 370 and Abell S1063.
**Injections.** `inject_shear.py` seed 50, 200 per mass at the five D-053 masses, same settings as the five
earlier fields, then a seven-field combination.
**What changes the plan.** A peak passing the rule goes to `/vet-candidate` (streak check on the image first;
S1063 showed streaks mimic the sign). R outside (0, 1.5] or < 20 calibrating rows is an error and the field
is dropped (D-053). If both fields are null, the seven-field limits replace the five-field ones.
**Known-case check.** R must fall near 0.41–0.48 (the cluster shear is seen); the El Gordo model frame offset
(0.224″, −0.016″) is applied by the model loader.
