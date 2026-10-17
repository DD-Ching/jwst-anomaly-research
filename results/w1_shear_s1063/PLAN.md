# Plan: W1 shear screen (D-053) on Abell S1063 (fixed before any screen run; 2026-10-17)

**Why.** TASKS.md "JWST focus" 1: the shear screen has run on Abell 2744, MACS0416, MACS1149 and Abell 370 only;
Abell S1063 has a pinned lens model (CATS) and DJA v7.5 photo-z (`docs/fields/abells1063.md`), so it can be added.

**Hypothesis (W1, hypothesis only).** A negative-mass lens imprints radial (S > 0) alignment of background shapes
around a dark centre. **Ordinary explanations:** residual cluster shear (wrong R or model), diffraction spikes,
PSF anisotropy and blends (the B-mode floor), edge apertures.

**Method.** Unchanged D-053 settings (`exotic_screens.py shear`, all thresholds ASSUMPTIONs as in
docs/exotic_limits.md): F200W GLIMPSE pipeline catalogue, CATS model, DJA zout photo-z, spike veto with Gaia stars,
R fitted per field (error if outside (0, 1.5] or < 3σ), 200 rotations, E/B rule (p < 0.05 and S > max |S_×|).

**What changes the plan.** A peak passing the E/B rule goes to `/vet-candidate` and stamps before any reading. A
null adds S1063 area to the W1 limits only if injections are run with the same `inject_shear.py` settings (else the
field is recorded as screened without a limit).
