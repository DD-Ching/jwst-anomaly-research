# Plan: `arcs` on Abell S1063 (written before the run, 2026-10-16)

**Hypothesis tested.** Background arcs in strong shear (g_min ≥ 0.2) are aligned with the stretching direction
predicted by the CATS v4.1 maps (`abells1063-cats`); the number of `anti` (≥ 60° off at every tested z) background
sources is consistent with the shape-noise expectation (`anti_expected`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments, member/ICL photo-z blends (13 of 46 DJA-matched images get z_phot < 0.6 here), diffraction
spikes, model error (the CATS image list fails its gate, D-035, because `arcs.txt` redshifts are placeholders).

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″ (extended to 100″ only as
a labelled sensitivity run), aligned ≤ 30°, anti ≥ 60°, z-margin 0.1. Known-case check first: elongated catalogued
images must be aligned (median offset ≲ 15°) or the model is not used.

**What changes the plan:** a known-case failure stops the run (model not trusted here). An anti excess with
P(≥ n) < 0.05 not explained by g_min ≥ 1 sources sends each anti source to cutouts and `/vet-candidate`.

**Inputs:** GLIMPSE 3293 o001 F200W pipeline catalogue (manifest `data/manifests/abells1063.ecsv`), DJA v7.5
eazy zout (sha256 `03bc4c28…`, docs/fields/abells1063.md).
