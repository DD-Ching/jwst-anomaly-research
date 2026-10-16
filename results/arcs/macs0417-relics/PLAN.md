# Plan: `arcs` on MACS J0417.5-1154 (written before the run, 2026-10-16)

**Why.** Abell S1063 has 4 background `anti` arcs against 0.52 expected (#273); under the clump-aware null with the
five-field clump rate it stays at P 0.031 (#279, exploratory). The S1063 excess was found post hoc in one of seven
fields. A field never used to tune or choose any `arcs` rule is the cheapest independent test of whether the excess
is a property of strong-lensing fields (a method or population effect) or of S1063 alone.

**Hypothesis tested (H0).** Background arcs in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
RELICS Lenstool v1 maps of MACS0417; the number of `anti` background sources (g_min < 1) is consistent with the
shape-noise expectation (`anti_expected`, σ_ε 0.26) and, in a second row, with the clump-aware null at the
five-field rate f̂ 0.274 (upper 0.364, #279).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak
shear, deblending fragments and clump pairs, member/ICL photo-z blends, diffraction spikes, lens-model error
(the RELICS v1 model is HST-era, Mahler et al. 2019).

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, aligned ≤ 30°, anti ≥ 60°,
z-margin 0.1, σ_ε 0.26. No image list is used (map model, as MACS1149 and S1063): the control is the background
alignment fraction.

**Outcome rule (frozen):**
- control: ≥ 70 % of strong-shear background rows aligned, else the model is not trusted here and no anti count
  is interpreted;
- **replicates** the S1063 pattern if P(≥ n_anti | g_min < 1) < 0.05 under the shape-noise null; then every anti
  row goes to cutouts and the clump-aware null;
- **does not replicate** if P ≥ 0.05; the S1063 excess then stays a single-field, exploratory result.

**What changes the plan:** a control failure stops the run (no field result). Replication sends the anti rows to
`/vet-candidate`-style ordinary tests (bands, ICL subtraction, second model: RELICS GLAFIC) before any reading.

**Inputs:** JWST program 1208 (CANUCS) MACS0417 cluster-pointing F200W pipeline catalogue (MAST level 3), CANUCS
DR1 `macs0417-clu` photometry with EAzY photo-z, RELICS Lenstool v1 maps (HLSP DOI 10.17909/T9SP45). Byte ranges
or small files only; nothing over 200 MB.
