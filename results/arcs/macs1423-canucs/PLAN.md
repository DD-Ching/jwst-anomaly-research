# Plan: `arcs` on MACS J1423.8+2404 (written before any `arcs` count; committed 2026-10-16 18:13Z, the run started 18:13:10–18:14Z)

**Why.** The Abell S1063 background anti excess (4 vs 0.52; P 0.031 clump-aware, #273, #279) was found post hoc in
one field. MACS0417 (#280; CANUCS v1 maps) did not replicate it (anti 3 vs 1.52, P 0.195); MACS0647 (#281) gave no
result (both HST-era CLASH models failed the control). The MACS0647 handoff said MACS1423 lacked a public model; in
fact CANUCS DR1 carries a JWST-era Lenstool model for it (found 2026-10-16 in the HLSP directory). MACS1423 has never
been used to tune or choose any `arcs` rule.

**Hypothesis tested (H0).** Background arcs in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak
shear, deblending fragments and clump pairs, member/ICL photo-z blends, diffraction spikes, detector artefacts,
lens-model error (the model is preliminary and constrained by only 5 spectroscopic systems, readme).

**Inputs (sha256 in `lensmodel.py` and RESULT.md):**
- JWST program 1208 (CANUCS) F200W pipeline catalogue `jw01208-o010_t005_nircam_clear-f200w_cat.ecsv` (6,399 s;
  MAST level 3). The MAST query (215.9495, +24.0786, ±0.06°) found no other public F200W level-3 pointing.
- CANUCS DR1 `macs1423-clu` photometry, converted with `inject_radial.canucs_zout` (as for MACS0417).
- Model `macs1423-canucs-maps`: CANUCS DR1 Lenstool maps (Desprez, preliminary v1), best fit, 100 mas,
  D_LS/D_S = 1, z_lens 0.545, flat ΛCDM H0 70, Ωm 0.3. `validate` must reproduce the published κ from the
  deflections (median relative difference ≤ 10⁻³), else no `arcs` run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2,
aligned ≤ 30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background
alignment fraction.

**Outcome rule (frozen; `outcome()` in `scripts/arcs_macs1423.py`, unit-tested before the run):** one model, no
fallback.
1. control: passes if ≥ 70 % of the strong-shear background rows are aligned (n = 0 fails); otherwise **no result**;
2. control passes: **replicates** if P(≥ n_anti) < 0.05 under the Poisson shape-noise expectation, else
   **does not replicate**.

**Pre-registered secondary checks (do not change the outcome):** background counts by g_min < 1 and ≥ 1; every
anti row gets an F200W cutout looked at by eye.

**What changes the plan:** replication sends the anti rows to the S1063 vetting chain (bands, ICL subtraction,
clump-pair census) before any reading; non-replication adds a second independent non-replication; no result leaves
the tally unchanged.

**Data:** three maps (87–89 MB each, gzip) and the photometry (25 MB) are cached in the gitignored data root and
deleted after the run. Nothing over 200 MB per file.
