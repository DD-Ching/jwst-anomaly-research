# Plan: `arcs` on MACS J0647.7+7015 (written before any `arcs` count, 2026-10-16 17:14Z)

**Why.** The Abell S1063 background anti excess (4 vs 0.52; P 0.031 clump-aware, #273, #279) was found post hoc in
one field and did not replicate in MACS0417 (#280; CANUCS maps, anti 3 vs 1.52, P 0.195). MACS0647 is a second
field never used to tune or choose any `arcs` rule, with JWST F200W imaging, a DJA v7 photo-z catalogue and two
public lens models. It is the cheapest further independent test of whether the S1063 excess is a property of
strong-lensing fields (a method or population effect) or of S1063 alone.

**Hypothesis tested (H0).** Background arcs in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak
shear, deblending fragments and clump pairs, member/ICL photo-z blends, diffraction spikes, detector artefacts,
lens-model error (both models are HST-era CLASH models, Zitrin et al. 2009, 2013; rms 1.4″ for LTM-Gauss).

**Inputs (pinned by sha256 in RESULT.md):**
- JWST program 1433 (PI Coe) F200W pipeline catalogues, MAST level 3: epoch 1 `jw01433-o010_t001` (primary) and
  epoch 2 `jw01433-o020_t001` (repeat), 2,104 s each. The MAST query (101.9574, +70.249, r 0.05°) found no other
  public F200W level-3 pointing.
- DJA v7.0 `macs0647-grizli-v7.0-fix.eazypy.zout.fits` (7,956 rows, 147 z_spec > 0), extracted from the photo-z
  tarball while streaming (the tarball itself is not stored).
- CLASH Zitrin v2 maps (HLSP, `missions/hlsp/clash/macs0647/models/zitrin/`), deflection in arcsec at
  D_LS/D_S = 1, 0.065″ pixels, z_lens 0.591: **primary `macs0647-clash-nfw`** (PIEMD + eNFW, parametric like the
  Lenstool CATS model used for S1063); **second `macs0647-clash-ltm`** (LTM-Gauss). Convention check (done before
  this plan, scratch): κ, γ1 and γ2 from the deflection gradients reproduce the published maps (median relative
  difference ≤ 6 × 10⁻⁴ for both models); the LTM-Gauss x-deflection carries a constant offset of about −47″, which
  does not enter κ or γ. `validate` repeats the κ and z = 2 μ check as part of the run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2,
aligned ≤ 30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background
alignment fraction.

**Outcome rule (frozen; applied by `outcome()` in `scripts/arcs_macs0647.py`, unit-tested before the run):**
1. control: the primary run passes if ≥ 70 % of its strong-shear background rows are aligned; otherwise the
   field gives no result from the primary model and the second model's run is reported with the same rule, labelled
   a deviation-free fallback (it is pre-registered here);
2. on the first run whose control passes: **replicates** if P(≥ n_anti) < 0.05 under the Poisson shape-noise
   expectation, else **does not replicate**;
3. both controls fail → **no result**.

**Pre-registered secondary checks (do not change the outcome):** the anti rows of the deciding run are matched to
epoch 2 (≤ 0.3″) and their epoch-2 class is reported; their class under the other model is reported; every anti
row gets an F200W cutout looked at by eye.

**What changes the plan:** replication sends the anti rows to the S1063 vetting chain (bands, ICL subtraction,
clump-pair census) before any reading; non-replication leaves the S1063 excess a single-field exploratory result
with two failed replications.

**Data:** the maps are 70.6 MB each (uncompressed FITS); 14 files are cached in the gitignored data root and
deleted after the run. Nothing over 200 MB per file.
