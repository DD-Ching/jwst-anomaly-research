# Plan: `arcs` on WHL0137−08 (Sunrise) (written before any `arcs` count, 2026-10-17 01:27Z)

**Why.** The Abell S1063 background anti excess (4 vs 0.52; P 0.032 clump-aware at the five-field clump rate, #273,
#279) is a single-field exploratory result. Independent fields so far: MACS0417 does not replicate, MACS0647 and
RXJ2129 give no result (controls fail), MACS1423 is explained by clump orientation. WHL0137 has never been run
through `arcs` and was never used to tune or choose any `arcs` rule. It has a JWST F200W pipeline catalogue
(program 2282 o010, pinned in `data/manifests/sunrise.ecsv`), a DJA v7.5 photo-z catalogue and four public RELICS
lens models, of which `whl0137-relics-lenstool` is already registered and used (D-033 addenda 3–4).

**Hypothesis tested (H0).** Background sources in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments and clump pairs, member/ICL photo-z blends, z_phot at the eazy grid edge (> 15), diffraction
spikes, detector artefacts, lens-model error (RELICS models are HST-era; WHL0137 is a lower-mass cluster with few
constraints, so the control may fail).

**Inputs (pinned by sha256 in RESULT.md):**
- `jw02282-o010_t001_nircam_clear-f200w_cat.ecsv` (manifest `sunrise`, sha256 `b0693065…`, jwst 2.0.1).
- DJA v7.5 `sunrise-grizli-v7.5-fix.eazypy.zout.fits` from
  `https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/sunrise-grizli-v7.5-fix.photoz.tar.gz` (51.2 MB, streamed; the
  tarball is not stored).
- Primary model `whl0137-relics-lenstool` (RELICS v1 maps, D_LS/D_S = 1, z_lens 0.566). Second model: RELICS v1
  GLAFIC κ/γ maps, used only as the pre-registered fallback and for the anti rows' class (same procedure as
  `scripts/arcs_second_models.py`), if it can be pinned in this run; otherwise reported as not run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2, aligned ≤
30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background alignment fraction.

**Outcome rule (frozen; `arcs_macs0647.outcome`, unchanged and unit-tested):**
1. control: a run passes if ≥ 70 % of its strong-shear background rows are aligned (n = 0 fails);
2. on the first passing run (primary, then second): **replicates** if P(≥ n_anti) < 0.05 under the Poisson
   shape-noise expectation, else **does not replicate**;
3. no passing run → **no result**.

**Pre-registered secondary checks (do not change the outcome):** z_phot > 15 rows among the anti are listed; every
anti row gets an F200W cutout looked at by eye; the F444W class of the anti rows if `arcs_bands.py` runs in time.

**What changes the plan:** "replicates" sends each anti row to the clump-pair census (#273 code) and
`/vet-candidate` before any reading. "Does not replicate" or "no result" adds one entry to the S1063 replication
tally and nothing else.
