# Plan: `arcs` on RXC J0600−2007 (MACS J0600.1−2008) (written before any `arcs` count, 2026-10-18 17:30Z)

**Why.** The Abell S1063 background anti excess (4 vs 0.52; P 0.032 clump-aware, #273, #279) is a single-field
exploratory result. Independent fields so far: MACS0417 and WHL0137 do not replicate, MACS0647 and RXJ2129 give no
result (controls fail), MACS1423 is explained by clump orientation. RXC J0600−2007 has never been run through any
screen in this repository and was never used to tune or choose an `arcs` rule. It has public JWST F200W imaging, a
DJA photo-z catalogue and two public RELICS lens models (map availability checked 2026-10-18).

**Hypothesis tested (H0).** Background sources in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments and clump pairs, member/ICL photo-z blends, z_phot at the eazy grid edge (> 15), diffraction
spikes, detector artefacts, lens-model error (both RELICS models are HST-era; Fujimoto et al. 2024,
arXiv:2404.03286, later found extended mass structure in this cluster, so the control may fail).

**Inputs (pinned by sha256 in RESULT.md):**
- JWST pipeline catalogue `jw06796-o010_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC). Choice rule, fixed here:
  of the three public F200W level-3 catalogues within 2′ (programs 4573, 6796, 6882), the one with the longest F200W
  exposure whose footprint centre lies within 30″ of the cluster centre: 6796 o010 (2,920 s, centre ≈ 13″ from the
  BCG). The other two are not used for the outcome.
- DJA v7.2 `rxcj0600-grizli-v7.2-fix.eazypy.zout.fits` streamed out of
  `https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rxcj0600-grizli-v7.2-fix.photoz.tar.gz` (33.8 MB; tarball not
  stored).
- Primary model: RELICS GLAFIC v3 maps (`rxc0600-relics-glafic`; the newest version in the HLSP). Second model:
  RELICS Lenstool v1 maps (`rxc0600-relics-lenstool`). Both at D_LS/D_S = 1; z_lens 0.43 (Fujimoto et al. 2021,
  arXiv:2101.01937); cosmology H0 70, Ω_m 0.3, as for the other RELICS fields. Each must pass `validate` (κ from the
  deflection maps against the published κ map) or it gets no run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2, aligned ≤
30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background alignment fraction.

**Outcome rule (frozen; `arcs_macs0647.outcome` through `arcs_rxj2129.decide`, unchanged and unit-tested):**
1. control: a run passes if ≥ 70 % of its strong-shear background rows are aligned (n = 0 fails);
2. on the first passing run (primary, then second): **replicates** if P(≥ n_anti) < 0.05 under the Poisson
   shape-noise expectation, else **does not replicate**;
3. no passing run → **no result**.

**Pre-registered secondary checks (do not change the outcome):** z_phot > 15 rows among the anti are listed; every
anti row gets an F200W cutout looked at by eye.

**What changes the plan:** "replicates" sends each anti row to the clump-pair census (#273 code) and
`/vet-candidate` before any reading. "Does not replicate" or "no result" adds one entry to the S1063 replication
tally and nothing else.
