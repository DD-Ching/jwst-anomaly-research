# Plan: `arcs` on MACS J0308.9+2645 (written before any `arcs` count, 2026-10-18 19:20Z)

**Why.** One more independent field for the Abell S1063 background anti-excess replication tally (4 vs 0.52; P
0.032 clump-aware, #273, #279). Tally so far: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no
result, WHL0137 no, RXC J0600 no, SPT0615 no result. MACS0308 has never been run through any screen here and was
never used to tune or choose an `arcs` rule. It has public JWST F200W imaging, two public RELICS map models and a
public RELICS HST BPZ catalogue, but no DJA photo-z catalogue, so the result is labelled exploratory (as SPT0615).

**Hypothesis tested (H0).** Background sources in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments and clump pairs, member/ICL photo-z blends, BPZ catastrophic outliers (HST-only photometry),
diffraction spikes, detector artefacts, lens-model error (both models are HST-era, light-traces-mass for Zitrin).

**Inputs (pinned by sha256 in RESULT.md / `lensmodel.py`):**
- JWST pipeline catalogue `jw05293-o005_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC, 1,578 s). Two public
  F200W level-3 catalogues lie within 2′ (5293 o005; 6882 o015, 1,245 s). Chosen by a coverage count only (no
  `arcs` quantity): S/N ≥ 10 rows within 50″ of the GLAFIC map centre are 705 (5293) vs 640 (6882).
- RELICS BPZ catalogue `hlsp_relics_hst_acs-wfc3ir_macs0308p26_multi_v1_cat.txt` (HLSP DOI 10.17909/T9SP45),
  converted exactly as for SPT0615 (`arcs_spt0615.bpz_zout`): `z_phot` = `zb`, `z160` = `zbmin`, `z840` = `zbmax`,
  `stel` ≥ 0.8 dropped. **ASSUMPTION:** the 95 % interval stands in for the 68 % one (conservative).
- Astrometry: the BPZ–JWST median offset (pairs within 1″, S/N ≥ 10) is removed if it exceeds 0.1″ in either axis
  (`arcs_spt0615.shift`), then the standard 0.3″ match.
- Primary model: RELICS GLAFIC v2 maps (`macs0308-relics-glafic`; the only GLAFIC version in the HLSP). Second:
  RELICS Zitrin LTM-Gauss v1 maps (`macs0308-relics-zitrin`). D_LS/D_S = 1; z_lens 0.356 (Acebron et al. 2018,
  ApJ 858, 42, Table 1); cosmology H0 70, Ω_m 0.3. Each must pass `validate` or it gets no run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2, aligned ≤
30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background alignment fraction.

**Outcome rule (frozen; `arcs_rxj2129.decide`, unchanged and unit-tested):**
1. control: a run passes if ≥ 70 % of its strong-shear background rows are aligned (n = 0 fails);
2. on the first passing run (primary, then second): **replicates** if P(≥ n_anti) < 0.05 under the Poisson
   shape-noise expectation, else **does not replicate**;
3. no passing run → **no result**.

**Pre-registered secondary checks (do not change the outcome):** every strong-shear background row of each run gets
F200W / F444W cutouts looked at by eye (program 5293 o005 i2d, S3 byte ranges).

**What changes the plan:** "replicates" sends each anti row to the clump-pair census (#273 code) and
`/vet-candidate` before any reading. "Does not replicate" or "no result" adds one entry to the S1063 replication
tally and nothing else.
