# Plan: `arcs` on SPT-CL J0615−5746 (written before any `arcs` count, 2026-10-18 18:25Z)

**Why.** One more independent field for the Abell S1063 background anti-excess replication tally (4 vs 0.52; P
0.032 clump-aware, #273, #279). Tally so far: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no
result, WHL0137 no, RXC J0600 no. SPT0615 has never been run through any screen here and was never used to tune or
choose an `arcs` rule. It has public JWST F200W imaging (program 4212), two public RELICS lens models and a public
RELICS HST photo-z catalogue. It has no DJA photo-z catalogue (DJA v7 index checked 2026-10-18), so the photo-z
source differs from the earlier fields; the result is labelled exploratory.

**Field survey for this choice (2026-10-18, recorded so it is not repeated).** DJA v7 photo-z exists for only two
lensing clusters not yet run (`cluster-j032728m1326`, `mrg0138`), neither with public maps found. The LAM Strong
Lensing Cluster Atlas (data.lam.fr/sl-cluster-atlas; ANIS API `server/search/sl-cluster-atlas/model`) lists 17
Lenstool models; of those clusters only Bullet (JWST 4598) has deep F200W, and VENUS (6882) visits are 1,245 s.
RELICS fields with JWST F200W: SPT0615 (4212, 2,920 s), MACS0308 (5293/6882, ≤ 1,578 s).

**Hypothesis tested (H0).** Background sources in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments and clump pairs, member/ICL photo-z blends, BPZ catastrophic outliers (HST-only photometry,
seven bands), diffraction spikes, detector artefacts, lens-model error (both models are HST-era).

**Inputs (pinned by sha256 in RESULT.md):**
- JWST pipeline catalogue `jw04212-o001_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC; the only public F200W
  level-3 catalogue within 2′).
- RELICS BPZ catalogue `hlsp_relics_hst_acs-wfc3ir_spt0615m57_multi_v1_cat.txt` (HLSP DOI 10.17909/T9SP45),
  converted to a zout-like table: `z_phot` = `zb`, `z160` = `zbmin`, `z840` = `zbmax`. **ASSUMPTION:** BPZ gives
  only the 95 % interval, so it stands in for the 68 % interval used elsewhere; this is conservative (fewer rows
  are called background). Rows with `stel` ≥ 0.8 are dropped (stars).
- Astrometry: the median RA/Dec offset of BPZ–JWST pairs within 1″ (S/N ≥ 10 JWST rows) is measured first; if it
  exceeds 0.1″ in either axis it is removed from the BPZ positions before the standard 0.3″ match.
- Primary model: RELICS GLAFIC v2 maps (`spt0615-relics-glafic`; the newest GLAFIC version in the HLSP). Second:
  RELICS Lenstool v1 maps (`spt0615-relics-lenstool`). D_LS/D_S = 1; z_lens 0.972 (Bleem et al. 2015; Paterno-Mahler
  et al. 2018, ApJ 863, 154); cosmology H0 70, Ω_m 0.3. Each must pass `validate` or it gets no run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2, aligned ≤
30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background alignment fraction.

**Outcome rule (frozen; `arcs_macs0647.outcome` through `arcs_rxj2129.decide`, unchanged and unit-tested):**
1. control: a run passes if ≥ 70 % of its strong-shear background rows are aligned (n = 0 fails);
2. on the first passing run (primary, then second): **replicates** if P(≥ n_anti) < 0.05 under the Poisson
   shape-noise expectation, else **does not replicate**;
3. no passing run → **no result**.

**Pre-registered secondary checks (do not change the outcome):** every anti row gets an F200W cutout looked at by
eye; the number of strong-shear background rows is reported (the high lens redshift limits it).

**What changes the plan:** "replicates" sends each anti row to the clump-pair census (#273 code) and
`/vet-candidate` before any reading. "Does not replicate" or "no result" adds one entry to the S1063 replication
tally and nothing else.
