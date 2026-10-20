# Plan: `arcs` on RXC J0018.5+1626 (MACS J0018.5+1626, Cl 0016+16) (written before any `arcs` count, 2026-10-20 14:15Z)

**Why.** One more independent field for the Abell S1063 background anti-excess replication tally (4 vs 0.52; P
0.032 clump-aware, #273, #279). Tally so far: MACS0417, WHL0137, RXC J0600 and MACS0308 do not replicate; MACS0647,
RXJ2129 and SPT0615 give no result; MACS1423 is explained by clump orientation. RXC J0018 has never been run through
any screen here and was never used to tune or choose an `arcs` rule. It has public JWST F200W imaging (VENUS 6882),
a DJA v7.5 photo-z catalogue (`rj0018p1626`) and one public RELICS lens model (Lenstool v1; the HLSP `models/`
directory holds no other, listing 2026-10-20). The 2026-10-18 field survey (`results/arcs/spt0615-relics/PLAN.md`)
missed this field because its DJA name is `rj0018p1626`.

**Hypothesis tested (H0).** Background sources in strong shear (g_min ≥ 0.2) align with the stretch predicted by the
model; the number of `anti` background sources is consistent with the shape-noise expectation
(`anti_expected_shape_noise_background`, σ_ε 0.26).

**Ordinary explanations for any `anti` source:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments and clump pairs, member/ICL photo-z blends, z_phot at the eazy grid edge (> 15), diffraction
spikes, detector artefacts, lens-model error (the model is HST-era).

**Inputs (pinned by sha256 in `lensmodel.RXC0018_RELICS_LENSTOOL` and `summary.json`):**
- JWST pipeline catalogue `jw06882-o001_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC, 1,245 s). Choice rule,
  fixed here: of the two public F200W level-3 catalogues within 3′ (6882 o001, 1,245 s; 12782 o001, 773 s), the one
  with the longer F200W exposure. The other is not used.
- DJA v7.5 `rj0018p1626-grizli-v7.5-fix.eazypy.zout.fits` streamed out of
  `https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rj0018p1626-grizli-v7.5-fix.photoz.tar.gz` (31.1 MB; tarball
  not stored).
- Model: RELICS Lenstool v1 maps (`rxc0018-relics-lenstool`), D_LS/D_S = 1; z_lens 0.546 (Zitrin et al. 2011,
  arXiv:1002.0521, as MACS J0018.5+1626); cosmology H0 70, Ω_m 0.3, as for the other RELICS fields. It must pass
  `validate` or it gets no run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2, aligned ≤
30°, anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map model): the control is the background alignment fraction.

**Outcome rule (frozen; `arcs_macs0647.outcome` through `arcs_rxj2129.decide`, unchanged and unit-tested).** With one
model the same run is passed as primary and second, so the rule reduces to:
1. control: the run passes if ≥ 70 % of its strong-shear background rows are aligned (n = 0 fails);
2. if it passes: **replicates** if P(≥ n_anti) < 0.05 under the Poisson shape-noise expectation, else **does not
   replicate**;
3. otherwise **no result**.

**Pre-registered secondary checks (do not change the outcome):** z_phot > 15 rows among the anti are listed; every
anti row gets an F200W cutout looked at by eye.

**What changes the plan:** "replicates" sends each anti row to the clump-pair census (#273 code) and
`/vet-candidate` before any reading. "Does not replicate" or "no result" adds one entry to the S1063 replication
tally and nothing else.

**Access note.** archive.stsci.edu served an incomplete TLS chain on 2026-10-20 (leaf issued by "InCommon
Intermediate CA - OVG2C", intermediate not sent). The maps were fetched with the intermediate from the leaf's own
AIA URL (`http://repository.emsign.com/certs/EEEMIncommonOVG2C.crt`, chaining to the bundled emSign root) added to
the CA file; verification stayed on. The sha256 pins make the route irrelevant to the result.
