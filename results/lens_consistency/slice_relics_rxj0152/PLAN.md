# SLICE JWST-era images vs three RELICS HST-era maps of RX J0152.7−1357 — frozen plan

Frozen and pushed before any model map is read (2026-10-25; the 12 map files were downloaded only to pin their
sha256). Rules are those of `../slice_relics_heldout/PLAN.md` (D-034 addendum bc, rules 1–6) and
`../slice_relics_glafic/PLAN.md` (addendum bd, rules 7–8); the differences are listed. Script
`scripts/slice_relics_rxj0152.py`; rule functions tested in `tests/test_slice_relics_rxj0152.py`.

## Question and hypotheses
- **Q1 (amplitude):** addendum bd derived "a RELICS Lenstool v1 map gives no reading alone" from two clusters
  (SMACS 0723, MACS J0553) where Lenstool v1 needed a D_LS/D_S factor ≈ 0.64–0.80 against GLAFIC, and two
  (RXC J0232, RXC J2211) where it did not. On a third cluster with three RELICS models, is the Lenstool v1
  factor (rule 7) `lenstool_like` while GLAFIC v2 and Zitrin LTM-Gauss v1 are `near_1`?
  - H_sys (systematic Lenstool v1 scaling): Lenstool v1 `lenstool_like`, the others `near_1`.
  - H_none (cluster-specific): all three `near_1` (or all in the same class).
  - Either outcome only refines the derived practical rule; neither is a lensing anomaly.
- **Q2 (positions):** does any image of the SLICE list miss every usable model by > 3″?
  - Ordinary expectation: no. All six HST-era families (1–6) were inputs of the RELICS models
    (Acebron et al. 2019 used the same six; SLICE text), so K0 is a self-consistency gate, not a held-out test.
- **What would change the plan:** an image missed by every usable model (rule 4 joint) → a
  `lensing_residual_input`, checked against the SLICE Lenstool model (LAM atlas, "restricted" at 10:10Z today)
  before any interpretation. A model failing P, K0 or F gives no reading.

## Inputs (observed)
SLICE (Cerny et al. 2026) Table 8 block "RX J0152.7−1357", parsed with `slice_relics_heldout.parse_block` from the
same e-print tarball (sha256 f7120d94…a9a7) into `slice_images.ecsv`; clump 1 of each source (20 rows: 7 sources).
- **Held-out (reading A):** family 7 only, the one source the SLICE text calls new in the JWST imaging (a
  galaxy-scale feature near a member overdensity in the SW; families 1–6 are "the six multiply-imaged galaxies used
  as constraints in Acebron et al. 2019").
- **Reading B:** family 1 at z_spec 3.93 (Acebron et al. 2019); the only spectroscopic family.
- **Known weakness, stated in advance:** family 7's two images are 1.45″ apart, inside the 3″ match radius, so a
  model predicting one image near them matches both. Reading A "no_flag" on family 7 is therefore uninformative;
  its predicted multiplicity (rule 6) is listed, not read.
- Text vs table: the SLICE text says nine source galaxies and names a Source #11; Table 8 lists seven. Table 8 is
  used as is.

| Model | Map set (`cl0152m13`, HLSP DOI 10.17909/T9SP45) | Registry name |
|---|---|---|
| Lenstool v1 | `lenstool/…_lenstool_v1` | `cl0152-relics-lenstool` |
| GLAFIC v2 | `glafic/v2/…_glafic_v2` | `cl0152-relics-glafic` |
| Zitrin LTM-Gauss v1 | `zitrin-ltm-gauss/v1/…_zitrin-ltm-gauss_v1` | `cl0152-relics-zitrin` |

x / y arcsec deflection at D_LS/D_S = 1 (RELICS readme v3), sha256-pinned in `jwst_anomaly.lensmodel`;
z_lens 0.8269 (SLICE Table 1); cosmology (70, 0.3) as for the other RELICS models.

## Rules (ASSUMPTIONs; per model unless stated)
1. **Redshift rule:** fitted z on the 60-step log grid z_lens + 0.2 … 12 (as addendum bc).
2. **Gate K0** (families 1–6, fitted z): rms ≤ 2.0″ **and** matched fraction (≤ 3″) ≥ 0.8.
3. **Frame check F** (families 1–6): |median| observed − nearest predicted ≤ 0.5″ in RA·cos δ **and** in Dec.
4. **Gate P** (addendum bd rule 8): median ratio of deflection-derived to published κ **and** of μ(z = 2) both in
   [0.95, 1.05] (`lens_consistency.py validate`). A model is **usable** only if it passes P **and** K0 **and** F.
5. **Reading A** (family 7, fitted z) and **reading B** (family 1 at z 3.93, no fit): per model, an image > 3″ from
   its nearest predicted image is listed; the **joint flag** is `ho.readings` over all usable models: an image
   missed by **every** usable model gets `lensing_residual_input` (A) / `redshift_residual_input` (B). With one
   usable model, its own misses are the joint flag (as addendum bc).
6. Families with more or fewer predicted than catalogued images are listed (not read).
7. **Amplitude class** (addendum bd rule 7, `slice_relics_glafic.amplitude_class`): median over families (not at the
   lower z-grid edge) of D_LS/D_S(z_fit) / D_LS/D_S(z_SLICE): `near_1` in [0.85, 1.15], `lenstool_like` ≤ 0.75,
   else `intermediate`; only for usable models. **Q1 answer** (`q1_answer`): `H_sys` if Lenstool v1 is
   `lenstool_like` and every other usable model `near_1`; `H_none` if every usable model is `near_1`;
   otherwise `undecided` (also when Lenstool v1 or both others are not usable).

## Not done here
No imaging, photometry or extra-image search. The other SLICE clusters without RELICS or CLASH maps on MAST
(Abell 68, MACS J0027, MACS J1621, RCS 0327, SPT J0358, SPT J0516, SPT J2011, PSZ1 G091, SMACS J2031) are not tested.
