# SLICE JWST-era images vs RELICS GLAFIC v2 maps (MACS J0553, RXC J2211) — frozen plan

Frozen and pushed before any GLAFIC map is read (2026-10-25). Second HST-era model for the two SLICE clusters of
D-034 addendum bc that have one on MAST (RXC J0232 has only Lenstool). Rules, image list and families are those of
`results/lens_consistency/slice_relics_heldout/PLAN.md` unchanged (same `slice_images.ecsv`, same rule functions);
only the model differs.

## Question and hypotheses
- **Q1 (MACS J0553 11.1):** does an independent HST-era model also miss 11.1 (3.49″ in Lenstool v1)?
- **Q2 (MACS J0553 amplitude):** Lenstool v1 needs one D_LS/D_S factor 0.637–0.708 for all 21 families (addendum bc,
  post hoc). Is that the Lenstool v1 product or model (H_L), or the SLICE list / redshift side (H_S)?
- **Ordinary expectations:** HST-era models differ in mass normalisation when built mostly on photometric z;
  RELICS Lenstool v1 SMACS 0723 maps also read about 1.25 × GLAFIC (addendum aw, post hoc). The RELICS lens-model
  readme v3 (MAST `hlsps/relics/docs/`, sha256 318f4c03…652d074) documents every deflection map as scaled to
  D_LS/D_S = 1 for all methods; addendum bc found the Lenstool κ and μ(z = 2) products consistent with that.
- **What would change the plan:** GLAFIC also missing 11.1 makes it a two-model lensing-residual input for the
  JWST-era SLICE model (still not a candidate); GLAFIC matching it closes 11.1 as a Lenstool-v1-only miss.
  A GLAFIC factor near 1 on MACS J0553 points to H_L; a factor near 0.66 points to H_S (the SLICE z or list).

## Inputs
- Images (observed): `results/lens_consistency/slice_relics_heldout/slice_images.ecsv`, read as is.
- Maps (model_prediction): `macs0553m33_glafic_v2` and `rxc2211m03_glafic_v2` best-fit x / y arcsec deflection,
  D_LS/D_S = 1 (readme v3), κ and μ(z = 2) for the product check; sha256-pinned in `jwst_anomaly.lensmodel`;
  z_lens 0.412 / 0.397 (SLICE Table 1) and cosmology (70, 0.3) as for the Lenstool runs.

## Rules (ASSUMPTIONs)
1–6. Rules 1–6 of the addendum-bc plan, per cluster, one model each: fitted-z grid, gate K0 (rms ≤ 2.0″, matched
   ≥ 0.8 within 3″ on families not held out), frame check F (|median| ≤ 0.5″), reading A (held-out families 11–21 /
   4–10 at fitted z, > 3″ → `lensing_residual_input`), reading B (spectroscopic families at z_spec, > 3″ →
   `redshift_residual_input`), image-count mismatches listed. Readings only for a model passing K0 and F.
7. **Amplitude factor (frozen here):** per family not at the grid's lower edge, f = (D_LS/D_S)(z_fit) / (D_LS/D_S)(z_SLICE)
   (the addendum-bc post hoc definition). Median over families: 0.85 ≤ median ≤ 1.15 → "near 1"; median ≤ 0.75 →
   "Lenstool-like"; otherwise "intermediate". Read only for a model passing K0 and F.
8. **Product check (gate P):** κ and μ(z = 2) from the deflection maps vs the published maps (as `lc.validate`);
   median ratio within 0.95–1.05 for both, else the model gives no reading.
9. **Two-model statement for 11.1:** "missed by both" only if GLAFIC passes K0, F and P and puts 11.1 > 3″.

## Not done here
No imaging; no other clusters; no re-run of the Lenstool numbers. The SLICE Lenstool models (LAM atlas, still
redirected to "restricted" at 09:20Z on 2026-10-25) stay the JWST-era reference.
