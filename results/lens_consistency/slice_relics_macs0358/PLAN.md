# SLICE JWST-era images of MACS J0358.8−2955 vs two RELICS HST-era maps of Abell 3192 — frozen plan

Frozen and pushed before any model map is read (2026-10-25; the 8 map files were downloaded only to pin their
sha256 and read their FITS headers). Rules are those of `../slice_relics_rxj0152/PLAN.md` (D-034 addendum be, itself
addenda bc rules 1–6 and bd rules 7–8); the differences are listed. Script `scripts/slice_relics_macs0358.py`; rule
tests in `tests/test_slice_relics_macs0358.py`.

## Question and hypotheses
- **Q1 (amplitude):** a fourth cluster for the RELICS Lenstool v1 amplitude question (addenda bd, be: needed a
  factor ≈ 0.64–0.80 on SMACS 0723 and MACS J0553, none on RXC J0232 and RX J0152). Rule 7 class of Lenstool v1 vs
  GLAFIC v2: `H_sys` (Lenstool v1 `lenstool_like`, GLAFIC `near_1`), `H_none` (both `near_1`), else `undecided`.
  Either outcome only refines the practical rule; neither is a lensing anomaly.
- **Q2 (positions):** does any SLICE image of families 4–6 miss every usable model by > 3″, or any image of the
  spectroscopic families 1–3 at its z_spec?
  - Ordinary expectation: no.
- **What would change the plan:** an image missed by every usable model → a `lensing_residual_input` (A) or
  `redshift_residual_input` (B), checked first against the input lists, the model's own image list and the SLICE
  Lenstool model (LAM atlas, "restricted" today) before any interpretation. A model failing P, K0 or F gives no reading.

## Inputs (observed)
SLICE (Cerny et al. 2026) Table 8 block "MACS J0358.8−2955" from the same e-print tarball (sha256 f7120d94…a9a7),
parsed with `slice_relics_heldout.parse_block` into `slice_images.ecsv`: 16 images, 6 sources (all clump 1).
SLICE Table 1: z_lens 0.425, centre (59.7201370, −29.9298850).
- **Families 1–3** (7 images): the three sources of Hsu et al. 2013, z_spec 1.261, 3.07, 3.07 (SLICE Table 8).
  Gate K0 and frame check F use them at the fitted z; reading B uses them at z_spec.
- **Families 4–6** (9 images; column `new_in_jwst` = True, meaning here "added by SLICE"): SLICE "add[s] three new
  sources" to Hsu et al. 2013; z are SLICE model values (2.44, 2.11, 1.57). Reading A, at the fitted z.
- **Known weaknesses, stated in advance:**
  - The RELICS model directories hold no readme or constraint list (checked 2026-10-25), so it is unknown whether
    the RELICS models used families 4–6. Reading A is therefore **exploratory, not held-out** (trial rule T5). A miss
    still says the RELICS map does not reproduce a SLICE image; a match does not test prediction.
  - Families 2 and 3 (same z_spec 3.07) lie 0.8–1.1″ apart, inside the 3″ match radius: a model can match one
    family's images with the other's predictions. K0 and reading B on them are weak; listed, not re-weighted.
  - Hsu et al. 2013 study the cluster's three-dimensional (multi-component) geometry; the maps are single-plane at
    D_LS/D_S = 1 and the fit uses z_lens 0.425 (ASSUMPTION, as RELICS and SLICE do).

| Model | Map set (`abell3192`, HLSP DOI 10.17909/T9SP45) | Registry name |
|---|---|---|
| Lenstool v1 | `lenstool/…_lenstool_v1` (1000² px, 0.2″) | `a3192-relics-lenstool` |
| GLAFIC v2 | `glafic/v2/…_glafic_v2` (2400² px, 0.1″) | `a3192-relics-glafic` |

## Rules (ASSUMPTIONs; per model unless stated; as addendum be)
1. **Redshift rule:** fitted z on the 60-step log grid z_lens + 0.2 … 12.
2. **Gate K0** (families 1–3, fitted z): rms ≤ 2.0″ **and** matched fraction (≤ 3″) ≥ 0.8.
3. **Frame check F** (families 1–3): |median| observed − nearest predicted ≤ 0.5″ in RA·cos δ **and** in Dec.
4. **Gate P:** median ratio of deflection-derived to published κ **and** of μ(z = 2) both in [0.95, 1.05]. A model is
   **usable** only if it passes P, K0 and F.
5. **Reading A** (families 4–6, fitted z) and **reading B** (families 1–3 at z_spec, no fit): per-model misses
   (> 3″) listed; the **joint flag** is an image missed by **every** usable model (one usable model: its own misses).
6. Families with more or fewer predicted than catalogued images are listed (not read).
7. **Amplitude class:** median over families (not at the z-grid lower edge) of D_LS/D_S(z_fit) / D_LS/D_S(z_SLICE):
   `near_1` in [0.85, 1.15], `lenstool_like` ≤ 0.75, else `intermediate`; usable models only; Q1 as above
   (`undecided` when either model is not usable).

## Not done here
No imaging, photometry or extra-image search. Post hoc listing of unmatched predicted images (`--posthoc`) only if
labelled post hoc.
