# MACS J2129: CLASH HST-era maps on the SLICE JWST-era image list (held-out test) — frozen plan

Frozen and pushed before any model map is read (2026-10-25). Inputs and rules mirror
`results/lens_consistency/macs0647_clash_heldout/PLAN.md` (D-034 addendum v); the differences are listed.
Script `scripts/macs2129_clash_heldout.py`, rule functions unit-tested in `tests/test_macs2129_clash_heldout.py`
(next stage; not written yet).

## Question and hypotheses
- **Question:** do the HST-era CLASH Zitrin v2 maps of MACS J2129.4-0741 (PIEMD + eNFW and LTM-Gauss, readme of
  2015-02-18, z_lens 0.590; 32 SL constraints, 6 free source redshifts, quoted rms 1.86″ / 2.45″) reproduce
  (a) the three source galaxies first identified in SLICE JWST imaging (sources 16–18) and (b) the families
  whose spectroscopic redshifts (Caminha+2019, Jauzac+2021) postdate the maps, at those redshifts?
- **Ordinary expectations:** a 2015 HST-era model misses some images by more than its fitted rms (model
  incompleteness, substructure, family assignment, a z that the model left free); a miss is a lensing-residual
  *input*, never a candidate.
- **What would change the plan:** both models fail gate K0 → the CLASH maps are no usable reference here and the
  MACS2129 route waits for the SLICE Lenstool maps (LAM atlas, offline on 2026-10-25). A held-out miss by every
  usable model → a lensing-residual input to check against the SLICE model before any reading.

## Inputs (observed)
`slice_images.ecsv`: SLICE Table 8 MACS J2129 block, 53 images, 17 families (no family 15), 13 families
(41 images) with spectroscopic z, 4 with SLICE model z (5, 16, 17, 18); `new_in_jwst` = 16, 17, 18 (9 images).
No uncertain flags in the table. Which images the CLASH model used is not documented (readme: 32 constraints);
families 1–14 are treated as possibly seen.

## Rules (ASSUMPTIONs; per model unless stated)
1. **Redshift rule (fitted z).** As MACS0647 rule 1: 60-step log grid from z_lens + 0.2 to 12, smallest rms among
   the z that reproduce every image within `MATCH_ARCSEC` = 3″, else smallest rms (`m0417_images.best_z`).
2. **Gate K0** (families 1–14, fitted z): rms ≤ 2.0″ **and** matched fraction (≤ 3″) ≥ 0.8.
3. **Frame check F** (same set, after rule 1): |median| observed − nearest predicted ≤ 0.5″ in RA·cos δ and Dec.
4. **Held-out reading A** (families 16–18, fitted z): reproduced when the nearest predicted image is ≤ 3″ away.
   An image not reproduced by **every** model passing K0 and F gets `lensing_residual_input`.
5. **Held-out reading B** (new vs MACS0647): the 13 spectroscopic families at their spectroscopic z (no fit),
   same 3″ match. An image not reproduced by every model passing K0 and F gets `redshift_residual_input`. Rms
   and matched fraction at the spectroscopic z are reported beside the fitted-z ones.
6. Families with more or fewer predicted than catalogued images are listed (stage-2 inputs; not read).

## Not done here
No photometry, no imaging, no extra-image search. Maps (x / y arcsec deflection) are sha256-pinned at the first
read. The SLICE Lenstool model itself is not available (LAM atlas closed); when it is, rerun reading A / B
against it as the second, JWST-era reference.
