# MACS0647: CLASH HST-era maps on the JWST-first systems (held-out test) — frozen plan

Frozen and pushed before any model read (2026-10-23). Named in D-034 addendum v. Script
`scripts/macs0647_clash_heldout.py`; rule functions unit-tested in `tests/test_macs0647_clash_heldout.py`.

## Question and hypotheses
- **Question:** do the HST-era CLASH Zitrin v2 maps (PIEMD + eNFW `macs0647-clash-nfw`, LTM-Gauss
  `macs0647-clash-ltm`; pinned sha256 in `lensmodel.MACS0647_CLASH_*`) reproduce the multiple-image
  systems first identified in JWST data (Meena+2023 systems 13–23, `new_in_jwst`), which were not among
  their constraints?
- **Ordinary expectations:** an HST-era model extrapolated to new systems typically misses some images by
  more than its fitted rms (model incompleteness, substructure, wrong family assignment, wrong photo-z);
  a missed image is a lensing-residual *input*, never a candidate, and needs a JWST-era model before any
  interpretation.
- **What would change the plan:** if both models pass gate K0 and reproduce every held-out image, the
  CLASH maps are adequate predictors here and the extra-image route could run on them; if they fail
  gate K0, the maps are not a usable reference on this field (route stays closed, D-034 addendum v).

## Inputs (observed)
`results/lens_consistency/macs0647_inputs/meena23_images.ecsv` (72 images, 25 families). Images marked
`uncertain` (8.3, 16.3, 16.4) are left out of every rule. A family enters only with ≥ 2 certain images.
No spectroscopic redshifts exist in the list.

## Rules (ASSUMPTIONs; per model unless stated)
1. **Redshift rule.** Each family takes a fitted z: the z on a 60-step log grid from z_lens + 0.2
   (0.791) to 12 (deviation from `m0417_images`' upper end 10: MACS0647-JD lies at z ≈ 10–11) that
   reproduces every image of the family within `MATCH_ARCSEC` = 3″, smallest rms among those, else
   smallest rms (`m0417_images.best_z`). The photometric z is only compared afterwards (reported, no
   rule).
2. **Gate K0** (validation set: certain images of the pre-JWST families 1–12, 6-1 / 6-2 / 6-3 as three
   families): image-plane rms ≤ 2.0″ **and** matched fraction (≤ 3″) ≥ 0.8. The 2.0″ (not the 1.0″ of the
   spectroscopic gates of D-034) because these families have fitted, not spectroscopic, z and the
   model's constraint set is not documented image by image in the readme.
3. **Frame check F** (validation set, after rule 1): median signed offset, observed minus nearest
   predicted image, along RA·cos δ and along Dec, each |median| ≤ 0.5″. Fails → no reading for that
   model.
4. **Held-out reading** (certain images of families 13–23, fitted z): an image is reproduced when its
   nearest predicted image lies ≤ 3″ away. The held-out rms and matched fraction are reported beside
   gate K0's. An image not reproduced by **every** model that passes both K0 **and** F gets the reading
   `lensing_residual_input`; with no such model, no image gets a reading.
5. Families with more or fewer predicted than catalogued images are listed (stage-2 inputs; not read).

## Not done here
No photometry, no reference catalogue, no extra-image search, no flux ratio. Kappa / μ map checks of
these maps were done when they were pinned (2026-10-16) and are not repeated.
