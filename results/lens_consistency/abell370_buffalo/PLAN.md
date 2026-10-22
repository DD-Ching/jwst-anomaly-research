# Abell 370: BUFFALO Niemiec+2023 model A maps, gate K (plan, frozen)

Frozen 2026-10-22 (commit time of this file), before the model was evaluated at any image position. Follows issue #465
and `../abell370_canucsref/RESULT.md` (handoff: pin BUFFALO model A, run gate K). Script: `scripts/a370_buffalo.py`;
tests `tests/test_a370_buffalo.py`. Read before freezing: the HLSP directory listing, the readme, the `input.par`, the
map FITS headers and pixel percentiles (loading facts only), and the CDS `table1.dat` class / z counts.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the published BUFFALO Abell 370 model A deflection maps (Niemiec et al. 2023, Lenstool,
  image-plane fit; `main` field, 0.12″ pixels, D_LS/D_S = 1) reproduce a published spectroscopic multiple-image list
  in the image plane (gate K)? If yes, Abell 370 gets an image-plane-validated model for the extra-image test.
- Expected: **pass or borderline**. Niemiec+2023 quote rms 0.90″ for model A on their own gold constraints (recorded
  in `../abell370_canucsref/RESULT.md`); our list is Lagattuta+2019's (overlapping, not identical, see below), and
  `m0417_images.gate_k` is a different solver (0.25″ grid, best predicted image within 3″).
- Ordinary explanations for a fail: images in our list that are not BUFFALO constraints, different image
  identifications or z between the two papers, our map interpolation near critical curves, a deflection-sign or
  frame error (checked by gate L). None is a lensing anomaly; a fail is a model-route result.
- What would change the plan: gate K passes → the next cycle pins the CANUCS DR1 Abell 370 catalogue and mosaic and
  runs stages 2–3 of `../macs0416_canucsref/PLAN.md` with this model. Gate K fails → no counter-image prediction or
  flux is read; the Abell 370 line of the extra-image route stays stopped.

## Inputs
- model_prediction: `abell370-buffalo-a-maps`, the `main` x / y deflection maps of
  `https://archive.stsci.edu/hlsps/buffalo/abell370/lens-models/lenstool_niemiec/a/sl/` (pinned by sha256 in
  `lensmodel.ABELL370_BUFFALO_A_MAPS`; 36.0 MB each), z_lens 0.375 and flat H0 70, Ωm 0.3 from the HLSP `input.par`.
  Model frame centre: the map reference pixel (= the par file's `reference 3 39.971340 -1.582260`).
- observed: Lagattuta et al. 2019, CDS `J/MNRAS/485/3738` `table1.dat` (138 images; sha256 pinned in the script).
  **Selection (frozen):** Class `gold`, `zspec` present, and Note not containing `1` (Note 1 marks images that
  Lagattuta's model does not use: 10.1 / 10.2 merged, 16.2 predicted but not seen). Expected 102 images. **ASSUMPTION:**
  z = zspec, fixed. No frame offset (the list and the model are both in the HST frame; no JWST pixel is read).
- Not used: the `full` / `par` map variants, the 100 MCMC sample maps, the magnification maps (their z is not
  documented in the readme beyond the file name).

## Frozen rules
1. **Gate L** (loading): `lens_consistency.map_check` of the model's κ (from the deflection-map Hessian) against the
   published `main` κ map: median relative difference ≤ 0.05 (ASSUMPTION, as `m1149_sharon.G1_MAX_MEDIAN_REL`). Fails →
   stop, a loading error, no gate K reading.
2. **Gate K** (`m0417_images.gate_k` via `m0417_images.run_field`, 0.25″ grid): image-plane rms of the selected images
   ≤ 1.0″ and ≥ 80 % reproduced within 3.0″. Fails → stop (no prediction read).
3. Reported, not rules: per-family rms; `run_field`'s `post_hoc` drop-unreproduced check; `extra_predicted` counts
   (written by `run_field`, not read for any conclusion in this cycle).
