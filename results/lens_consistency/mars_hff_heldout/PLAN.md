# MARS HFF maps vs JWST-era images absent from their catalogue (held-out stage; plan, frozen)

Frozen 2026-10-24 (commit time of this file), before the held-out residuals were computed. Follows D-034 addendum
ao (`../mars_hff_inputs/RESULT.md`: C, U and K pass on Abell 370, MACS0416, MACS0717, MACS1149). Script
`scripts/mars_hff_heldout.py`; tests `tests/test_mars_hff_heldout.py` (written before the run).

## Inputs
- model_prediction: MARS free-form deflection maps (Cha & Jee 2023; Zenodo 10.5281/zenodo.7804575, md5 checked;
  arcsec reading, as U chose), loaded exactly as in `mars_hff_inputs.py`; cosmology flat H0 70, Ωm 0.27.
- observed: MARS catalogue rows (all classes, all z) of the same field, used only to decide which images MARS saw.
- observed: JWST-era CANUCS DR1 Lenstool multiple-image files already stored here (positions and z as listed;
  same ASSUMPTION as their stages: listed z are fixed): Abell 370 `../abell370_canucsref/canucs_a370_images.ecsv`
  (115 images; Gledhill et al. 2025), MACS0416 `../macs0416_canucsref/canucs_m0416_images.ecsv` (303 images;
  Rihtarsic et al. 2025).
- Not tested: MACS0717 and MACS1149 (no JWST-era image list is stored or pinned here; the MACS1149 CANUCS readme
  names a gold catalogue but the file used before, `../macs1149_sharon/`, is the HST-era Sharon v4 list).

## Question, ordinary explanations, what changes the plan
- Question: does a model fitted before JWST reproduce JWST-era spectroscopic multiple images it never saw? This is
  the first test in this repository where the image positions were not used to fit the model under test.
- A pass makes MARS a validated second (free-form) model on these fields for extra-image and flux-ratio stages.
- A fail is expected from ordinary causes first: model incompleteness away from the HST-era constraints, coarse maps
  (1.26–1.52″ pixels, bilinear), a frame offset, the listed z, or the mean-source estimator. A held-out image with a
  large residual is **not** a flag in this stage; no flag, flux or extra-image reading is made.

## Frozen rules (functions with unit tests before the run)
Per field, independently; family = `lensmodel.image_family` of the CANUCS id.
1. **Frame offset O** `frame_offset`: median (ΔRA cos δ, ΔDec) from each CANUCS image to its nearest MARS row
   within 1.0″ (ASSUMPTION). CANUCS positions are shifted by O into the MARS frame. Reported; with < 5 pairs O = 0.
2. **Seen / held-out** `split_seen`: a CANUCS image is *seen* when a MARS row lies within 0.5″ (ASSUMPTION) of
   its shifted position, else *held-out*. Only images with listed z > 0 enter.
3. **Anchors** `heldout_residuals`: per family, the source is the mean back-traced position (arcsec reading,
   listed z) of its *seen* images when there are ≥ 2 (type E: new image of a known family); otherwise of all its
   images when there are ≥ 2 (type N: no seen image; type M: one seen image). Each held-out image is matched to the
   nearest predicted image (`lensmodel.find_images`, grid step 0.25″); no prediction counts as not reproduced.
   Images outside the map are dropped and counted.
4. **Gate H** `gate_h`, on all tested held-out images of a field: rms ≤ 1.0″ **and** ≥ 80 % within 3.0″
   (ASSUMPTION; same thresholds as gate K). Reported separately: type E only (the cleanest prediction), and the
   field's training gate K rms from addendum ao for comparison.
5. **Outcome** per field: "MARS reproduces the held-out JWST-era images" (H passes) / "does not" (H fails) /
   "not testable" (< 5 tested held-out images; ASSUMPTION).

## Stop rules
- No rule change after the first held-out residual is computed; any departure is labelled deviation or post hoc.
- Data: the 2.7 MB zip lives in the session scratchpad and is deleted after the run; only `summary.json` and
  per-field residual tables are stored.
