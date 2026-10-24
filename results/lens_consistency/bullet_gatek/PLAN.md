# Bullet Cluster: gate K of the CANUCS DR1 Lenstool best-fit model on its gold list (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any value of the Bullet deflection maps was read (the two map
files were downloaded and hashed only). Follows `../bullet_inputs/RESULT.md` (D-034 addendum x). Script
`scripts/bullet_gatek.py`; tests `tests/test_bullet_gatek.py`.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the published best-fit deflection maps, loaded by our map-model code, reproduce the gold
  multiple images the model was fitted to? A pass makes the maps usable for the next stages of the reference route
  (gate C2 and the extra-image test with the CANUCS Bullet catalogue); a fail stops the route on this field.
- Ordinary explanations of a fail: a loading or frame error on our side (map orientation, units, D_LS / D_S
  scaling, WCS), map resolution or extent, or the model itself. A gold image that is not reproduced is a model
  residual of a fit to these very images, so it is **not** a lensing-residual input and no flag is read here
  (exotic readings are not considered).
- What would change the plan: K passes → the next cycle freezes the gate C2 / extra-image plan with the CANUCS
  Bullet catalogue. K fails → check the loading against Lenstool's own `chires.dat` before anything else; if our
  numbers then match Lenstool's, the model is not used for the extra-image test.

## Inputs
- model_prediction: CANUCS Bullet DR1 Lenstool fiducial best-fit maps `bulletclu-dpl1-best-50mas.fits.gz` and
  `bulletclu-dpl2-best-50mas.fits.gz` (deflection in arcsec, D_LS / D_S = 1 per the README; sha256 pinned in
  `src/jwst_anomaly/lensmodel.py`), z_lens 0.296 (`best.par` potentials), flat ΛCDM H0 70, Ωm 0.3 (README).
- observed: `../bullet_inputs/bullet_gold_images.ecsv` (135 images, 49 systems, each at its catalogue system z;
  the model's constraints).
- model output, reported only: Lenstool `best_fit_model/chires.dat` (its own image-plane residuals of the fit).
- ASSUMPTION: frame offset 0. The model was fitted to the FK5 version of the list and we use the ICRS one; the two
  differ by far less than the gate tolerances (not a gate).
- **Download reason (> 200 MB per file, CLAUDE.md):** the best-fit deflection maps are published only at 50 mas,
  219 MB and 205 MB (gzip). The 200 mas maps are posterior samples, not the fiducial model, so they cannot test the
  fitted model as published. Cloud session: kept only in the session cache and deleted at the end of the run.

## Frozen rules (computable; `m0417_images.gate_k`, same code as MACS0417 / MACS1423 / MACS1149)
1. Back-trace each gold image to its source plane at its system z, average per system, solve the image plane on a
   0.25″ grid (`m0417_images` default), and take each image's distance to its nearest predicted image.
2. **Gate K** (combined over all 135 images): image-plane rms ≤ 1.0″ **and** ≥ 80 % of images within 3.0″.
   The function decides; both conditions must hold.
3. Reported, not gates: the per-system rms; predicted image counts above or below the catalogued count; Lenstool's
   own rms from `chires.dat` next to ours.

## Stop rules
- No change of any rule after the first map value is read; any departure is labelled "deviation" or "post hoc".
- No extra-image or flag reading in this stage.
