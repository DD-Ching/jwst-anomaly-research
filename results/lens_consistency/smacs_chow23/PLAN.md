# SMACS 0723 Chow+2023 (MARS) maps: input check and fixed-z gate K (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any value of the deflection maps was read. Handoff of
CHANGELOG 2026-10-24 (S1063 extra-image test): "another field with public JWST-era maps and spec-z systems". A Zenodo
search (2026-10-24) found no new cluster-map record beyond those already used except Zenodo 10.5281/zenodo.10158248,
the Chow et al. 2023 (arXiv:2310.09790) free-form SMACS 0723 model fitted to JWST ERO images, not yet used here. The
reference route on SMACS is closed with ICLv2 + DJA v7.4 (D-034 addendum s); this stage asks only whether a second,
JWST-era model loads and reproduces the fixed-z images, and records where the two models' predicted image counts
differ. Script `scripts/smacs_chow23.py`; tests `tests/test_smacs_chow23.py`.

## Inputs
- Record 10158248 (published 2023-11-20; listing, readme and the `alpha_x` header only seen so far): readme states
  `alpha_x` / `alpha_y` "in arcsec, normalised to Dls/Ds = 1"; 0.06″ / px; FOV 3.26′. Header of `alpha_x`:
  3260 × 3260, BITPIX −32, `RA---CAR` / `DEC--CAR`, CDELT ∓1.66667e-5° (north up, east left), CRPIX 1630, CRVAL
  (110.827212, −73.454661), no BUNIT. Each α file 42,514,560 B; only the two α files are fetched (85 MB, below the
  200 MB rule), sha256 recorded, deleted after the run.
- Image list: `results/lens_consistency/smacs_djaref/mahler22_images.ecsv` (Mahler+22 ICLv2 `arcs.dat`, observed);
  the **16 images of systems 1, 2, 3, 5, 7** with `z_spec` True (D-034 addendum r rule 0, ASSUMPTION), each at its z.
- Chow+2023 fitted JWST-identified multiple images, so **this is a reproduction check of a fit (method), not an
  independent test** of the model.

## Hypothesis, ordinary explanations, expected value
- Question (method): do the Chow maps, loaded by `MapLensModel`, reproduce the fixed-z multiple images? A pass gives
  SMACS a second, independent-team JWST-era model for cross-model readings (multiple images predicted by one model and
  not the other); a fail stops this model on the route. Nothing exotic is read in this stage.
- Ordinary explanations of a fail: loading (unit, sign of α_x, the CAR projection, WCS), map resolution, the mean
  source near caustics, different image identifications (Chow+23 used their own list), or the model itself.

## Frozen rules (computable; functions with unit tests written before the run)
Cosmology flat ΛCDM H0 70, Ωm 0.3 (ASSUMPTION; only D_ds / D_s enters). z_lens 0.39 (as `smacs_djaref`, from
ICLv2 `best.par`; ASSUMPTION for this model). Map value × reading scale = deflection at D_ds / D_s = 1.
1. **C (coverage)** `coverage` (`spt0546_inputs`): all 16 images inside the maps. Fails → stop.
2. **U (unit and x sign)** four readings: `A` (arcsec, as given), `A_neg_x` (arcsec, α_x negated), `P` (pixels,
   × 0.06″), `P_neg_x`. Per system at its z, the rms of `backtrace_images` `dtheta_arcsec`; a reading is accepted
   when the median over the 5 systems is ≤ 1.0″ (ASSUMPTION). U passes when exactly one reading is accepted
   (`spt0546_inputs.unit_check`); else stop.
3. **K (fixed-z gate)** `m0417_images.gate_k` on the 16 images with the accepted reading, exact image-plane solutions
   (`lensmodel.imageplane_residuals`, mean back-traced source, grid step 0.25″): passes when rms ≤ 1.0″ **and**
   ≥ 80 % of the images lie within 3.0″ of a predicted image.
4. **K-WLS** (reported, not a gate): the same with the WLS family source (`a2744_uncover.imageplane`).

Reported, not gates (fixed now): per-system rms and catalogued vs predicted image counts (`bullet_gatek.per_system`)
for Chow and, computed the same way in the same run, for ICLv2 (`smacs0723-iclv2`), with the list of systems whose
predicted counts differ between the two models; ICLv2 gate K rms of D-034 addendum r (0.321″) for comparison.

## Stop rules
- No change of any rule after the first map value is read; any departure is labelled "deviation" or "post hoc".
- No flux, flag or extra-image reading in this stage. A pass hands the cross-model comparison to a separately frozen
  plan.
