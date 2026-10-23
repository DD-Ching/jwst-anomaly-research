# W2: J2244−0550 public HST imaging test — frozen plan (v1: design, pools, split)

Owner brief 2026-10-23 after #473 (start 06:23 UTC).
- **Budget:** 4 h wall clock, 8 CPU-h, 4 cores, 3 GB download.
- **Scope:** one target and its controls; not a confirmation test for anything exotic.
- **Order:**
  - v1 (this file) is committed before any control or target is run with the pipeline. The only earlier fit was one
    tool test on SDSSJ0818+0601 F160W (§4).
  - v2 adds the per-band thresholds from the development controls and is committed before any validation control,
    injection or target is run.
  - Target residuals are computed only for a band whose gate passes.

## 0. Prior work and tool (step 1)
- **D-060 addendum (2026-10-21):** the Moffat N-image chain recovered 9/10 known lenses, but its negative control
  gave S/N 13.0 ≥ 5σ. That control was a no-lens stamp with empirical-PSF images fitted with Moffat profiles, so
  PSF mismatch alone mimicked lens light. It also missed the ~1″ lens SDSSJ0246−0825.
  - Revisit condition: a joint lens-light + point-source model (lenstronomy) keeping no-lens negative controls below
    threshold while recovering ≥ 0.9 of controls.
- **Earlier analyses of these HST data:** none found this time (MAST, arXiv / web searches for J2244−0550 and HST
  GO 17308 / 17455 results, 2026-10-23). "Not found this time" only.
- **Tool:** lenstronomy 1.14.2 (`FittingSequence`: PSO, then `psf_iteration`, then PSO; linear amplitudes;
  `UNLENSED` point sources, `SERSIC_ELLIPSE`, `UNIFORM`) and photutils 3.1.1 `EPSFBuilder` for the initial PSF.
  - They are installed in an isolated uv environment (Python 3.12, astropy 8.0.2, numpy 2.3.5, scipy 1.17.1); the
    project `.venv` is unchanged.
  - Code: `scripts/w2_j2244_hst.py` (data interface, set-up, result handling) and `scripts/w2_j2244_hst_rules.py`
    (frozen rules as functions, `Params`; tests in the project environment).
  - No optimiser, PSF solver or second method is written here.

## 1. Data (observed; byte-range cutouts from the public S3 copies of the pipeline products)
| band | product | exposures | date | set-up |
|---|---|---|---|---|
| F160W | `if8a5k030_drz.fits` (calwf3 3.7.3; 0.12825″/px, pixfrac 1, square kernel) | 4 FLT (`if8a5ktm/n/p/q`), total 422.07 s | 2024-07-05 | GO 17455, WFC3/IR IRSUB512 |
| F814W | `jf56c1010_drc.fits` (calacs, CTE-corrected; 0.05″/px, pixfrac 1, square kernel) | 2 FLC (`jf56c1z4q`, `jf56c1z6q`), 2 × 337 s | 2024-08-31 | GO 17308, ACS/WFC WFC1 |

- **One product per band.** The HAP SVM/MVM files resample the same exposures again, so they are not used and not
  counted as independent data.
- **Per-exposure products:** used only for the reproduction check of a detection (§6).
- **Data quality** (checked before any fit, `satcheck`): in F814W the brighter image's peak pixel is flagged saturated
  (DQ 256) in both 337 s exposures (peak 81 000–84 000 e⁻). F160W has no saturation flag at the target.
- **Pixel size is not resolution.** The PSF sets the resolution, and the PSF used is empirical.
- **Stored:** product names, header values, the stamp WCS origin, and the per-system fit parameters and statistics
  (JSON). Small PNGs show data, model and residual.

## 2. Model and statistic (identical for controls, injections and the target)
- **Stamp:**
  - 49 × 49 px (6.3″) in F160W, 101 × 101 px (5.05″) in F814W, centred on the catalogue position;
  - SCI in e⁻/s; WHT as the exposure map;
  - pixels with WHT < 0.5 × the stamp median are masked, which removes saturated cores dropped by the drizzle.
- **Noise:** background rms from 3σ-clipped stamp outskirts, plus Poisson noise via the exposure map, plus the PSF
  error map from lenstronomy's PSF iteration.
  - Drizzle correlated noise is not modelled per pixel. The threshold is set empirically from same-programme negative
    controls (§3) instead.
- **Initial PSF:** photutils `EPSFBuilder` (oversampling 1) on Gaia DR3 stars from fields of the same programme and
  filter (21 × 21 px F160W, 41 × 41 px F814W).
  - Star rules (`STAR_RULES`, `is_star`): RUWE ≤ 1.4; parallax/σ ≥ 3 or proper-motion S/N ≥ 5; G 16–21.5 (F160W)
    or 17.5–21.5 (F814W); no other Gaia source within 2.5″; no second peak above 5 % of the star's peak; ACS peak
    × per-exposure time ≤ 50 000 e⁻.
  - The two stars of a negative pair and any star within 3″ of a system are excluded from that system's PSF.
- **Models:**
  - M0: two point sources plus a uniform background;
  - M1: M0 plus one Sérsic, with R 0.04–2″, n 0.5–6, |e₁|, |e₂| ≤ 0.5 and centre within ±0.8″ of the image
    midpoint.
  - The two point sources are not forced to be lensed images of one source.
- **Fit:**
  - lenstronomy PSO (60 particles × 120 iterations, seed fixed), then `psf_iteration` (30 iterations, factor 0.5,
    median stack, symmetry 1, neighbour block 0.2″, PSF variance map kept, error-map radius 0.6″), then PSO again
    (σ scale 0.5).
  - Done separately for M0 and M1, so the PSF iteration can absorb extended light in either; injections test this
    (§5).
- **Statistic:** Δχ² = χ²(M0) − χ²(M1).
  - Both χ² use one noise model: M1's model and PSF error map.
  - Pixels within the core radius of either fitted image are excluded: 0.18″ in F160W, 0.10″ in F814W.
- **Status** (`rules.status`):
  - "fit unreliable": no two peaks, a fitted image > 0.15″ from its initial peak, or reduced χ²(M1) > 3;
  - Δχ² < T: "none";
  - Sérsic centre > 0.6″ from the image midpoint: "extended light offset";
  - otherwise "extended light".
  - Only "extended light" counts as a detection.

## 3. Controls (pools are computable; split fixed here)
- **Positive controls:** close lenses with independent evidence of a lens galaxy.
  - Taken from the D-056 / D-060 rebuilt quasar systems that LS DR10 could not decide (`results/w2_hst_blend/
    classes.ecsv`).
  - Each needs a catalogued lens redshift (`lens_z_known`). Caveat: the catalogue does not say whether every z_l
    comes from the lens galaxy itself.
  - Each uses the PSF library of its own programme: the stars from that programme's fields, in obs_id order, up to 40
    fields. Where a programme has several associations per target, the first obs_id is used.
  - Each must lie in an image of that band from:
    - F160W: GO 17455 (the target's programme) or GO 17199 (WFC3/IR F160W, 1000 s);
    - F814W: GO 17308 (the target's programme, 674 s).
  - Geometry on the stamp: exactly two point sources 0.7–1.8″ apart, from the model fit. A system with fewer than
    two peaks or a third peak is excluded and listed.
  - **Development** (fixed): SDSSJ0818+0601 and ULASJ0743+2457 in both bands. J0818 F160W was already fitted once in
    the tool test.
  - **Validation:** every other eligible pool member in system_id order, at most 6 per band:
    - F160W: SDSSJ0806+2006, SDSSJ1620+1203, VDESJ2325−5229;
    - F814W: SDSSJ0114+0722, J0140−1152, ULASJ0820+0812, SDSSJ0946+1835, SDSSJ1640+1932, WGD2021−4115 (then
      SDSSJ2146−0047, …, if some are ineligible).
- **Negative controls:** two real Gaia stars of the band's programme library, summed with their real noise.
  - Built by `make_pair`: integer-pixel shifts, no resampling. Separation 1.02″ (the target's), flux ratio as found,
    position angle 37° × k.
  - Fitted with a PSF built from the other library stars, so real star-to-star PSF mismatch is included. No galaxy
    light by construction.
  - Pairs use library stars in brightness order: (1, 2), (3, 4), …
  - The first 3 pairs per band are development; the next 6 are validation.
  - No system called a binary because no lens galaxy was seen is used as a negative.
- **Threshold, set from development only:** T_band = max(25, 2 × the largest Δχ² of the band's development negatives)
  (`rules.threshold`). It is written into v2 before any validation run.

## 4. Gate per band (validation controls only; `rules.gate`)
- **Pass requires all of:**
  - ≥ 5 decided validation negatives with 0 detections;
  - ≥ 3 decided validation positives with ≥ 80 % detected;
  - the PSF-absorption check (§5).
- Wilson intervals are reported. Zero false detections in a few negatives is not a zero false-positive rate.
- **If a band fails:** its target residual is not interpreted, no second method is tried, and the threshold is not
  changed. The failure is the result for that band.
- **Implementation errors:** fixed with before / after records kept. A method failure is never "fixed" by a new
  threshold.

## 5. Injections (simulated; after v2, before the target)
- **Into:** the first 2 validation negative pairs per band (real stars, real noise).
- **Profile:** Sérsic n = 4, q = 0.8, r_eff 0.2″ and 0.5″, random PA (seed = 1000 + index).
  - Placed on the line between the two images, 0.4 of the way from the fainter one.
  - Convolved with the pair's initial PSF.
  - Magnitudes F160W 19–23 and F814W 20–24 (AB, step 1).
- **Recovered:** status "extended light". Limit = the faintest magnitude of the unbroken run from the bright end
  (`rules.contiguous_limit`), per r_eff and stamp.
- **PSF absorption check:** at the brightest injection, the fitted Sérsic flux must be ≥ 0.5 × the injected flux, for
  both r_eff in both stamps.

## 6. Target (only bands that pass)
- **Run:** the same frozen run on J2244−0550.
- **If "extended light":**
  - check that it recurs in the independent per-exposure products (F160W: the four single-exposure HAP drizzles)
    and in the other band if that band passed;
  - consider a foreground galaxy, the quasar host, a neighbour and PSF residuals before any reading;
  - a single-band detection is reported as such.
- **If "none":** inject into the target stamp (§5 grid, same position rule) and report the limit only over the
  injected magnitude, size and position range. No claim of excluding all ordinary lenses, and no dark-lens or
  exotic reading.

## 7. Stop
- Stop when a band's gate fails, or when fewer than 3 positives or 5 negatives are decidable.
- Stop when the 4 h, 8 CPU-h or 3 GB budget would be exceeded.
- Stop when the work turns into tool building.

## v2 (2026-10-23 07:18 UTC): development settings and thresholds, frozen before any validation control, injection or target run
- **Development changes** (development controls only; each changed run is recorded):
  - **Iteration 1:** the v1 settings. Recorded in the RESULT log (not stored).
    - F160W negatives: Δχ² 1041, −3759, 27.5; reduced χ² 14–56.
    - F814W negatives: Δχ² 645, 176, 948.
    - Sérsic amplitudes came out negative.
    - Star stamps were drawn from the brightest library stars: in F160W 2–5× brighter than the target, in F814W
      ≈ 3.5× fainter.
  - **Before any control run (implementation):** stars are re-centred on their smoothed peak within ±6 px. The WHT
    rule now applies to the central 7 × 7 px, with ≤ 5 % low-weight pixels overall. The v1 rule (centre within
    1 px; any low-weight pixel) accepted 0 of 40 ACS fields' stars.
  - **v2 changes after the development runs:**
    1. Non-negative amplitudes (lenstronomy `check_positive_flux`). A negative "extended light" is subtracted PSF
       mismatch, not a galaxy.
    2. A detection needs a resolved Sérsic: R ≥ 0.15″ in F160W, R ≥ 0.075″ in F814W (≈ 1.2 and 1.5 px); smaller is
       "unresolved".
    3. Reduced χ²(M1) > 10, not > 3, marks a fit unreliable. All development positives have 1.9–4.0.
    4. Negative pairs are matched to the target's brightness: Gaia G within 0.75 mag of 18.0 (the images' 17.80 and
       18.21), paired in order of |G − 18.0|.
       - F814W stars with saturated cores are allowed up to G 19, because the target's brighter image is saturated
         too. They are never used for a PSF.
       - The F160W library uses all 144 GO 17455 fields (91 stars), the F814W library 80 GO 17308 fields
         (375 stars).
    5. Core radius 0.30″ in F160W and 0.15″ in F814W (v1: 0.18″ and 0.10″). Iteration 3 showed core dipoles from the
       undersampled PSF beyond 0.18″.
- **Iteration 4 = frozen settings** (`runs/dev.json`, `png/dev_*`):

| band | development negatives Δχ² (reduced χ²) | T = max(25, 2 × max) | development positives Δχ² → status |
|---|---|---|---|
| F160W | −838 (26.5), −899 (18.7), 3887 (3.7) | **7774.3** | J0743 7707 → none; J0818 5750 → none |
| F814W | 350 (8.6), 4592 (5.0), 930 (1.9) | **9184.5** | J0743 46151 → extended light; J0818 3681 → none |

- **Expectation, stated before validation:** the development numbers already put both F160W positives below T, so
  F160W is expected to fail its gate. The rule (§3, factor 2) is not changed; validation and the gate decide.
- **Validation:** 3 F160W and 10 F814W positives as listed in §3, plus 6 negative pairs per band (pairs 3–8). Then
  the injections of §5 into validation pairs 3 and 4 of each band. Then the target, only in a band whose gate passes.

## v3 (2026-10-23 ~07:55 UTC): implementation fixes after the independent review, before any rerun
The first validation (v2 settings) failed both gates; that result is kept (`runs/*_v2run*`, RESULT "Before"). The
review then found two implementation errors that make the v2 run differ from the plan:
1. **PSF error term zero.** `keep_psf_variance_map: True` made lenstronomy put back the absent initial variance map,
   so the noise had no PSF error term, contrary to §2. Fix: `keep_psf_variance_map: False`, which keeps the variance
   map the iteration estimates.
2. **F814W per-exposure time.** NCOMBINE sits in the SCI header, so the saturation cut used 674 s instead of 337 s
   (an effective 25 000 e⁻ per exposure instead of 50 000). Every star at the target's brightness was flagged
   saturated, and the F814W PSFs came from G 19–21 stars only.
   - Fix: read NCOMBINE from the SCI header.
   - The stored library's flags are re-derived from the stored stamps without re-download
     (`scripts/w2_j2244_hst_reflag.py`, `library_17308_F814W_v3flags.json`): 40 stars become unsaturated and 11 are
     dropped by the corrected rule. Stars the wrong cut had rejected (G > 19, peak 74–148 e⁻/s) are missing.
- **Negative pairs regenerated** with the same rule (|G − 18.0| order) on the corrected library.
  `scripts/w2_j2244_hst_jobs.py` now holds every job-building step.
- **Unchanged:** models, statistic, the threshold formula T = max(25, 2 × max development-negative Δχ²), status rules,
  gate, pools and split. Development is rerun first, T recomputed with the same formula and committed, then
  validation and injections. The target is run only in a band whose gate passes.
- The F160W library and pairs are unaffected by fix 2, but all runs change through fix 1.
- **v3 development rerun** (`runs/dev_v3.json`, `png/dev_*`), with thresholds by the unchanged formula, committed
  before validation:

| band | development negatives Δχ² (reduced χ²) | T | development positives Δχ² |
|---|---|---|---|
| F160W | 13601.5 (7.8), 649.3 (9.6), 1750.7 (2.3) | **27203.0** | J0743 2742.9, J0818 957.6 |
| F814W | 64.9 (7.9), 5465.4 (2.6), 805.3 (1.0) | **10930.8** | J0743 2249.0, J0818 3119.8 |

  Both bands' development positives are below T, so both are expected to fail validation. Nothing else changes.
- **File names:** the v2-run outputs were renamed after the review, so the references above to `runs/dev.json`,
  `runs/val.json`, `runs/inj.json`, `gate.json` and `runs/dev_iter3.json` now point to `runs/*_v2run*.json`,
  `gate_v2run.json` and `runs/dev_v2_iter3.json`. The v3 outputs carry the original names.
- **Errata (after the v3 runs, from the re-check):** the stored F160W T is 27202.92 (`runs/dev_v3.json`), shown as
  27203.0 above. Fix 1 keeps the estimated variance map only when a PSF-iteration step is accepted, and the outputs
  do not record which fits that was (RESULT §8).
