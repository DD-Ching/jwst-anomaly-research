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
