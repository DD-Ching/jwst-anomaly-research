# W3 ZTF v18 — PSF photometry check of the two v17 sources: validation passes in 1 of 9 NDs (1624) and 1 of 3 NDs (1625); both validated NDs are normal; both sources stay undetermined

- **Post hoc:** these are known items, so this is not a blind test.
  - `PLAN.md` (3d67081) fixed the method, controls, error terms and pass conditions after the data fetch and before
    any target fit.
  - Departures found in review are listed under "Deviations". No setting was chosen by the target outcome.
- **Labels:** IBE `_sciimg`, `_mskimg`, `_psfcat`, `_sciimgdao.psf` and `_sciimgdaopsfcent.fits` are `observed`;
  fits, ratios and classes are `derived`; injections are `simulated` on real images.
- **Downloads:** 48.4 MB (`manifest.json`).
- **Scope:** no search, no change to the v13 rule, no veto, no limit, no exotic interpretation.

## 1. Tools and product check (`manifest.json`, `daopsf.py.txt`)
- **PSF products:** the official per-exposure DAOPHOT PSF exists for all 20 exposures. Its form is a GAUSSIAN plus
  47×47 lookup tables with linear spatial terms.
- **Reader adapter** (`daopsf.py.txt`): it evaluates the model at the source's own quadrant position.
  - Rendered at the centre, it matches the official centre image to ≤ 5.24e-4 in all 20 exposures; the centre-image
    peaks are 0.038–0.212.
  - A 0.1 px shift gives ≈ 0.013 in one exposure (not stored).
  - The order of the spatial terms was also checked on real stars in review (not stored): swapping ΔX and ΔY raises
    the median reduced χ² from 3.7–6.0 to 17–20.
- **Fitter:** photutils 3.1.1 `PSFPhotometry` with an `ImagePSF` (oversampling 4), the official mask, a pixel error
  array and fixed positions. Each exposure is registered on ≥ 319 psfcat stars. No new solver.

## 2. Validation (`measure.json` → `validation`; `measure_v1a_trials.json`; `summary.json`)

### 1624, bad columns (per ND)
- **V1a, real stars:** 6 class-A stars (target magnitude ± 0.5) and 5 class-B stars (1 mag fainter), 1,287 trials.
  The controls are fitted with a *different* exposure's PSF, and the ND's bad-pixel geometry is moved by whole
  pixels.
- **V1b, injection on the real bad columns:** 6 positions × 2 levels per ND. The star is injected with the other
  exposure's PSF and fitted with the ND's own PSF under the real mask at the injection site.

| ND | info | V1a A: median r / σ / within 3σ | V1a B | V1b recovered/true (1.0, 0.398) | passes |
|---|---|---|---|---|---|
| 206151353 | 0.603 | 1.008 / 0.032 / **0.94** | 1.020 / 0.034 / **0.94** | **0.927**, **0.929** | no |
| 206250648 | 0.410 | 0.960 / 0.083 / 0.99 | **0.930** / 0.121 / 0.98 | 1.026, 0.958 | no |
| 206251129 | 0.887 | 1.000 / 0.002 / 0.99 | 1.000 / 0.002 / 1.00 | 1.039, 1.023 | **yes** |
| 206251369 | 0.888 | 1.000 / 0.025 / 0.96 | 0.996 / 0.027 / 0.95 | **0.927**, **0.905**; 1.0 level classed normal in **4/6** | no |
| 206350236 | 0.311 | 0.990 / 0.088 / **0.91** | 1.007 / 0.101 / 0.97 | **0.893**, **0.834** | no |
| 206350476 | 0.398 | 0.990 / 0.088 / **0.91** | 1.007 / 0.101 / 0.97 | 1.040, 1.022 | no |
| 206350717 | 0.303 | 0.990 / 0.088 / **0.91** | 1.007 / 0.101 / 0.97 | **1.072**, **1.050** | no |
| 206350958 | 0.323 | 0.969 / 0.083 / 1.00 | **0.923** / 0.119 / 0.98 | **0.932**, **0.839** | no |
| 206351197 | 0.678 | 1.000 / 0.025 / 0.96 | 0.996 / 0.027 / 0.95 | 0.977, **0.944** | no |

- Bold marks the condition that fails.
- **Same whole-pixel geometry:** 206350236, 476 and 717 share it, so their V1a results are identical.
  - V1a also tests 206251129 at milder masking than the target's (control info median 0.956 against 0.887, not
    stored).
- **Classification in V1b:** apart from 206251369 (4 of 6 at the 1.0 level), every 1.0-level injection was classed
  normal and every 0.398-level one dim. Within V1b, a 1-mag drop was not "repaired" into a normal flux.
- **What fails is the flux scale.**
  - In V1b the injected PSF (other exposure) and the fitted PSF (own) differ. With the core masked, the recovered
    flux moves by −17 % to +7 % (medians 0.834–1.072).
  - In V1a the formal errors under-cover in 4 NDs (2 distinct geometries).
- The 1624 seeing pairs differ by ≤ 0.25″ (`summary.json`), so the V1b offsets are not explained by large seeing
  differences.

### 1625, blend (V2, per exposure)
- **Method:** real faint stamps (g 17.7–18.7, 3–5 stars), injected at the target's offset from the real 6.6″
  neighbour, rotated 90°, 180° and 270°. Fitted jointly with a different exposure's PSF (the nearest in time).

| exposure | seeing (other PSF's) | recovered/true (1.0, 0.398) | robust σ | passes |
|---|---|---|---|---|
| 213941192 (detection) | 2.27″ (2.05″) | 1.016, 1.050 | 0.057, 0.154 | **yes** |
| 214131946 (detection) | 2.05″ (**4.05″**) | **−1.15, −4.64** | 0.53, 1.46 | no |
| 214237995 (ND) | 4.05″ (**2.05″**) | **1.34, 1.92** | 0.13, **0.37** | no |
| 214543142 (ND) | 4.64″ (**2.82″**) | **1.75, 2.95** | 0.06, 0.18 | no |
| 214647008 (ND) | 4.24″ (own; only field-1748 exposure) | 1.007, 1.008 | 0.054, 0.089 | yes, but **own PSF** (not independent) |
| 214648674 (detection) | 2.82″ (1.86″) | 1.058, **1.230** | 0.047, 0.110 | no |
| 214746563 (detection) | 1.86″ (2.82″) | 0.965, **0.605** | 0.076, **0.279** | no |

- **Degeneracy:** none. The neighbour's fitted flux changes ≤ 0.14 % between the two levels (median per exposure).
- **Failures:** they occur where the independent PSF's seeing differs by 0.96–2.0″, and the faint-source recovery
  ranges from −4.6 to 3.0× the truth. That is consistent with the faint source absorbing the PSF mismatch in the
  neighbour's wing; this reading is not separately tested.

## 3. Target measurements (`measure.json`, `summary.json`)
- **Reference fluxes:** F_ref is the median of the 4 detections. 1624: 3,232 µJy (scatter 39); 1625: 183 µJy
  (scatter 7). For 1625, 3 of those 4 detections failed V2 individually (deviation 2).
- **σ_R terms:**
  - fit error;
  - validation scatter;
  - bracket scatter;
  - PSF mismatch |F_own − F_other|;
  - a position term (deviation 1): the target shifted ±σ_pos, which is 0.095 px (1624) and 0.037 px (1625).

| source | judgeable ND | R ± σ_R | class | 1-mag drop (R 0.398) |
|---|---|---|---|---|
| 1624 | 206251129 (HMJD 59816.51) | 0.975 ± 0.057 | normal | excluded at > 3σ in this exposure |
| 1625 | 214647008 (HMJD 59900.47) | 0.874 ± 0.063 | normal | excluded at > 3σ, with caveats (below) |

- **Caveats for 1625 / 214647008:**
  - its V2 pass used its own PSF, so it is not independent;
  - σ_R has no PSF-mismatch term for it (zero by construction).
- **Detections** classed normal:
  - 1624: all 4 (R 0.98–1.01, σ 0.05–0.06);
  - 1625: 213941192 only (1.02 ± 0.08); the other 3 are unknown (validation failed).
- **Unknown NDs** (validation failed; their fluxes are **not** results): 8 of 9 in 1624 and 2 of 3 in 1625.
  - The stored, unvalidated ratios are 0.93–1.12 (1624) and 0.92–1.02 (1625). They are not used for any class.
  - 206350717 sits at the information threshold (0.303 against ≥ 0.30); it is unknown because its validation
    failed.
- **Residuals:**
  - 1624: target-box reduced χ² is 4.4–36.0 in the detections. A PSF mismatch at high S/N is an interpretation,
    not tested.
  - 1625: the neighbour-ring reduced χ² is 1.07–1.99 in all exposures; no systematic wing residual.
- **Per source (v17 PLAN section 4):** both are **undetermined**. No ND is dim, but most NDs are unknown.
  - Excluding a 1-mag drop in one exposure does not show that the source never varied.

## 4. Correction to #387 (v17) wording
- v17 D3's 214543142 value (2.01 ± 0.31) is not a case of an error too large to decide. It sits 3.2σ *above* the
  normal level and is unvalidated.
- v17 D3 was aperture photometry with rotated-aperture wing subtraction, not a PSF fit. v18 V2 (PSF fits, a
  different method) shows that wide-seeing blend fits can be biased high.
- So a bias in the v17 value is possible but untested there.

## 5. End of this follow-up
- **Done:** the method was validated, the passing exposures measured, the failing ones kept unknown.
- **Not continued:** no further method for these two sources in this round.
- **What would reopen each:**
  - **1624:** an independent PSF model accurate enough on the wings that fitting without the core reaches ≤ 5 %
    bias, or difference-image forced photometry. The latter needs a ZTF forced-photometry account (not registered).
  - **1625:** an independent PSF from an exposure of matching seeing (the same-quadrant exposures differ by 0.6–2.0″
    here), a good-seeing exposure inside HMJD 59896.38–59900.47, or a difference-image measurement.
- **Kept for reuse:**
  - `manifest.json` (every URL and size);
  - `stamps.npz` (61×61 science and mask stamps with WCS and headers, plus the parsed PSF models of all 20
    exposures; 0.92 MB, committed although `*.npz` is ignored, at the owner's request);
  - control lists and all trials (`measure.json`, `measure_v1a_trials.json`).
  - The full 600″ cutouts stay in `data/cache/w3_ztf/v18` (not committed; cloud containers are temporary).

## Deviations (post hoc, from review)
1. **Registration sign bug fixed.** Positions were off by twice the registration offset, up to 0.27 px. Fixing it
   changed no pass/fail, class or verdict, but changed stored numbers: 1625 214647008 R 1.000 → 0.874; the neighbour
   ring χ² in the wide-seeing NDs 3.5 → 1.07.
   - **Added a position term to σ_R** (not in PLAN section 3).
2. **References built partly from unvalidated detections.**
   - **1624:** the detections were taken as validated without a V1 test, though their fit boxes contain bad pixels
     (info 0.94–1.00). Their validation term is the repeatability of the unmasked control fits (0.041), which is not
     defined in PLAN.
   - **1625:** F_ref uses 3 detections that failed V2. One of them (214131946) has an other-PSF flux of −1.28× its
     own.
3. **V1a geometry and V1b setup.**
   - V1a moves the ND geometry by whole pixels (the sub-pixel phase is lost).
   - V1b fits under the real mask at the injection site, not geometry g, and its error array has no Poisson term
     for the injected star.
4. **Minimum counts not in PLAN:** V1a n ≥ 10, V1b ≥ 4 (PLAN: 6 positions; 6 were used), V2 ≥ 6.
   - The V2 neighbour test uses the median per exposure.
   - V2 classes use an absolute σ (max(s1, 0.398·s2)), while targets use a relative σ times |F|.
5. **Diagnostics only:** the 1625 target-box χ² was first computed with a target-only model, and the neighbour-ring
   residual was added. No class changed.
