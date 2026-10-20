# W3 ZTF v18 — PSF photometry check of the two v17 sources: validation passes in 1 of 9 NDs (1624) and 1 of 3 NDs (1625); both validated NDs are normal; both sources stay undetermined

- **Post hoc:** these are known items, so this is not a blind test.
  - `PLAN.md` (3d67081) fixed the method, controls, error terms and pass conditions after the data fetch and
    before any target fit.
  - The validation and target results below are those of the committed plan. No setting was changed after a
    result.
- **Labels:** IBE `_sciimg`, `_mskimg`, `_psfcat`, `_sciimgdao.psf` and `_sciimgdaopsfcent.fits` are `observed`;
  fits, ratios and classes are `derived`; injections are `simulated` on real images.
- **Downloads:** 48.4 MB (`manifest.json`).
- **Scope:** no search, no change to the v13 rule, no veto, no limit, no exotic interpretation.

## 1. Tools and product check (`manifest.json`, `daopsf.py.txt`)
- **PSF products:** the official per-exposure DAOPHOT PSF exists for all 20 exposures. Its form is a GAUSSIAN plus
  47×47 lookup tables with linear spatial terms.
- **Reader adapter** (`daopsf.py.txt`): it evaluates the model at the source's own quadrant position.
  - Rendered at the centre, it matches the official centre image to ≤ 5.24e-4 (peak ≈ 0.19) in all 20 exposures.
  - A 0.1 px shift gives ≈ 0.013.
- **Fitter:** photutils 3.1.1 `PSFPhotometry` with an `ImagePSF` (oversampling 4), the official mask, a pixel error
  array and fixed positions. Each exposure is registered on ≥ 319 psfcat stars. No new solver.

## 2. Validation (`measure.json` → `validation`; `measure_v1a_trials.json`)
Controls are fitted with a *different* exposure's PSF model (the nearest in time, same quadrant), so no control is
fitted with a PSF built from its own pixels.

### 1624, bad columns (per ND geometry)
- **V1a, real stars:** 6 class-A stars (target magnitude ± 0.5) and 5 class-B stars (1 mag fainter), 1,287 trials.
- **V1b, real-column injection:** stars at the target flux and at 0.398× it, injected on the same columns at the
  target's offset, 6 positions each.

| ND | info | V1a A: median r / σ / within 3σ | V1a B | V1b recovered/true (1.0, 0.398) | passes |
|---|---|---|---|---|---|
| 206151353 | 0.60 | 1.015 / 0.031 / **0.91** | 1.023 / 0.031 / **0.92** | **0.927**, **0.929** | no |
| 206250648 | 0.41 | 0.953 / 0.082 / 0.96 | **0.914** / 0.105 / 0.98 | 1.025, 0.957 | no |
| 206251129 | 0.89 | 1.000 / 0.003 / 1.00 | 1.000 / 0.002 / 1.00 | 1.037, 1.022 | **yes** |
| 206251369 | 0.89 | 0.996 / 0.025 / 0.96 | 0.990 / 0.033 / 0.95 | **0.928**, **0.906** | no |
| 206350236 | 0.31 | 1.000 / 0.086 / **0.92** | 1.022 / 0.116 / 0.97 | **0.892**, **0.831** | no |
| 206350476 | 0.39 | 1.000 / 0.086 / **0.92** | 1.022 / 0.116 / 0.97 | 1.041, 1.023 | no |
| 206350717 | 0.30 | 1.000 / 0.086 / **0.92** | 1.022 / 0.116 / 0.97 | **1.079**, **1.054** | no |
| 206350958 | 0.32 | 0.959 / 0.087 / 0.99 | **0.903** / 0.105 / 0.97 | **0.929**, **0.835** | no |
| 206351197 | 0.69 | 0.996 / 0.025 / 0.96 | 0.990 / 0.033 / 0.95 | 0.978, **0.945** | no |

- Bold marks the condition that fails.
- **Classification is not the failure:** in V1b every 1.0-level injection was classed normal and every 0.398-level
  one dim (`frac_normal_1` = `frac_dim_0398` = 1.0). The method does not "repair" a 1-mag drop.
- **What fails:** the flux scale. With the core masked, the other exposure's PSF biases the recovered flux by −17 %
  to +8 % (V1b medians 0.831–1.079), and the formal errors under-cover in 4 geometries.
- **Why, read off the data:** the bias appears even where the two seeings differ by ≤ 0.25″ (`summary.json`).

### 1625, blend (V2, per exposure)
- **Method:** real faint stamps (g 17.7–18.7, 3–5 stars per exposure), injected at the target's offset from the
  real 6.6″ neighbour, rotated 90°, 180° and 270°.

| exposure | seeing (other PSF's) | recovered/true (1.0, 0.398) | robust σ | passes |
|---|---|---|---|---|
| 213941192 (detection) | 2.27″ (2.05″) | 1.021, 1.063 | 0.052, 0.139 | **yes** |
| 214131946 (detection) | 2.05″ (**4.05″**) | **−1.19, −4.75** | 0.53, 1.41 | no |
| 214237995 (ND) | 4.05″ (**2.05″**) | **1.44, 2.16** | 0.10, 0.20 | no |
| 214543142 (ND) | 4.64″ (**2.82″**) | **1.77, 3.00** | 0.06, 0.16 | no |
| 214647008 (ND) | 4.24″ (own; only exposure in field 1748) | 0.991, 0.984 | 0.029, 0.077 | yes, but **own PSF** (not independent) |
| 214648674 (detection) | 2.82″ (1.86″) | 1.062, **1.235** | 0.045, 0.113 | no |
| 214746563 (detection) | 1.86″ (2.82″) | 0.965, **0.605** | 0.076, **0.279** | no |

- **Degeneracy:** none. The neighbour's fitted flux changes ≤ 0.1 % between the two levels everywhere.
- **What fails:** the flux scale, when the independent PSF comes from an exposure of different seeing (up to 2″
  apart). The faint source then absorbs the neighbour's wing mismatch.

## 3. Target measurements (`measure.json`, `summary.json`)
- **Reference fluxes:** F_ref is the median of the 4 detections. 1624: 3,231 µJy (scatter 48); 1625: 181 µJy
  (scatter 8).
  - For 1625, 3 of those 4 detections failed V2 individually, so its reference is only partly validated.

| source | judgeable ND | R ± σ_R | class | 1-mag drop (R 0.398) |
|---|---|---|---|---|
| 1624 | 206251129 (HMJD 59816.51) | 0.968 ± 0.052 | normal | excluded (> 3σ) |
| 1625 | 214647008 (HMJD 59900.47) | 1.000 ± 0.058 | normal | excluded, with the own-PSF caveat |

- **Detections** classed normal:
  - 1624: all 4 (R 0.98–1.02, σ 0.05);
  - 1625: 213941192 only (1.03 ± 0.08); the other 3 are unknown (validation failed).
- **Unknown NDs** (validation failed; their fluxes are **not** results): 8 of 9 in 1624 and 2 of 3 in 1625.
  - The stored, unvalidated ratios are 0.95–1.10 (1624) and 0.74–0.77 (1625). They are not used for any class.
- **Residuals:**
  - 1624: target-box reduced χ² is 2.5–36.6 in the detections (bright source, PSF mismatch at high S/N).
  - 1625: the neighbour-wing ring χ² is 3.5–3.6 in the two unvalidated wide-seeing NDs, against 1.3–2.2 elsewhere
    (systematic wing residual).
- **Per source (v17 PLAN section 4):** both are **undetermined**. No ND is dim, but most NDs are unknown.
  - Excluding a 1-mag drop in one exposure does not show that the source never varied.

## 4. Correction to #387 (v17) wording
- v17 D3's 214543142 value (2.01 ± 0.31) is not a case of an error too large to decide. It sits 3σ *above* the
  normal level.
- v18 V2 shows the same kind of wide-seeing blend measurement over-recovers by 1.4–3.0× with a mismatched PSF.
- So the v17 D3 values are better read as **unvalidated and probably systematically biased**, not as noisy.

## 5. End of this follow-up
- **Done:** the method was validated, the passing exposures measured, the failing ones kept unknown.
- **Not continued:** no further method for these two sources in this round.
- **What would reopen each:**
  - **1624:** a PSF model that is both independent of the target exposure and accurate on the wings, so that
    fitting without the core reaches ≤ 5 % bias (the V1 failures). Or difference-image forced photometry, which
    needs a ZTF forced-photometry account (not registered).
  - **1625:** an independent PSF from an exposure of matching seeing (the same-quadrant exposures differ by
    0.6–2.0″ here). Or a good-seeing exposure inside HMJD 59896.38–59900.47, or a difference-image measurement.
- **Kept for reuse:**
  - `manifest.json` (every URL and size);
  - `stamps.npz` (61×61 science and mask stamps with WCS and headers, plus the parsed PSF models of all 20
    exposures; 0.92 MB, committed although `*.npz` is ignored, at the owner's request);
  - the control lists and all trials (`measure.json`, `measure_v1a_trials.json`).
  - The full 600″ cutouts stay in `data/cache/w3_ztf/v18` (not committed; cloud containers are temporary).

## Deviations
- **Diagnostics only:** the residual diagnostics were fixed before recording. The first run computed the 1625
  target-box χ² with a target-only model, and the neighbour-ring residual promised in PLAN 2 was missing. No class
  or validation value changed (validation output identical, checked).
