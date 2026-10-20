# W3 ZTF v19 — validation design check: v18 mostly tested cross-exposure PSF mismatch; after one correction, 1 of 9 NDs (1624) and 2 of 3 NDs (1625) are reliably measured, all normal; both sources stay undetermined

- **Post hoc:** a method correction on known data, not a blind test.
  - `PLAN.md` (3599220) holds the audit and the single correction. Both were committed before the recomputation,
    which then ran once.
  - v18 results are kept and are now labelled a **stress test** (cross-exposure PSF).
- **Inputs:** cached v18 data only (0 MB downloaded). Same tools: photutils 3.1.1 `PSFPhotometry` and the v18
  DAOPHOT reader.
- **Labels:** images, masks, psfcat and PSF files are `observed`; fits and classes are `derived`; injections are
  `simulated` with real stamps on real images.
- **Scope:** no search, no rule change, no veto, no limit, no exotic interpretation.

## 1. Audit result (PLAN section 1)
- **V2 (1625) mainly measured cross-exposure PSF mismatch.**
  - It fitted real stamps with another exposure's model (up to 2.0″ different seeing). The targets use their own
    exposure's model.
  - Its truth fit used the same mismatched model.
  - The injected stamp noise was missing from σ.
- **V1a (1624) mixed three effects:** masking loss, cross-exposure mismatch amplified by the masked core, and error
  coverage. It also lost the sub-pixel phase.
- **V1b (1624)** injected an ideal model of another exposure, fitted with the own model. That is model-to-model
  mismatch.
- **The reference** came from detections that were not validated (1624) or failed V2 (1625).
- **v18's "own PSF, so not independent"** conflated three things:
  - **same exposure:** the normal measurement condition;
  - **same PSF-star set:** cannot be excluded, because the official PSF-star list is not in the products;
  - **same ideal model generated and fitted:** a code check only.
- **Not all failures were design:** 1624's columns remove 11–70 % of the PSF in the NDs (info 0.30–0.89), and v19
  still finds real effects there (section 2).

## 2. Corrected validation (`v19.json` → `validation`)
Own-exposure PSF everywhere.

### 1624, V1′: phase-matched real controls (6 A + 5 B stars, x phase within ±0.20 px), masked vs unmasked
- 9 of 13 exposures pass.
- **Failures:** 206151353 (B within-3σ 0.92), 206250648 (A median 0.962, within-3σ 0.90), 206350476 (A within-3σ
  0.92) and 206350958 (A within-3σ 0.90).
- **Repeatability:** unmasked control fluxes repeat across the 13 exposures with robust σ 1.4 %. That enters every
  R as a term.

### 1624, V1b′: real stamps from the same exposure, injected on its real columns at the target's column offset and phase

| ND | V1′ | V1b′ recovered/true (1.0, 0.398) | 1.0 classed normal | passes |
|---|---|---|---|---|
| 206151353 | fail | 0.997, 0.998 | 1.0 | no (V1′) |
| 206250648 | fail | **0.941, 0.875** | 1.0 | no |
| 206251129 | pass | 0.975, 0.951 | **0.5** | no |
| 206251369 | pass | 0.987, 0.974 | 1.0 | **yes** |
| 206350236 | pass | **0.931, 0.892** | 1.0 | no |
| 206350476 | fail | 0.976, 0.966 | 1.0 | no (V1′) |
| 206350717 | pass | **0.932, 0.892** | 1.0 | no |
| 206350958 | fail | **0.922, 0.830** | 1.0 | no |
| 206351197 | pass | 0.968, **0.934** | 1.0 | no |

- Every 0.398-level injection was classed dim (`frac_dim_0398` = 1.0 everywhere).
- **The detections all pass**, so the 1624 reference is validated. 206550548 has no bad pixel in its box; the other
  three pass V1′ and V1b′.
- **Reading, not tested further:**
  - Real stamps on the real columns come back 2–17 % low in several NDs, while phase-matched masking of controls
    does not (V1′ medians 0.96–1.03).
  - This is consistent with the unflagged pixels beside the flagged columns also being affected, which a mask
    cannot remove.
  - In 206251129 the bias (2.5 %) exceeds the validation σ (≈ 1 %), so half of the normal-level injections were not
    classed normal: the errors under-cover.

### 1625, V2′: real stamps with the own PSF and the injected noise counted

| exposure | seeing | recovered/true (1.0, 0.398) | robust σ | neighbour change / faint change | ρ(faint err, neighbour) | passes |
|---|---|---|---|---|---|---|
| 213941192 (det) | 2.27″ | 0.983, 0.956 | 0.040, 0.091 | 0.032 | 0.71 | yes |
| 214131946 (det) | 2.05″ | 1.032, 1.075 | 0.035, 0.070 | 0.041 | −0.27 | yes |
| 214237995 (ND) | 4.05″ | 0.988, 0.954 | 0.033, 0.061 | 0.034 | −0.73 | **yes** |
| 214543142 (ND) | 4.64″ | **1.134, 1.332** | 0.030, 0.045 | 0.041 | 0.45 | no |
| 214647008 (ND) | 4.24″ | 1.005, 1.008 | 0.071, 0.098 | 0.066 | 0.80 | **yes** |
| 214648674 (det) | 2.82″ | 1.022, 1.053 | 0.013, 0.026 | 0.012 | −0.37 | yes |
| 214746563 (det) | 1.86″ | 0.990, 0.968 | 0.005, 0.004 | 0.019 | −0.13 | yes |

- **Degeneracy:**
  - the neighbour's fitted flux moves by 1–7 % of the faint source's change between the two levels;
  - |ρ| reaches 0.80 (214647008), under the 0.9 failure limit but a strong coupling;
  - every injected dim case was classed dim.
- **214543142** (the worst seeing) over-recovers by 13–33 % even with its own PSF.
- **Repeatability (1625):** robust σ 6.2 % across exposures, entered in every R.

## 3. Measurements and reference (`v19.json` → `exposures`)
- **1624:** reference from the 4 validated detections, F_ref = 3,232 ± 41 µJy.
- **1625:** reference from all 4 detections, now validated, F_ref = 183 ± 7 µJy.

| source | ND | R ± σ_R | class | 1-mag drop (R 0.398) |
|---|---|---|---|---|
| 1624 | 206251369 (HMJD 59816.51) | 0.985 ± 0.027 | normal | excluded in this exposure |
| 1625 | 214237995 (HMJD 59896.38) | 0.918 ± 0.091 | normal | excluded in this exposure |
| 1625 | 214647008 (HMJD 59900.47) | 0.874 ± 0.091 | normal | excluded in this exposure |

- **Unknown NDs** (validation failed), whose fluxes are not results: 8 of 9 in 1624 and 1 of 3 in 1625 (214543142).
- **σ_R components:**
  - the fit error;
  - the validation scatter;
  - the cross-exposure repeatability of controls;
  - the position term;
  - the reference error.
- **PSF systematics:** represented only by real-control scatter with the own PSF. A control may have been a PSF star,
  so the term is a **lower bound**.
- **Per source (v17 PLAN section 4):** both are **undetermined**. No ND is dim, but not every ND is judgeable.
  - Excluding a 1-mag drop in an exposure does not show that the source never varied.

## 4. Claims that change
- **v18, 1624 / 206251129 "normal, 1-mag drop excluded": withdrawn.** In v19 the real-stamp injection shows its
  errors under-cover. The judgeable 1624 ND is now 206251369.
- **v18, 1625 / 214647008 "excluded, with caveats":** now validated with its own PSF and a validated reference. It
  stands with limitation 4b (possible PSF-star overlap of the controls).
- **New:** 1625 / 214237995 is judgeable and normal.
- **v18 validation tables:** a stress test of cross-exposure PSF mismatch, not of the measurement as made.

## 5. End
- **The one correction has run;** no further photometric attempt on these two sources follows.
- **What would reopen each:**
  - **1624:** a measurement that removes or models the influence of the unflagged pixels next to the bad columns.
    Difference-image forced photometry would do; it needs a ZTF account, not registered.
  - **1625 / 214543142:** a seeing-matched validation that passes, or a good-seeing exposure inside the run.
- **Common limit:** the official PSF-star lists would be needed to make the controls demonstrably independent.
