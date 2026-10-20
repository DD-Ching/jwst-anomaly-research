# W3 ZTF v19 — validation design check: v18 mostly tested cross-exposure PSF mismatch; after one correction, 1 of 9 NDs (1624) and 2 of 3 NDs (1625) pass the corrected validation, all normal; both sources stay undetermined

- **Post hoc:** a method correction on known data, not a blind test.
  - `PLAN.md` (3599220) holds the audit and the single correction. Both were committed before the recomputation,
    which ran once.
  - v18 results are kept and are now labelled a **stress test** (cross-exposure PSF).
- **Inputs and tools:** cached v18 data only (no download; not stored as a number). photutils 3.1.1 `PSFPhotometry`
  (version not stored) and the v18 DAOPHOT reader.
- **Labels:** images, masks, psfcat and PSF files are `observed`; fits and classes are `derived`; injections are
  `simulated` with real stamps on real images.
- **Scope:** no search, no rule change, no veto, no limit, no exotic interpretation.

## 1. Audit result (PLAN section 1)
- **V2 (1625) mainly measured cross-exposure PSF mismatch.**
  - It fitted real stamps with another exposure's model (up to 2.0″ different seeing). The targets use their own
    exposure's model.
  - Its truth fit used the same mismatched model.
  - The injected stamp noise was missing from σ.
- **V1a (1624) mixed three effects:**
  - masking loss;
  - cross-exposure mismatch amplified by the masked core (seeing pairs ≤ 0.24″ apart);
  - error coverage.
  - It also lost the sub-pixel phase.
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

### 1624, V1b′: real stamps from the same exposure, injected on its real columns at the target's column offset and phase

| ND | V1′ | V1b′ recovered/true (1.0, 0.398) | 1.0 classed normal (stored rule) | passes |
|---|---|---|---|---|
| 206151353 | fail | 0.997, 0.998 | 1.0 | no (V1′) |
| 206250648 | fail | **0.941, 0.875** | 1.0 | no |
| 206251129 | pass | 0.975, 0.951 | **0.5** | no (see review note 1) |
| 206251369 | pass | 0.987, 0.974 | 1.0 | **yes** |
| 206350236 | pass | **0.931, 0.892** | 1.0 | no |
| 206350476 | fail | 0.976, 0.966 | 1.0 | no (V1′) |
| 206350717 | pass | **0.932, 0.892** | 1.0 | no |
| 206350958 | fail | **0.922, 0.830** | 1.0 | no |
| 206351197 | pass | 0.968, **0.934** | 1.0 | no |

- Every 0.398-level injection was classed dim.
- **Detections:** the 1624 reference is built from the 4 detections that pass V1′ and V1b′. 206550548 has no bad
  pixel in its box.
- **Reading, consistent with but not tested:**
  - **The gap:** real stamps on the real columns come back 2–17 % low in several NDs, while phase-matched masking of
    controls does not (V1′ medians 0.96–1.03).
  - **One explanation:** the unflagged pixels beside the flagged columns are affected too.
  - **Alternatives not excluded:**
    - V1′ pools controls from 13 exposures, while V1b′ uses only the same exposure;
    - the 6 V1b′ rows lie close together (rows 30–156 of the cutout), so they are not independent.

### 1625, V2′: real stamps with the own PSF and the injected noise counted

| exposure | seeing | recovered/true (1.0, 0.398) | robust σ | neighbour change / faint change | ρ pooled (per level: 1.0, 0.398) | passes |
|---|---|---|---|---|---|---|
| 213941192 (det) | 2.27″ | 0.983, 0.956 | 0.040, 0.091 | 0.032 | 0.71 (0.69, 0.76) | yes |
| 214131946 (det) | 2.05″ | 1.032, 1.075 | 0.035, 0.070 | 0.041 | −0.27 (−0.22, −0.33) | yes |
| 214237995 (ND) | 4.05″ | 0.988, 0.954 | 0.033, 0.061 | 0.034 | −0.73 (−0.72, −0.77) | **yes** |
| 214543142 (ND) | 4.64″ | **1.134, 1.332** | 0.030, 0.045 | 0.041 | 0.45 (0.52, 0.50) | no |
| 214647008 (ND) | 4.24″ | 1.005, 1.008 | 0.071, 0.098 | 0.066 | 0.80 (0.76, **0.92**) | yes as implemented (review note 2) |
| 214648674 (det) | 2.82″ | 1.022, 1.053 | 0.013, 0.026 | 0.012 | −0.37 (−0.27, −0.60) | yes |
| 214746563 (det) | 1.86″ | 0.990, 0.968 | 0.005, 0.004 | 0.019 | −0.13 (−0.28, −0.04) | yes |

- **Sources:** the per-level ρ come from `posthoc.json`.
- **Every injected dim case** was classed dim.
- **214543142** (the worst seeing) over-recovers by 13–33 % even with its own PSF.

## 3. Measurements and reference (`v19.json` → `exposures`)
- **1624:** F_ref = 3,232 ± 41 µJy from 4 detections.
- **1625:** F_ref = 183 ± 7 µJy from all 4 detections, which pass V2′.

| source | ND | R ± σ_R | class | 1-mag drop (R 0.398) |
|---|---|---|---|---|
| 1624 | 206251369 (HMJD 59816.51) | 0.985 ± 0.027 | normal | excluded at ≥ 3σ in this exposure |
| 1625 | 214237995 (HMJD 59896.38) | 0.918 ± 0.091 | normal | excluded at ≥ 3σ in this exposure |
| 1625 | 214647008 (HMJD 59900.47) | 0.874 ± 0.091 | normal | excluded at ≥ 3σ, conditional on review note 2 |

- **Unknown NDs** (validation failed under the stored rule), whose fluxes are not results: 8 of 9 in 1624 and 1 of
  3 in 1625 (214543142).
- **σ_R components:**
  - the fit error;
  - the validation scatter;
  - the position term;
  - the reference error;
  - a cross-exposure repeatability term (deviation 1).
- **PSF systematics:** represented only by real-control scatter with the own PSF. A control may have been a PSF star,
  so this is a **lower bound**.
- **Per source (v17 PLAN section 4):** both are **undetermined**. No ND is dim, but not every ND is judgeable.
  - Excluding a 1-mag drop in an exposure does not show that the source never varied.

## 4. Claims that change
- **v18, 1624 / 206251129 "normal, 1-mag drop excluded":** not confirmed. It is not validated under v19's stored
  V1b′ criterion. That criterion's injection σ is narrower than the measurement's σ_R (review note 1), so this is a
  method defect, not evidence against the exposure. The claim is **suspended**, not reversed.
- **v18, 1625 / 214647008:** passes v19's V2′ as implemented, with a reference that passes V2′. It stays
  **conditional**:
  - the dim-level ρ is 0.92;
  - the neighbour change sits at the limit under PLAN's wording;
  - limitation 4b (possible PSF-star overlap of the controls).
- **New:** 1625 / 214237995 passes and is normal.
- **v18 validation tables:** a stress test of cross-exposure PSF mismatch, not of the measurement as made.

## 5. End
- **The one correction has run;** no further photometric attempt on these two sources follows.
- **What would reopen each:**
  - **1624:**
    - a measurement that removes or models the influence of the pixels next to the bad columns, such as
      difference-image forced photometry (needs a ZTF account, not registered);
    - or a validation rerun with the defects of review note 1 fixed.
  - **1625:**
    - for 214543142, a seeing-matched validation that passes, or a good-seeing exposure inside the run;
    - for 214647008, a per-level ρ that passes.
- **Common limit:** the official PSF-star lists would make the controls demonstrably independent.

## Review notes and deviations (post hoc; `posthoc.json` is computed from the stored trials, not a rerun)
1. **V1b′ classification σ.**
   - **The defect:** V1b′ classed injections with the v18-inherited σ, the fit error combined with V1′ scatter, a
     median 1.0 % for 206251129. The measurement's σ_R also carries the repeatability, position and reference
     terms (2.7 % for 206251129).
   - **With the measurement's composition:** all normal-level injections in every 1624 exposure would be classed
     normal (`posthoc.json`). 206251129 would then pass, since its medians are within ±0.05 and V1′ passes.
   - **Not adopted:** the one-run rule holds, so the stored result stands.
2. **V2′ formulas not defined in PLAN.**
   - ρ is pooled over both levels against the pooled neighbour median. Per level, 214647008 has ρ 0.92 at the dim
     level, which fails |ρ| ≥ 0.9.
   - The neighbour change is divided by the faint *change* (0.602 F_true). Normalised by the dim-level faint truth,
     214647008 gives 0.0998 against a limit of 0.10.
   - Neither changes the stored pass. Both make 214647008 marginal.
3. **The repeatability term is an addition to PLAN** (deviation 1).
   - It is conservative and includes control noise. For 1625, 11 of 14 controls are 0.5–1.5 mag fainter than the
     target, so it partly double-counts photon noise rather than measuring a pure PSF systematic.
   - Without it, 1625 / 214647008 has σ_R ≈ 0.07 and is still normal (not stored).
4. **Phase matching has no wrap-around** (pairs such as 0.95 vs 0.05 are skipped). That shrinks n for targets near
   phase 0 or 1. V1b′ requires n ≥ 6 trials, which is stricter than v18's 4 and never binds (minimum 18).
5. **Geometry:** it is taken from floor(target) and the fit box is centred on round(x). All target geometry columns
   lie inside the fit box. Only detection 206150633 (columns at +4 and +5, phase 0.66) can lose an edge column for
   controls on the other side of phase 0.5.
