# W3 ZTF v19 — validation design check and one correction of the v18 PSF test (post hoc)

Plan fixed 2026-10-20 after reading the v18 code and results, before any recomputation. Owner authorization: after
#389, 2 h, ≤ 100 MB. Only 708205400009336 (1624) and 709212400017577 (1625).

**Post hoc:** this is a method correction on known data, not a blind test. No new solver, no search, no change to the
v13 rule or to the v18 pass thresholds. The v18 results stay in place, relabelled a **stress test** (cross-exposure
PSF).

## 1. Audit of v18 (`../v18/psfphot_script.py.txt`)

| step | image measured | source (injected or real) | PSF used for the fit | the v18 target measurement uses |
|---|---|---|---|---|
| V1a (1624) | every exposure e (13) | real control stars in e | **another** exposure's model (nearest in time) | own-exposure model |
| V1b (1624) | ND e (real columns) | **ideal model** of another exposure, plus Poisson noise | own model of e | own-exposure model |
| V2 (1625) | exposure e | real faint stamps from e | **another** exposure's model (seeing up to 2.0″ different) | own-exposure model |
| targets | e | the target | own model of e; the other-exposure fit only as a "mismatch" term | — |
| F_ref | detections | the target | own model | detections were **not** validated (1624) or failed V2 (1625) |

### Findings
1. **V2 mainly measured cross-exposure PSF mismatch** (e.g. a 4.05″ image fitted with a 2.05″ PSF), not the target
   measurement, which uses the exposure's own model.
   - Its "truth" (the stamp's isolated fit) was also made with the mismatched PSF.
   - The injected stamp carries its own sky noise, scaled. The error array was not updated for it, so that noise
     was left out of σ while still being present in the image.
2. **V1a mixed three effects:** loss of the masked pixels; cross-exposure PSF mismatch (amplified when the core is
   masked); and the error estimate (the within-3σ test).
   - The ND geometry was moved by whole pixels, so the sub-pixel phase against the column was lost.
3. **V1b** injected an *ideal* model of another exposure and fitted it with the own model.
   - That tests model-to-model mismatch, not real stars.
   - The error array lacked the injected star's Poisson term.
4. **"Own PSF, so not independent" (v18, 214647008) was too coarse.** Three different things are involved:
   - (a) **same exposure:** the model was built from that image's stars. This is the normal measurement condition
     and is what the target uses;
   - (b) **same PSF-star set:** a control that was itself a PSF star gets an optimistic fit. The official PSF-star
     list is not in the products, so this cannot be excluded; it is a stated limitation;
   - (c) **same ideal model** generated and fitted. That is a code check only.
5. **The reference:** F_ref was built from detections that were not validated, and σ_ref was their scatter only.
6. **PSF systematics:** with no second PSF, v18's mismatch term for 214647008 was 0 by construction. That is not a
   complete error.
- **Not all failures come from the design:** in 1624 the columns remove 10–70 % of the PSF (info 0.30–0.89), and
  that loss is real.

## 2. The one correction (`v19_script.py.txt` → `v19.json`)
- **Same tools:** photutils `PSFPhotometry`, the v18 reader, registration (sign-fixed), background, mask bits, boxes
  and information condition.
- **Data:** only the cached v18 data (no download).

### Model
- Every fit (controls, injections, targets, references) uses **the exposure's own PSF** at the source position, as
  the real measurement does.

### V1' (1624), per exposure e, ND and detection alike
- **Geometry:** the bad pixels in e's target box, relative to the target.
- **Controls:** real class-A and class-B stars (the v18 lists), in any of the 13 exposures, fitted in their own
  exposure. A star is used for e's geometry only when its x sub-pixel phase matches the target's within ±0.20 px;
  the columns run in y, so the x phase sets the geometry.
- **Per control:** r = masked fit / unmasked fit.
- **Passes when** (v18 thresholds, A and B separately): n ≥ 10, |median r − 1| ≤ 0.05, robust σ ≤ 0.20, and ≥ 95 %
  of |r − 1| ≤ 3σ_tot (σ_tot = robust σ combined with each fit's formal error).

### V1b' (1624), real-stamp injection on the real columns of e
- **Stamps:** real isolated A or B stars from e itself (the background-subtracted 21×21 stamp), with x phase within
  ±0.20 px of the target's.
- **Placement:** an integer shift onto the same columns at the target's column offset, at 6 rows ≥ 15″ from any
  psfcat source.
- **Scale:** to the target's reference flux and to 0.398× it.
- **Truth:** the scale times the stamp's unmasked own-PSF fit.
- **Errors:** the error array adds the injected noise, (k · err_stamp)². The stamp noise is real and is counted once.
- **Passes when:** |median recovered/true − 1| ≤ 0.05 at both levels; the 1.0 level is classed normal and the 0.398
  level dim in ≥ 90 % of trials.

### V2' (1625), per exposure
- **Method:** as v18 V2 (real stamps from e, rotated 90°, 180° and 270° about the real neighbour), but with the own
  PSF for the truth fit and the joint fit, and with the injected noise added to the error array.
- **Passes when:** the v18 thresholds hold, |median recovered/true − 1| ≤ 0.10, robust σ ≤ 0.20, and class fractions
  ≥ 0.9. **New and stricter:**
  - |median neighbour-flux change| between the two levels ≤ 0.10 × the faint truth flux;
  - the correlation between the faint-flux error and the neighbour-flux error across trials is reported, and
    |ρ| ≥ 0.9 fails.

### Representativeness
- **1624:** the target's column offset and x phase are kept. Seeing and noise are those of the same exposure.
- **1625:** the neighbour, separation, contrast and seeing are the real ones of the same exposure.
- **Positions:** fixed at the injected truth. Registration uncertainty enters σ as the v18 position term.

## 3. Reference and errors
- **F_ref** = median own-PSF flux of the **validated** detections only (≥ 2 required; otherwise every ratio is
  unknown).
- **σ_ref** = max(std, median per-detection σ_F) / √n.
- **σ_F of each exposure combines:**
  - the fit error;
  - the validation scatter of that exposure (V1' or V2', relative);
  - the position term.
- **PSF systematics:**
  - represented only by the validation scatter of real controls fitted with their own PSF;
  - since a control may have been a PSF star (limitation 4b), this term is a **lower bound** and is reported so;
  - no separate PSF term is set to zero and called complete.
- **R** = F_e / F_ref, with σ_R from both. The v18 class rules and 1-mag tests are unchanged.

## 4. One run, then stop
- **One run:** this correction runs once.
- **Passing exposures:** both the event and the reference exposures are measured and reported.
- **Failing exposures:** they stay unknown. No further photometric attempt on these two sources follows in this
  round, and no other method is tried.
- **Wording:**
  - the v18 1625 "1-mag drop excluded" is narrowed to "conditional on an unvalidated reference and an incomplete
    PSF error" until v19 settles it;
  - v18 results are kept and labelled a stress test.
