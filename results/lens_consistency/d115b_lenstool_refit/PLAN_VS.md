# D115b stage VS: is Lenstool's source-plane χ² continuous at `best.par`? (plan, frozen 2026-10-26, before any VS run)

## Question and why this stage
Stage V3′ (RESULT_V3P.md) found the image-plane χ² (`forme -1`) discontinuous at `best.par` itself: single offsets of
3 × 10⁻⁵ of a parameter's range switch family D52b's χ² by +12.1. Every local re-fit validated against that objective
is uninformative. Lenstool's `image` keyword `forme` ≥ 0 computes the χ² in the source plane instead ("<0 image plane
chi2. >=0 source plane chi2"; `forme 0` uses "an amplification averaged over all the images of a system"; Lenstool
wiki, Image section, accessed 2026-10-26, SOURCES.md). A source-plane χ² needs no image finding, so it has no
lost-image or branch-switch mechanism. VS asks only whether that objective is continuous at `best.par`; it fits nothing.

## Hypothesis, ordinary expectations and what would change the plan
- H: with `forme 0`, small single-parameter offsets change χ² smoothly (no jumps).
- Against it: the amplification average can switch when a critical curve crosses an image (μ changes sign or diverges),
  which is a source-plane discontinuity of its own; Lenstool may also print no χ² in this mode.
- Plan change: `pass` → a frozen F-S plan (the PLAN_F2L.md Nelder–Mead, with the source-plane χ² as objective and the
  image-plane χ² reported at its end point), validated from a perturbed start first. `fail` → the source-plane route is
  closed as well; the next handling is image-plane χ² with D52b's matching fixed. `undetermined` → no reading; the
  missing evaluations are listed.

## Run file and evaluations (tool `scripts/d115b_lenstool_refit.py --stage vs`)
- The stage V′ run file (`build_vprime`: `best.par`, `input.par`'s full-precision reference, `inverse 0`, `chi2 1`,
  `input.par`'s image section) with the single line `forme -1` replaced by `forme 0` (`build_vs`; nothing else
  changes; 8515 at its `best.par` σ).
- Objective: the last `INFO: chi2 <value>` Lenstool prints (`read_chi2pos`). A run that fails, exceeds 60 s
  (PLAN_F2L.md timeout) or prints no such line has no value.
- VS0: `best.par` itself. Without a value there, the gate is `undetermined` and nothing else runs.
- Offsets: the 28 O1–O6 parameters of PLAN_F2L.md (`free_params`), each sign s ∈ {+1, −1}, steps
  δ ∈ (3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2) of the parameter's range (ASSUMPTION): 336 evaluations. A (parameter,
  sign) whose largest offset leaves [0, 1] in u is skipped (not clipped), so every step is the stated one.
- Recorded per (parameter, sign): the six Δ values (below); `chires.dat` is not required in this mode.

## Gate VS (function `gate_vs`, unit-tested before the run)
Δ_k = χ²(u_best + s δ_k e_j) − χ²(VS0), for k = 1…6 in increasing δ.
- A (parameter, sign) series has a **jump** iff for some k ≤ 4 (δ_k ≤ 1e-3) both |Δ_k| > 2 and
  |Δ_k| ≥ 0.6 |Δ_{k+1}| (ASSUMPTIONs: between steps that differ by 3–3.3× a smooth series changes by a factor of
  about 3 (linear) to 11 (quadratic), so a ratio ≥ 0.6 marks a step; 2 is an absolute floor for numerical noise, not
  measured; V3′'s jumps were +10 to +14).
- A series is **incomplete** iff any of its six values is missing.
- `pass` iff VS0 has a value, ≥ 1 sign per parameter is evaluated, and no series has a jump or is incomplete;
  `fail` iff any series has a jump; otherwise `undetermined`.

## Not part of this stage
No optimiser, no D115b Δχ², no change to TASKS priorities. Source-plane and image-plane χ² are not compared to each
other.

## Provenance
All χ² are `model_prediction` of the pinned `macs0416-canucs` model (D-044) and its single-parameter variants; the
image list is `observed`; the steps, the jump rule and its thresholds are ASSUMPTIONs.
