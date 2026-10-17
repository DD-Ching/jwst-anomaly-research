# Per-stamp leave-one-out kernel fits: plan (written before any comparison run)

Owner instruction 2026-10-17 (follow-up to #288). Budget 02:23Z → 05:23Z. One method only.

## 1. Diagnosis (B2 F150W, first 25 stamp fits profiled with cProfile; session, not stored)
- **Per stamp:** about 1.5 s.
  - `numpy.linalg.lstsq` on a freshly built weighted design of about 97,000 × 82 runs once per clipping iteration
    (3 per fit), at 0.365 s per call: about 73 %.
  - The rest is mostly rebuilding the same design matrix: the same pixel gather from `ref`, and the same finite
    footprint (`minimum_filter`) over the whole box.
- **What all stamp fits of one band share:** the stamp pixel set, its design rows, the weights `1/sqrt(var)` and
  the targets.
- **What each fit differs in:** fit i leaves out the rows of stamp i (`mask = stamp & (lab != i)`), plus its own
  clipped rows in each iteration.
- B1 has 217 stamps and B2 has 1,199 (#288 `batch/timing.json`), so the cost grows as stamps × stamp pixels.

## 2. The one method: shared design with Gram-matrix downdating (`--kernel-stamp-fit shared`)
- **Built once per band:**
  - the row set R = all stamp pixels that pass the existing `good & fin` tests;
  - the weighted design A (R × 82) and targets;
  - G = AᵀA and g = Aᵀt.
- **Each fit i keeps the existing algorithm step by step:**
  - The rows are R minus stamp i. The tested source never takes part in its own fit.
  - Each iteration solves the normal equations of its current kept rows:
    G − A_iᵀA_i − A_cᵀA_c, where c are this fit's clipped rows. When the kept rows are fewer than half, they are
    summed directly.
  - Residuals z are computed for all of this fit's rows. The robust σ, the 4σ clip, the convergence test and
    `n_iter = 3` are unchanged.
  - Clipping is per fit and may re-admit rows, exactly as now. No fit shares another's solution.
- **Exact fallbacks:**
  - If a fit would exceed `max_pix` (200,000 rows), the old `kernel_fit` runs, because its random subsample depends
    on the row list.
  - A fit with < 10 × 82 rows is skipped as now (keeps the whole-box kernel).
- **Numerics:** normal equations instead of SVD least squares. The solve is in float64 with column equilibration
  (`scipy.linalg.solve`, `assume_a="pos"`; on failure `lstsq` on the small 82 × 82 system).
  - Iterations, precision and sample are not reduced.
  - The condition number of A is measured on B2 and reported.
- **Kept as is:**
  - the default `refit` path (unchanged code);
  - the whole-box kernel;
  - stamp selection;
  - convolution.

## 3. Fixed tolerances (shared vs refit, same inputs)
1. **Kernels:**
   - per stamp, max|ΔK| ≤ 1e-6 × max|K|;
   - every clipped-row count identical, or else listed.
2. **Maps:** difference map and S, max|ΔS| ≤ 0.01 over finite pixels, per band.
3. **Detections:**
   - identical hit set, `near_bright_core` flags, morphology verdicts and review list;
   - every hit with 4.9 ≤ |S| ≤ 5.1 in any band, and every changed classification, listed individually.
4. **Injections:** every injected source's stage outcome identical (hit, sign, flag, morph, review), or listed
   individually.
- Any failure of 1–4 that is not explained by a listed borderline item fails the validation. The `refit` default
  then stays, and the new path is kept only behind its flag.

## 4. Test data (fixed now)
- **Small sample:**
  - (a) a synthetic unit test, shared vs refit on a random image with 12 compact stamps;
  - (b) B1 F150W: all stamp kernels compared at the kernel level (`--kernel-stamp-fit` both ways on the same
    arrays).
- **Dense region B2** (#288 box: 15.717097, −49.262488; F150W + F444W; the #288 final-chain settings):
  - one screen with no injections, run with `refit` and with `shared`;
  - one injection run each with the 60 unchanged / partial / vanish rows of `sample_b2_s40.csv` (`appear` rows
    dropped), through `--inject-pre`, so injections pass the matching step under test.

## 5. Timing (same hardware, same settings)
- **Concurrency:** the four B2 runs (refit screen, refit injection, shared screen, shared injection) run at once,
  one process each on the 4 cores, so old and new share the same load.
- **Reported:**
  - matching step (`resample_register_match`, plus the shared design build reported separately as its
    preparation cost);
  - total wall;
  - peak RSS (`getrusage` of the child);
  - with no warm cache. The cache directory is not used, so first-use costs are included.

## 6. Out of scope
Algorithm swap, GPU, new pipeline, new sky regions, SMACS, re-running history.
