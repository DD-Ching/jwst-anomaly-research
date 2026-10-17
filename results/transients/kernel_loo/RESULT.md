# Per-stamp leave-one-out kernel fits: result

Plan: `PLAN.md` (diagnosis, the one method, tolerances and test data), committed before any comparison. Addendum 1
records the post-review changes (deviations). Owner instruction 2026-10-17; budget 02:23Z → 05:23Z. Every number is
`derived`.

## Which repeated work was removed
- **Before:** for every compact-source stamp i, `kernel_match --kernel-stamps` called `kernel_fit` with
  `mask = stamp & (lab != i)`.
  - Each call rebuilt the whole weighted design of all *other* stamp pixels: 86,892 rows in all, about 86,800 per
    fit, × 82, in B2 F150W.
  - It then ran an SVD least-squares on that design once per clipping iteration (3 per fit), about 1.4 s per stamp.
- **Now (`kernel_fit_stamps`, `--kernel-stamp-fit shared`):** the design rows, weights, targets and their 82 × 82
  Gram matrix are built **once per band** (0.15–0.29 s).
  - Each fit and each iteration subtracts stamp i's own rows and its clipped rows. When those carry more than half
    of the Gram's trace, the Gram is summed from the kept rows instead.
  - The small system is solved with Cholesky and refined twice with the exact residual of the kept rows.
- **Unchanged:** leave-one-out, per-fit clipping (rows can be re-admitted), `n_iter = 3`, the 4σ clip and float64.
  Fits over 200,000 rows use the original routine.
- **Libraries:** numpy matrix products and the `scipy.linalg` Cholesky solve (`lstsq` fallback); no new dependency.

## Validation (tolerances fixed in PLAN §3 before comparing; test changes in addendum 1)
| Test | Data | Result | Tolerance |
|---|---|---|---|
| 1. Kernels, max\|ΔK\| / max\|K\|, and clipped χ² | unit tests (incl. 1e5 and 1e6 flux stars); **every** stamp of B1 (175 F150W, 43 F444W) and B2 (1,076 F150W, 33 F444W) (`check2/`) | ΔK ≤ 1.2e-13; Δχ² ≤ 1.1e-15 relative (same clipping) | 1e-6 |
| 2. Maps (B2 screen and injection run) | F150W, F444W | max\|ΔS\| ≤ 3.0e-8; difference maps ≤ 1.9e-9; 0 finite-pixel mismatches (`compare_b2.json`) | 0.01 |
| 3. Hits, flags, morphology, review list | B2 screen (2 hits); injection run (36 hits, 29 review entries) | identical | identical |
| 4. Injections (60 unchanged / partial / vanish rows through `--inject-pre`) | B2 | every stage outcome identical | identical |

- **Near-threshold items:**
  - Pixels with 4.9 ≤ |S| ≤ 5.1: 25 (F150W) and 8 (F444W) in the screen; 51 and 23 in the injection run.
  - None crosses |S| = 5 between versions.
  - No hit has 4.9 ≤ |S| ≤ 5.1 in any band.
- **No value or class changed beyond float rounding,** so **no extra misses.**
- Two repeats of each shared run are byte-identical.
- **Stamp counts:** the counts here (1,076 F150W stamps and 86,892 rows in B2, from the chain's registered images)
  replace the #288 session estimates (1,199 stamps, 97,883 pixels, counted before registration).
- **Scope:** these tests show that shared reproduces refit on these data. They do not re-measure the search's
  sensitivity, which stays as in #288.

## Speed (B2, same machine and settings, no cache; `timing_b2.json`)
| Run | Whole run | Matching step | Per-stamp kernels | Preparation | Peak memory |
|---|---|---|---|---|---|
| screen, refit → shared | 1,508 → 98–100 s (**×15**) | 1,492 → 83 s (×18) | 1,456 → 47 s (×31) | 0.28–0.29 s | 860 → 860 MB |
| injection run, refit → shared | 1,341 → 167–176 s (**×8**) | 1,242 → 76–77 s (×16) | 1,206 → 40 s (×30) | 0.15–0.19 s | 1,467 → 1,463–1,465 MB |

- **Load:**
  - Refit was timed with 4 runs started together. Its runs spent their last ~20 min on a half-loaded machine,
    which favours refit.
  - Shared (final code) was timed as 4 shared runs together.
  - The loads are similar, not identical. Refit has one run per variant; shared has two.
- **What still costs time** (shared, F150W):
  - the rest of matching, about 28 s (screen, 74.1 − 46.5 s): resampling, registration, `psf_match`, the
    whole-box kernel and the stamp median filters;
  - in injection runs, the empirical-PSF reads (`pre_injections`): 53 s in F150W under this load.
  - Neither is changed here.
- **Not measured:** B3.
- **Default:** `--kernel-stamp-fit` (and `kernel_match`) now default to `shared` (D-094); `refit` stays available.

## Not done (owner scope)
- algorithm swap, GPU, new pipeline;
- new regions or a batch;
- SMACS registration check;
- re-running history.
