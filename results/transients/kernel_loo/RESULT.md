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
  Gram matrix are built **once per band** (≤ 0.38 s per run).
  - Each fit and each iteration subtracts stamp i's own rows and its clipped rows. When those carry more than half
    of the Gram's trace, the Gram is summed from the kept rows instead.
  - The small system is solved with Cholesky and refined with the exact residual of the kept rows until the step
    is ≤ 1e-12 of the solution. If that fails, the Gram is re-summed directly.
- **Unchanged:** leave-one-out, per-fit clipping (rows can be re-admitted), `n_iter = 3`, the 4σ clip and float64.
  Fits over 200,000 rows use the original routine.
- **Libraries:** numpy matrix products and the `scipy.linalg` Cholesky solve (`lstsq` fallback); no new dependency.

## Validation (tolerances fixed in PLAN §3 before comparing; test changes in addenda 1–2)
| Test | Data | Result | Tolerance |
|---|---|---|---|
| 1. Kernels, max\|ΔK\| / max\|K\|, and clipped χ² | unit tests (incl. 1e5 and 1e6 flux stars and an outlier inside a stamp; `max_pix` fallbacks); **every** stamp of B1 (175 F150W, 43 F444W) and B2 (1,076 F150W, 33 F444W) (`check3/`) | ΔK ≤ 1.2e-13; Δχ² ≤ 1.1e-15 relative (same clipping) | 1e-6 |
| 2. Maps (B2 screen and injection run) | F150W, F444W | max\|ΔS\| ≤ 3.0e-8; difference maps ≤ 1.9e-9; 0 finite-pixel mismatches (`compare_b2.json`) | 0.01 |
| 3. Hits, flags, morphology, review list | B2 screen (2 hits); injection run (36 hits, 29 review entries) | identical | identical |
| 4. Injections (60 unchanged / partial / vanish rows through `--inject-pre`) | B2 | every stage outcome identical | identical |

- **Near-threshold items:**
  - Pixels with 4.9 ≤ |S| ≤ 5.1: 25 (F150W) and 8 (F444W) in the screen; 51 and 23 in the injection run.
  - None crosses |S| = 5 between versions.
  - No hit has 4.9 ≤ |S| ≤ 5.1 in any band, in either run's maps.
- **No value or class changed beyond float rounding,** so **no extra misses.**
- Two repeats of each shared run are identical: ΔS = 0, same hits and outcomes (`compare_b2_repeat.json`).
- No fit needed the refinement fallback (`kernel_stats.stamp_unconverged` = 0).
- **Stamp counts:** the counts here (1,076 F150W stamps and 86,892 rows in B2, from the chain's registered images)
  replace the #288 session estimates (1,199 stamps, 97,883 pixels, counted before registration).
- **Scope:** these tests show that shared reproduces refit on these data. They do not re-measure the search's
  sensitivity, which stays as in #288.
- **The first version failed:** plain Gram downdating, without the cancellation guard or refinement, missed the
  tolerance for a dominant bright stamp in a synthetic scene (1.1e-6 at flux 1e5 and 1.2e-4 at 1e6; review
  measurement, not stored). It was replaced (PLAN addendum 1).

## Speed (B2, same machine and settings, no cache; `timing_b2.json`)
| Run | Whole run | Matching step | Per-stamp kernels | Preparation | Peak memory |
|---|---|---|---|---|---|
| screen, refit → shared | 1,508 → 106 s (**×14**) | 1,492 → 89–90 s (×17) | 1,456 → 52–53 s (×28) | 0.38 s | 860 → 860–861 MB |
| injection run, refit → shared | 1,341 → 157–176 s (**×8–9**) | 1,242 → 77–78 s (×16) | 1,206 → 40–42 s (×29–30) | 0.22–0.26 s | 1,467 → 1,438–1,465 MB |

- **Load (deviation from PLAN §5, addendum 2):**
  - Refit and the final shared code were timed in different batches, each of 4 processes on 4 cores.
  - The refit runs spent their last ~20 min on a half-loaded machine, which favours refit.
  - Refit has one run per variant; shared has two.
- **Preparation** (building the shared design per band): 0.22–0.38 s per run. In the kernel checks it was
  0.006–0.367 s per band (`check3/`).
- **What still costs time** (shared, F150W):
  - the rest of matching: 27–29 s (resampling, registration, `psf_match`, the whole-box kernel, the stamp median
    filters);
  - in injection runs, the empirical-PSF reads (`pre_injections`): 42–59 s.
  - Neither is changed here.
- **Not measured:** B3.
- **Default:** `--kernel-stamp-fit` (and `kernel_match`) default to `shared` (`STAMP_FIT_DEFAULT`, D-094). `refit`
  stays available; any other value raises.

## Not done (owner scope)
- algorithm swap, GPU, new pipeline;
- new regions or a batch;
- SMACS registration check;
- re-running history.
