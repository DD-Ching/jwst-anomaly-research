# Per-stamp leave-one-out kernel fits: result

Plan: `PLAN.md` (diagnosis, the one method, tolerances and test data), committed before any comparison. Owner
instruction 2026-10-17; budget 02:23Z → 05:23Z. Every number is `derived`.

## Which repeated work was removed
- **Before:** for every compact-source stamp i, `kernel_match --kernel-stamps` called `kernel_fit` with
  `mask = stamp & (lab != i)`. Each call rebuilt the whole weighted design of all *other* stamp pixels (about 87,000
  rows × 82 in B2 F150W) and ran an SVD least-squares on it once per clipping iteration (3 per fit), about 1.4 s
  per stamp.
- **Now (`kernel_fit_stamps`, `--kernel-stamp-fit shared`):** the design rows, weights, targets and their 82 × 82
  Gram matrix are built **once per band** (0.15–0.18 s).
  - Each fit subtracts stamp i's own rows and its own clipped rows, then solves the small system.
  - Leave-one-out, per-fit clipping (rows can be re-admitted), `n_iter = 3`, the 4σ clip and the float64 precision
    are unchanged. Fits over 200,000 rows fall back to the original routine.
- **Libraries:** numpy matrix products and the `scipy.linalg` Cholesky solve (with an `lstsq` fallback); no new
  dependency.

## Validation (tolerances fixed in PLAN §3 before comparing)
| Test | Data | Result | Tolerance |
|---|---|---|---|
| 1. Kernels, max\|ΔK\| / max\|K\| | synthetic (unit test); B1 F150W all 175 stamps, F444W all 43; B2 F150W first 40 of 1,076 (`check/`) | ≤ 2.7e-11 | 1e-6 |
| 2. Maps, max\|ΔS\| (B2 screen and injection run) | F150W, F444W | ≤ 4.8e-7; difference maps ≤ 3.7e-9; 0 finite-pixel mismatches (`compare_b2.json`) | 0.01 |
| 3. Hits, flags, morphology, review list | B2 screen (2 hits); B2 injection run (36 hits, 29 review entries) | identical; max\|Δ max S\| 3.5e-9 | identical |
| 4. Injections (60 unchanged / partial / vanish rows through `--inject-pre`) | B2 | every stage outcome identical | identical |

- **Near-threshold items:**
  - Pixels with 4.9 ≤ |S| ≤ 5.1: 25 (F150W) and 8 (F444W) in the screen; 51 and 23 in the injection run.
  - None crosses |S| = 5 between the two versions.
  - No hit has 4.9 ≤ max|S| ≤ 5.1.
- **No value or class changed beyond float rounding,** so **no extra misses.**
- The design's condition number is 320–1,520 (`check/*.json`), so normal equations lose at most about 1e-10
  relative.
- These tests show that shared reproduces refit. They do not re-measure the search's sensitivity, which stays as in
  #288.

## Speed (B2, same machine and settings, no cache; `timing_b2.json`)
| Run | Whole run | Matching step | Per-stamp kernels | Preparation | Peak memory |
|---|---|---|---|---|---|
| screen, refit → shared | 1,508 → 67 s (**×22**) | 1,492 → 53 s (×28) | 1,456 → 16 s (×92) | 0.18 s | 860 → 859 MB |
| injection run, refit → shared | 1,341 → 123 s (**×11**) | 1,242 → 50 s (×25) | 1,206 → 14 s (×86) | 0.15 s | 1,467 → 1,467 MB |

- **Load:** the four runs started together, one process each on the 4 cores.
  - The shared runs ended first, so the refit runs spent their last ~20 min on a half-loaded machine. This favours
    refit, so the ratios are not inflated by load.
  - One run per variant: no repeat timing, no spread.
- **What still costs time** (`timing_b2.json`):
  - The rest of F150W matching takes about 29 s in the shared screen: resampling, registration, `psf_match`, the
    whole-box kernel and the stamp median filters.
  - In injection runs, the empirical-PSF reads (`pre_injections`) take about 42 s in F150W.
  - Neither is changed here.
- **Not measured:** the B3 run.
- **Default:** `--kernel-stamp-fit` now defaults to `shared` (D-094). `refit` stays available, and job dicts without
  the key, from older callers, use `refit`.

## Not done (owner scope)
- algorithm swap, GPU, new pipeline;
- new regions or a batch;
- SMACS registration check;
- re-running history.
