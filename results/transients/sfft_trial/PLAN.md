# SFFT vs `arc_diff` (#291 baseline): plan (fixed before any comparison run)

Owner instruction 2026-10-17. Budget 04:02Z → 06:02Z. Development comparison on old data; a small sample supports
only a trial choice, not "generally better".

## 1. Integration check (done 04:02–04:15Z)
- **Paper:** arXiv:2309.09143 (Hu & Wang 2024, AJ 167, 5) is the JWST method: B-spline spatially varying kernels,
  separate photometric scaling and Tikhonov regularization. Code: `sfft` (MIT, PyPI 1.7.3, github.com/thomasvrussell/sfft).
- **Blocker for the paper method:** `sfft.BSplineSFFT` (`BSpline_Packet`, `SingleSFFTConfigure.SSC`,
  `ElementalSFFTSubtract.ESS`) supports only the CuPy GPU backend (double-precision GPU). Its Numpy path prints
  "Numpy backend is not supported for B-Spline SFFT in current version" and returns `None`.
  - No GPU is available, and a GPU port or purchase is out of scope.
  - **The paper method is not tested here.**
- **Usable without changes:** ordinary SFFT, `Customized_Packet.CP` with polynomial spatial variation and the
  Numpy CPU backend.
  - Its automatic masking packets (Easy/Auto) call the external SExtractor binary, which is not installed. They
    are not used.
  - Test: 22 s on a synthetic 770 × 770 image, 1 thread, including numba compilation.
  - Install: `pip install sfft==1.7.3` adds sfft, sep, pyfftw and fastremap, with no version changes. Session
    only; not added to `pyproject.toml`.
- **Everything below is ordinary SFFT, not the JWST B-spline method.**

## 2. What is compared
- **The only step swapped is PSF/kernel matching.** `arc_diff` (#291: pixel kernel plus per-stamp leave-one-out
  kernels, `--kernel-px 4 --kernel-stamps`, shared) against `sfft.Customized_Packet.CP`.
- **Shared with both:** reading, resampling, tile-median registration, local shifts, the `psf_match` blur sign
  (which side is convolved), significance, detection, flags, morphology and the review rule. They are identical:
  `scripts/sfft_trial.py` replaces `arc_diff.kernel_match` and nothing else.
- **SFFT settings (fixed, no tuning):**
  - REF = the side `arc_diff` convolves, SCI = the other side, `ForceConv='REF'`;
  - `GKerHW=4` (9 × 9, as `--kernel-px 4`), `KerPolyOrder=2`, `BGPolyOrder=0`, `ConstPhotRatio=True`;
  - Numpy backend, 1 thread.
  - Masked images = the images with NaN set to 0, and no source masking (SFFT's automatic masking needs
    SExtractor). This is SFFT without its usual transient and variable masking: a stated handicap.
- **Noise (stated difference):** `arc_diff` propagates variance through its kernels (var·K² + var). For SFFT the
  variance is var_ref + var_sci, without kernel propagation. Both then rescale S by the robust std of S over the
  box (`significance`), so the false-alarm threshold is the same empirical 5σ in both.
- **Same physical injections, both through `--inject-pre`:**
  - B1: held-out samples `fade_calib/sample_s30.csv` and `sample_s31.csv`, with the empirical PSF;
  - B2: `fade_calib/sample_b2_s40.csv`;
  - classes unchanged / appear (brighten) / partial (fade) / vanish.
- **Data:** B1 (15.7300, −49.2530) and B2 (15.717097, −49.262488); F150W + F444W; #288 / #291 settings.

## 3. Metrics (computed by fixed code)
1. **Cost:** total wall time per run, including reads, preparation and SFFT compilation; peak RSS; matching-step
   time. Both run under the same load (4 processes at once, 1 per core).
2. **Injections** (`fade_calib.py eval` stages):
   - unchanged false hits (any hit);
   - per class × magnitude review-with-sign recovery;
   - per-source outcome, where every source recovered by one tool and not the other is listed.
3. **Screens without injections:**
   - hits, unflagged hits, the review list;
   - robust std and max|S| of S in three regions:
     - bright cores: within 0.8″ of epoch-1 pixels above the 99.9th percentile, the `near_bright_core` definition;
     - compact-source stamps: `compact_stamps` of both epochs, the crowded proxy;
     - the B1 spike annulus: 1.5–4.0″ around the bright star (15.727222, −49.253939).

## 4. Decision rule (fixed now)
- **"Try SFFT in the next batch"** only if all of these hold on B1 + B2 together:
  - (a) unchanged false hits ≤ `arc_diff`'s;
  - (b) review-with-sign recoveries of fades (partial + vanish) ≥ `arc_diff`'s, and brightenings not lower by more
    than 2;
  - (c) screen review-list entries ≤ `arc_diff`'s;
  - (d) wall time ≤ 2 × `arc_diff`'s.
- **"Keep `arc_diff`"** if (a) or (b) is worse, or if SFFT fails to run.
- **Otherwise "evidence insufficient".**
- No parameter changes after the first run. A failure is recorded, not tuned.
