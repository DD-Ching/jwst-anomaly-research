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

## Addendum 1 (2026-10-17T04:20Z, after `/code-review` of the first results; corrections, not tuning)
- **The review found three adapter or accounting faults.** They penalised SFFT, or both tools, and none is an SFFT
  parameter:
  1. **Masking:** each masked image was zero-filled only at its own NaNs. Now both are zeroed on the union NaN
     mask, grown by `GKerHW`, as SFFT's own packets do.
  2. **Variance and kernel:**
     - Before, the variance was var_ref + var_sci; with the clamp of the S rescale at 1, this raised SFFT's bar.
     - Now an effective 9 × 9 kernel is fitted to SFFT's own matched image (`kernel_fit`; its constant term absorbs
       SFFT's background).
     - That kernel propagates var_ref, gives the kernel sum and convolves post-matching injections.
  3. **Accounting:** `summarize` judged "mixed" signs without per-band S. It now uses `fade_calib.stage_rows`.
- **What was re-run:** the five SFFT runs, once. The first results stay in `runs_v1/`.
- **Unchanged:** the `arc_diff` runs, the settings, the decision rule.
- **Load difference:** the corrected SFFT runs ran as 4 at once plus 1 alone, where the `arc_diff` runs were 2 + 2.
- **Deviation:** the compact-stamp region in `regions` uses epoch-1 stamps only, because the maps store no epoch-2
  image. PLAN §3.3 said both epochs.
- **Second review (04:31Z):** the union mask had been grown by a diamond, not the kernel footprint. It is now the
  full 9 × 9 footprint, and the five SFFT runs were re-run once more.
  - Results by version (`summary*.json`), `arc_diff` unchanged throughout:
    - first SFFT runs: `runs_v1/`;
    - diamond mask: `runs_v2/`;
    - final: `runs/`.
  - Also stated as a limitation: SFFT's variance uses one mean kernel, so it is biased low where SFFT's local kernel
    departs from the mean.
- **Load (deviation from PLAN §3.1):**
  - `arc_diff`: batches of 2 + 2 (screens, B1 seeds) and 1 (B2 seed 40), each alongside the same number of first
    SFFT runs.
  - Final SFFT: 4 at once, then 1 alone.
  - The wall-time ratio is therefore indicative only.
