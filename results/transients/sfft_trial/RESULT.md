# SFFT vs `arc_diff` (#291 baseline): result

Plan: `PLAN.md`, written before any run. Owner instruction 2026-10-17; budget 04:02Z → 06:02Z. A development
comparison on old data (B1, B2); every number is `derived`; injected sources are simulations. A small sample
supports only a trial choice.

## Verdict: keep `arc_diff` (#291) as the baseline
- **The PLAN §4 rule (`summary.json`) gives "keep arc_diff".** This holds after the post-review corrections of PLAN
  addendum 1.
- **Ordinary SFFT:**
  - fails (a), unchanged false hits: 7/60 vs 1/60;
  - fails (b), fades: 83/120 vs 92/120;
  - passes (c), screen review entries: 1 vs 1;
  - passes (d), wall time: 441 s vs 580 s.
- **The JWST B-spline method of the paper (arXiv:2309.09143) was not tested.** In sfft 1.7.3 it is GPU-only (CuPy);
  the Numpy path refuses: "Numpy backend is not supported for B-Spline SFFT in current version". No GPU is
  available, and a port or purchase is out of scope.
- **Not tested either:** SFFT with its own source masking (the Easy/Auto packets need the SExtractor binary).
- **What was tested:** ordinary SFFT (`Customized_Packet`, polynomial kernel, Numpy CPU) without source masking. It
  replaced only the kernel-matching step of an otherwise identical chain.
- **Conclusion:** the difference is modest on this small sample (2 boxes, 180 injected sources). It does not show
  that SFFT is generally worse; it shows no gain worth a switch. Keep the #291 tool; no second tool port.

## Integration (≤ 30 min gate)
- **Licence and package:** MIT; PyPI `sfft` 1.7.3.
- **Dependencies added:** sep, pyfftw, fastremap. No version changes; session install only.
- **CPU:** the polynomial mode runs on CPU: 22 s on a synthetic 770 × 770 image with 1 thread, including numba
  compilation.
- **Blockers for the paper method:**
  - GPU (double-precision CuPy) only;
  - automatic masking needs SExtractor;
  - the paper's PSF-rotation handling sits in the B-spline path.

## Comparison (B1 seeds 30 + 31, B2 seed 40; same physical pre-matching injections; corrected SFFT in `runs/`, first SFFT runs in `runs_v1/`)
Review list with the correct sign, B1 / B2 (`runs/*_eval.json`):

| Class | `arc_diff` B1 | SFFT B1 | `arc_diff` B2 | SFFT B2 |
|---|---|---|---|---|
| vanish 26.0 / 27.0 / 27.5 | 14/14, 12/14, 10/12 | 14/14, 12/14, 10/12 | 5/7, 6/7, 4/6 | 5/7, 6/7, 4/6 |
| partial 26.0 / 27.0 / 27.5 | 14/14, 8/12, 4/14 | 11/14, 6/12, 1/14 | 7/7, 5/6, 3/7 | 7/7, 5/6, 2/7 |
| appear 26.0 / 27.0 / 27.5 | 9/12, 8/14, 12/14 | 9/12, 6/14, 11/14 | 6/6, 5/7, 5/7 | 6/6, 5/7, 5/7 |
| unchanged: any hit | 1/40 | 7/40 | 0/20 | 0/20 |

- **Totals (`summary.json`):**
  - fades 92/120 (`arc_diff`) vs 83/120 (SFFT);
  - brightenings 45/60 vs 42/60;
  - unchanged false hits 1/60 vs 7/60, all of SFFT's in B1;
  - screen review entries 1 vs 1.
- **Per source:** 27 injected sources differ in hit or review between the tools (listed in `summary.json`).
  - 15 reach the review list with the right sign only with `arc_diff`: 11 partial, 3 appear, 1 vanish.
  - 5 only with SFFT: 2 unchanged (false), 2 partial, 1 vanish.
  - 7 differ only in being hit.
- **Signal:** vanishing 26 AB sources in B1 have median S of −32.7 (F150W) and −33.2 (F444W) with SFFT, against −46.3
  and −39.0 with `arc_diff` (`diagnostics.json`). SFFT's S is lower. This is consistent with more faint partial fades
  falling below threshold, but that is not checked per source.
- **First SFFT runs (superseded, `runs_v1/`):**
  - They had a non-union NaN mask and unpropagated variance.
  - Results: fades 58/120, brightenings 29/60, false hits 11/60 (`runs_v1/summary_v1_recounted.json`).
  - In B1 F150W the vanishing sources had S of +35 (`diagnostics.json`).
  - The first `summary.json` also miscounted mixed signs for both tools (`runs_v1/summary_v1_miscounted.json`).
- **Residual statistics** (screens, robust std / max|S| of S; `runs/*_regions.json`):
  - B1 F150W:
    - SFFT 2.68 / 254 on bright cores vs `arc_diff` 2.95 / 202;
    - stamps 1.68 vs 1.57;
    - spike annulus 1.11 / 16.0 vs 1.19 / 18.4.
  - B1 F444W:
    - cores 1.89 vs 1.85;
    - stamps 1.73 vs 1.34;
    - spike annulus 1.08 / 17.1 vs 1.18 / 11.3.
  - B2 F444W: cores 1.06 vs 1.28; stamps 0.99 vs 1.01.
  - Screens: B1 47 hits (16 unflagged) vs 40 (6); B2 1 vs 2.
  - Stamp regions use epoch-1 stamps only (deviation, addendum 1).
  - Neither tool is cleaner everywhere. Recovery and false hits, above, decide.
- **Noise handling:**
  - Both now propagate variance through a kernel; for SFFT this is the effective whole-box kernel, since SFFT's own
    kernel varies over the box.
  - Both rescale S by its robust std, clamped at 1.
  - The SFFT B1 F444W screen has robust std 0.93 after rescale, so its bar there is slightly above an empirical 5σ.
- **Cost:**
  - Total wall time for the 5 runs: 580 s (`arc_diff`) vs 441 s (SFFT), ×1.31.
  - Peak memory: 0.86–1.50 GB vs 2.16–2.23 GB.
  - SFFT's numba compilation is included in each run.
  - Loads differ slightly (addendum 1).

## Batch-search denominator (for the next budgeted batch; `venus_pairs.json`)
- **Pairs:**
  - VENUS (6882) has 56 NIRCam F150W level-3 targets.
  - 18 have at least one other public NIRCam F150W level-3 program within ±0.05° (MAST query 2026-10-17), about 31
    candidate epoch pairs.
  - This is an **upper bound**: F444W, footprint overlap and pipeline versions are not checked.
- **Boxes** (`diagnostics.json`): from the catalogued image lists, non-overlapping 24″ boxes with ≥ 5 images:
  - El Gordo 5 (3 done);
  - MACS0416 12;
  - SMACS 4 (1 done).
  - Other fields need their image lists.
- **Cost per box** (#291 chain, 1 core; `runs/arc_*/rss.json`): screen 1–2 min plus one injection run 2–3 min, so about
  4 min per core, or about 1 min per box on 4 cores, plus range reads.
- **Review burden:**
  - B1–B3 gave 1, 0, 0 review entries per box (`fade_calib/` and `kernel_loo/`). MACS0416 tiling gave about 0.3 per
    tile (TASKS: 9 open passes in 28 tiles).
  - Expect about 0.3–1 entry per box. Each takes about 10–20 min to vet (an estimate, not stored).
