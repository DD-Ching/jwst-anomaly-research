# SFFT vs `arc_diff` (#291 baseline): result

Plan: `PLAN.md`, written before any run. Owner instruction 2026-10-17; budget 04:02Z → 06:02Z. A development
comparison on old data (B1, B2); every number is `derived`; injected sources are simulations. A small sample
supports only a trial choice.

## Verdict: keep `arc_diff` (#291) as the baseline
- **The PLAN §4 rule (`summary.json`) gives "keep arc_diff".**
  - Ordinary SFFT fails (a) false hits, (b) fades and (c) screen review entries.
  - It passes (d) wall time.
- **The JWST B-spline method of the paper (arXiv:2309.09143) was not tested.** In sfft 1.7.3 it is GPU-only
  (CuPy); the Numpy path refuses: "Numpy backend is not supported for B-Spline SFFT in current version". No GPU is
  available, and a port or purchase is out of scope.
- **Not tested either:** SFFT with its own source masking. That is the Easy/Auto packets, which call the external
  SExtractor binary, which is not installed.
- **What was tested:** ordinary SFFT (`Customized_Packet`, polynomial kernel, Numpy CPU) **without** source masking,
  a stated handicap. It replaced only the kernel-matching step of an otherwise identical chain.
- **Conclusion:** keep the #291 tool. No second tool port is started.

## Integration (≤ 30 min gate)
- **Licence and package:** MIT; PyPI `sfft` 1.7.3.
- **Dependencies added:** sep, pyfftw, fastremap. No version changes; session install only.
- **CPU:** the polynomial mode runs on CPU: 22 s on a synthetic 770 × 770 image with 1 thread, including numba
  compilation.
- **Blockers for the paper method:**
  - GPU (double-precision CuPy) only;
  - automatic masking needs SExtractor;
  - the paper's PSF-rotation handling sits in the B-spline path.

## Comparison (B1 seeds 30 + 31, B2 seed 40; same physical pre-matching injections; `runs/`)
Review list with the correct sign, B1 / B2 (`runs/*_eval.json`):

| Class | `arc_diff` B1 | SFFT B1 | `arc_diff` B2 | SFFT B2 |
|---|---|---|---|---|
| vanish 26.0 / 27.0 / 27.5 | 14/14, 12/14, 10/12 | 8/14, 5/14, 3/12 | 5/7, 6/7, 4/6 | 5/7, 6/7, 5/6 |
| partial 26.0 / 27.0 / 27.5 | 14/14, 8/12, 4/14 | 7/14, 3/12, 2/14 | 7/7, 5/6, 3/7 | 7/7, 5/6, 2/7 |
| appear 26.0 / 27.0 / 27.5 | 9/12, 8/14, 12/14 | 9/12, 3/14, 1/14 | 6/6, 5/7, 5/7 | 6/6, 5/7, 5/7 |
| unchanged: any hit | 1/40 | 8/40 | 0/20 | 3/20 |

- **Totals (`summary.json`):**
  - fades 85/120 (`arc_diff`) vs 30/120 (SFFT);
  - brightenings 39/60 vs 27/60;
  - unchanged false hits 1/60 vs 11/60;
  - screen review entries 1 vs 2.
- **Per source:** 71 injected sources differ in hit or review between the tools.
  - 47 are on the review list only with `arc_diff`: 16 appear, 16 vanish, 15 partial.
  - 7 are on it only with SFFT: 3 unchanged (false review entries), 3 partial, 1 vanish.
  - The remaining 17 differ only in whether they were hit, with neither on the review list.
  - Listed in `summary.json`.
- **Where SFFT fails:** in B1 F150W.
  - Its fitted photometric ratio is 0.895, against `arc_diff`'s 0.991.
  - The raw S scatter is 6.6 × the noise model, against 1.16 for `arc_diff` (`noise_scale` in `runs/*/summary.json`).
  - Every source leaves a positive residual: the median S of vanishing 26 AB sources is +35, against −46 for
    `arc_diff` (session check, not stored).
  - F444W behaves like `arc_diff`: vanish −33.5 vs −39.0; noise scale 1.04.
  - **Likely cause, not verified:** the unmasked fit includes the box's saturated star and bright cores, which SFFT
    normally masks. In B2 (whether it holds a saturated star was not checked), SFFT is close to `arc_diff` on
    vanishing sources (16/20 vs 15/20) and partial fades (14/20 vs 15/20) but has more false hits (3/20 vs 0/20).
- **Residual statistics** (screens, robust std / max|S| of S; `runs/*_regions.json`):
  - B1 F150W:
    - SFFT 1.72 / 139 on bright cores against `arc_diff` 2.95 / 202;
    - spike annulus 0.61 / 7.6 against 1.19 / 18.4.
  - B1 F444W: SFFT is worse on cores (2.28 vs 1.85) and stamps (2.60 vs 1.34).
  - B2 F444W: SFFT is worse on cores (1.68 vs 1.28) and stamps (1.94 vs 1.01).
  - Lower residuals on bright cores in one band do not mean better detection. Recovery and false hits, above, are
    the deciding measures.
- **Noise handling (stated difference):**
  - `arc_diff` propagates variance through its kernels; the SFFT adapter uses var_ref + var_sci.
  - Both rescale S by its robust std, so the 5σ threshold is the same empirical one.
- **Cost (same load: each batch held equal numbers of `arc_diff` and SFFT runs, 1 process per core):**
  - total wall time for the 5 runs: 580 s (`arc_diff`) vs 427 s (SFFT), ×1.36;
  - peak memory 0.86–1.50 GB vs 2.22–2.27 GB;
  - SFFT's numba compilation is included in each run.

## Batch-search denominator (for the next budgeted batch; `venus_pairs.json`)
- **Pairs:**
  - VENUS (6882) has 56 NIRCam F150W level-3 targets.
  - 18 have at least one other public NIRCam F150W level-3 program within ±0.05° (MAST query 2026-10-17), about 31
    candidate epoch pairs.
  - This is an **upper bound**: F444W, footprint overlap and pipeline versions are not checked.
- **Boxes:** from the catalogued image lists, non-overlapping 24″ boxes with ≥ 5 images:
  - El Gordo 5 (3 done);
  - MACS0416 12;
  - SMACS 4 (1 done).
  - Other fields need their image lists.
- **Cost per box (#291 chain, 1 core, from these runs):** screen 1–2 min plus one injection run 2–3 min, so about
  4 min per core, or about 1 min per box on 4 cores, plus range reads.
- **Review burden:**
  - B1–B3 gave 1, 0, 0 review entries per box. MACS0416 tiling gave about 0.3 per tile.
  - Expect about 0.3–1 entry per box; each takes about 10–20 min to vet.
