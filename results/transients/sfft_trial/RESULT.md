# SFFT vs `arc_diff` (#291 baseline): result

Plan: `PLAN.md`, written before any run. Owner instruction 2026-10-17; budget 04:02Z → 06:02Z. A development
comparison on old data (B1, B2); every number is `derived`; injected sources are simulations. A small sample
supports only a trial choice.

## Verdict: keep `arc_diff` (#291) as the baseline
- **The PLAN §4 rule (`summary.json`) gives "keep arc_diff".** This holds after the two rounds of post-review
  corrections in PLAN addendum 1, and holds in every version.
- **Final ordinary SFFT:**
  - fails (a), unchanged false hits: 5/60 vs 1/60;
  - fails (b), fades: 83/120 vs 92/120;
  - passes (c), screen review entries: 0 vs 1;
  - passes (d), wall time: 418 s vs 580 s, indicative only (loads differ, addendum 1).
- **The JWST B-spline method of the paper (arXiv:2309.09143) was not tested.** In sfft 1.7.3 it is GPU-only (CuPy);
  the Numpy path refuses: "Numpy backend is not supported for B-Spline SFFT in current version". No GPU is
  available, and a port or purchase is out of scope.
- **Not tested either:** SFFT with its own source masking (the Easy/Auto packets need the SExtractor binary).
- **What was tested:** ordinary SFFT (`Customized_Packet`, polynomial kernel, Numpy CPU) without source masking. It
  replaced only the kernel-matching step of an otherwise identical chain.
- **Conclusion:** the difference is modest on this small sample (2 boxes; 180 injected changes, 60 each of
  appear, partial and vanish, plus 60 unchanged controls), and it shows no gain
  worth a switch. Keep the #291 tool; no second tool port.

## Integration (≤ 30 min gate)
- **Licence and package:** MIT; PyPI `sfft` 1.7.3.
- **Dependencies added:** sep, pyfftw, fastremap. No version changes; session install only.
- **CPU:** the polynomial mode runs on CPU: 22 s on a synthetic 770 × 770 image with 1 thread, including numba
  compilation.
- **Blockers for the paper method:**
  - GPU (double-precision CuPy) only;
  - automatic masking needs SExtractor;
  - the paper's PSF-rotation handling sits in the B-spline path.

## Comparison (B1 seeds 30 + 31, B2 seed 40; same physical pre-matching injections; final SFFT in `runs/`)
Review list with the correct sign, B1 / B2 (`runs/*_eval.json`):

| Class | `arc_diff` B1 | SFFT B1 | `arc_diff` B2 | SFFT B2 |
|---|---|---|---|---|
| vanish 26.0 / 27.0 / 27.5 | 14/14, 12/14, 10/12 | 14/14, 12/14, 10/12 | 5/7, 6/7, 4/6 | 5/7, 6/7, 4/6 |
| partial 26.0 / 27.0 / 27.5 | 14/14, 8/12, 4/14 | 12/14, 5/12, 1/14 | 7/7, 5/6, 3/7 | 7/7, 5/6, 2/7 |
| appear 26.0 / 27.0 / 27.5 | 9/12, 8/14, 12/14 | 9/12, 6/14, 11/14 | 6/6, 5/7, 5/7 | 6/6, 5/7, 5/7 |
| unchanged: any hit | 1/40 | 5/40 | 0/20 | 0/20 |

- **Totals (`summary.json`):**
  - fades 92/120 (`arc_diff`) vs 83/120 (SFFT): vanish 51/60 in both, partial 41/60 vs 32/60;
  - the equal vanish totals are not the same sources: one vanishing source reaches the review list only with
    `arc_diff` (b1_s30 pre019), one only with SFFT (b1_s31 pre007), and one is hit only by SFFT (b1_s30 pre011)
    (`summary.json`);
  - brightenings 45/60 vs 42/60;
  - unchanged false hits 1/60 vs 5/60, all of SFFT's in B1;
  - screen review entries 1 vs 0.
- **Per source:** 26 injected sources differ in hit or review between the tools (listed in `summary.json`).
  - 15 reach the review list with the right sign only with `arc_diff`: 11 partial, 3 appear, 1 vanish.
  - 5 only with SFFT: 2 unchanged (false), 2 partial, 1 vanish.
  - 6 differ only in being hit.
- **Signal** (median S at injected positions, B1 seeds 30 / 31; `diagnostics_median_s.json`):
  - partial 26 AB, F150W: SFFT −27.1 / −21.2 vs `arc_diff` −32.6 / −25.1;
  - unchanged 26 AB, F150W: residual +29.9 / +31.6 vs +19.8 / +19.8.
  - The weaker partial-fade signal and the larger residual of unchanged sources match where SFFT loses recoveries and
    gains false hits. This is not checked per source.
  - Vanishing sources are recovered identically.
- **Earlier SFFT versions (superseded):**
  - First runs (`runs_v1/`, non-union NaN mask, unpropagated variance): fades 58/120, brightenings 29/60, false hits
    11/60 (`runs_v1/summary_v1_recounted.json`). B1 seed-30 vanishing 26 AB sources had F150W S +35.4
    (`diagnostics_median_s.json`).
  - Diamond mask (`runs_v2/summary_v2.json`): 83/120, 42/60, 7/60.
  - The first `summary.json` also miscounted mixed signs for both tools (`runs_v1/summary_v1_miscounted.json`).
- **Residual statistics** (screens, robust std / max|S| of S; `runs/*_regions.json`):
  - B1 F150W:
    - bright cores: SFFT 2.77 / 268 vs `arc_diff` 2.95 / 202;
    - stamps 1.63 vs 1.57;
    - spike annulus 1.11 / 16.9 vs 1.19 / 18.4.
  - B1 F444W:
    - cores 2.01 vs 1.85;
    - stamps 1.83 vs 1.34;
    - spike annulus 1.14 / 17.9 vs 1.18 / 11.3.
  - B2 F444W: cores 1.06 vs 1.28; stamps 0.97 vs 1.01.
  - Screens: B1 47 hits (13 unflagged) vs 40 (6); B2 1 vs 2.
  - Stamp regions use epoch-1 stamps only (deviation, addendum 1).
  - Neither tool is cleaner everywhere. Recovery and false hits, above, decide.
- **Noise handling:**
  - `arc_diff` propagates variance through its per-stamp and whole-box kernels.
  - SFFT's variance uses one effective whole-box kernel fitted to its own matched image. That is biased low where
    SFFT's local kernel departs from the mean: a remaining limitation.
  - Both rescale S by its robust std, clamped at 1.
- **Cost:**
  - Total wall time for the 5 runs: 580 s (`arc_diff`) vs 418 s (SFFT), ×1.39, indicative only (loads differ,
    addendum 1).
  - Peak memory 0.86–1.50 GB vs 2.16–2.23 GB, a ratio of 1.49–2.51 (`summary.json`).
  - SFFT's numba compilation is included in each run.

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
  4 min per core for one box. Four boxes in parallel on 4 cores give a throughput of about 1 box per minute; one
  box still takes about 4 min. Range reads come on top.
- **Review burden:**
  - B1–B3 gave 1, 0, 0 review entries per box (`fade_calib/` and `kernel_loo/`). MACS0416 tiling gave about 0.3 per
    tile (TASKS: 9 open passes in 28 tiles).
  - Expect about 0.3–1 entry per box. Each takes about 10–20 min to vet (an estimate, not stored).

## Corrections (2026-10-17, owner; wording only, no re-run)
- Vanish recovery was printed as 51/51. Over all three magnitudes and both boxes it is 51/60 for each tool; per
  magnitude 19/21, 18/21 and 14/18 at 26.0 / 27.0 / 27.5 AB (`runs/*_eval.json`).
- The sample is 180 injected changes (60 appear, 60 partial, 60 vanish) plus 60 unchanged controls, not "180 injected
  sources".
- Equal vanish totals do not mean the same sources were recovered (see the per-source list).
- "About 1 min per box" was a parallel-throughput estimate; one box takes about 4 min on one core.
