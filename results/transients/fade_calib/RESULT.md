# Fading / vanishing calibration of `arc_diff`, measured speed-up, small batch search: result

Plan: `PLAN.md` §0–6, committed before any run (2225589); addenda 3–5 are dated and were each committed before the
runs they govern (35748f4, d9399f6, d0c1ea3). Owner instruction 2026-10-16; budget 23:32Z → 03:32Z. Every number is
`derived` unless labelled; injected sources are simulations; thresholds are ASSUMPTIONs.

## Outcome in short
- **The first held-out run (seeds 20 + 21) failed, and it was invalid.**
  - Unchanged false hits were 0.25–0.28.
  - `/code-review` then found that the injection geometry was wrong: epoch-2 sources were not registered like real
    sources.
- **After the correction (addendum 4), the fresh held-out seeds 30 + 31 pass the pre-registered rule (§4).** Pass
  rule ≤ 0.10:

  | Injection model | Vanish at 26.0 AB, review list with sign | Unchanged false hits |
  |---|---|---|
  | Empirical PSF (primary) | 14 / 14 | 1 / 40 = 0.025 |
  | Gaussian | 14 / 14 | 4 / 40 = 0.10, exactly at the limit |

- **The pass is post hoc in one respect:** the injection code was corrected after a held-out run had been seen.
  The pass rule, thresholds and samples are unchanged, and seeds 30 + 31 had never been run.
- **The margin is small.** The corrected dev rerun (seeds 10 + 11) gives 6 / 40 = 0.15 for both models. Over
  all four corrected seeds, false hits are 7 / 80 (empirical) and 10 / 80 (Gaussian).
- **Batch search:** see section "Batch search" below.

## Chain and tool changes
- **`--inject-pre`:**
  - Sources are added to both raw epochs before resampling, registration and matching.
  - Epoch-2 positions follow the chain's own injection-free global registration (`pair_shift`, `epoch2_pixel`), so
    an unchanged source lands where a real sky-fixed one would.
  - Local tile shifts (2 px cap) are not emulated.
- **The one chain fix: `--reg-tile-arcsec 6`.** The global shift is the median over 6″ tiles. Tiles < 90 %
  finite, or pinned at the ±3 px window, are dropped.
  - Before: the whole-box correlation of this pair was pinned at that edge, F150W (1.71, −3.00) px and F444W
    (0.53, −0.95) px.
  - Now: F150W (0.64, −0.52) px from 12 of 16 tiles, and F444W (0.11, −0.18) px from 16 tiles
    (`box1_v2/compare_286.json`).
  - The default stays 0, so historical runs reproduce.
- **Injection model:**
  - Primary: per-epoch empirical PSF. The 5 most concentrated isolated S/N ≥ 50 sources are chosen on epoch 1 and
    median-stacked per epoch.
  - The planned band-FWHM Gaussian is run beside it.
  - Choosing the empirical model as primary was post hoc (addendum 3).
- **`--cache-dir`:** the raw boxes (as header text, never a pickled WCS), the empirical PSFs and the injection-free
  registration. Every function argument is part of the key.

## Calibration (box 1, F150W + F444W, PEARLS o241 → VENUS o051)
- **Sample:**
  - 80 positions per seed on epoch-1 F150W pixels above the 90th percentile, 1.5″ apart (`sample_s*.csv`). Epoch 2
    played no part.
  - Classes: unchanged, appear, partial (m → m + 0.75), vanish; at 26.0 / 27.0 / 27.5 AB.
  - Each seed has 6–7 per cell; two seeds combined give n = 12–14.
- **Stages counted:** hit (|S| ≥ 5 in ≥ 2 bands) → not `near_bright_core` → morph `pass` → sign. The review list is
  the first three.
  - Counts per cell are in the `*_eval.json` files.
  - Intervals are binomial 68 % Wilson intervals.

**Held-out seeds 30 + 31 (`heldout2/`), review list with the correct sign:**

| Class / AB | Empirical PSF (primary) | Gaussian |
|---|---|---|
| vanish 26.0 | **14 / 14** [0.93, 1.00] | 14 / 14 [0.93, 1.00] |
| vanish 27.0 | 12 / 14 [0.74, 0.93] | 12 / 14 [0.74, 0.93] |
| vanish 27.5 | 10 / 12 [0.70, 0.91] | 11 / 12 [0.80, 0.97] |
| partial 26.0 | 14 / 14 [0.93, 1.00] | 13 / 14 [0.83, 0.97] |
| partial 27.0 | 8 / 12 [0.52, 0.79] | 8 / 12 [0.52, 0.79] |
| partial 27.5 | 4 / 14 [0.18, 0.42] | 6 / 14 [0.31, 0.56] |
| appear 26.0 / 27.0 / 27.5 | 9 / 12, 8 / 14, 12 / 14 | 9 / 12, 8 / 14, 13 / 14 |
| **unchanged: any hit** | **1 / 40** [0.01, 0.06] | **4 / 40** [0.06, 0.16] |
| unchanged: on the review list | 0 / 40 | 0 / 40 |

- **Missed fades** (`heldout2/*_eval.json` stage counts, empirical / Gaussian):
  - **Partial fades at 27.5 AB (a flux drop equal to a 28.25 AB source):** 9 / 14 and 6 / 14 are not detected at
    all. This is the detection limit, not a later stage.
  - **Partial fades at 27.0:** 3 / 12 and 1 / 12 are not detected; the rest of the losses are flags and morphology.
  - **Vanish at 27.0:** 1 / 14 is flagged `near_bright_core` in each model; one more is lost to morphology.
  - **Vanish at 27.5:** the empirical model misses 1 / 12 at detection and loses 1 to morphology; the Gaussian loses
    1 to morphology.
  - No review-list entry of a fade has the wrong sign.
- **Unchanged sources:** all their hits are removed before the review list, by morphology in every case.
- **Corrected dev rerun (`dev2/`, diagnosis only):**
  - unchanged 6 / 40 for both models, every one removed by morphology;
  - vanish 26.0 AB 12 / 14 for both.
- **Invalid runs, kept for the record (wrong epoch-2 geometry):**
  - `dev/`: the plan chain, empirical PSF alone, registration fix with the Gaussian, registration fix with the
    empirical PSF;
  - `heldout/`: seeds 20 + 21, unchanged false hits 11 / 40 and 10 / 40.
  - Their failure measured the injection artefact (a 0.9 px dipole for each unchanged source), not the chain.
- **Blind spots (not calibrated):**
  - extended fades;
  - positions within 0.8″ of bright cores (flagged, not reviewed);
  - blank-sky positions;
  - fainter than 27.5 AB;
  - the Poisson noise of injected sources;
  - PSF wings beyond 0.6″;
  - local registration residuals;
  - the density of 80 injections per run, which can shift the 99.9th-percentile core threshold that sets
    `near_bright_core`.

## Effect on #286 (box 1 re-screened, no injections; `box1_v2/`)
- **Screen:** 50 → 40 hits; 18 → 6 unflagged; review list 4 → 1.
- **The remaining entry** (h021) is at the #286 h024 position (0.0″). It was already vetted: a diffraction-spike
  residual with no forced-photometry change ≥ 2σ.
- **Nothing new; the #286 conclusion (nothing open) stands.**
- **#286's brightening efficiencies (0.85 / 0.78 / 0.54)** were measured after the injection point with the pinned
  registration. They are superseded for this box by the table above.
- **Other fields:** a scan of every stored `results/transients/**/summary.json` `shift_px` finds the ±3 px edge only
  in El Gordo runs. The others are:
  - MACS0416 |shift| ≤ 1.2 px;
  - Sunrise ≤ 0.3 px;
  - Abell 2744 grid ≤ 1.8 px;
  - SMACS GO-4043 ≤ 0.5 px.
  - SMACS ERO → VENUS has an unpinned (2.11, −1.82) px in F150W, not compared with a tile median. This is a
    pending check; that search recovered SN Atalanta.

## Speed-up (PLAN §5)
SPEEDUP_PLACEHOLDER

## Batch search (PLAN §6, addendum 5)
BATCH_PLACEHOLDER

## Regenerate
- **Samples:**
  `fade_calib.py make --uri <o241 F150W i2d> --ra <box> --dec <box> --seed S --out sample_*.csv`.
- **Boxes:**
  `fade_calib.py boxes --uri <4 i2d> --exclude 15.7300 -49.2530 --out boxes.json`.
- **`arc_diff.py`:** the command is in each `summary.json`. Maps go to the session scratchpad and are not stored.
- **Then:**
  - `arc_diff_morph.py --hits … --maps … --all-hits`;
  - `fade_calib.py eval --run DIR MORPH [--run …] --out eval.json`.
