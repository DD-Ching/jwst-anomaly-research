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
- **What the numbers are:**
  - 14 / 14 is a small-sample result for this box and these positions (68 % interval 0.93–1.00), not a general
    reliability.
  - The unchanged-source false-hit fraction counts residuals of injected steady point sources on galaxy and arc
    pixels. It is **not** the false-alarm rate of a real-sky search, which depends on real residuals (cores, spikes,
    edges) and is handled by the morphology and vetting stages.
- **Known limitation (final review, not corrected, no rerun per the owner's 2026-10-17 instruction):**
  - With injections present, the chain's own registration moves away from the injection-free shift used to place
    epoch-2 sources, by 0.02–0.15 px in F150W and ≤ 0.03 px in F444W (`inject_shift_residual.json`). The
    addendum-4 error was 0.9 px.
  - Injected sources therefore keep a small misregistration. This can add false hits for unchanged sources and
    lower fade recovery slightly.
  - The injection-free screens of B1–B3 are not affected.
- **The margin is small.** The corrected dev rerun (seeds 10 + 11) gives 6 / 40 = 0.15 for both models. Over
  all four corrected seeds, false hits are 7 / 80 (empirical) and 10 / 80 (Gaussian).
- **Batch search (B1–B3):** 0 review entries in B2 and B3. B1's single entry is the already-vetted #286 h024
  position. Nothing is pending.

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
  - SMACS ERO → VENUS had an unpinned (2.11, −1.82) px in F150W; the tile median is (0.20, −0.40) px, and the
    re-screen keeps the outcome (Atalanta recovered, 0 `vet`; `../smacs_arcdiff/o001_o057_tile/RESULT.md`).

## Speed-up (PLAN §5)
- **What is cached (`--cache-dir`):** raw-data products only: the raw boxes, the empirical PSFs and the
  injection-free registration. Matching depends on the injected sources and is never cached.
- **Fair comparison** (`speedup.json`, final code):
  - Box 1, seed 30, empirical and Gaussian injections, the same sample and settings, 2 runs concurrently.

  | | Empirical injections | Gaussian injections | Pair of runs (set wall) |
  |---|---|---|---|
  | No cache | 161.8 s | 111.7 s | 162.6 s |
  | First use (cache built) | 149.4 s | 103.2 s | 150.1 s |
  | Later use (cache filled) | **84.1 s (×1.9)** | **85.0 s (×1.3)** | **85.6 s (×1.9)** |

  - All 4 cached runs are **byte-identical** to the uncached run: same hits, same stage outcome for every injected
    source, |ΔS| = 0. **Extra misses: 0.**
- **Scope of the gain:**
  - It comes from not re-reading the 80″ PSF boxes and the analysis boxes, so it applies to repeated injection runs
    on the same box.
  - Measured once, on one box with one seed; treat it as this box's figure.
  - The cache takes 25 MB of local scratch per box pair (`speedup.json`), deleted after use. A reprocessed archive file needs a
    cleared cache (the key is URI + box).
  - **B2 and B3 do not gain:** their F150W matching takes 1,466 s and 1,948 s against 71 s in box 1, and the cache
    saves under 4 % there.
- **Main cost found:** the per-stamp leave-one-out kernel fit (`kernel_match --kernel-stamps`).
  - B2 and B3 hold 1,199 and 1,065 compact-source stamps, against 217 in box 1 (`batch/timing.json`).
  - Speeding up this step is not attempted in this unit (owner, 2026-10-17); it is left for a decision.
- **Superseded first attempt:** `speedup_v1.json` (cache pickled the WCS, not exact; old injection geometry).

## Batch search (PLAN §6, addendum 5)
- **Regions (all F150W + F444W, PEARLS o241 → VENUS o051, 24″ boxes, the final chain):**
  - **B1:** the #286 box, re-screened (section "Effect on #286").
  - **B2** (15.717097, −49.262488): 8 catalogued images.
  - **B3** (15.731568, −49.270822): 7 catalogued images.
  - B2 and B3 were chosen from the image list and coverage only (`boxes.json`, committed before any difference),
    100 % finite.
- **Quality:**
  - Registration (`batch/*/summary.json`):
    - B2: F150W (0.54, −0.61) px over 15 tiles, F444W (0.09, −0.18) px over 16;
    - B3: F150W (0.75, −0.34) px over 16 tiles, F444W (0.14, −0.14) px over 16.
  - None is pinned.
- **Local sensitivity:** one injection seed per box (B2 seed 40, B3 seed 41); n = 6–7 per cell; review list with the
  correct sign, 68 % intervals (`batch/b*_inj/eval.json`). No pass rule is applied to these boxes; the numbers
  describe them.

  | Class / AB | B2 | B3 |
  |---|---|---|
  | vanish 26.0 | 5 / 7 [0.53, 0.85] | 6 / 7 [0.68, 0.94] |
  | vanish 27.0 | 6 / 7 [0.68, 0.94] | 6 / 7 [0.68, 0.94] |
  | vanish 27.5 | 4 / 6 [0.46, 0.82] | 6 / 6 [0.86, 1.00] |
  | partial 26.0 | 7 / 7 [0.88, 1.00] | 5 / 7 [0.53, 0.85] |
  | partial 27.0 | 5 / 6 [0.64, 0.93] | 2 / 6 [0.18, 0.54] |
  | partial 27.5 | 3 / 7 [0.26, 0.61] | 2 / 7 [0.15, 0.47] |
  | unchanged: any hit | 0 / 20 [0, 0.05] | 0 / 20 [0, 0.05] |

  - Most vanish losses in B2 at 26.0 AB (2 of 7) come from the `near_bright_core` flag.
  - **B1's efficiencies do not transfer:** each box has its own, and with n = 6–7 the intervals are wide.
- **Screen:**
  - B2: 2 hits, both `near_bright_core`; morphology dipole and sign_flip.
  - B3: 6 hits, 1 unflagged (h005, sign_flip); morphology 3 edge, 2 sign_flip, 1 pass (h001, flagged).
  - **Review lists: B2 0, B3 0** (cap 15, not reached).
- **Spot-check of rejected hits (sheets `batch/b2/sheet.png`, `batch/b3/sheet.png`, looked at):**
  - B2 has only 2 hits and B3's single unflagged one is h005, so all hits of both boxes were looked at.
  - B2:
    - h000: a sign-flipping residual on a bright compact core;
    - h001: a ± ridge along an edge-on galaxy.
  - B3:
    - h000, h002 and h003: ridge or dipole residuals of elongated galaxies at the box edge;
    - h001: a compact-core residual;
    - h004: an edge-on galaxy ridge;
    - h005: a two-band hit with mixed sign. F444W is a ± dipole (S −14.6 / +11.1); F150W has only a marginal
      +5.5 peak and no compact feature on the sheet (S extremes within 0.3″, session check, not stored). A dipole
      is a residual, not a point-like fade.
  - All 8 rejections are consistent with their verdicts. None is a missed point-like fade.
- **Pending signals:** none. Box 1's only entry is the #286 h024 position, already vetted.
- **Zero entries is not an exclusion.** Fades, simulations and passing CI are not evidence of new physics.

## Regenerate
- **Samples:**
  `fade_calib.py make --uri <o241 F150W i2d> --ra <box> --dec <box> --seed S --out sample_*.csv`.
- **Boxes:**
  `fade_calib.py boxes --uri <4 i2d> --exclude 15.7300 -49.2530 --out boxes.json`.
- **`arc_diff.py`:** the command is in each `summary.json`. Maps go to the session scratchpad and are not stored.
- **Then:**
  - `arc_diff_morph.py --hits … --maps … --all-hits`;
  - `fade_calib.py eval --run DIR MORPH [--run …] --out eval.json`.
