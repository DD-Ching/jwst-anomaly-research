# Fading / vanishing calibration of `arc_diff`, measured speed-up, small batch search: result

Plan: `PLAN.md`, sections 0–6 committed before any run (2225589). Addendum 3 was written after the development seeds
and before any held-out run (35748f4). Owner instruction 2026-10-16; budget 23:32Z → 03:32Z. Every number is
`derived` unless labelled; injected sources are simulations; thresholds are ASSUMPTIONs.

## Outcome in one line
**The calibration fails its pre-registered pass rule, so the batch search (section 6) was not run.**
- Vanishing sources are recovered: at 26.0 AB, 14 / 14 reach the review list with the correct sign.
- But unchanged sources produce hits at 0.25–0.28, against the rule's ≤ 0.10. Almost all of these hits are removed
  later by the morphology stage.

## What changed in the chain
- **New `--inject-pre`:** injects sources into the raw boxes of both epochs before resampling, registration and
  matching (PLAN §1). Each epoch gets its own grid and `PIXAR_SR`.
- **The one fix:** `--reg-tile-arcsec 6` (addendum 3).
  - **What was wrong:** the whole-box FFT registration of this pair is pinned at its ±3 px window edge in F150W:
    (1.71, −3.00) px, against a 6″ tile median of (0.70, −0.57) px. F444W gives (0.53, −0.95) against
    (0.11, −0.18). Values: `box1_regfix/compare_286.json`.
  - **Why it mattered:** the local tile shift is capped at 2 px, so parts of the F150W box stayed misregistered.
  - **Effect on #286:** this applied to the #286 screen too. Its injections were added after this step, so they
    could not show it.
  - **The fix:** the global shift is the median over tiles. The default stays 0, so historical runs reproduce.
- **Injection model (post-hoc deviation, addendum 3):**
  - Primary: per-epoch empirical PSFs. The 5 most concentrated isolated S/N ≥ 50 sources are chosen on epoch 1, in
    a 40″-half box, and cut at the same sky positions in both epochs.
  - The planned band-FWHM Gaussian (the same in both epochs) is run beside it on the held-out seeds.
  - Peak fractions:
    - F150W 0.243 / 0.214 (epoch 1 / epoch 2);
    - F444W 0.104 / 0.107 (`heldout/hoE_s20/summary.json`).
- `--cache-dir` (section 3 below).

## Calibration (PLAN §2–§4)
- **Sample:**
  - 80 positions per seed on epoch-1 F150W pixels above the 90th percentile, 1.5″ apart (`sample_s*.csv`). Epoch 2
    played no part.
  - Classes (unchanged, appear, partial m → m + 0.75, vanish) × 26.0 / 27.0 / 27.5 AB.
  - Dev seeds 10, 11; held-out seeds 20, 21. Each seed has n = 12–14 per cell.
- **Stages counted:** hit (≥ 2 bands) → not `near_bright_core` → morph `pass` → correct sign. The review list is
  the first three.
  - Per-cell stage counts: `heldout/ho{E,G}_eval.json`, `dev/*_eval.json`.
  - Intervals below are binomial 68 % Wilson intervals.

**Held-out (seeds 20 + 21), review list with the correct sign:**

| Class / AB | Empirical PSF (primary) | Gaussian (plan §1) |
|---|---|---|
| vanish 26.0 | **14 / 14** [0.93, 1.00] | 14 / 14 [0.93, 1.00] |
| vanish 27.0 | 12 / 14 [0.74, 0.93] | 12 / 14 [0.74, 0.93] |
| vanish 27.5 | 10 / 12 [0.70, 0.91] | 10 / 12 [0.70, 0.91] |
| partial 26.0 | 13 / 14 [0.83, 0.97] | 13 / 14 [0.83, 0.97] |
| partial 27.0 | 7 / 12 [0.44, 0.71] | 10 / 12 [0.70, 0.91] |
| partial 27.5 | 5 / 14 [0.24, 0.49] | 8 / 14 [0.44, 0.70] |
| appear 26.0 / 27.0 / 27.5 | 8 / 12, 9 / 14, 11 / 14 | 10 / 12, 10 / 14, 13 / 14 |
| **unchanged: any hit (false hit)** | **11 / 40 = 0.275** [0.21, 0.35] | **10 / 40 = 0.25** [0.19, 0.32] |
| unchanged: on the review list | 1 / 40 (27.5 AB) | 0 / 40 |

- **Pass rule (§4):**
  - vanish 26.0 ≥ 0.80: passed;
  - unchanged false hits ≤ 0.10: **failed** under both injection models.
  - **Calibration: failed.** No threshold was tuned, and the one fix is used.
- **Development seeds (10 + 11), unchanged false hits** (`dev/`):
  - plan chain 12 / 40;
  - empirical PSF alone 12 / 40;
  - registration fix with Gaussian 8 / 40;
  - registration fix with empirical PSF 4 / 40.
  - Development and held-out disagree (4 / 40 against 11 / 40), so with n = 40 the dev value was optimistic.
- **Where the false hits are:**
  - Of the 11 held-out (empirical) false hits, 3 lie within 0.16″ of a hit in the un-injected run. The other 8 are
    caused by the injected source itself (session match against `box1_regfix/hits.ecsv`; not stored).
  - Most unchanged hits are F150W-dominated residuals of a compact 26 AB source on a galaxy or arc pixel.
  - Every one except one is removed by morphology (rosette / dipole / edge). The review list therefore holds
    0–1 / 40.
  - The frozen rule counts hits, not review entries, so it fails. A rule on review entries would be a post-hoc
    criterion change and needs the owner.
- **Blind spots (not calibrated):**
  - extended fades;
  - positions within 0.8″ of bright cores (not reviewed);
  - sky and blank-field positions;
  - fainter than 27.5 AB;
  - injected-source Poisson noise (ASSUMPTION: negligible);
  - PSF wings beyond 0.6″ (truncated);
  - only one box and one pair.

## Effect on #286 (box 1 re-screened with the fixed registration, no injections; `box1_regfix/`)
- **Screen:** 50 → 39 hits; 18 → 7 unflagged; review list 4 → 1.
- **The remaining entry**, h020, is #286 h024 at the same position (0.0″). It was already vetted as a spike residual
  with no forced-photometry change ≥ 2σ.
- **The other #286 entries:**
  - h031 and h049 are no longer separate hits;
  - h043's nearest hit (1.3″) is morph `extended`.
- **No new entry; the #286 conclusion (nothing open) stands.**
- **Its sensitivity numbers do not transfer.**
  - The #286 brightening efficiencies were measured after the injection point and with the old registration.
  - Fading sensitivity with the pre-matching injections is in the table above.
- **Other fields:** they used the whole-box registration too. A scan of every stored `results/transients/**/summary.json`
  `shift_px` shows the ±3 px edge is reached **only in the El Gordo runs**. The other runs are:
  - MACS0416 tiles |shift| ≤ 1.2 px;
  - Sunrise ≤ 0.3 px;
  - Abell 2744 grid ≤ 1.8 px;
  - SMACS GO-4043 ≤ 0.5 px.
  - SMACS ERO → VENUS (`smacs_arcdiff/o001_o057`, `arc_chain/runs/smacs*`) has a large but unpinned shift of
    (2.11, −1.82) px in F150W. It is **not** compared with a tile median, so it stays a pending check. That search
    recovered SN Atalanta.

## Speed-up (PLAN §5; `speedup.json`)
- **Profile (one F150W band, 78.5 s):**
  - kernel matching 55.1 s (`lstsq` 30.7 s, median filters 21.0 s);
  - significance 6.8 s;
  - S3 reads 8.5 s.
  - With empirical-PSF injections, a second 40″-half read and its stacking add about 40 s per band.
- **What was changed:** only raw-data products are cached (`--cache-dir`): the analysis boxes and the per-epoch
  empirical PSFs. Matching depends on the injected sources and is never cached.
- **Fair comparison:**
  - Same 4 held-out runs, the same samples and settings, 4 concurrent processes.
  - Per-run wall time:
    - baseline 168–171 s (empirical) and 96–98 s (Gaussian);
    - warm cache 86–88 s for all four, which is about ×1.9 for the empirical runs and ×1.1 for the Gaussian runs.
  - The cold cache matches the baseline: set wall 163 s.
- **Agreement (post hoc wording):**
  - Plan §5 required byte-identical outputs. They are not byte-identical, but the un-cached chain itself is not
    byte-reproducible run to run: in the cold G21 run, S differs by 2.4 × 10⁻⁸ without any cache effect.
  - With the cache, S differs by at most 4.8 × 10⁻⁷ relative.
  - Hit sets, flags and signs are identical, and every injected source's stage outcome is identical in all 8
    comparisons: **0 extra misses**.
  - Under the plan's wording the factor is therefore not a validated speed-up; it is reported as measured.
- **Cost:**
  - 25 MB of local scratch per box pair; it must be deleted after use.
  - The cache key is URI + box, so a reprocessed archive file needs a cleared cache.
  - No change for single un-injected screens: reads are about 8 s of about 90 s.

## Batch search (PLAN §6)
**Not run:** the calibration failed. No extra regions were chosen, and there is no review list or vetting beyond
the box-1 re-screen above.

## Regenerate
- Samples:
  `fade_calib.py make --uri <o241 F150W i2d> --ra 15.7300 --dec -49.2530 --seed S --out sample_sS.csv`.
- `arc_diff.py`: the commands are in each `summary.json`. Maps went to the session scratchpad and are not stored.
- Then `arc_diff_morph.py --hits … --maps … --all-hits` and
  `fade_calib.py eval --run DIR MORPH [--run …] --out eval.json`.
