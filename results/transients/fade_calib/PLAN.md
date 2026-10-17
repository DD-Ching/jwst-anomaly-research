# Fading / vanishing calibration of `arc_diff`, measured speed-up, small batch search: plan (before any run)

Owner instruction 2026-10-16; 4 h total from 23:32Z. Cross-session state: `STATE.md`. Long-term aim: testable
evidence for special lensing (W3: a compact lensed image that dims or vanishes). An ordinary fading is not wormhole
evidence.

## 0. What earlier recoveries mean
`arc_diff.py --inject` and `--inject-grid` add stars to epoch 2 *after* resampling, registration and PSF / kernel
matching. Their recoveries (#275, #276, #277, #286) are therefore **efficiencies after the injection point**. They
do not test whether registration or matching erase or create a signal. Historical results are not re-run.

## 1. Pre-matching injections (new option `arc_diff.py --inject-pre CSV`)
- Each row is `name, ra, dec, mag1, mag2`. An empty magnitude means "no source in that epoch".
- Sources are added to the **raw epoch-1 and epoch-2 boxes**, each on its own WCS and pixel grid with its own
  `PIXAR_SR`, **before** resampling, registration, PSF and kernel matching. Everything after that runs unchanged.
- **Point source shape:** a Gaussian with the band's FWHM (`PSF_FWHM`) in both epochs. ASSUMPTION: same instrument,
  band and pixel scale; wings, spikes and the 0.84° roll are not modelled.
- **Noise:** none is added for the injected source. ASSUMPTION: its Poisson noise is negligible next to the
  background at ≥ 26 AB in these exposures.
- **Not a calibration:** subtracting real image blocks, or flipping the sign of a difference.

## 2. Fixed sample (El Gordo box of #286: 24″ at (15.7300, −49.2530), F150W + F444W, PEARLS o241 → VENUS o051)
- **Positions:** `grid_positions` on the **epoch-1** F150W box only. These are pixels above the epoch-1 90th
  percentile (galaxies, arcs, ICL), 1.5″ apart, 1″ margin, up to 80 per seed. Epoch 2 plays no part in choosing
  them.
- **Classes and magnitudes:** classes cycle in order, and magnitudes cycle independently over (26.0, 27.0, 27.5)
  AB.

  | Class | Epoch 1 | Epoch 2 |
  |---|---|---|
  | unchanged | m | m |
  | brighten (appear) | none | m |
  | partial fade | m | m + 0.75 (flux × 0.5) |
  | vanish | m | none |

- **Seeds:**
  - development 10 and 11, used for diagnosis and for at most one fix;
  - **held-out 20 and 21**, never used for any adjustment. The pass rule is evaluated on them only.
- **Settings:** `--kernel-px 4 --kernel-stamps --min-bands 2 --half-arcsec 12`, as #286.

## 3. Stages counted per injected source (`scripts/fade_calib.py`; match radius 0.15″ = `arc_diff` default)
1. Hit: |S| ≥ 5 in ≥ 2 bands, so the multi-band requirement is part of this stage.
2. Not `near_bright_core`.
3. `arc_diff_morph.py --all-hits` verdict `pass`.
4. Correct sign: `-` or `mixed` with a negative maximum for fades; `+` for appear.
5. Review list = 1 ∧ 2 ∧ 3. This is the programmatic end of the chain.

- **False hits:** for the unchanged class, any hit within 0.15″.
- **Not programmatic and not counted:**
  - by-eye sheet review;
  - forced-photometry vetting;
  - positions within 0.8″ of bright cores, which are reported but never reviewed.

## 4. Pass rule (held-out seeds 20 + 21 only)
- **Vanish at 26.0 AB:** review-list fraction with the correct sign ≥ 0.80.
- **Unchanged (all magnitudes):** false-hit fraction ≤ 0.10.
- Everything else (partial fades, 27.0 / 27.5, appear) is reported with binomial 68 % intervals and no threshold.
- **One fix allowed:** at most one, for a cause diagnosed on the development seeds. If made, the held-out seeds are
  run only after it. No threshold is tuned.

## 5. Speed-up (≤ 1 h; same sample and settings)
- Time each stage (S3 box reads, resample / registration / matching, detection, injections, morphology).
- Improve only the largest one: candidates are a local cache of the raw boxes so that each box is read once, and
  running several seeds or classes in one process.
- Never cache anything that depends on the injected sources (registration, kernels, maps): pre-matching
  injections must go through matching every time.
- **Fair comparison:** the same seeds give the same injected rows. Report wall time, stage counts and false hits
  with and without the change. Outputs must agree exactly, otherwise no speed-up factor is reported.

## 6. Small batch search (≤ 1 h; only if section 4 passes)
- **Box 1:** the #286 box was already screened with the unchanged detection chain. If no fix changes the chain, its
  4 fading entries (all explained, #286) stand, and only the calibration is new; the screen is not redone.
- **Up to 2 more boxes** in the same El Gordo pair, chosen **before any difference is seen**, from coverage and
  data quality only: inside both epochs' F150W and F444W footprints, not overlapping box 1, ≥ 5 catalogued images,
  epoch-1 finite fraction ≥ 0.95.
- **Each new box:**
  - one pre-matching injection seed (vanish / partial / unchanged as section 2);
  - the review cap is **15 entries per box**. Above it, save all flags and stop for that box; no new rule is
    added.
  - Vetting as #286: sheet, forced photometry, spike and on-spike checks, bright-source controls where they apply.
  - Spot-check 3 rejected (non-pass or flagged) hits per box on the sheet.
- **Outcome:** an entry with no ordinary or instrumental explanation goes to `/vet-candidate`, and the owner is
  told. Zero entries gives no exclusion of wormholes or of their abundance.

## Addendum 3 (2026-10-16T23:55Z, after the development seeds, before any held-out run)
Development runs are stored in `dev/` (`<variant>_s10`, `<variant>_s11`, `<variant>_eval.json`).

- **Diagnosis (dev seeds 10 + 11, plan settings):**
  - Unchanged-class false hits are 12 / 40 (`dev/dev_eval.json`); 26.0 AB unchanged sources give F150W S of about +40
    to +97.
  - The F150W whole-box registration is (1.71, −3.00) px, pinned at the ±3 px window edge. This value is the same
    without any injection (#286 `o241_o051/summary.json`).
  - On 16 non-overlapping 6″ tiles the WCS-resampled epoch 2 needs a median of (0.70, −0.57) px, consistent with the
    ≤ 0.02″ catalogue offsets (`reg6882/`). A few tiles near the saturated star pull the whole-box peak.
  - The local tile shift is capped at 2 px. After it, several tiles stay 1–2 px off in x (session diagnostic, not
    stored).
  - Earlier injections were made after this step and could not see it.
- **Diagnostic that is not the fix:**
  - Per-epoch empirical PSF injections alone (`dev/fix_*`) leave false hits at 12 / 40.
- **The one fix (chain change):**
  - New `arc_diff.py --reg-tile-arcsec 6`: the global shift is the median of per-tile registrations over
    non-overlapping 6″ tiles.
  - It is the default for this unit, and from here on for every box in section 6. Default 0 keeps historical runs
    reproducible.
  - On dev: false hits 8 / 40 with Gaussian injections (`dev/reg_eval.json`) and 4 / 40 with empirical PSFs
    (`dev/regE_eval.json`); vanish 26.0 AB review with sign 12 / 14 for both.
- **Injection-model deviation (post hoc, chosen after seeing the dev results):**
  - Primary model: per-epoch empirical PSF (`--inject-pre-psf empirical`). The 5 most concentrated isolated sources
    with S/N ≥ 50 are chosen on epoch 1 in a 40″-half box, cut at the same sky positions in both epochs and
    median-stacked per epoch. This implements the owner's "per-epoch point-source shape"; the §1 band-FWHM Gaussian
    is identical in both epochs.
  - Because the choice was made after the Gaussian variant failed the dev false-hit rule, the held-out seeds are run
    with **both** models and both are reported. The pass rule is applied to the empirical model.
  - The Gaussian result is reported beside it, and its failure, if it fails, is stated.
- **Consequence for section 6:** the chain changed, so box 1 (#286) is re-screened with `--reg-tile-arcsec 6`.
- **Unchanged:** thresholds, the pass rule, the samples and the held-out seeds.

## Addendum 4 (2026-10-17T00:25Z, after `/code-review` of the held-out run, before any further run)
- **Review finding (correctness):**
  - `--inject-pre` placed each epoch's source at the WCS position of (ra, dec).
  - The chain then shifts WCS-resampled epoch 2 by the measured registration (F150W ≈ (0.70, −0.57) px), so every
    injected epoch-2 source ended up displaced by that shift. Real sources are brought into line by the same shift.
  - An "unchanged" injection was therefore a 0.9 px dipole. The held-out result of seeds 20 + 21 and every dev
    variant measured that injection artifact, not the chain.
  - Those results stay in `heldout/` and `dev/` and are reported as **invalid injection geometry**. They are not
    used for the pass rule.
- **Corrections (none changes a threshold, the samples or the pass rule):**
  1. Epoch-2 sources are placed where a real sky-fixed source lands (`epoch2_pixel`). The chain's global
     registration is measured once on the injection-free boxes (`pair_shift`, same method as the chain).
     - Local tile shifts (≤ 2 px cap; El Gordo residual median ≈ 0.2 px) are not emulated.
  2. `register_median` drops tiles < 90 % finite in either image or pinned at the window edge. This is a review
     robustness patch of the addendum 3 fix, not a new fix. If it changes box 1's shift, box 1 is re-screened.
  3. The cache no longer pickles a WCS; it keeps the header text, so a cached box is exact. Every function argument
     is in the cache key.
- **New held-out seeds 30 and 31** (`sample_s30.csv`, `sample_s31.csv`, written now, never run).
  - Seeds 20 and 21 were seen, so they are no longer held out.
  - Dev seeds 10 + 11 are re-run first for diagnosis only. **Nothing is changed after the dev rerun, whatever it
    shows.**
  - Then 30 + 31 run once with both injection models. The pass rule (§4) is applied to the empirical model; the
    Gaussian is reported beside it.
- **Kept (review point noted):**
  - `near_bright_core` and the morphology epoch-1 map include the injected epoch-1 sources. A real fading source is
    also in epoch 1, so this is realistic for one source.
  - The density effect of about 60 injected epoch-1 sources on the 99.9th-percentile threshold is a stated
    limitation. Per-stage counts are reported.
