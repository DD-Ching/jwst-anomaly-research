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
