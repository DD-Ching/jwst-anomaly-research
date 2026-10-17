# Bounded VENUS batch search: plan (fixed before any box selection run or difference)

Owner instruction 2026-10-17 (after #293 / D-095). Budget 06:01Z → 10:01Z. One reproducible search of new data with
the fixed #291 tool. An ordinary fading is not wormhole evidence; zero entries is not an exclusion.

## 1. Already done (not repeated)
- **Searched VENUS-pair boxes** (excluded by position):
  - El Gordo B1 (15.7300, −49.2530), B2 (15.717097, −49.262488) and B3 (15.731568, −49.270822): #286, #288;
  - Sunrise (24.346822, −8.464511), o010 → o052: `sunrise_arcdiff`;
  - SMACS system-5 box (110.824442, −73.452577), ERO → VENUS, including the tile-median registration re-run (#292).
- **Other fields:** the MACS0416 and Abell 2744 work used non-VENUS pairs. A VENUS pair there is new data.
- **Open PRs:** none claims this unit (checked 06:01Z).

## 2. Box selection (`scripts/batch_pairs.py`; metadata only)
- **Input:** `results/transients/sfft_trial/venus_pairs.json`, the 18 VENUS targets with another public NIRCam F150W
  level-3 program.
- **Pair, per target:**
  - one VENUS observation and one other-program observation, each with F150W and F444W level-3 mosaics whose MAST
    footprints contain the target centre;
  - if several qualify, the pair with the longest shallower F150W exposure (expected sensitivity).
- **Ordering:** pairs by that depth, deepest first.
- **Boxes:**
  - a 3 × 3 grid of 24″ boxes centred on the target (the cluster core, where lensed images concentrate), ordered by
    distance from the centre;
  - kept if the centre and corners lie in all four footprints, it does not overlap a searched box, and `read_box`
    finds ≥ 95 % finite pixels in all four images;
  - up to 2 per pair, in pair order, up to 20 boxes.
- **Recorded:** headers (DATE-OBS, CAL_VER, CRDS_CTX, PA_V3), reasons for every rejection, and the effective area
  (`boxes.json`).
- No difference image is made before `boxes.json` is committed.

## 3. Fixed program and settings
- **Code:** `main` at 7c975ac, with `arc_diff.py` unchanged (#291, D-094 shared stamp kernels, D-095).
- **Settings:** `--half-arcsec 12 --kernel-px 4 --kernel-stamps --min-bands 2 --reg-tile-arcsec 6 --procs 1`,
  bands F150W + F444W.
- **Chain:** `arc_diff_morph.py --all-hits`. Review list = hit (|S| ≥ 5 in both bands) ∧ not `near_bright_core` ∧
  morphology `pass`.
- **Execution:** `scripts/batch_run.py` runs the same commands over the boxes in a 4-worker pool (`OMP_NUM_THREADS=1`)
  with one shared raw-data cache. No per-box changes.

## 4. Per-box calibration and controls (fixed now)
- **Injections:** one pre-matching injection run per box (`--inject-pre`, `--inject-pre-psf empirical`).
  - The sample comes from `fade_calib.py make` on the box's epoch-1 F150W image: 80 positions above the 90th
    percentile, 1.5″ apart.
  - Classes cycle unchanged / appear / partial (m → m + 0.75) / vanish at 26.0 / 27.0 / 27.5 AB.
  - Seed = 100 + box index.
  - This gives about 20 unchanged controls, 20 vanish, 20 partial and 20 appear per box.
- **Pass condition per box (computed by `batch_run.py`):**
  - vanish at 26.0 + 27.0 AB on the review list with the correct sign ≥ 0.70 (about 10/14);
  - unchanged-control false hits (any hit) ≤ 2.
- **Failing boxes are "undeterminable".** Their screen entries are still saved, but the box is not reported as
  searched with zero candidates.
- **Thresholds:** none is changed for any box.
- **Valid range for a passing box:** fades that remove ≥ 26–27 AB of flux, i.e. a source that vanishes. Partial
  fades (a 50 % drop) and 27.5 AB are reported with their per-box efficiencies, without a pass claim. The calibration
  is valid only on galaxy and arc pixels (the injection positions), not within 0.8″ of bright cores.

## 5. Review order and capacity (fixed now)
- **Order:**
  - review-list entries of passing boxes first, by max|S| descending;
  - then those of undeterminable boxes, by max|S|.
- **Capacity:** at most 10 detailed reviews in the batch, within the time budget.
- **A detailed review consists of:**
  - the contact sheet, looked at;
  - forced photometry in every covering mosaic (`hit_lightcurve.py`);
  - the class rule.
- **Classes:**
  - `artifact`: a sheet veto, or edge, spike, persistence or registration;
  - `no significant change`: no band ≥ 3σ in forced photometry;
  - `real change`: same sign in both bands, ≥ 3σ in at least one, and sky-fixed;
  - `real change` without an ordinary explanation is `open`: it goes to `/vet-candidate` and the owner is told.
- **Entries beyond capacity or time:** saved in `pending.json` as pending, never dropped or called non-candidates.
- No follow-on batch.

## 6. Hard stops
- 09:30Z: no new computation; records only.
- 10:01Z: stop.
