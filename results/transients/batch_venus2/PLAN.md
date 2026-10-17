# VENUS batch 2: plan (fixed before any box selection run or difference)

Owner instruction 2026-10-17 (after #300). Budget 10:55Z → **14:55Z**, covering calibration, search, review and
records. Same fixed rules as batch 1 (`../batch_venus/PLAN.md` §2–§5); only what differs is written here. An
ordinary fading is not wormhole evidence; zero entries is not an exclusion.

## 1. Already done (not repeated)
- Batch 1 (#296): 20 boxes (`../batch_venus/boxes.json`), including the 8 undeterminable ones, which are **not
  re-run**. Earlier searched boxes: El Gordo B1–B3, Sunrise, SMACS system 5 (batch-1 `--exclude` list).
- Open PRs and `claude/*` branches checked 10:56Z: no other claim on VENUS boxes.
- Tools unchanged since #296 (`arc_diff.py`, `arc_diff_morph.py`, `fade_calib.py`, `hit_lightcurve.py`,
  `batch_run.py`, `batch_review.py`), except the selection fix below.

## 2. Box selection (`scripts/batch_pairs.py`, metadata only; batch-1 rule unchanged)
- Same input (`results/transients/sfft_trial/venus_pairs.json`), same pair rule, depth ranking, 3 × 3 grid of 24″
  boxes, footprint and ≥ 95 % finite-pixel checks, `--per-pair 2`, `--max-boxes 20`.
- `--exclude`: the batch-1 list plus the 20 batch-1 box centres.
- **Deviation (selection bug fix, labelled):** the overlap test treated boxes that only touch a searched box (centres
  exactly 24″ apart, stored rounded to 1e-6°) as overlapping at random; 13 of 40 RA-neighbours of batch-1 boxes were
  wrongly rejected by it. `overlaps()` now allows a 0.01″ tolerance (unit test in `tests/test_batch_rules.py`). It
  uses no difference information and changes no search threshold.
- Box ids are batch-local (`<target>_<k>` counted within this batch); positions identify boxes across batches.
- `boxes.json` is committed before any difference image.

## 3. Search, calibration, review (batch-1 §3–§5, unchanged)
- `arc_diff.py --half-arcsec 12 --kernel-px 4 --kernel-stamps --min-bands 2 --reg-tile-arcsec 6 --procs 1`,
  F150W + F444W; `arc_diff_morph.py --all-hits`; review list = hit ∧ not near a bright core ∧ morphology pass.
- Per-box injections (`fade_calib.py make`, 80 positions, seed 100 + box index) and the pass rule: vanish 26 + 27 AB
  on the review list with sign ≥ 0.70 and ≤ 2 unchanged-control false hits. Failing boxes are undeterminable.
- Review: passing boxes first, by max|S|; capacity 10 detailed reviews (sheet looked at, `hit_lightcurve.py` forced
  photometry, full class rule including sky-fixed offset and bright neighbour as in the #296 post-review checks);
  the rest saved in `pending.json`.
- A `real change` without an ordinary explanation goes to `/vet-candidate` and the owner is told.
- Wording: a source absent later is "not detected in later observations", not "vanished".

## 4. Hard stops
- 13:55Z: no new computation; records only.
- 14:55Z: stop. No follow-on batch.
