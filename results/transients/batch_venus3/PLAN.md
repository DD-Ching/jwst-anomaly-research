# VENUS batch 3: plan (fixed before box selection or any difference)

Owner instruction 2026-10-17 (after #309). Budget 16:06Z → **20:06Z** (selection, calibration, search, review,
records). Rules of batch 1 (`../batch_venus/PLAN.md` §2–§5) and batch 2 (`../batch_venus2/PLAN.md`) apply unless changed
here. An ordinary transient is not wormhole evidence; zero entries is not an exclusion.

## 1. Already done / not repeated
- Batches 1–2: 40 boxes (`../batch_venus/boxes.json`, `../batch_venus2/boxes.json`), excluded by position, including
  their 16 undeterminable boxes (kept on hold, not re-run). Earlier exclusions: the batch-1 `--exclude` list.
- Paused units (S6a, X4/R1, others) not restarted; MACS0417 delay work ends at #309.
- Open PRs and branches checked 16:06Z: no other VENUS claim.

## 2. Selection (`scripts/batch_pairs.py`, metadata only)
- Same input (`results/transients/sfft_trial/venus_pairs.json`), pair rule, 3 × 3 grid of 24″ boxes, footprint check,
  ≥ 95 % finite pixels in all four images, `--per-pair 2`, `--max-boxes 20`, touching-box tolerance (#304).
- **New (owner instruction), `--prefer-extra-epochs`:** candidate boxes from all pairs are ranked by depth tier
  (factor-2 bins of the shallower F150W exposure, deeper first), then by the number of extra public F150W+F444W visits
  whose F150W footprint contains the box centre (more first), then pair rank and distance. Extra epochs are counted
  among the programmes listed in `venus_pairs.json` and other VENUS visits only (newer programmes are not counted).
  No difference image or source detection is used; a second-epoch detection is never a selection criterion.
- Fewer than 20 qualifying boxes → fewer boxes. `boxes.json` is committed before any difference image.

## 3. Search, calibration (unchanged rules)
- `arc_diff.py --half-arcsec 12 --kernel-px 4 --kernel-stamps --min-bands 2 --reg-tile-arcsec 6 --procs 1`
  (F150W + F444W), `arc_diff_morph.py --all-hits`; review list = hit ∧ not near a bright core ∧ morphology pass.
- Per box: `fade_calib.py make` (80 positions, unchanged / appear / partial / vanish at 26.0 / 27.0 / 27.5 AB),
  **seed = 300 + box index**; pass = vanish 26 + 27 AB on the review list with sign ≥ 0.70 and ≤ 2 control false hits.
- Status per box: completed + passed / completed + calibration failed / timed out / not started. Only passed boxes
  count as reliable coverage.
- **Time caps (fixed now, from batches 1–2: normal screens ≤ 384 s, injections ≤ 723 s; the one per-stamp-fallback box
  took 2,294 s and > 5,000 s):** screen 900 s, injection 1,200 s per task (`--timeout-screen 900 --timeout-inj 1200`).
  A task over its cap is killed, saved and its box marked undeterminable (timed out). Caps are not extended after
  results are seen.
- **No new search task starts after 18:20Z** (`--no-start-after`). 4 workers, `OMP_NUM_THREADS=1`, shared cache.
- Expected blind spot of the caps: boxes whose kernel stamps exceed `max_pix` (dense, bright cluster cores), which
  take the slow per-stamp path; they are recorded, not searched.
- Tool versions: `main` at d190eeb plus this PR's `batch_run.py` (timeout, stop time) and `batch_pairs.py`
  (`--prefer-extra-epochs`); `arc_diff.py`, `fade_calib.py`, `arc_diff_morph.py`, `hit_lightcurve.py`,
  `batch_review.py` unchanged.

## 4. Review (unchanged rules; capacity 10)
- Order: review-list entries of passed boxes by max|S| descending, then those of undeterminable boxes by max|S|.
- Each: sheet looked at; `hit_lightcurve.py` forced photometry in every covering mosaic; coded class; the stored #296
  post-review checks and the #304 full-rule function (sky-fixed ≤ 0.1″, all-mosaic |z| ≥ 3 ⇒ undetermined).
- Entries beyond capacity are saved in `pending.json`, not judged. Names: `B3-<target>-<box>-<hit> (ra, dec)`.
- A real change gets catalogue/literature checks (a subagent may help; TNS tried once) and, if a concrete question
  remains, a `/vet-candidate`-style record naming the cheapest discriminating next check. No new unit is opened.
- Review computation ends by 19:00Z; 19:00Z–20:06Z records, review, PR. Stop at 20:06Z; no next batch.
