# VENUS batch 4: plan (fixed before box selection or any difference)

Owner instruction 2026-10-18 (after #323), branch B: the same-source pairing feasibility test failed its frozen pass
rule (`../../same_source/a2744_feasibility/`), so this round returns to the fixed transient tools. Budget of the whole
round 03:02Z → **07:02Z**; at least 45 min kept for review, records and PR. Rules of batch 3 (`../batch_venus3/PLAN.md`
§2–§4) apply unchanged unless stated here. An ordinary transient is not wormhole evidence; zero entries is not an
exclusion; calibration-failed boxes are not reliable zero-candidate regions.

## 1. Already done / not repeated
- Batches 1–3: 60 boxes, excluded by position (batch-3 `--exclude` list plus the 20 batch-3 box centres), including
  the 26 held undeterminable boxes (not re-run). Paused units are not restarted.
- Open PRs and `claude/*` branches checked 03:03Z: no other VENUS claim.

## 2. Selection (unchanged tool and options)
- `batch_pairs.py --pairs results/transients/sfft_trial/venus_pairs.json --exclude <batch-3 list + 20 batch-3 centres>
  --prefer-extra-epochs --per-pair 2 --max-boxes 20 --min-finite 0.95` (metadata only). `boxes.json` is committed
  before any difference image. Fewer qualifying boxes → fewer boxes.

## 3. Search and calibration (unchanged)
- `arc_diff.py`, `arc_diff_morph.py`, `fade_calib.py` as batch 3; seed = 400 + box index; pass = vanish 26 + 27 AB
  recovered on the review list with sign ≥ 0.70 and ≤ 2 control false hits.
- Caps unchanged: screen 900 s, injection 1,200 s per task. **No new search task starts after 04:50Z.** 4 workers,
  `OMP_NUM_THREADS=1` (4 cores, 15 GB measured 03:03Z).

## 4. Review (unchanged; capacity 10)
- Order, sheet look, forced photometry, coded class, #296 post-review checks and #304 full-rule function as batch 3.
  Names `B4-<target>-<box>-<hit> (ra, dec)`.
- If a review entry is a real change worth independent checking, no further search volume is added; the remaining
  time goes to its evidence, ordinary explanations and the cheapest discriminating check.
- Review computation ends by 06:00Z; records, review and PR after that. Stop at 07:02Z; no next batch.
