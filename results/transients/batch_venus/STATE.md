# Bounded VENUS batch search with fixed `arc_diff` (#291, D-095): cross-session state

Owner instruction 2026-10-17. Budget 4 h, from 06:01Z to **10:01Z** (hard stop), including data checks,
calibration, search, review and records. Branch `claude/batch-venus`.

Fixed rules (from the instruction): ≤ 20 unsearched boxes chosen from pair metadata only (common bands, coverage,
dates, versions, expected sensitivity), never from differences; program version, screen rules, per-box injection
calibration and unchanged controls fixed before any search; boxes failing calibration are "undeterminable"; one
fixed program over independent boxes in parallel, shared cache; ≤ 10 detailed reviews, ordered by a pre-fixed rule;
the rest saved as pending. No follow-on batch, no tool work, no SFFT.

## Progress
- 06:02Z #293 report wording corrected (no re-runs).
- 06:04Z PLAN committed; 06:09Z boxes.json committed (20 boxes, 10 pairs, 3.2 arcmin², metadata only).
- 06:10Z program fixed (d9640f8; `arc_diff.py` unchanged from main 7c975ac); batch started 06:10Z
  (`scripts/batch_run.py run`, 4 workers, maps and cache in the session scratchpad).

## Next action
When the batch finishes: `batch_run.py summarize` (also run by `run`), look at sheets, `batch_review.py`
(capacity 10), RESULT and state files. Hard stop for computation 09:30Z, for everything 10:01Z.
