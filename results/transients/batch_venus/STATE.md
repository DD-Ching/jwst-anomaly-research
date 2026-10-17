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

- 06:32Z batch done (40 runs, all rc 0): 12 boxes passed, 8 undeterminable; review queue 5.
- 06:35Z review done: 5 / 5 reviewed (capacity 10), 0 pending; PLCK G165 h002 real fade (ordinary transient).
- 06:40–06:50Z review fixes: full-rule post-review checks (h003 → undetermined), h002 vetted
  (`docs/candidates/batch_venus_PLCKG165_h002.md`), SCI PA_V3, rule tests, deviations labelled.
- Results: `RESULT.md`.

## Next action
`/code-review`, PR ready, CI, merge under the merge policy; final report. No follow-on batch.
