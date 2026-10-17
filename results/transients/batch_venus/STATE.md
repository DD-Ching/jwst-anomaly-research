# Bounded VENUS batch search with fixed `arc_diff` (#291, D-095): cross-session state

Owner instruction 2026-10-17. Budget 4 h, from 06:01Z to **10:01Z** (hard stop), including data checks,
calibration, search, review and records. Branch `claude/batch-venus`.

Fixed rules (from the instruction): ≤ 20 unsearched boxes chosen from pair metadata only (common bands, coverage,
dates, versions, expected sensitivity), never from differences; program version, screen rules, per-box injection
calibration and unchanged controls fixed before any search; boxes failing calibration are "undeterminable"; one
fixed program over independent boxes in parallel, shared cache; ≤ 10 detailed reviews, ordered by a pre-fixed rule;
the rest saved as pending. No follow-on batch, no tool work, no SFFT.

## Next action
Correct the #293 report wording (no re-runs), then the pair metadata check.
