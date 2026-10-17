# Fade calibration + batch speed-up + small search: cross-session state

Owner instruction 2026-10-16 (replaces the separate "fading calibration" and "batch speed-up" instructions).
Total budget 4 h, from 2026-10-16T23:32Z to **2026-10-17T03:32Z**. Any session taking over must read this file
first and continue from "Next action". Do not re-claim; continue on branch `claude/fade-calib` (PR #288).

| Phase | Budget | Status |
|---|---|---|
| 1. Orient, claim | — | done 23:32–23:40Z |
| 2. Fading calibration (≤ 2 h) | to 01:40Z | done 23:40–00:00Z: **failed** the pre-registered rule (unchanged false hits 0.25–0.28 > 0.10); one fix used (`--reg-tile-arcsec`) |
| 3. Measured speed-up (≤ 1 h) | — | done 00:00–00:05Z: raw-data cache, warm ×1.9 (empirical) / ×1.1 (Gaussian), 0 extra misses, not byte-identical (baseline is not either) |
| 4. Small batch search (≤ 1 h, only if phase 2 passes) | — | **not run** (phase 2 failed); box 1 re-screened as a check of #286 only |

Results: `RESULT.md`. Plan and addendum 3: `PLAN.md`.

## Next action
State files, `/code-review`, PR ready, CI, merge under the merge policy; final report to the owner. After that the
unit is closed. Re-opening needs an owner decision: either a criterion on review-list entries instead of raw hits
(a science-criterion change), or a different fix for unchanged-source residuals.
