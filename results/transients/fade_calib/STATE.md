# Fade calibration + batch speed-up + small search: cross-session state

Owner instruction 2026-10-16 (replaces the separate "fading calibration" and "batch speed-up" instructions).
Total budget 4 h, from 2026-10-16T23:32Z to **2026-10-17T03:32Z**. Any session taking over must read this file
first and continue from "Next action". Do not re-claim; continue on branch `claude/fade-calib`.

| Phase | Budget | Status |
|---|---|---|
| 1. Orient, claim | — | done 23:32–23:40Z |
| 2. Fading calibration (≤ 2 h) | to 01:40Z | in progress |
| 3. Measured speed-up (≤ 1 h) | — | not started |
| 4. Small batch search (≤ 1 h, only if phase 2 passes) | — | not started |

## Done
- Existing injections in `scripts/arc_diff.py` (`--inject`, `--inject-grid`) are added to the *matched* epoch 2,
  after resampling, registration and PSF / kernel matching. Their recoveries are efficiencies *after the
  injection point* (#275, #276, #277, #286). This is noted here; historical results are not re-run.

## Next action
Write `PLAN.md` (sample, pass rules) before any run; then implement pre-matching injection.
