# Per-stamp kernel-fit cost (`arc_diff.kernel_match --kernel-stamps`): cross-session state

Owner instruction 2026-10-17 (follow-up to #288): one unit, one main session, no `/batch`. Budget 3 h, from
2026-10-17T02:23Z to **05:23Z**. Read this file first and continue from "Next action"; branch `claude/kernel-loo`.

Scope: only the measured per-stamp leave-one-out kernel-fit cost (#288 `batch/timing.json`). At most one
implementation. No algorithm swap, GPU or new pipeline, no new sky regions, no SMACS check, no re-run of history.

## Next action
Write PLAN.md (diagnosis, the one method, fixed tolerances and test data) before any comparison run.
