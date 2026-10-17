# Per-stamp kernel-fit cost: cross-session state

Owner instruction 2026-10-17. Budget 02:23Z → **05:23Z**; branch `claude/kernel-loo` (PR #291).

Done: plan; shared-design implementation (`--kernel-stamp-fit shared`, now the default, D-094); validation (all
PLAN §3 tests pass on B1 and B2); post-review precision fix (guard + refinement, PLAN addendum 1); every B1/B2 stamp re-checked; second review fixes (addendum 2); final timing (B2 whole run ×14 screen, ×8–9 injection run). Results: `RESULT.md`.

## Next action
`/code-review`, PR ready, CI, merge under the merge policy. This optimisation line ends here.
