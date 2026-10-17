# Fade calibration + batch speed-up + small search: cross-session state

Owner instruction 2026-10-16. The owner's follow-up on 2026-10-17 keeps the hard stop at **2026-10-17T03:32Z**: no
new research items, and no restarts of finished computations. Continue on branch `claude/fade-calib` (PR #288).

| Phase | Status |
|---|---|
| 2. Fading calibration | done. Seeds 20 + 21 are invalid (injection geometry); fresh held-out seeds 30 + 31 pass (addendum 4) |
| 3. Measured speed-up | done. Exact cache; box 1 warm ×1.9 / ×1.3; byte-identical; no gain in B2 / B3 |
| 4. Small batch search | done. B1 review 1 (vetted), B2 0, B3 0; spot-checks looked at |

Results: `RESULT.md`. Plan and addenda 3–5: `PLAN.md`.

## Next action
Final `/code-review`, PR ready, CI, merge under the merge policy, then the final report to the owner. If time runs
out: leave the PR open with this file as the handoff.
