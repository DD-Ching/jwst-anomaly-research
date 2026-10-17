# SFFT vs `arc_diff` (#291 baseline): bounded tool comparison, cross-session state

Owner instruction 2026-10-17. Budget 2 h, from 04:02Z to **06:02Z** (hard stop; save and stop). Branch
`claude/sfft-trial`. Purpose: choose the tool for a future batch search; no long-term tool development.

Step budget: ≤ 30 min to decide whether SFFT can be used as is (paper implementation, licence, CPU, dependencies).
If integration needs a large rewrite, a GPU port or new infrastructure, record the blockers and stop. Only B1 / B2
small boxes and the existing injection settings are reused; comparison rules are fixed before any run.

Prior repository state: DECISIONS addendum 2026-10-15c ("HOTPANTS/SFFT not adopted ... revisit if a field with a
spatially varying PSF (box > 1′) leaves ridges"); no earlier SFFT evaluation.

## Status
Done: integration check (paper method GPU-only), comparison with two rounds of post-review corrections (keep
`arc_diff`, D-095), VENUS pair denominator.
Results: `RESULT.md`.

## Next action
`/code-review`, PR ready, CI, merge under the merge policy. No further SFFT work.
