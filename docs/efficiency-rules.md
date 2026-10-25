# Efficiency rules (owner rule 2026-10-26)

These rules apply to every round. They change the workflow only: no research method, threshold, review
count or quality requirement changes. Items marked **[auto]** are checked by a script or pre-commit
hook; **[manual]** items rely on this list.

## 1. Self-check before review
- Before the independent review, work through `docs/pre-review-checklist.md` and fix what it finds.
- The review and the targeted re-check use the same checklist as their baseline.
- **[auto, advisory]** `python scripts/check_results.py --review results/<dir> --md CHANGELOG.md DECISIONS.md`
  covers checklist items 1, 2, 3, 5 and 6.
- **[manual]** Item 4 (by-construction controls) and the TASKS numbers.

## 2. Outputs and file size
1. Any per-row or per-window output that may exceed 500 kB is written as `.json.gz` (or `.ecsv.gz`)
   from the start, never first as plain JSON.
   **[auto]** The pre-commit hook `results-output-size` (`scripts/check_results.py --sizes`) fails on a
   tracked plain-text file under `results/` above 500 kB. Older files are listed in
   `scripts/check_results_grandfathered.txt`; never add new entries.
2. A raw download-ledger JSONL is gitignored from the moment it is created; only a gzip snapshot is
   committed.
   **[auto]** `.gitignore` covers `net_ledger*.jsonl`, and the same hook fails on a tracked raw
   ledger.
3. Run pre-commit before every commit, on tracked files only (`pre-commit run --all-files`, or
   `--files` listing tracked paths only), and confirm that no file exceeds 1 MB.
   **[auto]** `check-added-large-files` (1 MB), and the 500 kB hook above.

## 3. State files and merging
1. **[manual]** CHANGELOG, TASKS, DECISIONS, SOURCES and `docs/state/trial.md` are edited last, in the
   same commit that goes to review.
2. Before every push, fetch `main`. If it has moved, merge it first, so that the PR does not turn
   "dirty" and stop triggering CI.
   **[auto]** `scripts/push_branch.sh` fetches `main`, merges it when the branch is behind (it stops on
   a conflict), runs pre-commit and pushes.
3. **[manual]** If CI does not appear after a push, check the PR's `mergeable_state` (GitHub MCP
   `pull_request_read get`) instead of waiting.

## 4. Data reuse
1. **[manual]** When a request is sent again and the earlier result is merged and still current, cite
   it; do not re-query the archives. Run one minimal query only when currency must be confirmed, and
   say why.
2. **[manual]** If the next round will probably reuse a light-curve cache, ask the owner whether to keep
   it. Do not delete it by reflex. Large caches are still deleted only after vetting and re-runs are
   done.
3. **[manual]** Skip products that are clearly unused, e.g. Kepler Q0 files. A change to what is
   downloaded goes into that round's frozen plan.

## 5. Test scope
1. **[manual]** Locally, run only the tests and lint of the affected modules; CI runs the full offline
   suite. Run the full suite locally only when a shared module changes (e.g. `scripts/w3_kepler.py`,
   `scripts/w3_ztf.py`, `src/jwst_anomaly/`).
2. **[manual]** If a local dependency that this round does not use is missing and test collection
   fails, exclude that file with `--ignore` and say so in the PR body.

## 6. Not covered here (needs a separate owner authorization)
- Vectorising computations (e.g. screening loops in numpy): a tool change; it needs its own round and
  a re-run comparison.
- A full-pool streaming search, or switching to target-driven star selection: it changes the budget or
  the research question; the owner decides.
