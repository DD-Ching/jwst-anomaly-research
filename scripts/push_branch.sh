#!/usr/bin/env bash
# Push the current branch after syncing with main (docs/efficiency-rules.md section 3.2).
# Fetches origin/main; if the branch is behind, merges it (stops on a conflict, never rebases);
# runs pre-commit on tracked files; then pushes with -u. Retries the push on network errors.
set -euo pipefail
branch=$(git rev-parse --abbrev-ref HEAD)
if [ "$branch" = "main" ]; then echo "refusing to push main" >&2; exit 1; fi
git fetch -q origin main
if ! git merge-base --is-ancestor origin/main HEAD; then
  echo "main has moved: merging origin/main into $branch"
  git merge --no-edit origin/main || { echo "merge conflict: resolve, commit, re-run" >&2; exit 1; }
fi
pre=".venv/bin/pre-commit"; [ -x "$pre" ] || pre="pre-commit"
"$pre" run --all-files
for wait in 0 2 4 8 16; do
  sleep "$wait"
  if git push -u origin "$branch"; then exit 0; fi
done
echo "push failed after retries" >&2
exit 1
