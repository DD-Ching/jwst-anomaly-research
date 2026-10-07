<!--
Every change goes through a PR with green CI. External PRs need the owner's approving review; the
owner's agents merge their own PRs under CLAUDE.md's merge policy. AI-agent PRs carry the `agent` label;
add `needs-human` when a decision is scientific, irreversible, costly, legal or credential-related.
See CONTRIBUTING.md.
-->

## Summary

<!-- What changed and why. Link the issue or task, e.g. "Closes #12". -->

## Decision

<!-- The D-NNN entry in DECISIONS.md that this PR adds or relies on, or "None". -->

## Evidence

<!-- Trimmed test/lint output, run IDs, figures. Label claims: observed / derived / simulated /
model_prediction / assumption / hypothesis (docs/methodology.md). -->

## Limitations

<!-- Known gaps, untested paths, caveats. -->

## Follow-ups for TASKS.md

<!-- New tasks, owner actions, or state-file updates for the coordinator to fold in. -->

## Checklist

- [ ] Offline tests pass: `python -m pytest -q`
- [ ] Network tests pass, or N/A (needed when archive or catalog access changed): `python -m pytest -q --run-network`
- [ ] `ruff check src tests scripts` and `ruff format --check src tests scripts` pass, as in CI. CI is authoritative, because pre-commit's pinned ruff can lag behind CI's latest.
- [ ] New tables set `meta["provenance"]` and `meta["source"]`; new claims and figures carry a provenance label
- [ ] Reuse decisions recorded in DECISIONS.md (`D-NNN`); new external sources in SOURCES.md (URL, version/DOI, access date)
- [ ] No data files and no file over 1 MB; no secrets, tokens or personal data
- [ ] State files updated where applicable (TASKS.md, CHANGELOG.md, or listed under Follow-ups)
- [ ] If a public stage signature changed: docs/architecture.md and `pipeline.py` updated in this PR
- [ ] No anomaly score is presented as evidence of new physics; every citation was checked at its source
