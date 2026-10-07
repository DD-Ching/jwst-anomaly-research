# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.

## Now: M0 bootstrap batch (2026-10-07), one PR per unit
1. Archive access: `query.py`, `acquire.py`, manifests, `scripts/fetch_reference_sample.py`
2. Catalog ingestion: `catalog.py` (pipeline `_cat.ecsv` → cross-band source table)
3. Features + baseline ranking: `features.py`, `rank.py` (+ injection-recovery)
4. Cutouts + visualization: `cutouts.py`, `viz.py` (S3 byte-range cutouts)
5. External cross-check: `crossmatch.py` (SIMBAD/NED/Gaia/XMatch)
6. Candidate store + runner + CLI: `candidates.py`, `provenance.py`, `pipeline.py`, `cli.py`
7. OSS hygiene: CITATION.cff, CONTRIBUTING, templates, pre-commit, Dependabot
8. Agent operations: `.claude/skills/*`, `.claude/settings.json`, `claude.yml`, `docs/operations.md`
9. Reuse landscape survey: `docs/landscape.md`

## Next
- Integration: merge the bootstrap units, run `configs/reference_sample.yaml` on real data, and record results (M1).
- Owner setup: GitHub secret `CLAUDE_CODE_OAUTH_TOKEN`; cloud routine with a network allowlist (docs/operations.md).

## Later
- Forced photometry or HLSP catalogs for consistent colors (see D-001 "Revisit if").
- Zenodo DOI on first tagged release.
