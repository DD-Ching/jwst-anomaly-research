# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: CHANGELOG 2026-10-07 (runs `20261007T033832Z-cfcf6032`, `20261007T035338Z-8ea21f89`).

## Now (M1 follow-through)
1. **First vetting with `/vet-candidate`:**
   - SMACS NIRCam `jw02736-o001_t001_nircam_f200w_2925`, `_2559`, `_1096` (galaxy stratum, run
     `20261007T035338Z-8ea21f89`): elongated, arc-like; compare them with the published SMACS 0723
     lens models (docs/landscape.md).
   - MIRI rank 17 (`jw02736-o002_t001_miri_f770w_94`): the known NED `G_Lens` object, which exercises the
     known-object path.
2. **Faint stars:** about 4 PSF-like sources without Gaia counterparts remain in the SMACS galaxy top 20.
   Add a size–magnitude stellar locus to `classify` (D-012 "Revisit if").
3. Align the `/vet-candidate` skill with the real `jwst-anomaly candidates` subcommands and confirm
   `docs/candidates/` as the place for vetting records.
4. Report the top-k artifact and star fractions per run as a tracked metric.
5. Quality gate follow-ups (D-011 "Revisit if"): per-band gating for colours; rank shallow regions as a
   separate stratum instead of excluding them (55% of MIRI sources are excluded now).

## Next (M2)
- Ingest the DJA v7.4 `smacs0723` and `ceers-full` catalogs and photo-z for consistent matched-aperture
  photometry. Validate CEERS colours against CEERS DR1.0. The `ceers-full` files are 250–400 MB, so state
  the reason before downloading them.
- Re-fetch program 2736 after MAST reprocesses it with jwst ≥ 3.0 (expected around mid to late October 2026).
  Until then, compare only photometric columns across programs (D-010).
- Image embeddings (Zoobot via `timm`, DINOv2) on cutouts, evaluated against the baseline with
  injection-recovery. Keep the evaluation fields out of any training data.
- Crossmatch: propagate Gaia positions to the JWST epoch; identify lenses by name or literature, not
  only by otype.
- Cutouts: north-up panels and an option for shared brightness scaling across panels.

## Owner setup (needs-human)
- Run `claude.exe` interactively in the repo once and accept workspace trust.
- `claude.exe setup-token`, then `gh secret set CLAUDE_CODE_OAUTH_TOKEN`; install the Claude GitHub App;
  test with an `@claude` issue (docs/operations.md §5).
- Cloud routine: `/schedule` with the Custom network allowlist (docs/operations.md §3).
- Zenodo DOI at the first tagged release; optionally a separate Code of Conduct contact.

## Later
- M3: published SMACS 0723 lens models (Mahler+2022, RELICS, Caminha+2022); AnomalyMatch for lens finding.
- A weekly scheduled link check of the state files (`scripts/check_links.py`).
- Switch Dependabot from pip to uv once a `uv.lock` exists. Keep the version in `pyproject.toml`,
  `__init__.py` and `CITATION.cff` in sync on release.
- NIRISS imaging of program 2736 (F115W, F200W).
