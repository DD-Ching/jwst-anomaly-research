# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: CHANGELOG 2026-10-07 (run `20261007T020124Z-4bfabaaa`).

## Now (M1 follow-through)
1. **Quality gating before ranking.** Exclude or separately rank sources on low weight or the image edge
   (i2d `WHT` via S3 byte-range reads, or the MAST `s_region` footprint) and hot-pixel-like sharpness
   extremes. Metric: the artifact fraction of the top 20 per sample, compared with the run above.
2. **Star/galaxy separation.** Use Gaia plus `CI`/`is_extended` and rank the two populations separately,
   because bright stars dominate the SMACS NIRCam top 20.
3. **First vetting with `/vet-candidate`:**
   - SMACS NIRCam ranks 6 and 19 (`jw02736-o001_t001_nircam_f200w_2925`, `_2559`): elongated and arc-like;
     compare them with the published SMACS 0723 lens models (docs/landscape.md).
   - MIRI rank 17 (`jw02736-o002_t001_miri_f770w_94`): the known NED `G_Lens` object, which exercises the
     known-object path.
4. Align the `/vet-candidate` skill with the real `jwst-anomaly candidates` subcommands and confirm
   `docs/candidates/` as the place for vetting records.
5. Report the top-k artifact and star fractions per run as a tracked metric.

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
