# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: CHANGELOG 2026-10-07 (run `20261007T120528Z-378af5a1`; size bias via
`scripts/feature_size_bias.py`).

## Now (M1 follow-through)
1. **S/N floor for ranking:** with size-robust colours (D-013), faint near-noise sources reach the SMACS
   galaxy top 20 (run `20261007T120528Z-378af5a1`). Rank only sources above a reference-band S/N (rank the rest separately or
   not at all), then re-check the top 20.
2. **Faint stars:** about three PSF-like sources without Gaia counterparts remain in the galaxy top 20. Add
   a size–magnitude stellar locus to `classify` (D-012 "Revisit if").
3. **Vet the new top candidates** with `/vet-candidate` and `scripts/vet_evidence.py`:
   - thin arcs `jw02736-o001_t001_nircam_f200w_1032` (NED `[MJR2023] 028.2`, G_Lens) and `_438`;
   - `jw02736-o001_t001_nircam_f277w_829` (SIMBAD `[YML2023] F150DB-C-4`);
   - MIRI rank 17 (`jw02736-o002_t001_miri_f770w_94`), the known NED `G_Lens` object, which exercises the
     known-object path.
4. Teach `/vet-candidate` to use `scripts/vet_evidence.py` and the real `jwst-anomaly candidates`
   subcommands (`add-vetting --outcome pass|fail|inconclusive`).
5. Report the top-k artifact and star fractions per run as a tracked metric.
6. Quality gate follow-ups (D-011 "Revisit if"): gate each band for colours; rank shallow regions as a
   separate stratum instead of excluding them.

## Next (M2)
- Matched-aperture colours for other samples (D-013 is done for SMACS NIRCam):
  - CEERS: DJA `ceers-full` (250–400 MB, so state the reason before downloading), validated against
    CEERS DR1.0.
  - MIRI: evaluate fixed-aperture colours.
- DJA photo-z (`smacs0723-grizli-v7.4-fix.photoz.tar.gz`, 59 MB) as a feature or for vetting.
- Re-fetch program 2736 after MAST reprocesses it with jwst ≥ 3.0 (expected around mid to late October 2026).
  Until then, compare only photometric columns across programs (D-010).
- Image embeddings (Zoobot via `timm`, DINOv2) on cutouts, evaluated against the baseline with
  injection-recovery. Keep the evaluation fields out of any training data.
- Crossmatch and classify: propagate Gaia positions to the JWST epoch, so fast-moving stars and brown
  dwarfs land in the star stratum (D-012). Identify lenses by name or literature, not only by otype.
- Cutouts: north-up panels and an option for shared brightness scaling across panels.

## Owner setup (needs-human)
- Run `claude.exe` interactively in the repo once and accept workspace trust.
- `claude.exe setup-token`, then `gh secret set CLAUDE_CODE_OAUTH_TOKEN`; install the Claude GitHub App;
  test with an `@claude` issue (docs/operations.md §5).
- Cloud routine: `/schedule` with the Custom network allowlist (docs/operations.md §3).
- Zenodo DOI at the first tagged release; optionally a separate Code of Conduct contact.

## Later
- Reuse `classify`'s full-catalog SIMBAD/Gaia matches for the top-k cross-match (re-query only NED).
- M3: published SMACS 0723 lens models (Mahler+2022, RELICS, Caminha+2022); AnomalyMatch for lens finding.
  Evaluate model shear and magnification at candidate positions (open question in
  `docs/candidates/jw02736-o001_t001_nircam_f200w_1096.md`).
- A weekly scheduled link check of the state files (`scripts/check_links.py`).
- Switch Dependabot from pip to uv once a `uv.lock` exists. Keep the version in `pyproject.toml`,
  `__init__.py` and `CITATION.cff` in sync on release.
- NIRISS imaging of program 2736 (F115W, F200W).
