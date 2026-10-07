# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: CHANGELOG 2026-10-08 (run `20261007T130207Z-c8ffc470`).

## Now (M1 follow-through)
1. **Vet the new candidates** with `/vet-candidate` and `scripts/vet_evidence.py` (run `20261007T130207Z-c8ffc470`):
   - catalogued lens-related `jw02736-o001_t001_nircam_f277w_829` (`[YML2023] F150DB-C-4`, a red F150W
     dropout, rank 3) and `jw02736-o001_t001_nircam_f200w_1032` (`[MJR2023] 028.2`, rank 8);
   - arc-like `..._f200w_438` (rank 11);
   - red faint `..._f356w_1243` (`[RBI2023] 18`, rank 12);
   - MIRI's known `G_Lens` object.
2. **Blended stars:** three PSF-like blends (`940`, `2242`, `1571`) remain in the SMACS galaxy top 20
   (D-015 "Revisit if"). Add image-based PSF-spike detection on cutouts.
3. Teach `/vet-candidate` to use `scripts/vet_evidence.py` and the real `jwst-anomaly candidates`
   subcommands.
4. Report the top-k artifact and star fractions per run as a tracked metric.
5. **Per-sample confirmation** (D-014 limitation): only SMACS NIRCam can confirm single-band detections
   and apply the stellar locus. MIRI passes 103 of 530, and CEERS drops every single-band source. Add
   matched photometry or confirmation for each sample (DJA for CEERS, NIRCam counterparts for MIRI).
6. **Control comparability:** CEERS differs from SMACS in colours (aper50), gate (no confirmation) and
   classification (no locus). Science-vs-control comparisons stay invalid until it gets DJA photometry.
7. Diffraction-spike mask for multi-band spike detections, only if they appear (D-014 "Revisit if").

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
