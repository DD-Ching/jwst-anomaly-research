# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: CHANGELOG 2026-10-08 (run `20261007T162133Z-2511977c`).

## Now (M1 follow-through)
1. **Vet `jw02736-o001_t001_nircam_f200w_2804`** (galaxy rank 6 in run `20261007T144011Z-11342d69`). It is point-like (r50 2.13 px,
   `mag_auto` 23.4) with F150W−F444W +3.7 and very red F356W−F444W, a brown-dwarf-like signature (D-016).
   Check proper motion, published brown-dwarf lists and the F356W methane band.
2. **Spike-flagged sources in the galaxy ranking** (D-018 "Revisit if"): `940`, `2242` and `1571` are now flagged
   `spikes` but keep their ranks. Decide whether to backfill the galaxy top k without them, or move them to the star
   stratum. Then track the top-k contamination fraction.
3. **Star-stratum sources with NED galaxy matches** (run `20261007T162133Z-2511977c`). Check their classification and photometry.
   - `f090w_1719`: NED `[NDA2023] 00713`, type G, 0.03″.
   - `f200w_2054`: a stellar-locus member 0.23″ from the lensed image NED `[MJR2023] 002.3`; it has spikes, and an
     arc lies next to it.
4. **Vet MIRI's known `G_Lens` object** with `/vet-candidate` and `scripts/vet_evidence.py`.
5. **Track metrics per run:** the fractions of top k that are artifacts, stars and catalogued objects.
   Recovering catalogued lensed images and high-z candidates validates the ranking; uncatalogued high-rank
   sources (like `438`) form the discovery set.
6. **Per-sample confirmation** (D-014 limitation): only SMACS NIRCam can confirm single-band detections
   and apply the stellar locus. MIRI passes 103 of 530, and CEERS drops every single-band source. Add
   matched photometry or confirmation for each sample (DJA for CEERS, NIRCam counterparts for MIRI).
7. **Control comparability:** CEERS differs from SMACS in colours (aper50), gate (no confirmation) and
   classification (no locus). Science-vs-control comparisons stay invalid until it gets DJA photometry.
8. Diffraction-spike mask for multi-band spike detections, only if they appear (D-014 "Revisit if").

## Next (M2)
- Matched-aperture colours for other samples (D-013 is done for SMACS NIRCam):
  - CEERS: DJA `ceers-full` (250–400 MB, so state the reason before downloading), validated against
    CEERS DR1.0.
  - MIRI: evaluate fixed-aperture colours.
- DJA photo-z as a ranking feature or stratifier (it is used for vetting since cycle 7).
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
