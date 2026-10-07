# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: CHANGELOG 2026-10-08 (run `20261007T202419Z-65728e80`).

## Now (M1 follow-through)
1. **Vet `jw02736-o001_t001_nircam_f200w_2804`** (galaxy rank 6 in run `20261007T144011Z-11342d69`). It is point-like (r50 2.13 px,
   `mag_auto` 23.4) with F150W−F444W +3.7 and very red F356W−F444W, a brown-dwarf-like signature (D-016).
   Check proper motion, published brown-dwarf lists and the F356W methane band.
2. **Star-stratum sources with NED galaxy matches** (run `20261007T162133Z-2511977c`). Check their classification and photometry.
   - `f090w_1719`: NED `[NDA2023] 00713`, type G, 0.03″.
   - `f200w_2054`: a stellar-locus member 0.23″ from the lensed image NED `[MJR2023] 002.3`; it has spikes, and an
     arc lies next to it.
3. **Vet MIRI's known `G_Lens` object** with `/vet-candidate` and `scripts/vet_evidence.py`.
4. **Spike streaks confirmed by a coincident DJA object** (D-014 "Revisit if"; run `20261007T202419Z-65728e80`): CEERS
   `f200w_2842` (#4) is detected in F200W only, with ellipticity 0.95 and area 376 px. A bright source lies 4.6″
   away, and a DJA object sits 0.04″ off. Test single-band, highly elongated detections near bright sources
   (a spike-geometry check), without dropping lensed arcs.
5. **MIRI confirmation** (D-014 limitation): MIRI passes 103 of 530. Confirm MIRI detections with NIRCam
   counterparts. SMACS and CEERS now both have DJA confirmation.
6. **Control comparability, what remains** (D-022): pipeline-catalog features are jwst 3.0.0 in CEERS and 2.0.1
   in SMACS, until 2736 is reprocessed. CEERS has no stellar locus (too few stars). DJA colours are still to be
   validated against CEERS DR1.0. Then compare the science and control feature distributions.
7. Diffraction-spike mask for multi-band spike detections, only if they appear (D-014 "Revisit if").

## Next (M2)
- Matched-aperture colours for other samples (D-013 is done for SMACS NIRCam):
  - CEERS: done (D-022). Validation against CEERS DR1.0 is still open.
  - MIRI: evaluate fixed-aperture colours.
- DJA photo-z as a ranking feature or stratifier (it is used for vetting since cycle 7).
- Re-fetch program 2736 after MAST reprocesses it with jwst ≥ 3.0 (expected around mid to late October 2026).
  Until then, compare only photometric columns across programs (D-010).
- Image embeddings (Zoobot via `timm`, DINOv2) on cutouts, evaluated against the baseline with
  injection-recovery. Keep the evaluation fields out of any training data.
- Crossmatch and classify: propagate Gaia positions to the JWST epoch, so fast-moving stars and brown
  dwarfs land in the star stratum (D-012). Identify lenses by name or literature, not only by otype.
- Cutouts: north-up panels and an option for shared brightness scaling across panels.
- NIRCam F444W grism spectra of program 4043 (SOURCES) for top-k sources they cover: redshifts and lines. Needs
  a reuse check (grizli or the jwst WFSS pipeline products).

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
