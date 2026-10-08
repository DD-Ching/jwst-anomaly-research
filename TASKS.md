# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: D-023 (owner, 2026-10-08), D-025/D-026 and the CHANGELOG entries of 2026-10-08.

## Now (M3 lensing-violation search, D-023)
1. **Lens-model consistency on SMACS** (D-024: model validated; arc orientations give a null result). Next:
   - counter-images: done for the catalogued systems (D-029: 0 of 11 testable uncatalogued images absent). Left:
     - bright single arcs;
     - a BCG/ICL-subtracted test for system 17;
     - reference fluxes for systems 11, 16 and 26;
   - parity and flux-ratio checks for the catalogued systems (DJA photometry against model magnification ratios);
   - critical curves at z_s = 1, 2, 4 and arc curvature against them (not yet produced);
   - shapes for the `arcs` test: reject blends, require S/N ≥ 50, compare F150W and F444W;
   - the exotic-specific screens of docs/exotic_lensing.md once these ordinary checks are done.
2. **Cluster fields** (done: #30–#32, D-026). Run the lens-model checks per field:
   - El Gordo: Caminha+2023 multiple images and magnification maps (CDS);
   - Abell 2744: UNCOVER v2.0 maps;
   - **Done (D-030):** `validate --model elgordo-caminha23 | abell2744-bergamini23` reproduces Lenstool's χ²
     (82.53 vs 80.22; 146.64 vs 146.60). No image-position anomaly. Next:
     - `lens_consistency.py images --forced-image` per field (El Gordo: F277W `jw01176-o241_t012`; Abell 2744
       needs DJA photometry, the 233 MB catalogue, with a DECISIONS entry, because pipeline catalogs miss core arcs);
     - search radii from `bayes.dat` instead of a fixed 1″;
     - Abell 2744 multiplicity residual (D-030): do `bayes.dat` samples split 3.2a/b, 34.1a/b, 700.1a/b?
     - El Gordo magnification-map check (CDS `magnification_best_fit_z2.fits`) as a `validate` map test.
   - Sunrise: RELICS or Scofield+2025.
   Field follow-ups (docs/fields/*.md):
   - vet Abell 2744 `5904`, `7987`, `4731`, `264`, `7298`;
   - vet El Gordo `2828` (a linear feature) and `1601`;
   - vet Sunrise `1756`, `1907`, `1578`;
   - Sunrise `1869` is a two-band detection on a star spike (D-014 spike mask);
   - Sunrise has a rel_weight plateau at 0.45–0.55 (1,909 sources), so the D-011/D-021 0.5 cut splits it on noise;
   - duplicate merged sources (Abell 2744 `3999`/`7179`, 0.25″ apart);
   - extended-veto fallback for sources without DJA r50 (Abell 2744 `7694`);
   - El Gordo module A is a flanking field;
   - reword `acquire._warn_if_reprocessed`.
3. **Two-epoch search** (D-027): null results for SMACS/VENUS and Sunrise o010/o120; Earendel is steady
   (docs/fields/sunrise.md). Next:
   - divide forced-photometry significances by the control std (1.2–1.5) before thresholding;
   - add a third epoch: Sunrise VENUS 6882 o052 (F150W, F444W, jwst 3.0.0; forced photometry only);
   - vet `c0049`'s epoch-2 streak, if it recurs elsewhere (a satellite or asteroid trail, or scattered light).

## Then (M1 follow-through, after the M3 items)
4. **Star-stratum sources with NED galaxy matches** (run `20261007T162133Z-2511977c`). Check their classification and photometry.
   - `f090w_1719`: NED `[NDA2023] 00713`, type G, 0.03″.
   - `f200w_2054`: a stellar-locus member 0.23″ from the lensed image NED `[MJR2023] 002.3`; it has spikes, and an
     arc lies next to it.
5. **Vet MIRI's known `G_Lens` object** with `/vet-candidate` and `scripts/vet_evidence.py`.
6. **Spike streaks confirmed by a coincident DJA object** (D-014 "Revisit if"; run `20261007T202419Z-65728e80`): CEERS
   `f200w_2842` (#4) is detected in F200W only, with ellipticity 0.95 and area 376 px. A bright source lies 4.6″
   away, and a DJA object sits 0.04″ off. Test single-band, highly elongated detections near bright sources
   (a spike-geometry check), without dropping lensed arcs.
7. **MIRI confirmation** (D-014 limitation): MIRI passes 103 of 530. Confirm MIRI detections with NIRCam
   counterparts. SMACS and CEERS now both have DJA confirmation.
8. **Control comparability, what remains** (D-022): pipeline-catalog features are jwst 3.0.0 in CEERS and 2.0.1
   in SMACS, until 2736 is reprocessed. CEERS has no stellar locus (too few stars). DJA colours are still to be
   validated against CEERS DR1.0. Then compare the science and control feature distributions.
9. Diffraction-spike mask for multi-band spike detections, only if they appear (D-014 "Revisit if").

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
- Second epochs and medium bands: check the galaxy top k against VENUS (6882) and other later programs
  with `scripts/epoch_compare.py` (D-017). Moving sources are stars; variable ones are AGN or transients.
- NIRCam F444W grism spectra of program 4043 (SOURCES) for top-k sources they cover: redshifts and lines. Needs
  a reuse check (grizli or the jwst WFSS pipeline products).

## Owner setup (needs-human)
- Run `claude.exe` interactively in the repo once and accept workspace trust.
- `claude.exe setup-token`, then `gh secret set CLAUDE_CODE_OAUTH_TOKEN`; install the Claude GitHub App;
  test with an `@claude` issue (docs/operations.md §5).
- Cloud routine: live since 2026-10-08 (docs/cloud-routine-prompt.md). Check that the environment's Custom network
  allowlist has every host in docs/operations.md §3. Merge any PR labelled `merge-ready` (cloud runs merge with the
  GitHub MCP tool when they can).
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
