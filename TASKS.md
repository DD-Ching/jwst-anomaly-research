# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: D-023 (owner, 2026-10-08), D-025/D-026 and the CHANGELOG entries of 2026-10-08.

## Now (owner direction 2026-10-09: System A/B hypothesis rounds, D-069; docs/hypotheses/)
0. **E1 causal event network** (D-074; owner idea 4): GBM × ICECAT-1 × GWTC × CHIME Cat 2 pair counts by lag and
   separation. The five requested channels are null (family p = 0.28). The CHIME–CHIME 1 h–1 d wide excess (z ≈ 3.9)
   is reproduced by a calibrated rate-modulated null (week-scale detection-rate modulation; its cause is a
   hypothesis). The 100 s–1 h wide cell drops to z = 2.5 under it (family-wise p = 0.26, D-074 addendum). Next:
   - GW sky maps, so GW channels get directions;
   - a time-resolved CHIME uptime series (only from the collaboration: owner decision); it would also decide
     the 100 s–1 h residual;
   - signed-lag ("which comes first", addendum 3) and antipodal (addendum 2) channels are null; redo both with GW sky maps;
   - IceTracks-DR2, Swift, Einstein Probe.
0. **D1 distance self-consistency** (D-073 + addenda; owner idea 1): H0LiCOW (R and D_dt LOO), 94 localized FRBs
   (NE2001 and YMW16) and TDCOSMO 2025 power-law D_dt LOO on 8 lenses (addendum 2: max 1.12σ with κ_ext): all null.
   Reach: R ×2–8 / ×0.2–0.4; D_dt ×1.36–1.85 / ×0.41–0.78. Next: composite-model D_dt chains as a model-choice check;
   new per-lens D_d needs a numpy-only reader for the hierArc `*_const_processed.pkl` likelihoods.
1. **S1 burst twins** (D-071): Fermi GBM null, < 7.2 × 10⁻³ twin pairs per eligible burst (ratio 1, 15 % band
   mismatch). Next: CHIME/FRB Catalog 2 (CANFAR downloads reset by the proxy; try another route), TTE for short and
   faint bursts, an s ≠ 1 chain, a generative pulse-model null.
2. **S2 flat-kernel SN residuals** (D-070 + addendum): Pantheon+ × LS DR9 null (γ < 0.025), DES-SN5YR × LS DR9 null
   (γ < 0.019 mag per unit T/⟨T⟩, α = 2, one-sided 95 %, dilution- and chain-corrected); α = 1 is kernel-degenerate. The
   lensing positive control is only ~1.7σ (DES) / ~2σ combined. Next: a deeper galaxy column (DES Y3 Gold or LS DR10
   z < 22) so the control detects, and the full Pantheon+ / DES-SN5YR covariances.
3. **S3 hybrid images** (D-072): COSMOGRAIL doubles null, sensitivity only r ≈ 0.4–2.4 (red quasar variability makes
   the copy collinear with the main image). Next, low priority: SN Refsdal / SN H0pe imaging at the model lags.
4. **System A rounds paused** (D-076: rounds 2–3 gave no testable survivor). Restart only with a new owner brief.
   Optional: the B0218+357 VLA polarisation sign test (R2-B).

## Continuing (owner direction 2026-10-08: evidence of wormholes / negative mass / warp in any public dataset; D-054)
1. **W3 in published microlensing samples** (D-057–D-059, D-061). Mróz 2019 bulge (5,790) and 2020 disk (460) fitted and
   vetted: 127 + 6 flags, **0 survive** (the D-059 chunk survivors BLG667.04.62161 and BLG624.18.69573 fail
   `feature_coverage`). **The published samples cannot limit W3:** 0 of 600 injected W3 events pass the emulated Mróz
   selection (PSPL controls 15–43 %), so the chunk re-fits only re-test a selection that excludes the signal. Next:
   light curves from **before** a PSPL selection. **MOA-II 9-year Cut-0 light curves keep W3 (D-062):** gb22
   null (30 flags, 0 survive), Γ₉₅ ≈ 1.2–4.4 × 10⁻⁶ per star per yr (t_E 10–300 d; 0.6–1.1 × 10⁻⁵ at 3 d). Calibrate the variable-baseline test (largest efficiency loss) and add injections per cell. Next: a field with published N_s
   (Nunota et al. 2024 Table 1; light-curve tar size first), then more fields to scale the exposure. Gaia DR3 `vari_microlensing` is ruled out (D-061: 0 / 240 W3
   injections both selected and flagged; one real flag, not a candidate). Next: Gaia DR3 variables with epoch
   photometry (vari_summary; check by injection whether W3 survives the variability classifier first), KMTNet
   public seasons, OGLE EWS (owner decision, terms). Optional:
   re-fit the bulge under D-058's bounded π_E (`merge-chunks`, then `vet` + `revet`; can only add flags); re-run the
   6 arXiv name queries that errored.
2. **W1/W2 in wide imaging** (D-056: published lens catalogues null; f_dark < 0.31 (typical, 3/25 after the 2026-10-08
   pair-match and one-entry-per-lens amendment) for quasar/radio-selected lenses, the only selections sensitive to a
   dark lens). Next: HSC PDR3 photometry or PSF-subtracted HST
   image models for the 307 blended or too-close lensed quasars (makes the test decisive; HST *catalogue* photometry
   is not: D-060, efficiency 0.46); the CHITAH lens models of the 3 open SuGOHI IX pairs; rejected lensed-quasar pairs: LS DR10 cannot decide them
   (D-064, control efficiency 0/5 at 1.9–2.6″); of the 11 colour-matched ones, J0130+0725 has no lens light in HST to
   F814W ≈ 23 (D-064 addendum: binary vs dark lens needs spectra of both images or two-epoch flux ratios), J0728+2607
   has none either with an empirical PSF (F814W ≈ 23), and 9 have no HST (Euclid DR1 when public); catalogued image separations for the 15
   pair-decided systems the pair check cannot test yet (incl. the 3 CHITAH pairs; HSC lens models); HST lens photometry (CASTLES-type) for the blended quasar systems; ALMA positions for submm systems;
   deeper imaging (HSC PDR, Euclid DR1) for the conservative variant; a lens list with image positions for W1.
3. **W5** (D-063): Legacy Surveys DR10, 340.5 deg², 40 flags, 0 survivors; n₉₅ ≈ 0.012–0.018 deg⁻² at θ_E = 8–32′.
   Euclid Q1 counts would not open θ_E < 6′ (D-065: clustering-limited, ×1.0–1.4 S/N on a fifth of the area).
   **Euclid Q1 radial-shear screen (D-066, D-067):** all three Deep Fields (60 deg²) null; R = 0.56 ± 0.16 from SZ
   clusters (R = 0.5 kept); n₉₅ ≈ 0.049 / 0.051 deg⁻² at θ_E = 2′ / 4′; 1′ not limited (efficiency 0.43–0.58). Next:
   per-tile star-ellipticity gradient test; Euclid DR1 when public; θ_E ≤ 1′ and ≈ 1° are still open.
4. **Warp:** monitor the literature (checked 2026-10-08: no distant-observer EM template — Lentz & Felton 2024 give fluxes only; Clough et al.
   2024 waveform not public; O3 superluminal-burst search already null, Kuwahara & Cannon 2023). Recheck monthly.

## JWST focus (D-047 screens; continues under the direction above)
1. **Injection-recovery, then limits** (Phase 2; `exotic_sim.inject_images` / `inject_light_curve`):
   - W1: the shear screen (D-053, `exotic_screens.py shear`, `inject_shear.py`) is built; four fields null; limits
     7.7 × 10³ / 1.3 × 10³ / 6.5 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉, still blind at ≤ 2 × 10¹¹ M☉. Next:
     lower the B-mode floor (Abell 2744 max S_× 4.21: PSF-anisotropy model from stars, blend rejection, drop edge
     apertures with < 50 % coverage); add SMACS 0723 and El Gordo when DJA photo-z is reachable (tarballs 404
     2026-10-08) and MACS0717 / Abell S1063 with photo-z (DJA v7.5, 75 MB); stack S around `radial` and orphan-pair
     centres; other lens redshifts → volume density vs Takahashi & Asada;
   - W2 / dark deflectors (D-048 clusters, D-051 deep fields: null, limits; D-055: the flanking-field orphan
     excess was a null (e) artefact, no excess under the fixed (e) or the companion-aware (f); D-048 clusters re-run, null): segmentation-map
     adjacency for the same-galaxy and visible-lens rules (the main efficiency losses); pixel-level injections;
     optionally CEERS / GOODS-S / PRIMER (~600 MB each);
   - W3 (D-052, null; headline limits from MACS0416 only): re-run `forced` for all three fields when S3 cutout jobs
     work, ≥ 300 multi-epoch controls so NEXUS and Abell 2744 become calibrated; replace the single-epoch veto for
     vanish flags with the D-039 persistence test (restores full-vanish sensitivity); add JADES and new NEXUS epochs;
     SN/TNS check for any survivor; leave the Sunrise transient track to its owner run;
   - W5: counts N(>S) around `radial` centres.
2. W3 inside caustic-crossing arcs needs a microlens with macro shear (Chang-Refsdal-type; reuse-check first).
3. Warp: recheck only when a paper gives an electromagnetic prediction for a distant observer.

## Supporting (M3 lensing-violation search, D-023; serves the focus above)
1. **Lens-model consistency on SMACS** (D-024: model validated; arc orientations give a null result). Next:
   - counter-images: done for the catalogued systems (D-029: 0 of 11 testable uncatalogued images absent). Left:
     - bright single arcs;
     - a BCG/ICL-subtracted test for system 17;
     - reference fluxes for systems 11, 16 and 26;
   - flux ratios and colours: done for SMACS and El Gordo (D-032: no anomaly; most images untested). Left:
     - Abell 2744 `fluxratios` (needs the 233 MB DJA catalogue, see below);
     - totals for arcs inside cluster-galaxy/ICL light (inspect DJA `_phot_apcorr.fits`);
     - `bayes.dat` μ errors, then recheck SMACS 6.3 (0.6–0.7 mag brighter than 6.1+6.2 predict);
     - an observational parity test (image orientation or resolved structure; not yet designed);
   - critical curves at z_s = 1, 2, 4 and arc curvature against them (not yet produced);
   - shapes for the `arcs` test: reject blends, require S/N ≥ 50, compare F150W and F444W;
   - exotic screens (D-031, D-034, D-036–D-038, D-040–D-043): SMACS, El Gordo, Sunrise (radial), Abell 2744,
     MACS1149, MACS0717, MACS0416, Abell 370 and Abell S1063 (radial) null. All six HFF clusters are screened. Next:
     - multi-plane checks are in the library (D-046, `LensModel.split_planes`); a `--plane NAME=Z[:SIGMA]` option for
       `lens_consistency.py images/validate` is next;
     - orphan pairs (D-048, null): injection-recovery of simulated dark-lens pairs (limit per deg²); CANUCS segmentation
       maps for the same-galaxy test; blending model for pairs < 1″; re-run on Abell 2744 (UNCOVER), MACS0417, MACS1423;
       MACS1149 image list when CANUCS v2 is out;
     - CANUCS models pinned (D-044: `macs0416-canucs`, `abell370-canucs`); next, the 100 MCMC sample maps for μ spreads;
     - `radial`: low-weight veto (relative WHT < 0.5) or aper50 S/N floor; pass `--spike-stars` (Gaia, D-043) on every
       field and re-run earlier fields' radial screens with it;
     - re-run MACS0717 radial when a photo-z catalogue exists;
     - `forced_check`: match and search radii scaled to the image-plane rms when it exceeds 1″; classify model copies of
       unpredicted catalogued images automatically (D-041);
     - second models as pinned `MapLensModel` entries: Sharon v4cor (MACS1149, MACS0717), Keeton v4 (MACS0717);
     - MACS1149: `macs1149-sharon` in `MODELS`; the D-037 rules in `forced_check`; re-run `images` with CANUCS DR1
       photo-z; fix system 16's z (spectroscopic z, or forced photometry on its far-image track); evaluate the CANUCS
       lens models (see the CANUCS item above);
     - `bayes.dat` position spreads (`posterior`, D-045) as `forced_check` search radii; the UNCOVER v2.0 cross-check
       for Abell 2744.
2. **Cluster fields** (done: #30–#32, D-026). Run the lens-model checks per field:
   - El Gordo: Caminha+2023 multiple images and magnification maps (CDS);
   - Abell 2744: UNCOVER v2.0 maps;
   - **Done (D-030):** `validate --model elgordo-caminha23 | abell2744-bergamini23` reproduces Lenstool's χ²
     (82.53 vs 80.22; 146.64 vs 146.60). No image-position anomaly. Next:
     - `images --forced-image`: done for El Gordo (D-034) and Abell 2744 (D-036); 700.1 needs a finer z scan;
     - search radii from `bayes.dat` instead of a fixed 1″;
     - ~~Abell 2744 multiplicity residual~~ — done (D-045): model resolution at folds. Next: `bayes.dat` μ errors in
       `fluxratios` (Abell 2744 and El Gordo chains validated; El Gordo's chain uses σ² = a·b, D-045);
     - ~~El Gordo magnification-map check~~ — done (issue #68): `validate` reproduces the z=2 and z=4 maps (median
       ratio 1.00001, parity 100 % at |μ|<10), so El Gordo μ and parities are validated. Its chain is validated too (`forme -10`, D-045).
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
   (docs/fields/sunrise.md). Third epoch (o052, calibrated with `--controls`): 0 of 57 candidates. Next:
   - `n0022` and `n0150` are detector persistence (D-039). Run `scripts/persistence_check.py` on every future
     single-epoch candidate before vetting. Next: `/vet-candidate n0153` (on sky in all 8 dithers, F150W ×1.9 over
     2.9 yr; host, AGN colours, a pipeline-version check on neighbours);
   - forced photometry on all catalogued sources, not only catalog-stage candidates (the two-band rule misses
     blue transients);
   - vet `c0049`'s epoch-2 streak, if it recurs elsewhere (a satellite or asteroid trail, or scattered light).

- **Speed** (owner focus, step 2): `find_images` is vectorised. Still open: evaluate published deflection maps
  (UNCOVER v2.0, RELICS, HFF) for fields without a Lenstool `best.par`.

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
