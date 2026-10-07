# Changelog

Newest first. Results, failed approaches and the handoff state; not a diary.

## 2026-10-08: Sunrise two-epoch search, a null result; Earendel steady (D-027 amendment)
- WHL0137 2282 o010 against o120 (same pipeline): 372 catalog candidates → 58 in ≥ 2 bands → 52 after the bright-Gaia
  mask → 0 by recentred forced photometry in all four bands. 2 pass in one or two bands, both ordinary on visual check
  (an epoch-2 streak; galaxy outskirts at the footprint edge). docs/fields/sunrise.md.
- Earendel: |Δm| ≤ 0.18 mag (≤ 1.5σ) in F090W/F115W/F277W/F356W over 164 days.
- Noise calibration on 150 ordinary sources: ERR-based significances are 1.2–1.5× too large.
- New: `scripts/transient_combine.py` (band coincidence and bright-Gaia mask) and
  `transient_forced.py --recentre-arcsec`.
- **Failed approach:** a fixed aperture at a literature position. Earendel seemed to brighten by 0.76 mag at 5.7σ,
  because the aperture sat 0.15″ off the source and the PSF rotates by about 180° between the epochs.
- **Handoff:** the lens-model PR #35 (another run was working on it), then lens-model checks per field.
  Time domain: rescale the thresholds by the control std, and add VENUS 6882 o052 (F150W, F444W) as a third epoch.
## 2026-10-08: Hourly cloud routine
- Routine `trig_01PNAmgcfqef8CvhPAY8ggbP` runs `/research-cycle` loops hourly at :07 UTC (Opus 5.5, no connectors)
  with the prompt in docs/cloud-routine-prompt.md.
- First test run: clone, venv, skills and MAST downloads work.
  - GitHub GraphQL is blocked in cloud sessions, so `gh pr` / `gh issue` fail. The REST equivalents and a single
    background CI-wait loop are in docs/operations.md §3, and the research-cycle skill points to them.
  - The merge API is denied by `.claude/settings.json`, so cloud runs stop at merge-ready with the label
    `merge-ready`. The owner or a local session merges.
  - Coordination: runs last about 40 minutes and skip `local-wip` PRs and branches committed to in the last 15
    minutes. An unmerged PR's handoff is read from its branch.
- **Handoff:** unchanged (TASKS "Now (M3)"); the cloud routine continues it.

## 2026-10-08: Two-epoch transient search, a null result; exotic-lens signatures (D-027)
- `scripts/transient_search.py` produces catalog-level candidates; `scripts/transient_forced.py` checks them with
  forced aperture photometry on byte-range cutouts of both epochs.
- SMACS against VENUS: catalog 441 (F444W) and 628 (F150W) → 144 in both bands → 4 by forced photometry (6 with
  WCS-centred apertures) → 0 after a visual check. All are bright stars or a source next to a bright neighbour, where
  fixed-aperture systematics dominate. No credible transient at |Δm| ≥ 0.3 mag, ≥ 5σ.
- **Failed approach:** catalog-only comparison across pipeline versions, which picks up deblending differences.
- docs/exotic_lensing.md records the verified exotic-lens signatures (wormhole demagnification, radial images for
  negative mass), the degeneracies with ordinary lensing, and existing limits. Warp drives have no imaging
  prediction to test.
- **Handoff:** the lens-model stage (TASKS "Now (M3)" 1), then Sunrise's same-pipeline epoch pair.
## 2026-10-08: Three cluster fields; robustness rules from them (D-025, D-026, cycle 17)
- Parallel worktree workers added Abell 2744 (#32), El Gordo (#31) and Sunrise/WHL0137 (#30). PR #29 fixed downloads
  for MAST products whose listed size is stale after the 2026-10-01 reprocessing. Proposals are folded into
  D-026 and SOURCES.
- **D-025** (found on Abell 2744):
  - catalogue "stars" with DJA r50 too large for a point source (bright cluster-galaxy cores in Gaia) go back to the
    galaxy ranking, before the locus calibrates;
  - very red spiky sources are exempt from the D-020 host test.
  - Run `20261007T211229Z-883a5066`: 8 vetoed, the locus recalibrated, the star top 10 is now mostly point sources with spikes, and
    `5904` is kept.
- **Handoff:** TASKS "Now (M3)". The lens-model stage (worker, `claude/lens-model`) is in flight; then per-field
  model checks and the transient forced photometry.


## 2026-10-08: Owner priorities, lensing-violation search first (D-023)
- The owner merged PR #20 (`2804`, inconclusive) and reprioritised: M3 lens-model consistency first; more clusters in
  parallel (Abell 2744, El Gordo, Sunrise; MACS J0416 has no DJA v7 catalog); a two-epoch transient search; no
  pauses; one review per PR.
- **Handoff:** TASKS "Now (M3)" items 1–3. Worker PRs with label `batch-clusters` are in flight.


## 2026-10-08: CEERS control on DJA matched photometry (D-022, cycle 16)
- `ceers_t021` now uses the DJA v7.4 `ceers-full` catalog, a 250.5 MB download with the reason stated in D-022.
  The control field is built the same way as SMACS for colours, confirmation, screening and the top-k metric.
- Run `20261007T202419Z-65728e80`:
  - 3,069 of 7,582 sources match DJA, and 3,156 pass the gate (2,757 before);
  - the top 20 is 85% known objects with 0% cutout-flagged, and 4 sources were screened;
  - the stellar locus does not apply (2 catalogued stars).
- **Found:** `f200w_2842` (#4) is a single-band, ellipticity-0.95 streak, probably a spike from the bright source
  4.6″ away. A coincident DJA object confirmed it (D-014 "Revisit if").
- **Handoff:** TASKS "Spike streaks confirmed by a coincident DJA object", then the star-stratum NED matches.


## 2026-10-08: Low-weight screening of the top-k pool (D-021, cycle 15)
- `stages.cutouts.screen_low_weight` re-applies D-011's 0.5 weight threshold with the cutout's own core weight. It
  runs on every non-star stratum, and `screened.ecsv` now records each source's `reason`.
- Run `20261007T194733Z-ac5f1b59`: the stripe artifact `3034` (#3) is removed: short-wavelength only, not in DJA, core weight 0.468.
  The galaxy top 20 is now 40% known objects, with 1 flagged source (`2915`, kept on purpose) and 4 screened
  (one by D-021, two by D-019, one by D-020).
- **Handoff:** TASKS "Star-stratum sources with NED galaxy matches", then MIRI's `G_Lens` vetting.


## 2026-10-08: Host test for spiky sources (D-020, cycle 14)
- `cutouts.host_ratio` measures annulus light against the peak around the spike peak. D-019 screening now also
  removes spiky sources with no host light (< 0.004), so the two-star segment `940` (0.0017) leaves the top 20.
  `2915` (0.0076, a galaxy nucleus) stays.
- Run `20261007T190843Z-db3c83fe`: the galaxy top 20 is 40% known objects, with 1 spiky source (`2915`, legitimate) and 3 screened.
- **Next contamination:** `3034` (#3) is a diagonal stripe artifact in a low-weight region (cutout flag `low_weight`).
  It passed the D-011 gate.
- **Handoff:** TASKS "Stripe artifact at galaxy rank 3".


## 2026-10-08: Top-k composition metric; `940` re-diagnosed (cycle 13)
- Every stratum now reports its top-k composition: known, lens-related, catalogued star, cutout-flagged,
  spikes and screened. It is stored in `run_record.json` as `samples[].topk` and defined in
  docs/methodology.md. The crossmatch rows now keep `is_lens_related`.
- Baseline, run `20261007T183007Z-4342c5c2`: the galaxy top 20 is 35% known objects, 10% lens-related, 0% stars and 15% flagged.
- **`940` is not saturated.** Its F200W cutout has no no-data pixels. The pipeline merged two stars into one
  5,788 px segment whose centroid lies 0.20″ from the brighter star's peak, so aper50 at the centroid is faint
  (25.7 against isophotal 20.1). It is a single-band detection confirmed by the DJA match (0.18″) to the bright
  star.
- **Failed approaches** (saturation tests):
  - aper50 − isophotal magnitude: cluster galaxies reach 10–21 mag, so `940` (5.6) is not an outlier;
  - a no-data core within 0.3″: `940` has none, and only `1345` does.
- **Handoff:** TASKS "Star-pair segments in the galaxy ranking" (a host-vs-PSF profile test or a centroid–peak
  offset rule) and "Star-stratum sources with NED galaxy matches".


## 2026-10-08: Vetted `f200w_2915`, a catalogued quiescent galaxy with a compact core (cycle 12)
- Run `20261007T174424Z-540be92d` (clean `a27d628`), galaxy rank 21, kept ranked by D-019 as a spiky source with non-stellar colours.
- It is SIMBAD `[VBG2023] SMACS 1060` (Valentino+2023 atlas of colour-selected quiescent galaxies at z > 3). DJA
  z_phot is 2.84 (95% 2.68–3.10), mildly below that selection. Its F090W−F150W break is 2.5 mag. The spikes come from its compact core. D-019 behaved as
  intended.
- **Two epochs:** the position residuals disagree between bands (7 mas F444W, 34 mas F150W), as expected for centroid
  shifts of an extended source, and Δm is F444W only (+0.11; F150W +0.01). So it shows neither motion nor
  variability. `epoch_compare.py` now documents that its significance assumes a point-like target.
- **New resource:** program 4043 has NIRCam F444W grism spectra over part of SMACS 0723 (SOURCES).
- **Handoff:** TASKS "Saturated stars in the galaxy ranking", then the star-stratum NED matches.


## 2026-10-08: Spike screening of the galaxy top k (D-019, cycle 11)
- The galaxy strata get cutouts for 2k. Spike-flagged sources with stellar colours (`1571`, `2242`) are removed
  and backfilled. Spiky sources with non-stellar colours stay ranked and are noted: `940` (a saturated star with
  corrupted colours) and `2915` (a bright nucleus, possibly an AGN).
- **Failed approach:** screening every spiky source removed `2915`. Spikes also mark bright galactic nuclei.
- Run `20261007T170539Z-b369bfa2`: the top 20 now runs to rank 22, with no star+galaxy blends left.
- **Handoff:** TASKS "Vet `2915`" (a bright red nucleus), "Saturated stars in the galaxy ranking" and "Star-stratum
  sources with NED galaxy matches".


## 2026-10-08: Diffraction-spike flag on cutouts (D-018, cycle 10)
- `cutouts.spike_statistic`: hexagonal-harmonic power around the brightest peak near the target. NIRCam cutouts with
  `spike_s6` ≥ 3 get the `spikes` flag, and the report lists them.
- Run `20261007T162133Z-2511977c`: it flags exactly the three PSF-like blends in the galaxy top 20 (`940`, `2242`, `1571`) and all of
  the star-stratum top 10. `940` is really a saturated star: its pipeline aper50 is 5 mag fainter than DJA's
  measurement, and the stellar locus missed it because its colours are corrupted.
- **Failed approach:** centring on the catalog centroid, which misses blends whose star is offset.
- PR #20 (`2804`, needs-human) was revised after review: F410M−F444W is now an observed colour, and the line
  is a hypothesis.
- **Handoff:** TASKS "Spike-flagged sources in the galaxy ranking" (do flagged sources leave the top k?) and
  "Star-stratum sources with NED galaxy matches".

## 2026-10-08: Vetted `f200w_2804`, inconclusive, a little-red-dot-like source (cycle 9b)
- Run `20261007T151316Z-5706ae30` (clean `6fade35`), galaxy rank 6. Uncatalogued, unresolved, F444W 22.9.
- **Second epoch:** no proper motion over 4.0 years (VENUS 6882, 2026-06-05; 1.9 mas = 0.95σ, so
  |μ| ≲ 1.5 mas/yr) and no F444W variability (0.004 mag).
- F444W is 0.725 mag brighter than F410M in the same epoch. That is an observed colour; a line at 4.25–4.98 µm
  and F410M absorption are both hypotheses.
- It meets Kokorev+2024 little-red-dot colour cuts. A brown dwarf is disfavoured on two counts: the 1.5–2.8 µm
  colours are red rather than blue, and there is no motion.
- Verdict: `inconclusive: needs NIRSpec spectroscopy`. Record: `docs/candidates/jw02736-o001_t001_nircam_f200w_2804.md`.
- **Handoff:** the owner decides whether to pursue it (`needs-human`).


## 2026-10-08: Two-epoch vetting tool (D-017, cycle 9a)
- `scripts/epoch_compare.py`: proper motion and variability of one target from two public level-3 catalogs of one
  filter, with a local median frame tie. Program 6882 (VENUS, 2026-06-05) gives SMACS 0723 a 4-year second epoch in
  F150W, F182M, F210M, F300M, F410M and F444W (SOURCES "Second epochs").
- **Handoff:** the vetting record for `f200w_2804` follows in its own PR (inconclusive, `needs-human`).


## 2026-10-08: One-sided stellar-locus size test (D-016, cycle 8)
- Bright calibration stars have inflated r50 (5–8 px under 20.5 mag), so D-015's two-sided band around 2.86 px
  missed the unsaturated stellar sequence at r50 ≈ 2.0 px. The size test is now 0.5–1.2 × r50_psf.
- Run `20261007T144011Z-11342d69`: the locus adds 90 stars, not 16. The star stratum has 109; its top 10 is all PSF-like.
  `f200w_1874` is now a star.
- Point-like sources left in the galaxy top 20 are faint (mag > 24), blends, or `2804`, whose red
  F356W−F444W makes it a brown-dwarf-like candidate to vet. Details are in D-016.
- **Handoff:** TASKS Now 1 (vet `2804`), Now 2 (blended stars).


## 2026-10-08: Vetted four high-ranked SMACS galaxies; DJA photo-z for vetting (cycle 7)
- PR #16 (D-015) review fixes merged. Run `20261007T132538Z-96e911bb` (clean `725eac3`) reproduces the D-015 numbers.
- `scripts/vet_evidence.py --photoz`: the nearest DJA eazy entry, labelled model_prediction (SOURCES "Vetting
  (cycle 7)"). `/vet-candidate` now names the script and the real `candidates` subcommands.
- All four are ordinary; none is an artifact. Records are in `docs/candidates/`.

  | Source | Rank | Verdict |
  |---|---|---|
  | `829` | 3 | red, dusty galaxy at spectroscopic z = 2.74 (NED), a published F150W-dropout "z ≈ 11–20" candidate |
  | `1032` | 8 | lensed arc `[MJR2023] 028.2` (Mahler et al. 2023), a tangential pair with 028.1; blue colours likely from emission lines (hypothesis) |
  | `438` | 11 | uncatalogued; Balmer-break galaxy at z_phot 4.8 (hypothesis, 6-band photo-z) |
  | `1243` | 12 | catalogued red point-like high-z candidate, photometric z ≈ 5.6–5.75 |

- **Lesson:** after D-011/D-014/D-015 the top of the galaxy ranking is dominated by real but rare populations
  (lensed arcs, dusty and high-z galaxies), most already catalogued. The pipeline recovers known lensing
  features. Uncatalogued sources such as `438` are the discovery set.
- **Side finding:** the D-015 stellar locus missed the PSF-like star `f200w_1874` (F200W 21.4 mag, spikes)
  because its r50 (2.10 px) is 27% below r50_psf (TASKS Now 1).
- **Handoff:** TASKS Now 1 (stellar-locus lower bound) and Now 2 (blended stars).


## 2026-10-08: Stellar locus for stars missing from Gaia (D-015, cycle 6)
- `classify.stellar_locus`: a source counts as a star when its DJA detection-image r50 is within 20% of the
  catalogued stars' median and its colours are stellar. The DJA join now carries `r50_pix` and `mag_auto`.
- Run `20261007T130207Z-c8ffc470`: the locus added 16 stars to SMACS, and 3 of the 6 PSF-like sources left the galaxy top 20.
  `829`, `1032` and `438` are ranks 3, 8 and 11. Details are in D-015.
- **Failed approaches:**
  - pipeline CI and `semimajor_sigma` do not separate stars from galaxies, because spikes and
    saturation inflate them;
  - size alone admits compact galaxies, because their colours are too broad.
- **Handoff:** TASKS Now 1 (vet `829`, `1032`, `438`, `1243`). Three blended stars remain in the galaxy top
  20 (TASKS Now 2).


## 2026-10-07: Rank only confirmed detections (D-014, cycle 5)
- Quality gate: best-band S/N floor (`low_snr`) and single-band confirmation (`single_band`, confirmed by a
  DJA match). Too few survivors now fall back to the D-011 tests instead of ranking ungated.
- S/N is now inverted exactly from the pipeline's `abmag_err` everywhere: the features S/N ≥ 3 gate,
  the quality gate and the scripts. The linear approximation overstated S/N at low S/N.
- Run `20261007T124944Z-bbad4ab3`: spike and stripe detections in the SMACS galaxy top 20 fell from 5–6 to 1, and near-noise
  sources from about 6 to 0. `[YML2023] F150DB-C-4` and `[MJR2023] 028.2` are ranks 4 and 9. Evidence and
  counts are in D-014.
- **Failed approach:** an S/N floor in the reference band only removed red dropouts, including `829`
  (D-014).
- **Handoff:** TASKS Now 1 (stellar locus: about 5–6 faint PSF-like stars remain), Now 2 (vet `829`,
  `1032`, `438`, `1243`).

## 2026-10-07: Size-robust colours from DJA matched apertures (D-013, cycle 4)
- New `photometry.py`. The DJA v7.4 catalog is verified against its sha256 before it is cached, and is
  joined one-to-one within 0.2″. Colours now come from 0.5″ matched apertures.
- `features`: `ref_mag` and the morphology features keep the pipeline aperture's S/N gate, so sources
  without a DJA match keep their morphology. With matched photometry, a band counts as detected when DJA
  measured it at S/N ≥ 3. DAOFind sharpness/roundness apply only when 0 < CI_50_30 ≤ 1.8.
- **The size bias is gone, measured on the same 1,423 rows** (run `20261007T122054Z-c52935ec`, F200W−F277W, F200W aper50
  S/N > 10). Spearman between log area and colour:
  - aper50: +0.181 (p = 6.7e-12);
  - isophotal: −0.312 (p = 1.5e-33);
  - DJA 0.5″: −0.012 (p = 0.66).

  Re-measure with `python scripts/feature_size_bias.py --run-dir <outputs>/runs/20261007T122054Z-c52935ec --sample
  smacs0723_nircam --compare dja05`.
- Join coverage: 2,729 of 5,254 sources matched one-to-one. Nearest-neighbour matching gave 2,920, of
  which 191 were fragments sharing one DJA object.
- **Ranking:** the sources vetted in #13 dropped (`2925` to rank 1297, `2559` to 482, `1096` to 296). Two
  catalogued lens-related objects are in the galaxy top 20 (observed cross-matches):
  - rank 7 `f277w_829`, SIMBAD `[YML2023] F150DB-C-4` (also NED G_Lens);
  - rank 13 `f200w_1032`, NED `SMACS J0723-73:[MJR2023] 028.2` (G_Lens).
- **Top-20 purity is still mixed** (visual, unvetted):
  - 4 faint PSF-like stars missing from Gaia;
  - 5–6 detections on diffraction spikes or parallel stripes (a new artifact class the D-011 gate misses);
  - about 6 faint, compact or near-noise sources, one of them hot-pixel-like;
  - 3 arc-like sources;
  - 1 catalogued elongated galaxy.
- Failed approaches, both fixed after `/code-review` (15 findings, all addressed):
  - nearest-neighbour joining (fragments shared photometry);
  - gating morphology on DJA S/N, which erased the morphology of 44% of sources.
  The first run, `20261007T120528Z-378af5a1`, used both.
- **Handoff:** TASKS Now 1 (mask detections on diffraction spikes and stripes), Now 2 (S/N floor for
  ranking), Now 3 (stellar locus), Now 4 (vet the catalogued lensed images and arcs).

## 2026-10-07: First vetting: three arc-like galaxy candidates (cycle 3)
- New `scripts/vet_evidence.py` gathers vetting evidence:
  - per-band catalog rows;
  - radius and tangential alignment from both the catalog orientation and image moments (they agree
    to within 1–7°);
  - nearest multiple image in Mahler+2022 `arcs.dat`;
  - nearest star;
  - SIMBAD/NED/Gaia cross-match at 1″ and 3″;
  - six-band S3 cutouts.
- Vetted SMACS NIRCam `f200w_2925`, `f200w_2559` and `f200w_1096` (records in `docs/candidates/`; store status
  `explained`):
  - none is among the 62 Mahler+2022 constraint images (nearest 24–45″); this does not show they are
    singly imaged, which needs the model's critical curves (M3);
  - 2925 and 2559 are not tangentially aligned (45–60°); 1096 is at 19°, consistent with weak shear
    (`hypothesis`);
  - NED knows 2925 (z 1.98) and 1096 (z 1.36, with catalogued clumps) as background galaxies.
- **Systematic finding (derived):** the baseline's colour and shape features are size-biased.
  - aper50 F200W−F277W reddens with isophotal area (median +0.77 → +0.99, Spearman 0.21,
    p ≈ 1e-16), while isophotal colour does the opposite (+0.71 → +0.02). EE apertures differ in angular
    size between SW and LW.
  - Reproduce with `python scripts/feature_size_bias.py --run-dir <outputs>/runs/20261007T040938Z-01527ace
    --sample smacs0723_nircam`: sources detected in F200W and F277W with F200W aper50 S/N > 10, F200W
    isophotal-area bins [0, 50, 200, 1000, ∞) px.
  - DAOFind sharpness reaches +145σ on clumpy extended galaxies.

  This explains these three candidates (verdict `catalog effect`) and sets the next priority (TASKS Now 1).
- **Handoff:** TASKS Now 1, size-robust features (DJA matched-aperture colours; sharpness only for compact
  sources).

## 2026-10-07: Stars ranked as their own stratum (D-012)
- New `classify.classify_sources`. One bulk Gaia DR3 + SIMBAD cross-match over the whole catalog
  (about 8 s per sample) labels stars. The runner ranks them as `<sample>-stars` with their own top k;
  galaxies are ranked without them.
- Failed approach: the pipeline's `is_extended` flag labels every bright, saturated star in the top 20
  as extended.
- Real data, run `20261007T040938Z-01527ace`:
  - Stars found: SMACS 52, MIRI 12, CEERS 7. Of these, 40, 8 and 7 passed the gate and were ranked as
    `-stars`.
  - Galaxy-stratum top 20s: 0 cross-matched stars (SMACS had 7), 1/60 image-quality flags.
- Visual check (unvetted): the SMACS galaxy top 20 is now mostly interacting, clumpy or elongated
  galaxies, including three arc-like sources (`f200w_2925`, `f200w_2559`, `f200w_1096`). About 4
  faint PSF-like sources without Gaia counterparts remain.
- `/code-review` found 15 issues, all fixed. The main ones:
  - a classification leaving fewer than 2 galaxies could abort the run;
  - stars were dropped entirely when there were fewer than `min_stars`;
  - the offline suite made live CDS calls;
  - the rule used any match instead of the nearest one.
- **Handoff:** TASKS Now 1 (vet the arc-like sources against published lens models), then Now 2
  (stellar locus for faint stars).

## 2026-10-07: Quality gate before ranking (D-011)
- New `quality.assess_sources` with `cutouts.sample_weight_map`. A coarse WHT map (1″ cells,
  concurrent S3 byte ranges, about 10 s per 1.8 GB mosaic) gives each source's relative weight and
  edge distance. A CI_50_30 test
  catches sources sharper than the PSF. Ranking now sees only sources that pass.
- Real data, run `20261007T033832Z-cfcf6032`: top-20 sources with a cutout image-quality flag fell from
  20/60 to 2/60.
  - SMACS NIRCam: 6→1. MIRI: 10→0. CEERS: 4→1.
  - Gated out: 28% of SMACS NIRCam, 55% of MIRI, 28% of CEERS.
- `/code-review` found 15 issues in the first version, all fixed:
  - edge distance was snapped to cells, so `edge` never fired for LW or MIRI;
  - a gate passing too few sources could abort the run;
  - the gate ran even when unconfigured;
  - noise was labelled as artifacts.
- Visual check (unvetted): the SMACS NIRCam top 20 is now about 12 bright stars plus 8 galaxies, which
  include interacting pairs and the two arc-like sources (`f200w_2925`, `f200w_2559`). Low-weight
  noise and streaks are gone. Stars dominate, so star/galaxy separation is next (TASKS Now 1).
- The agent self-merge policy was adopted (owner decision; CLAUDE.md "Merge policy and version
  control"), and PRs #1–#10 were merged.
- **Handoff:** TASKS Now 1 (star/galaxy separation).

## 2026-10-07: M0 bootstrap batch landed; first real-data run (M1 slice)
- Nine parallel units became PRs #1–#9: archive, catalog, rank, cutouts, crossmatch, runner/CLI, OSS,
  agent harness, landscape. The integration branch merges all of them and adds a contact sheet to the
  runner. 297 offline tests pass.
- First end-to-end run on real public data: run `20261007T020124Z-4bfabaaa`, config `reference_sample_v1`,
  about 3 minutes on the owner's laptop.
  - Merged sources: SMACS 0723 NIRCam (6 bands) 5,254; MIRI (4 bands) 530; CEERS t021 (7 bands) 7,582.
  - The top 20 per sample get cutouts through S3 byte-range reads (no full i2d download) and a
    SIMBAD/NED/Gaia cross-check.
- Results. These come from the integrating agent's visual classification of the contact sheets
  plus the cross-match; nothing is vetted.
  - SMACS NIRCam top 20: mostly bright stars with diffraction spikes (7 confirmed by Gaia/SIMBAD,
    about 11 PSF-like), 2 image artifacts (a streak and low-weight noise), and galaxies or blends.
    Ranks 6 and 19 are elongated and arc-like; they are hypotheses to vet against published lens models.
  - MIRI top 20: 10/20 carry image-quality flags (edge or low weight). One known lens-related object
    (NED `G_Lens`) was recovered; it is a known object, not a discovery.
  - CEERS control: with F200W only, about 16/20 were artifacts (edges, hot-pixel-like compact sources,
    linear streaks). With all 7 bands, about 5/20 are artifacts and the rest are structured galaxies, so
    cross-band features matter.
  - Interpretation: the baseline mostly surfaces stars and instrument or processing artifacts, as
    expected. Quality gating and star/galaxy separation come next (TASKS Now 1–2).
- Injection-recovery (unit 3, NIRCam): only LOF clearly beats random on every injection type. The
  combined score reaches precision@50 of 0.06–0.20. Nothing is tuned toward one injection type until
  vetted labels exist (D-004).
- Data caveats: program 2736 catalogs were made with jwst 2.0.1 and CEERS with jwst 3.0.0, and deblending
  changed in between (D-010). The CEERS `s_ra/s_dec` is the MIRI target; use `s_region` instead.
- Process lessons:
  - `setup-uv` has no floating major tag, so pin the full version.
  - Workers collided in the shared scratchpad; each now gets its own subdirectory (CLAUDE.md).
  - Unit 5's `gh pr create` was blocked, so the coordinator opened #9 with an inline body;
    `--body-file` is denied by project settings.
- Guardrails:
  - `main` is protected: a PR and 4 CI checks are required, admins included.
  - A ruleset blocks force-push on every branch.
  - Private vulnerability reporting and Dependabot security updates are on.
  - Direct commits to `main` happened only in Step 0, before protection existed.
- **Handoff:** the owner merges #1–#9 in any order, then this integration PR (updated by merging `main`
  into it). The next cycle starts at TASKS "Now" item 1.

## 2026-10-07: M0 bootstrap started
- Created repository skeleton: interfaces (`schema.py`, stage stubs with fixed signatures), state files,
  CI, BSD-3-Clause license, verbatim agent charter (`docs/agent-charter.md`) imported by `CLAUDE.md`.
- Verified on MAST: program 2736 level-3 NIRCam/MIRI imaging is public; catalogs ~3 MB vs NIRCam i2d
  ~1.8 GB, hence the catalog-first slice (D-001). The S3 mirror key pattern is confirmed.
- Owner decisions: public repo; every change via PR and the owner merges; chat in Chinese.
