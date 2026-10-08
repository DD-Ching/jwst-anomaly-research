# Changelog

Newest first. Results, failed approaches and the handoff state; not a diary.

## 2026-10-08: W3 in MOA-II Cut-0 light curves (gb22): the pre-screen keeps W3 events; chunk 1/8 fitted (D-062)
- Cloud run; continues the unpushed-PR branch `claude/w3-moa` of the previous run (adapter, script, tests), merged
  with main (its D-060 renumbered **D-062**). `gb22.tar` (3.5 GB, sha256 verified) and `metadata.ipac.tar.gz` fetched.
- `prescreen` on all 18,599 gb22 light curves: 0 errors, 137 s wall; **1,058 passes** (`derived`).
- `inject --prescreen-only` (seed 60; 300 W3, 100 PSPL on real quiet light curves): W3 survives emulated Cut-0 and the
  pre-screen at 23 / 43 / 48 / 38 / 28 % (t_E 3 / 10 / 30 / 100 / 300 d); PSPL 0 / 100. This is the first sample in
  this project whose selection keeps W3 (OGLE 0 / 600, Gaia 0 / 240).
- `fit --chunk 1/8`: 133 passes, 0 errors, 855 s wall; best ordinary PSPL 130 / FSPL 3; ΔBIC(min exotic) 5/50/95 % = −13,497 / −631 / −193; **131 / 133 flag** (`derived`). Every pass flags (dip-shaped variables prefer a negative-flux model), so `vet` is
  the discriminating step; no conclusion before it runs on all chunks and the contact sheet is inspected.
- New: `fit --chunk K/N` writes `results/w3_moa/fits_gb22_chunkKofN.ecsv.gz`; `merge-chunks --n N` refuses missing,
  stale-Params or wrong-membership chunks (the pre-screen is deterministic and recomputed each session, ~2 min).
- **Failed approach:** waiting with `until ! pgrep -f '<script>'` and `pkill -f` from the shell — both match the
  waiting shell itself (lost ~8 min; the D-059 lesson again). Wait on a log line or the output file instead.
- **Handoff / next:** chunks 2–8 (`OMP_NUM_THREADS=1 python scripts/w3_moa.py prescreen && ... fit --chunk K/8`,
  ~15 min each), then `merge-chunks --n 8`, `vet`, `sheet`, full `inject`, `limit`.

## 2026-10-08: W3 in the Gaia DR3 microlensing candidates: the published selection rejects W3 too (D-061)
- Cloud run. Hypothesis: the Gaia DR3 candidates (363; Wyrzykowski et al. 2023) are selected less PSPL-shaped than
  Mróz et al., so they could limit W3. Tested the selection before fitting, as D-057 requires.
- New `gaia_mulens.GaiaDR3Microlensing` adapter (TAP + DataLink epoch photometry, 8 parallel batches of 12 ids,
  ~3 s per source on the server; Table D.1 sample labels from the pinned arXiv source) and `scripts/w3_gaia.py`.
- Emulated Sample A cuts pass 143 / 163 real Sample A events. W3 injections: **2 / 240 selected, 0 / 240 selected
  and flagged**; PSPL controls 17 / 120. W3 dimming fails skewness < 0 and the skew–Abbe cut, as designed for
  brightenings. The fits of all 363 give one flag (4053892503992268288, ΔBIC −40.3). It is an event truncated at the
  window end on a variable baseline (light curve inspected), not a candidate.
- **Failed approaches (rules):** the Extractor cuts (n points, duration > 135 d, max σ > 50) cannot be emulated from
  the paper. The guessed definitions fail 126 / 163 real events, so they are left out. A single-id DataLink request
  answers with bare CSV, not a zip. Truncated chunked replies happen, so retry them.
- **Next:** W3 needs light curves taken before any microlensing selection. Gaia DR3 has epoch photometry only for
  its variable sources (vari_summary; ~11.7 M per the DR3 release, recheck), selected by variability, not shape. Next, check whether W3 survives that
  classifier (`vari_classifier_result`) by injection, and whether a sky-region subset is small enough to screen.
  Other options: KMTNet public seasons; OGLE EWS (owner decision, terms).

## 2026-10-08: W2 deflector test at HST resolution from the Hubble Source Catalog: not decisive (D-060)
- Cloud run. Hypothesis: HST resolution decides the lensed quasars that Legacy Surveys left blended or too close
  (D-056). `scripts/w12_hsc_probe.py`: HSC v3 summary sources within 4″ of each of the 444 galaxy-scale quasar/radio
  systems (MAST catalogs API, 58 s); ≥ 2 point sources (CI < 1.3) are the images, an extended source (CI ≥ 1.5)
  near their centroid and > 0.2″ from an image is the deflector; sources in < 2 HSC images are dropped as likely
  artifacts (MAST's recommendation; #90 review) (ASSUMPTIONs in `Params`).
- 71 of 444 systems have HSC sources in ≥ 2 HSC images (91 with any). Validation on systems with a published lens redshift (a lens galaxy is known
  to exist): 13 deflector, 15 none, 38 undecided → **efficiency 13/28 = 0.46** (`derived`; 15/35 = 0.43 without the
  artifact cut). Misses include quads and doubles (13 of 16 "none" have 2 point images). In four inspected misses (H1413+117, HE1104−1805, SBS0909+532,
  HE2149−2745) the HSC rows within 4″ are only the quasar images: the lens galaxy is absent from the catalogue, not
  mis-typed (likely lost in the quasar PSF; hypothesis, no cutouts inspected). Without a lens redshift: 1 none
  (HS0810+2554), 4 undecided. The 3 open SuGOHI IX pairs (D-056) have no HSC sources.
- **Failed approach (rule):** HST *catalogue* photometry cannot decide a dark deflector in lensed quasars — a
  "none" is more likely a missed lens (0.54) than a dark one. No limit, no candidate; HS0810+2554 is not flagged.
- **Next:** PSF-subtracted HST image modelling (e.g. drizzled frames from MAST, quasar PSF + Sérsic fit) validated
  on the same known-lens set, or HSC PDR3 photometry; until then W2 in wide imaging stays at D-056.

## 2026-10-08: Warp literature check: still nothing testable
- Subagent search (arXiv API 2023–2026, INSPIRE citations of Clough et al. 2024): no imaging or lensing prediction
  for a distant observer; Lentz & Felton 2024 give order-of-magnitude EM fluxes for a bubble 100 lyr away but no
  template that separates it from ordinary transients (found by review); the collapse-burst waveform is not public; an O3 search for superluminal-source GW
  bursts (Kuwahara & Cannon 2023) is already null. D-047 stands; recorded in SOURCES.md and docs/exotic_lensing.md.

## 2026-10-08: W3 OGLE bulge (all 5,790 events): no candidate; the published selection rejects every W3 event (D-057)
- Worktree worker, unbounded π_E (predates D-058; bounding can only add flags). Bulge: 0 fit failures; best ordinary
  PSPL 5,377 / PAR 401 / FSPL 12; ΔBIC(best exotic) 5/50/95 % = −3.6 / 5.7 / 12.0; 127 flags < −10 (`derived`).
- Vetting, cumulative: 127 → 113 refit all ordinary → 80 robust errors → 73 variable baseline → 14 season offsets →
  9 season drifts → 7 binary source/lens, VSX/Gaia, arXiv → 1 `feature_coverage` → 1 jackknife → **0** two unrelated
  events (new `revet` tests: ≥ 3 epochs where the models differ by > 3σ; drop up to 3 influential epochs keeping ≥ 3 in
  the feature; two independent PSPL bumps). BLG519.21.110304's exotic spikes sat on the 2011 event and a 1-day bump in
  2015 (two PSPL bumps better by ΔBIC 24.4). Disk: 6 → 0.
  Contact sheet of the 7 late survivors inspected. The D-059 chunk survivors fail `feature_coverage` here.
- Injections: 600 W3 events (n = 1, ε < 0; t_E 3–300 d, ρ 0.01/0.1) on real bulge cadences + 300 PSPL controls. The
  fitter flags 42–97 %, but **0 / 600 pass the emulated Mróz selection** (controls 15–43 %; cuts failed most: one
  bump, PSPL χ², χ₃₊). The emulation passes 63.9 % of the real selected events (somewhat stricter). No rate limit is
  derivable; ε_W3/ε_PSPL < 0.12–0.33 (95 %). Wall time ≈ 5 CPU h (bulge fit 10.4 ks, vet 6.6 ks, inject 2.0 ks).
- **Failed approaches (rules):** exotic fits started only from the PSPL solution miss the W3 geometry — start with
  caustic spikes on pairs of maxima and on absolute t_E; ΔBIC alone is not a candidate test (season blends and
  feature coverage remove 120 of 127); a jackknife must keep ≥ 3 epochs in the exotic feature or it kills real short
  events; an exotic fit whose spikes sit on two bumps years apart needs the two-unrelated-events test; pin BLAS threads (`OMP_NUM_THREADS=1`) with process pools; a published
  PSPL-selected sample cannot constrain a non-PSPL signal — inject through its selection before fitting it.

## 2026-10-08: W3 OGLE bulge, chunks 2–3/12: 20 flags, no candidate; chunk tables tracked (D-059)
- Cloud run. Chunk 1/12's fit table lived only in the ephemeral session and is lost. `fit --chunk K/N` now also
  writes a deterministic gzipped copy to `results/w3_ogle/` (~0.22 MB per chunk); `merge-chunks --n 12` joins
  chunks 1..12 into the table `vet` reads and marks it the whole sample only when every chunk is present, fitted
  with the current `Params` and holds exactly its own events (none skipped) (D-059). A chunk takes ~9 min on 4 cores (not ~17).
- Chunk 2/12 (483 events, 0 errors): best ordinary PSPL 454 / PAR 28 / FSPL 1; ΔBIC(min exotic) 5/25/50/75/95 % =
  −4.7 / 3.8 / 5.8 / 6.6 / 11.2; 49 below 0, 9 flags below −10 (`derived`). `vet`: 8 fail; **BLG624.18.69573**
  (no EWS name; t_E ≈ 180–240 d, best ordinary PAR) passes every automated test (N1neg ΔBIC −19.0 vs PAR; season
  offsets −15.0, drifts −15.5; binary source / lens −10.1; 0 VSX / Gaia matches). Contact sheet and residuals
  inspected: the N1neg model puts its first caustic spike inside a season gap (no data at t − t0 ∈ [−50, 0] d; the
  models differ by > 5 % over [−89, −7] d); its Δχ² comes from 2 peak points (−9.4) and 16 post-peak points (−6.7).
  An unsampled caustic plus a sparsely sampled peak is not evidence: not a candidate (ASSUMPTION-level judgement).
- Chunk 3/12 (483 events, 0 errors): PSPL 447 / PAR 36 / FSPL 0; ΔBIC(min) 5/50/95 % = −2.6 / 5.7 / 11.1; 39 below
  0, 11 flags; **0 survive** `vet` (season offsets/drifts remove 10, the refit of all ordinary models the 11th; one
  VSX match). Contact sheet inspected; in several flags the exotic and ordinary curves also differ mainly in gaps
  (e.g. BLG597.28.9837 has no peak data).
- **Failed approach:** chaining chunk runs with `while pgrep -f 'chunk 2/12'` — the waiting shell matches its own
  pattern and never starts the next chunk.
- **Next:** a `gap_coverage` vetting test (require data where the exotic and best ordinary models differ, else the
  flag fails); refit chunk 1 and fit chunks 4–12 (two or three per run); `merge-chunks`, `vet`, `sheet`; the
  empirical ΔBIC null and xallarap fit for BLG667.04.62161 (scratch null-simulation design: PSPL best fit plus white
  noise, and plus season-wise circularly shifted residuals); then `inject` / `limit`.

## 2026-10-08: W3 OGLE bulge, chunk 1/12 (483 events): one marginal flag survives automated vetting
- Cloud run. `w3_microlensing.py fit --chunk K/N` fits events K−1, K−1+N, … so sessions fit disjoint, field-balanced
  parts of the bulge sample (measured ~10 s CPU per bulge event, ~4 h for all 5,790 on 4 cores; one chunk of 12 is
  ~17 min wall). Chunk 1/12: 483 events, 0 errors; best ordinary PSPL 449 / PAR 32 / FSPL 2; ΔBIC(min exotic)
  5/25/50/75/95 % = −3.5 / 3.7 / 5.6 / 6.5 / 12.4; 41 below 0, 12 flags below −10 (`derived`).
- `vet`: 11 of 12 fail; **BLG667.04.62161 (OGLE-2015-BLG-1250) passes every automated test**: E2neg ΔBIC −12.1 vs
  PSPL; season offsets −12.7, season drifts −12.4, binary source / binary lens no better than PSPL, 0 outliers,
  baseline χ²/dof 0.87, 0 arXiv records, no VSX / Gaia variable within 1″. Contact sheet inspected: faint source
  (peak flux ~2× baseline), large scatter, a flattened peak and a few low points ~+20…+50 d; the E2neg plateau fits
  those. ASSUMPTION-level judgement: marginal, not a candidate — with 483 trials and 12 flags below −10 from a
  heavy-tailed ΔBIC distribution, one −12 survivor is expected without any exotic lens; untested ordinary
  explanations: xallarap, per-season error underestimation near the peak, blending/difference-imaging systematics
  of a faint source, and a calibrated null (the injection/limit stage).
- **Next:** chunks 2–12 (one per run: `fit --chunk K/12`, then `vet` / `sheet`); for BLG667.04.62161, an
  empirical ΔBIC null from the same chunk's PSPL-simulated light curves (does −12 occur at rate ≥ 1/483?) and an
  xallarap fit before any further attention; then `inject` / `limit`. `limit` refuses vetting of a single chunk
  (#86 Codex), so a `merge-chunks` step (concatenate `fits_bulge2019_chunk*of12.ecsv`, then `vet`) comes first.

## 2026-10-08: W3 OGLE disk re-fit with bounded parallax: still null (D-058)
- Cloud run. ASSUMPTION `Params.pie_max = 5`: PAR fits with |π_E| > 5 are rejected in every fit and vetting refit
  (the unbounded fits reached π_E ~ 30–1,400). Disk (460 events, 376 s): best ordinary PSPL 408 / PAR 52 / FSPL 0;
  36 of 368 PAR fits sit on the bound; ΔBIC(min exotic) 5/50/95 % = −4.0 / 3.6 / 6.9.
- 7 flags (was 6; new GD1217.10.8703 at ΔBIC −10.0, `E2pos`), **0 survive** `vet` (all tests complete). The new flag
  is four post-peak points 0.15–0.35 below baseline that neither model fits; per-season offsets remove it (ΔBIC 0.9).
  Contact sheet inspected.
- Absorption check (simulated: 40 E2pos/N1neg injections, 700 uniform epochs, white noise σ = 0.05–0.3): flags 25
  bounded vs 25 unbounded; ΔBIC shifts ≤ 0.46. On this cadence unbounded parallax does not absorb exotic signals,
  so the bound is a physical prior, not a sensitivity gain. Not tested: seasonal gaps and correlated systematics.
- Environment: install the `mulens` extra (`-e ".[dev,cloud,mulens]"`); without MulensModel `fit` silently skips PAR.
- Follow-up (#84 Codex): checkpoint rows carry a `Params` hash; `fit` refits rows from other Params, so a chunked
  bulge run never mixes bounded and unbounded parallax fits.
- **Next:** bulge `fit` in chunks or locally, then `vet` / `sheet`; `inject` / `limit` with season-drift vetting.

## 2026-10-08: W1/W2 in published lens catalogues: no dark deflector; weak limits (D-056)
- Worktree worker: lenscat (32,838), Euclid Q1 Discovery Engine (2,584) and SuGOHI (3,961) merged into 35,862 systems
  (`src/jwst_anomaly/lenscats.py`, `scripts/w12_lenscats.py`; pinned by sha256). Deflector test in Legacy Surveys
  DR10 Tractor; footprint and depth from the DR10 brick summary, not from detected sources.
- 17,555 galaxy-scale systems in the footprint. Galaxy-finder (17,102) and sub-mm (110) systems cannot show a dark lens
  and give no limit; only lensed-quasar (325) and radio-interferometric (18) systems are tested (pair test: two point
  images, nothing bright enough near the expected deflector). Blended, too-close and faint-galaxy cases are undecided.
- 29 decided (13 with a deflector, 16 without); 13 of the 16 have a literature lens galaxy, and 3 SuGOHI IX CHITAH
  pairs (090434−005328, 091517+040747, 104122−005618) are open only in the typical variant — a lens below the local
  LS depth explains them conservatively. Not candidates. Cutout sheets inspected.
- 95 % limits on the dark-deflector fraction, test completeness assumed (not measured): typical f_dark < 0.48
  (quasar, k = 3, N = 16), < 0.23 (radio, k = 0, N = 13), < 0.27 (all, k = 3, N = 29); conservative < 0.50 (quasar,
  N = 6). Earlier 1.5 × 10⁻⁴ / 8.2 × 10⁻³ / 0.13 are withdrawn (PR #81 reviews). No W1 geometry (no image positions).
- **Failed approaches (rules):** a limit is valid only over systems where the test could have found the signal —
  prove it by injection; a recovery factor from deleting deflectors and re-running the same code is 1 by
  construction; coverage must come from footprint/depth products, never from "a source nearby"; evaluate exclusion
  flags on every system; Data Lab TAP takes no table uploads or q3c (batch box ORs); lenscat types cluster-survey
  entries as "galaxy", has AGEL declination and SPT position errors and rounded positions, and keeps rejected
  candidates; "no lens redshift" ≠ "no lens".

## 2026-10-08: W3 OGLE-IV disk sample: null (D-057)
- Cloud run. `jwst_anomaly.ogle` (Mróz et al. 2019/2020 adapter) and `scripts/w3_microlensing.py`
  (`fit` / `vet` / `sheet` / `inject` / `audit` / `limit` / `manifest` / `summary`). Manifest
  `data/manifests/ogle_mroz.ecsv`. `www.astrouw.edu.pl` reachable from the cloud (2026-10-08).
- **Disk (Mróz 2020, all 460 Table B1 events; 525 s on 4 cores):** best ordinary PSPL 406 / PAR 54 / FSPL 0.
  ΔBIC(min exotic) quantiles 5/50/95 % = −3.9 / 3.6 / 6.9; 6 flags below −10, all `E2pos`. `vet`: 0 survivors;
  every flag loses the exotic preference under per-season baseline offsets and/or drifts (ΔBIC −4.6 … +8.0), and
  BLG568.12.9169 is already −5.1 on refit. arXiv 0 mentions, no VSX / Gaia DR3 variable within 1″. Contact sheet
  inspected: two flags have sparse peak coverage; two have post-peak points below baseline (the Ellis
  demagnification shape) that one season's drift absorbs.
- **Limits:** none yet (injections use bulge light curves); the disk null is a flag count, not a rate limit.
- Review follow-up (#82 Codex): binary-lens α starts were passed in radians to MulensModel (degrees), so only
  0.5–5.8° was searched; fixed and disk re-vetted: still 0 survivors (BL BICs move by ≤ 18). Injection vetting now
  refits PAR with season trends as `vet` does; failed XMatch queries leave a flag unvetted; `limit` refuses a
  zero-event limit unless the bulge vetting is a complete null.
- **Caveats:** unbounded parallax fits reach π_E ~ 30–1,400 (unphysical); they can absorb an exotic signal and cut
  sensitivity. The season-drift test may also absorb real W3 dips: calibrate both with injections.
- Timing: the bulge fit is ~1.1 s/event on 4 cores (~1.8 h for 5,790); `fit` checkpoints to
  `fits_<key>.partial.jsonl`, but cloud disks are ephemeral: fit in chunks per run or locally.
- **Next:** bound π_E (ASSUMPTION, e.g. |π_E| < 5) and re-fit disk; bulge `fit` (local or chunked), `vet`, `sheet`;
  `inject` / `limit` with season-drift vetting inside the injection loop.

## 2026-10-08: W1 shear (aperture-mass) screen: four clusters null; limits 5–8× stronger than radial (D-053)
- Cloud run. `exotic_screens.py shear` builds D-050: PSF-deconvolved catalogue ε, cluster shear removed, spike
  segments vetoed, Schirmer 10″ aperture-mass S/N map, rotation null, B-mode check. `scripts/inject_shear.py` reuses
  D-049's W1 painting.
- **Validation on known signals.** Measured ε along the cluster model's g rises with |g| (MACS0416, Abell 2744), with
  responsivity R = 0.41–0.48, not 1. The screen now removes R g; injected images keep R of their lens shear.
- **Real fields** (Abell 2744, MACS0416, MACS1149, Abell 370; 30.2 arcmin²): **null**. The first run (no spike veto)
  had Abell 370 E 4.51 / B 4.68 and Abell 2744 E 3.92, all p ≤ 0.01 against the rotation null. Diffraction spikes
  were the cause: the veto (77–178 segments per field) brought them to 3.48 / 3.46 and 3.39. B-mode extremes still
  beat the rotation null in three fields, so an E peak must also beat the field's max |S_×| (post hoc,
  conservative). None passes both.
- Injections (200 per field and mass; spike rows never painted): 33 / 197 / 405 of 800 at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (radial 8 / 84 /
  156 of 1,600); 1 at 2 × 10¹¹, 0 at 2 × 10¹⁰. 95 % limits 7.7 × 10³ / 1.3 × 10³ / 6.5 × 10² deg⁻² (radial 6.1 × 10⁴ /
  7.0 × 10³ / 4.0 × 10³). docs/exotic_limits.md.
- **Failed approaches (rules):** subtracting the full model g from catalogue moments (leaves a W1-signed radial
  residual); any radial-alignment statistic without the spike veto; trusting the rotation null without a B-mode
  check; point-mass 1/x² and top-hat filters (lost to Schirmer 10″ in injections).
- Data: the DJA SMACS v7.4 and El Gordo v7.0 photo-z tarballs return 404 (2026-10-08), so these fields are out. MAST
  catalogues and CANUCS DR1 catalogues (~30 MB each) downloaded fine.
- Wall time: 58–168 s per field for 1,000 injections; the base screen takes 6–13 s. Review rounds changed the injection model (spike and
  near-core and unresolved rows not painted, lens change from raw moments); every re-run kept the four real fields
  null. The counts moved from 59/264/459 (first model) to 33/197/405 (final); not painting unresolved rows was the
  largest step (from 50/238/449; docs/exotic_limits.md caveats).
- Final review: painted images now face the spike veto too (409 → 405 at 2 × 10¹³ M☉; limits unchanged at two
  digits). `psf_sigma_px` (now shared with `inject_radial.py`) raises when no S/N > 50 row has a size, where it
  returned NaN. Left as is (maintainability only): the E/B summary is computed in both `cmd_shear` and
  `inject_shear.run_field`; the null recomputes the |e|² noise term per draw.
- **Handoff:** TASKS "JWST focus" 1 W1: lower the B-mode floor (PSF anisotropy, blends), more fields, stacking.
## 2026-10-08: Orphan-pair null (e) fixed, companion-aware null (f): no deep-field excess (D-055)
- Cloud run. Hypothesis: the D-051 flanking-field orphan excess (246 vs 211, P = 0.010) comes from null (e)'s
  colour cell (member i only, non-finite colours in the 0–0.3 bin) or from physical companions the 10–30″
  reference misses.
- `pair_cells` uses both members' colour bins (unordered), a non-finite bin and a ≥ 5-pair cell floor (S/N × size
  fallback); null (f) conditions on 3–6″ pairs. Nine fields re-run (≈ 1 min each, in parallel), orphans unchanged.
  Deep fields: (e) P = 0.13, (f) P = 0.35; clusters P ≥ 0.13. **The excess was a null-model artefact.** Table:
  docs/orphan_pairs.md.
- Injections re-run (6 fields in parallel, ~4 min): no-candidate limits reproduce D-051; background-aware limits
  tighten by ×0.76 (docs/exotic_limits.md). Orphan set unchanged, so D-051's inspected contact sheets stand.
- Not separated: the shares of the cell changes.
- **Handoff:** W2 next is the segmentation-map same_galaxy rule (TASKS "Now" 1); PR #78 (W1 shear) in flight.
## 2026-10-08: Survey-agnostic signature layer; any public dataset in scope (D-054)
- Owner direction: find evidence of traversable wormholes / negative-mass objects or warp-drive spacetimes in any
  public dataset. `jwst_anomaly.signatures` registers W1/W2/W3/W5 (prediction, injection, screens, ordinary
  mimics, limits) and defines the `LightCurveSurvey` / `CatalogueSurvey` adapters.
- Reuse-check: MulensModel (extra `mulens`) for ordinary microlensing fits; OGLE-IV Mróz et al. 2019/2020 samples
  first for W3 (with published efficiencies), then Gaia DR3 `vari_microlensing`, then KMTNet; lens catalogues
  (lenscat, Euclid Q1, SuGOHI) for W1/W2. No published survey light-curve search for negative-mass or Ellis
  events exists, so a W3 limit would be new (owner first). OGLE EWS use waits for the owner (its terms).
- **Handoff:** W3 re-analysis of the Mróz OGLE-IV samples (TASKS "Now" 1).

## 2026-10-08: W3 multi-epoch dimming / inverted-microlensing screen: null, limits (D-052)
- Worktree worker: `scripts/dimming_screen.py` cross-matches per-epoch level-3 catalogues (F200W + F444W), flags
  vanishing, achromatic-dimming and rise-dip-rise sources, vetoes catalogue effects (incl. a `bright_neighbour`
  veto: a source ≥ 100× brighter within 1.5″), and confirms with S3 byte-range forced photometry. Fields:
  NEXUS-Center (8 epochs), MACS0416 (5; F200W baseline 0.35 yr, F444W 3.26 yr), Abell 2744 (3). El Gordo rejected
  (one shared band); JADES deferred.
- 5,177 catalogue flags → 1,441 after catalogue tests → 6 confirmed by forced photometry (MACS0416) → 0 after the
  bright-neighbour veto (all six sit on or beside saturated stars; cutout sheet inspected: the star's wings and
  spikes rotate through the aperture with each epoch's mosaic orientation). No candidate.
- Injection (`exotic_sim.inject_light_curve`, n = 1, ε < 0) with per-copy vetoes: efficiency 0.04–0.24. Headline
  95 % limits (MACS0416 only, the one field with calibrated forced errors, 34 sources): rate < 0.056 per source
  per yr at t_E = 1 yr (umbra fraction τ < 0.18); all fields indicative < 0.015 (τ < 0.049). The per-copy
  single-epoch veto makes full vanishes of two-epoch sources undetectable (main efficiency loss).
- **Failed approaches (rules):** a catalogue non-detection is zero flux only if the source would have been ≥ 10σ;
  vanish only when every testable band vanishes; `is_extended == False` is not a point-source cut (use CI_70_30);
  inject multiplicatively with matched noise; Gaia masks must drop the source's own match before looking for a
  neighbour; forced fluxes need each epoch's pixel solid angle; never credit exposure to sources flagged before
  injection; don't run parallel or multi-target S3 cutout jobs through the proxy (s3fs "bucket does not exist").
- Caveats: forced counts are from the 06:39–07:03 UTC run (later S3 failures), re-calibrated offline; NEXUS and
  Abell 2744 forced errors uncalibrated; PEARLS (MAST lists 2.0.1, headers 3.0.0) + CANUCS processing differ.
- Wall time: screen 4–76 s, injections 20–136 s, forced photometry 3–32 min per field.
- Final /code-review fixes (cloud run): the forced-stage SIMBAD/NED label read non-existent `*_otype` columns
  (always empty; now `best_match_*` of `crossmatch.XMATCH_COLUMNS`); `calibrated` is False when an epoch image was
  unreadable; the dead saturated-star branch of `bright_neighbour` removed (`near_star` covers it). Counts and
  limits unchanged (no survivor reached the cross-match). A failed SIMBAD/NED service is now named in
  `forced.ecsv` meta (`n_<service>` == -1), `inject` also treats `unread_images` as uncalibrated, and
  `check_params` guards the saturated-star / near_star coupling (#77).

## 2026-10-08: Orphan-pair cutout footprints: frame-token polygons, deterministic visit order
- Cloud run. Review follow-ups to #73 that its final squash did not carry: `_in_region` parses `POLYGON ICRS …`
  (a frame token used to match nothing, so every cutout read "outside"); MAST observations are sorted by obs_id
  before "first covering visit" is chosen (archive row order no longer picks the cutout). Counts and limits do
  not depend on either (cutouts only). TASKS: null (e) colour cell asymmetry.
- **Handoff:** TASKS "Now" 1, the D-050 W1 shear screen, is next; it needs the SMACS `_cat.ecsv` and DJA zout.

## 2026-10-08: Dark-deflector (orphan-pair) screen in deep fields, with injection-recovery: null (D-051)
- Worktree worker: `orphan_pairs.py` now reads CANUCS and DJA catalogues through one column layout (D-048 counts
  reproduced bit for bit) and runs on six deep fields: the five CANUCS NIRCam flanking fields and DJA GOODS-N
  (108 arcmin²). `scripts/inject_pairs.py` paints `exotic_sim` point-mass, W2 Ellis and W1 pairs into the real
  catalogues and runs the unchanged screen.
- 355 orphans against 315.4 expected (strictest null, P = 0.015); the excess sits in the flanking fields (246 vs
  211, P = 0.010), GOODS-N matches every null. Read as unmodelled physical companions (lensing would need ~140 dark
  galaxy-mass deflectors per arcmin²); open, not a candidate. 90 top orphans inspected: 50 ordinary (knots,
  satellites, groups, artefacts), 40 faint chance-like pairs. Nothing for `/vet-candidate`.
- Per-deflector efficiency ≤ 0.8 % (point mass, W2) and 0.1–18.6 % (W1, θ_E 0.3–1.5″). Background-aware 95 % limits:
  W1 < 1.3 × 10⁴ deg⁻² at θ_E 1.5″ (4.5 × 10¹¹ M☉); W2 < 2.9 × 10⁵ deg⁻² at 0.7″ (throat ≈ 10 pc). This covers
  θ_E ≤ 1.5″, where `radial` (D-049) is blind. docs/exotic_limits.md, docs/orphan_pairs.md.
- **Failed approaches (rules):** DJA Kron apertures for the same-galaxy rule (2.6× CANUCS; use 3.3 × flux_radius);
  the standalone DJA GOODS-N zout (older catalogue; stream the tarball member); a global null (c) in deep fields
  (under-predicts orphans; use the conditioned null (e)); unretried S3 reads (spurious NoSuchBucket via the proxy).
- Downloads > 200 MB (GOODS-N catalogue 224 MB, photo-z tarball 371 MB streamed) stated in D-051.
- Wall time: 3–20 s per field search; 190–266 s per field of injections.

## 2026-10-08: W1 limits re-run: mass-parametrised, blend-aware, independent nulls (D-049 update)
- #70 merged an intermediate version. This re-run fixes the review findings: one lens mass with θ_E per source
  redshift (photo-z, else z_s = 2), overlapping image pairs painted as one blend, an independent 200-draw null per
  batch of 10 trials, screen grid/null/defaults shared with `cmd_radial` (byte-identical SMACS output), a stated
  S/N ≥ 5 detection floor, measured footprint-border excess (3–16 %).
- Recovered of 1,600 per mass (2 × 10¹⁰ / 2 × 10¹¹ / 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉): 0 / 0 / 8 / 84 / 156.
  Headline 95 % limits (six photo-z fields, 38.0 arcmin²): none below 10¹² M☉; < 6.1 × 10⁴, 7.0 × 10³, 4.0 × 10³
  deg⁻² at 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉ (all eight fields, optimistic: 3.7 × 10⁴, 3.8 × 10³, 2.1 × 10³).
- **Failed approaches (rules):** one θ_E for every source; two lines for overlapping images; one fixed null for all
  trials; a `nanmin(S/N)` detection floor; recovery at the 3-line peak instead of p_random; catalogue shapes
  without PSF deconvolution; headline limits including fields without photo-z.
- Wall time: 371–993 s per field (1,000 injections, 100 independent nulls; 4 parallel).
- **Handoff:** a W1-specific screen benchmarked with this harness (TASKS).

## 2026-10-08: W1 shear screen reuse-check (D-050)
- `/reuse-check` for the W1-specific screen: build a catalogue aperture-mass map (−M_ap: a negative-mass lens gives
  negative tangential shear) with scipy cKDTree. TreeCorr NG only gives a stacked ⟨N M_ap⟩(R) and has no Windows
  wheels; lenspack works on pixelised maps. Rough S/N 1–2 at θ_E = 1″ in a 5″ aperture, so gains are expected
  mainly at θ_E ≥ 2–3″.
- **Handoff:** implement D-050 in `exotic_screens.py` and benchmark it with `scripts/inject_radial.py` on SMACS
  first (smallest field; needs its MAST `_cat.ecsv` and DJA zout in the cache).

## 2026-10-08: PR #70 merged after review fixes; Abell 2744 radial doc brought to the post-D-034 result
- Cloud run. PR #70 (D-049) was conflicted with main: merged `origin/main`, then `/code-review` on the final diff.
  Main finding: MACS0717 and Abell S1063 have no photo-z, so the injector painted their members and foreground
  galaxies as W1 images (Abell S1063 had the highest efficiency of all fields). The headline limits now use the six
  photo-z fields (38.1 arcmin²): < 7.0 × 10⁴ / 6.8 × 10³ / 3.2 × 10³ deg⁻² at θ_E = 3″ / 6″ / 10″, about 27× weaker
  than Takahashi & Asada; the all-field values stay in `limits.json` as optimistic. `inject_radial` now reads the
  screen defaults from `exotic_screens.radial_defaults()`.
- docs/fields/abell2744.md: radial result with the D-034 spike veto (134 `anti` arcs, max 5 lines, p 0.505; null).
- PR #61 (n0153, `needs-human`) brought up to date with main; it still waits for the owner.
- **Handoff:** the W1-specific screen (TASKS "Now" 1) is next; it fits one full cycle.

## 2026-10-08: W1 injection-recovery through `radial`: the screen is blind to negative-mass lenses (D-049)
- Worktree worker: `scripts/inject_radial.py` paints `exotic_sim` W1 lenses (n = 1, ε < 0) into the real catalogues
  of all eight null `radial` fields and runs the unchanged screen (`exotic_screens.radial_candidates`, a pure
  refactor of `cmd_radial`'s selection). 200 lenses per field and θ_E; 51.2 arcmin² screened.
- Recovered: 0 / 1,600 at θ_E = 0.3″ and 1″; 10 at 3″; 80 at 6″; 181 at 10″. 95 % limits on W1 lens surface
  density (six photo-z fields): none below 3″; < 7.0 × 10⁴ deg⁻² at 3″ (|M| ≈ 1.4–4.3 × 10¹² M☉), < 3.2 × 10³ deg⁻²
  at 10″. MACS0717 and Abell S1063 (no photo-z, members get lensed) are excluded as optimistic. About 27×
  weaker than Takahashi & Asada spread over 0 < z < 1. docs/exotic_limits.md.
- Why: an image reaches e ≥ 0.5 only for β ≲ 2.3 θ_E, and `anti` against the cluster keeps a third, so a lens puts
  1–3 arcs into a screen whose null needs 5–8 lines.
- **Failed approaches (rules):** recovery is the screen's p_random, not the 3-line peak; cache the null draws once
  per field (0.2–1 s per lens instead of 15–150 s); deconvolve the PSF before applying the lens Jacobian.
- Wall time: 322–877 s per field for 1,000 lenses (5 θ_E).
- **Handoff:** a W1-specific screen (collinear radial pairs flanking an empty centre, orientation against the
  candidate centre, local null), benchmarked with this harness.

## 2026-10-08: El Gordo lens model reproduces the published magnification maps (issue #68)
- `validate` now compares Lenstool-par models with published signed μ maps (`mag_map_files`, fetched only by
  `validate`); `map_check` reports `parity_agree` for signed maps. El Gordo (CDS J/A+A/678/A3, z=2 and z=4, |μ|<10,
  37-px sub-grid): median |μ| ratio 1.00001 / 1.00001, p95 relative difference 0.04 % / 0.07 %, parity 100 %;
  χ²pos unchanged (82.53). At full resolution (250,000 points), >20 % mismatches sit only on critical curves or
  within ~1″ of member cores (0.4″ map pixels).
- El Gordo `bayes.dat` `Chi2` explained (D-045 amendment): the sampling run uses `forme -10`, an image-plane χ² with
  σ² = a·b from the image list (19 of 56 images at 1.24″). It reproduces `Chi2` row by row to 0.1 % and the
  ln(Lhood) offset (75.904) exactly. The chain is validated. **Failed approach:** a source-plane χ² with free
  per-family weights (held-out ρ 0.81) fit only partly; the definition came from Lenstool's source.
- **Handoff:** El Gordo μ, parities and the MCMC chain are validated. Next: posterior μ spreads in `fluxratios`.

## 2026-10-08: Exotic-lens predictions first: wormhole / negative-mass searchable, warp not (D-047)
- Owner focus (2026-10-08): search only for signatures of traversable wormholes / negative-mass lenses and warp-drive
  spacetimes, predictions first. Research worker; every citation fetched from arXiv / Crossref.
- Searchable: W1 negative-mass dark lens (radial pair beside an empty centre), W2 Ellis pair without deflector, W3
  inverted microlensing (umbra between caustic spikes), W5 count deficit. Not searchable: the 4 % Ellis gutter, µas
  shifts, and every warp signature (Alcubierre exterior is flat; no published imaging/photometric prediction for a
  distant observer). The warp branch is stopped.
- `jwst_anomaly.exotic_sim`: Kitamura+2013 power-law lens family (either sign of ε), finite-source light curves,
  `inject_images` / `inject_light_curve` (`simulated`) for injection-recovery.
- Review fixes: an exact finite-source integral, checked against inverse ray shooting, puts the negative-mass spike
  peaks at ×7.0 / 2.35 / 1.53 (ρ = 0.01 / 0.1 / 0.3; first version 9.2 / 2.6 / 1.7). Exact demagnification onset
  added (2/(n+1) is KNA13's large-n estimate); KNA13's n = 3 "~10 %" is rounding of their Fig. 2c (13–14 %). NaN
  epochs no longer read as an umbra.
- **Handoff:** injection-recovery for `radial` (W1) and the dark-lens search (W2) to turn nulls into limits; a
  dimming class for the transient screen (W3); counts around `radial` centres (W5).

## 2026-10-08: Multi-plane lens models (D-046)
- `LensModel.split_planes` moves named potentials (e.g. a foreground galaxy fitted as a member) to their own
  redshift; `MultiPlaneLensModel` solves the multi-plane lens equation, and `find_images` / `backtrace_images` /
  `imageplane_residuals` accept it.
- Reproduces the D-042 system-51 result with library code: σ ≤ 70 km/s at z 0.268 leaves exactly 3 images.
- **Handoff:** use `split_planes` in `/vet-candidate` for any extra or missing image near a non-member galaxy.

## 2026-10-08: Orphan image pairs, a blind dark-deflector screen: null (D-048)
- Worktree worker: `scripts/orphan_pairs.py` looks for SED-matched close pairs with no published system and no visible
  galaxy between them in the CANUCS DR1 catalogues of MACS0416, MACS1149 and Abell 370.
- Pair excesses come from same-redshift groups (null (c)); 38 orphans against 34.5 expected. The 35 top orphans are
  knots, group members or chance matches, all at |μ| ≈ 1–2.5. No candidate.
- **Failed approaches:** a midpoint-only lens rule (missed galaxies between members); no per-band S/N cut (25 % false
  SED matches).
- Wall time: 5–10 s per field per search, 30–60 s with cutouts.
- **Handoff:** injection-recovery so the null becomes a limit; segmentation-map same-galaxy test; more clusters.

## 2026-10-08: Lenstool MCMC posteriors; Abell 2744 multiplicity residual explained (D-045)
- `lensmodel.read_lenstool_bayes` / `posterior_par` and `lens_consistency.py posterior` rebuild published models at
  MCMC samples (potfile rescaling, sampled family redshifts). Validated on Abell 2744: best.par is a chain row, and
  two random rows give χ²pos 173.5 / 173.0 against the chain's 178.1 / 174.1.
- Abell 2744 (12 samples plus best.par): 3.2a/b, 34.1a/b and 700.1a/b stay merged in every model. CATS v4.1 splits
  34.1a/b, and 3.2a/b sit on the caustic in both models. That is model resolution at folds; 0 surviving. 4.2c moves
  only 0.2–0.4″ (μ 8.7–10.0), well below its galaxy-scale systematics (D-036).
- **Failed approach:** the El Gordo CDS chain's `Chi2` column (54–77) does not track our χ²pos (93–106), not even in
  rank, so that chain is not validated.
- Wall time: about 70 s per Abell 2744 sample on a 0.25″ grid; 15 min for 13 models.
- **Handoff:** feed the posterior μ spread into `fluxratios` (then recheck SMACS 6.3 once a SMACS chain is pinned) and
  the position spread into `forced_check` search radii.

## 2026-10-08: CANUCS Lenstool models pinned (D-044)
- `macs0416-canucs` reproduces Lenstool's image-plane χ²pos (330.8 against 344.30; rms 0.51″ over 303 images).
  `abell370-canucs` is a source-plane fit (image-plane rms 2.3″).
- They are the independent second model for vetting (the D-042 system-27 and system-51 checks used scratch code).
- `fluxratios` now refuses gated image lists and map models, as `images` does; its `--offset-arcsec` defaults to
  the model's pinned frame offset.
- **Handoff:** a two-plane option in `LensModel`; use `macs0416-canucs` in `images` runs as a direct cross-check.

## 2026-10-08: Abell 370 and Abell S1063 radial screens: null (D-043)
- Worktree worker: CANUCS 1208 for Abell 370 (CANUCS DR1 photo-z; frame offset −0.121″, −0.015″ pinned) and GLIMPSE
  3293 for S1063 (DJA v7.5 photo-z). Only `radial` applies (image lists gated, D-035).
- Abell 370 raised two flags (15 and 8 lines, p ≤ 0.01). Both are diffraction-spike chains from Gaia stars that are
  off the mosaic or saturated and absent from the catalogue. The Gaia-seeded veto is now in the code
  (`radial --spike-stars`, `scripts/gaia_stars.py`); with it, p 0.945. S1063: p 0.435, unchanged by the veto.
- **Failed approach:** spike vetoes seeded from the pipeline catalogue miss saturated and off-mosaic stars, whose
  spikes reach 37″.
- Wall time: under 2 min of pipeline per field, plus cutout vetting.
- **Handoff:** a low-weight veto or an aper50 S/N floor in `radial`; use `--spike-stars` on every field.

## 2026-10-08: MACS0416 counter-images and flux ratios: null; CANUCS model as second model; two-plane check (D-042)
- With the D-040 solver, the MACS0416 image list is open. `images` and `fluxratio` raised 21 flags; 0 survive. The most
  persistent was system 51's fourth image: two independent models predict it, and the 51.1–51.3 photometry says it
  should appear at 9–25σ, but it is not seen. It comes from a z 0.268 foreground galaxy modelled as a member. A
  scratch two-plane model with that galaxy at σ ≤ 60 km/s gives exactly the 3 observed images.
- System 27's two bright `absent` predictions (S/N 247–341) are CATS-only galaxy-scale caustics. The JWST-era CANUCS
  Lenstool model (222 potentials, 111 spectroscopic systems) reproduces system 27 with exactly its 3 images, and the
  cutouts show empty sky there.
- Flux-ratio flags: 45.2 is blended with a bright galaxy 0.5″ away; 38.1 is marginal.
- Wall time: about 3 min of pipeline plus about 5 min of vetting.
- **Handoff:** pin the CANUCS models (MACS0416, MACS1149, Abell 370) as `MODELS` entries for second-model vetting.

## 2026-10-08: Sunrise transient candidates `n0022` and `n0150` are detector persistence (D-039)
- New `scripts/persistence_check.py`: per-exposure photometry on level-2 `_cal` files (S3 byte ranges), plus the
  same detector pixel in earlier exposures on that detector. Validated on a synthetic afterimage and on two real
  sources (`n0153`, Earendel: 0 suspect detections).
- `n0022`: afterimages of a bright galaxy (2.0 % and 0.8 % of its flux) in o010 dithers 3 and 4; the dither geometry
  stacks both on one sky position. `n0150`: in all three epochs an afterimage (0.04–0.07 %) of a saturated star, one
  exposure later; its "motion" is the dither vector. Cutouts inspected (dithers 1–4 at fixed detector pixels).
- `n0153` is on sky in all 8 dithers and F150W is ×1.9 brighter in 2025 at the `_cal` level: still a variable
  candidate.
- **Lesson:** a single-epoch source in a mosaic can be an afterimage that two dithers place on one sky position.
  The archived `_cal` files carry no DQ flag there.
- **Handoff:** `/vet-candidate n0153`; persistence-check new single-epoch candidates before vetting.

## 2026-10-08: MACS0717 screens: null (D-041)
- Worktree worker on VENUS 6882 o029 (10 bands, the only public NIRCam association) with `macs0717-cats` (rms 3.21″).
- 51 flags, 0 surviving. 29 are the model's own copies of catalogued images it does not reproduce (1.6–5.6″ off), 6 are
  CATS-only extra images (Sharon v4cor and Keeton v4 predict none), and 5 have μ more than 2× model-dependent.
  System 65's flux ratio is a 0.6″ catalogue offset. `radial` p ≥ 0.70.
- **Failed approach:** a fixed 1.5″ match radius for a 3.2″-rms model makes most catalogued images "unpredicted",
  and their model copies then flag as absent or confused.
- Wall time: about 7.5 min of pipeline plus about 20 min of vetting.
- **Handoff:** in `forced_check`, rms-scaled match and search radii and automatic copy classification; Sharon v4cor and
  Keeton v4 as pinned `MapLensModel` entries.

## 2026-10-08: `find_images` fold refinement; frame offsets for map models (D-040)
- Cells on a critical curve near the source are subdivided into ≤ 0.02″ sub-cells. MACS0416 system 26 is now solved
  (rms 1.57″ → 0.76″), and its image list is open with offset (0.208, −0.025). No other model changes beyond 0.03″
  rms; Abell 2744 `validate` takes 41 s instead of 29 s.
- **Bug fixed:** `apply_frame_offset` was a no-op for map models (the maps are looked up by sky position), and
  `radial` / `arcs` skipped it for them. No earlier result used a map model with an offset. MACS0416's radial screen
  was re-run in the JWST frame: still null.
- **Handoff:** `images` / `fluxratio` on MACS0416 (CANUCS photo-z); radial on Abell 370 and Abell S1063.

## 2026-10-08: MACS0416: system 26 is a solver-grid artefact; radial null (D-038)
- Worktree worker on CANUCS 1208 (`jw01208-o004_t002`, 8 bands). The CATS system-26 residual (11″) is the 0.25″
  `find_images` grid missing a merging pair near the critical curve; the source lies 0.001–0.005″ from the caustic. On a 0.1″ grid the rms is 0.811″ (1.13×
  quoted). The model stays map-only until the solver refines its grid near high |μ|.
- `radial`: 98 centres against 99.6 random; the 7-line centre (p 0.29) is low-S/N noise segments. With CANUCS DR1
  photo-z (background cut): 77 against 85.0, p 0.225. 0 flags.
- **Failed approach:** jwst 3.0.0 isophotal S/N admits noise segments (66 of 137 anti arcs have aper50 S/N < 3). The
  radial screen needs an aper50 S/N floor.
- Wall time: about 4.5 min (lens and exotic scripts plus vetting stamps; no pipeline `run`).
- **Handoff:** adaptive grid refinement in `find_images`; an aper50 S/N floor in `radial`.

## 2026-10-08: CANUCS DR1 photo-z for MACS1149 (D-037 addendum)
- **Failed approach:** the worker checked only DJA, but docs/landscape.md already listed CANUCS DR1 (PSF-matched
  EAzY photo-z for A370, MACS0416, MACS0417, MACS1149 and MACS1423). Check docs/landscape.md before reporting that a
  field has no photo-z.
- MACS1149 system 16.2: z_phot 2.25 (95 % 0.23–2.33). This disfavours the CATS-fitted 4.419.
- Radial re-run with the background cut: 12 peaks against a null mean of 12.8; null. 6 of 68 matched lensed images
  get a blended z < 0.6.

## 2026-10-08: MACS1149 screens: null (D-037)
- Worktree worker on CANUCS 1208 (8 bands) with `macs1149-cats`. The image-plane rms is 0.673″ (gate passes) and
  the frame offset is under 0.02″.
- 3 flags, 0 surviving: system 16 (fitted z 4.419; at z 2.5 in CATS or 3.0 in Sharon v4cor the image lands on 16.2), system 2 (CATS topology error,
  checked against Sharon v4cor), Refsdal-host knot 1192 (next to the BCG; μ differs >2× between models).
  `fluxratio` 0 flags; `radial` p = 0.91.
- **Failed approach:** forced photometry with an annulus on a BCG core gives negative fluxes; use high-pass.
- No DJA mosaic for MACS1149 (v7).
- Wall time: about 5 min of pipeline plus about 10 min of vetting.
- **Handoff:** `macs1149-sharon` (or the CANUCS models) in `MODELS`; the D-037 rules in `forced_check`; system 16's
  redshift (spectroscopic z, or forced photometry on its far-image track); `images` with CANUCS photo-z.

## 2026-10-08: Abell 2744 screens: null (D-036)
- A worktree worker ran `validate`, `images --forced-image`, `fluxratio` and `radial` on Bergamini+2023b, with
  cutouts of every flag. 15 of 16 flags were ordinary. The survivor, 4.2c, was re-run under the D-034 rules: it
  is `no_reference` (resolved-knot references, BCG halo), not a candidate.
- Radial: 35 peaks against a random mean of 50.4, p ≥ 0.575.
- **Lesson:** a "252σ absent" image can come from a resolved-knot reference. D-034's compact-reference rule now
  catches this.
- Family 4's c images were vetted (4.1c 4–8× underluminous after BCG subtraction; 4.2c undetected). They are
  explained by μ systematics next to member 34423 (3.9–28.7 under ±30 %; CATS 7.3). **Rule:** a μ that moves by
  more than 2× under member perturbation is untestable.
- **Handoff:** HFF field runs (D-035 models), `bayes.dat` position spreads, UNCOVER v2.0 cross-check.

## 2026-10-08: Sunrise third epoch with calibrated significances; two open transient candidates (D-027)
- `transient_forced.py --controls`: noise scale from ordinary sources, applied before thresholding.
- Sunrise o010 against VENUS o052 (2.9 yr; F150W, F444W): scales 1.30/1.18; **0 of 57 catalog candidates pass**;
  Earendel steady (Δm ≤ 0.12 mag, < 1σ).
- Among the 200 controls, `n0022` (gone after 2022-07 in four SW bands), `n0150` (a different position in each
  epoch) and `n0153` (+0.5–0.6 mag in both bands) change. Not vetted; ordinary explanations (supernova, moving
  object, AGN, edge artefact) are untested. docs/fields/sunrise.md has the numbers.
- **Failed approach:** the catalog stage plus the two-band rule misses single-pair, blue transients.
- Follow-up #53: an empty control selection is flagged uncalibrated instead of being skipped silently.
- **Handoff:** `/vet-candidate` for `n0022` and `n0150` (level-2 `_cal` exposures per filter: is the source in
  every dither? Does `n0150` move within one visit?); then a grid of forced photometry (all sources, not only
  catalog candidates) for every epoch pair.

## 2026-10-08: Six HFF clusters as CATS map models (D-035)
- MACS0416, MACS1149, Abell 370, MACS0717, Abell S1063 and Abell 2744 (CATS v4/v4.1) are pinned map models. All six
  reproduce their published z = 2 magnification maps (median 9e-5 to 6e-3).
- Image lists pass the rms gate for MACS1149 and MACS0717. MACS0416 has a `params.txt` but fails the gate (system 26).
  Abell 370, Abell S1063 and Abell 2744 have no `params.txt`, and their placeholder redshifts give 9–12″
  residuals. All four stay map-only.
- **Failed approach:** using CATS `arcs.txt` redshifts as published.
- **Handoff:**
  - per-cluster JWST field runs (configs, catalogs, `radial` on all six, `images` / `fluxratio` on MACS1149 and
    MACS0717);
  - `/vet-candidate` for Abell 2744 4.2c (the field worker's survivor).

## 2026-10-08: El Gordo counter-images and radial screen: no candidate (D-034)
- Forced photometry: 3 `absent` images on the first pass, all ordinary on vetting (cutouts plus numbers):
  - 23: model position error;
  - 6: reference on a galaxy wing; photo-z-consistent counterparts at 2.6–3.1″;
  - 7: HST→JWST frame offset of 0.22″.
- Radial screen: max 4 lines, p = 0.64. A 6-line "centre" was a star's diffraction spikes.
- Four new rules (D-034): frame offset, compact and consistent references, a residual-scaled search radius, and a
  spike-segment veto. SMACS re-run: unchanged, 0 absent.
- Tally: El Gordo screened 17 uncatalogued predicted images and 37 anti arcs; flags 3 + 1; **surviving 0**.
  Wall time: about 50 s for `images`, 21 s for `radial`.
- **Handoff:** Abell 2744 (worker running), then the HFF/RELICS map fields via `MapLensModel` (#47).

## 2026-10-08: Map-based lens models; Sunrise radial screen null (D-033)
- `lensmodel.MapLensModel` evaluates published deflection maps. The first one is `whl0137-relics-lenstool`, which
  reproduces RELICS κ to 3.5e-5 and μ(z = 6.2) to 2.0e-5 (medians).
- Sunrise `exotic_screens radial`: 29 usable anti arcs, 3 centres against a null mean of 2.2, max 3 lines, p 0.885.
  **No candidate.** Wall time 17 s, with cached maps and catalogs.
- Cycle tally:
  - SMACS + Sunrise screened: 60 + 0 images, 34 + 29 arcs;
  - flags: SMACS fluxratio 2 (removed by the two-band compactness rule), radial 0 significant;
  - surviving vetting: 0.
- **Handoff:** more map fields (HFF: Abell 2744, MACS0416, MACS1149, Abell 370; RELICS clusters), El Gordo and
  Abell 2744 `images --forced-image` and screens.

## 2026-10-08: Flux-ratio and colour test of catalogued images: no anomaly in SMACS or El Gordo (D-032)
- New: `lens_consistency.py fluxratios`.
  - It compares each image's DJA `mag_auto` + 2.5 log|μ| and its F150W−F444W colour with the other images of its
    system (leave-one-out).
  - It drops blends, segments larger than 20,000 px and counterparts whose photo-z excludes the system redshift.
- Results (`derived`):
  - SMACS: 8 images in 4 systems flux-tested, rms 0.60 mag.
  - El Gordo (0.3″ after a +0.221″ frame shift): 23 images in 11 systems, rms 0.49 mag.
  - No colour outlier (rms 0.11 and 0.07 mag).
  - The flagged pairs SMACS 6.1/6.3 and El Gordo 18b/18c (and 7b/7c, at 0.5″ without the shift) fail vetting.
    Forced photometry brings SMACS 6 down to 0.61–0.78 mag. The El Gordo pairs are low-S/N or chromatic, so they
    are measurement failures. **No anomaly.**
- **Failed approaches**, each checked on cutouts or SEDs:
  - DJA aperture × `tot_corr` is not a total flux for arcs.
  - SMACS 1.1 looked 1.4 mag too bright, from host-halo light in a 32,864 px segment.
  - El Gordo 9a/9c differed by 1.2 mag in colour, because 9a's DJA counterpart is a z_phot 0.89 interloper.
  - Bare eazy 95 % intervals exclude good images.
  - A median that includes the image itself halves pair differences.
- Limit: DJA misses most arcs inside cluster light (SMACS: 38 of 60 images unmatched within 0.3″).
- **Handoff:**
  - El Gordo `images --forced-image` (TASKS "Now" 2);
  - Abell 2744 `fluxratios` with the 233 MB DJA catalogue (needs a DECISIONS entry);
  - `bayes.dat` μ uncertainties, then recheck SMACS 6;
  - BCG/ICL-subtracted totals, to test the core images.
## 2026-10-08: Exotic-lens screens; SMACS null (D-031)
- New `scripts/exotic_screens.py`:
  - `fluxratio`: two-band forced photometry, luminosity ratio against sibling images, compactness and chromatic
    gates;
  - `radial`: anti-tangential arcs whose axes converge on a dark centre, with a false-alarm rate from randomised
    position angles.
- SMACS:
  - fluxratio: 6 compact images consistent, 0 flags. A first, single-band compactness gate had flagged system 7,
    an ordinary knot-vs-whole-arc mismatch;
  - radial (background sources only): 12 centres against a null mean of 9.0 (p95 15); max 4 lines, p 0.945.
  - **No exotic candidate.**
- **Failed approaches (now rules):**
  - fixed-aperture flux ratios on resolved arcs: surface brightness is conserved, so the ratios scale with 1/|μ|
    (systems 5 and 10);
  - a single-band compactness gate (knots of a clumpy arc pass in F150W);
  - a uniform-angle null for the radial screen (the selected arcs point at the mass centre).
- **Handoff:** fan out per cluster (El Gordo, Abell 2744, Sunrise, then HFF/RELICS). Each runs `validate`,
  `images --forced-image`, `exotic_screens fluxratio` and `radial`, with contact sheets of all flags.

## 2026-10-08: `find_images` 3–5× faster with identical images
- Seeds are pre-filtered with boolean sign tests on the mapped grid corners, and Newton steps run for every seed in
  one `fields_xy` call.
- Benchmark: all catalogued systems at a 0.25″ grid, images identical (max |Δ| 0 arcsec, same counts):
  - SMACS: 23.2 → 6.7 s;
  - El Gordo: 44.3 → 13.0 s;
  - Abell 2744: 118 → 23.5 s (at ±190″).
- **Handoff:** the deflection grid itself (one-time and cached) is now the main cost. Published deflection maps
  (UNCOVER, RELICS, HFF) could replace it for fields without a Lenstool model (TASKS).

## 2026-10-08: Image-plane χ² reproduces Lenstool for SMACS, El Gordo and Abell 2744 (D-030)
- Merged #40 (counter-images, D-029) after its last commit, which GitHub had not attached to the PR, was picked up
  by a follow-up commit.
- `validate` now computes the exact image-plane χ² for every model; El Gordo and Abell 2744 are `MODELS` entries.
  Issue #41's numbers are now reproducible in the repository: χ² 30.87/30.91 (SMACS), 82.53/80.22 (El Gordo),
  146.64/146.60 (Abell 2744). No catalogued image is off by more than 3σ in any field.
- Parser fixes: letter-suffixed image ids, several ids per `z_m_limit`, 6-decimal `_kpc` rounding.
- **Failed approach:** the per-image error column for El Gordo (χ² 52.0); a uniform 0.621″ matches Lenstool.
- CDS reset connections through the proxy in this run; the cache was seeded from a `curl` copy with the same sha256.
- **Handoff:** run `lens_consistency.py images` with forced photometry on El Gordo and Abell 2744 (TASKS "Now" 2);
  then the exotic screens (docs/exotic_lensing.md).

## 2026-10-08: SMACS counter-images: no predicted image is absent (D-029)
- New: `lensmodel.find_images` (an image-plane solver on a cached deflection grid), `lens_consistency.py images`, and
  `--forced-image`, which runs forced photometry on S3 byte-range stamps.
- Results:
  - ICLv2 reproduces all 60 catalogued images within 0.04–0.91″, with 4 demagnified central images.
  - Of the 11 testable uncatalogued images: 3 are recovered (systems 9 and 8, and system 17 marginally: flux ratio
    0.36 on the BCG gradient), 1 is
    confused, 1 is undetectable, 6 have no reference flux, and **0 are absent**.
  - System 8's model z = 11.76 is contradicted by its F090W/F150W detections.
- **Failed approach:** pipeline-catalog flux references near cluster galaxies. They falsely made system 9 `missing`.
- The Mahler κ tarball now comes from raw.githubusercontent.com, because github.com/raw returns 403 in cloud runs.
- **Handoff:** TASKS "Now" 1 is done for counter-images. Next:
  - parser fixes (issue #41), then the same test on El Gordo and Abell 2744;
  - the exotic screens of docs/exotic_lensing.md (demagnified images, radial images around a dark centre).

## 2026-10-08: El Gordo and Abell 2744 lens models, preliminary: no image-position anomaly (issue #41)
- Preliminary validation on published products. The numbers come from scratch code and are not yet reproducible
  in the repository; they need a re-run with #40's solver.
  - `lensmodel.py` reproduces the Caminha+2023 El Gordo magnification maps (median |Δμ|/μ 4e-5).
  - Exact image-plane solves match Lenstool's χ²: El Gordo 82.5 vs 80.22 (rms 0.754″, paper 0.75″); Abell 2744
    (Bergamini+2023b) 146.64 vs 146.60 (rms 0.427″).
- Results:
  - No image-position anomaly in either field.
  - Predicted uncatalogued counter-images: El Gordo 16, inconclusive, because the MUSE Lyα images are too faint in
    continuum. Abell 2744 30, not yet testable, because the pipeline catalogs miss most arcs in the core.
  - docs/fields/elgordo.md, docs/fields/abell2744.md.
- **Failed approach:** source-plane back-trace χ² for image-plane-optimised models (El Gordo 121.6 against 80.22).
- **Bugs found on `main`:** the parser fails on letter-suffixed image IDs and on 6-decimal `_kpc` rounding (#41).
- Also merged #38 after bringing it up to date with `main`.
- **Handoff:** #40 (image solver, local session); then reproduce these numbers in the repository, the parser fixes,
  and the counter-image flux test on DJA photometry (TASKS "Now" 2).

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
- **Handoff:** lens-model checks per field (#35 merged; El Gordo and Abell 2744 model validation in issue #41).
  Time domain: rescale the thresholds by the control std, and add VENUS 6882 o052 (F150W, F444W) as a third epoch.

## 2026-10-08: Cloud runs use the GitHub MCP tools; conflicting PRs get no CI (D-028)
- The first routine run opened #38 with the session's GitHub MCP tools. The prompt, the research-cycle skill and
  docs/operations.md §3 now prefer them, merge with `mcp__github__merge_pull_request` after checking the policy,
  and fall back to REST plus the `merge-ready` label.
- **Failed approach:** #37, a REST squash-merge allow rule, was closed. The glob `pulls/*/merge` also matches other
  `pulls/...` writes, so the Bash merge API stays denied.
- **Lesson:** a PR that conflicts with `main` gets no `pull_request` CI at all (#38). Merge `origin/main` in first.
- The cloud merge gate checks author, head repository (no forks), branch, labels, every page of files and the CI jobs
  (not the skipped `claude` runs), and pins `expectedHeadSha`. D-028 records the decision.
- **Handoff:** TASKS "Now" 1 (counter-images, then parity and flux ratios; see the lens-model entry below). #38
  needs `origin/main` merged in before its CI and merge.

## 2026-10-08: Lens-model stage validated; SMACS arc orientations agree with ICLv2 (D-024)
- `lensmodel.py` is a Lenstool dPIE port (from PyAutoGalaxy, MIT) that evaluates a published `best.par`. It
  reproduces Mahler+2022 ICLv2:
  - κ map: median |Δκ| 1.6e-5;
  - back-trace χ² of the 60 catalogued images: 31.18, against Lenstool's 30.91.
- `scripts/lens_consistency.py` provides `validate` and `arcs`.
- SMACS arcs: 21 of 25 elongated strong-shear background sources, each tested at its own photo-z range, are
  aligned with the predicted stretch (p = 2.6e-7). The 6 anti candidates are all ordinary:
  - 2 segmentation blends;
  - 1 galaxy at z ≈ 0.77 with an intrinsic shape;
  - 3 noisy low-surface-brightness shapes.
  This is a null result.
- **Lesson:** moment orientations from pipeline segments pick up blends and low-S/N shapes, so filter them before
  calling a source anti-tangential.
- **Handoff:** TASKS "Now" 1, which leaves counter-images (predicted but missing, or observed but not predicted)
  and parity/flux ratios. Then El Gordo and Abell 2744; check their model profiles first.


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
