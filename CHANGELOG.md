# Changelog

Newest first. Results, failed approaches and the handoff state; not a diary.

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
