# Exotic-lens limits from injection-recovery

What the null exotic screens (docs/exotic_lensing.md "Screens and results") rule out. We inject simulated signals
into the real catalogues and run the unchanged screens on them. Every injected object is **simulated**. Recovery
efficiencies and limits are **derived** from those simulations. θ_E for a given mass is a **model_prediction** of
the point-lens formula. A limit describes how sensitive a screen is to a hypothetical lens population. It is not
evidence that any such lens exists or does not exist beyond what the stated efficiency covers.

## W1 negative-mass lenses (radial screen)

D-049; `scripts/inject_radial.py`, tests in `tests/test_inject_radial.py`. The signature is D-047's W1: an n = 1 lens
with ε < 0. There are no images of sources within β < 2 θ_E (the umbra), and two radially stretched images on the
source's side of the lens for every source beyond that. The lens **mass |M|** is the parameter. θ_E at z_s = 2 is
given only as a label.

### Method

1. **Fields and inputs.** The eight fields where `exotic_screens.py radial` came back null. Each uses the
   catalogue, photo-z, `--max-radius` and spike veto from its field doc (`FIELDS` in the script).
   - The models are those of the field docs, in the JWST frame.
   - The screen thresholds are read from `exotic_screens.radial_defaults()`, the CLI defaults.
   - The base screen reproduces each documented result (table below). Abell 2744 is the exception: its doc
     predates the D-034 spike veto. The current code drops 42 spike segments and gets 134 arcs, max 5 lines,
     p 0.505 (null).
   - Abell S1063 runs without photo-z, the configuration of D-043 that needs no 75 MB tarball. The documented
     result with photo-z is also null.
2. **Lens.** A lens of mass |M| ∈ {2 × 10¹⁰, 2 × 10¹¹, 2 × 10¹², 8 × 10¹², 2 × 10¹³} M☉ sits at the cluster
   redshift (ASSUMPTION), at a uniform random position in the **screened footprint**. The footprint is the set of
   1″ grid points within `--max-radius` that have a catalogue source within 4″ (ASSUMPTION).
3. **θ_E per source.** Each lensed row gets its own θ_E = θ_E(|M|, z_l, z_s), with θ_E² ∝ |M| D_LS / (D_L D_S)
   (Planck18; `theta_e_arcsec`).
   - z_s is the row's photo-z (`z_phot`, nearest match within 0.3″).
   - Rows without a photo-z take z_s = 2 (ASSUMPTION; the screen's middle grid redshift). They are 30–54 % of the
     lensable rows in fields with photo-z, and all rows in MACS0717 and S1063. The median photo-z of the matched
     lensable rows is 1.6–2.3.
4. **Lensed sources.** The lenses act on the field's own catalogue rows behind the lens position.
   - A row counts as background under the screen's own rule: either it has no photo-z, or its photo-z z160 is
     above z_l + 0.1. Point sources brighter than AB 20 (stars) are never lensed.
   - Rows with β < 2 θ_E are removed (umbra).
   - Rows with 2 ≤ β ≤ 4 θ_E (the D-047 range, ASSUMPTION) are replaced by their images from
     `exotic_sim.inject_images(n = 1, sign = −1)`.
   - Rows beyond 4 θ_E are left unchanged.
5. **Painting** (ASSUMPTIONs):
   - **Shape.** Start from the source's catalogue second moments. Remove an isotropic PSF, map the result by A⁻¹,
     then add the PSF back. A⁻¹ has the signed eigenvalues 1/λ_r along the lens direction and 1/λ_t across it.
     The PSF σ is the 1st percentile of semimajor × (1 − ellipticity) among S/N > 50 rows: 1.0–1.44 px
     (**derived**).
   - **Brightness and noise.** Magnitude − 2.5 log₁₀|μ|; isophotal area × |μ|; S/N × √|μ|.
   - **Blending.** The two images of one source are painted as a single row when their separation is below the
     sum of their isophotal equivalent radii √(area |μ| / π). That stands in for the catalogue's deblending scale:
     overlapping isophotes form one segment. The merged row sits at the flux-weighted centroid. Its second moments
     are the flux-weighted moments plus the spread of the two centres, its flux is the sum, and its
     S/N = S/N_src √(Σ|μ|).
     - On average 0.01–0.56 merges happen per lens, rising with |M|. Merges happen near the caustic, where the
       stretch is largest.
   - **Detection.** Rows with S/N < 5 are dropped as undetected (ASSUMPTION). This only affects the image
     bookkeeping (`n_images`): the screen itself requires S/N ≥ 10. The catalogues' faintest rows have
     S/N 2.4–6.9, and the screen itself needs S/N ≥ 10.
   - The pixel scale, 0.0312″ in all eight catalogues, is measured from their pixel and sky centroids.
6. **Screen.** The injected catalogue goes through `exotic_screens.radial_candidates`, the code `cmd_radial`
   itself uses. Then come `radial_grid`, `line_counts` and `convergence_peaks`, with the same defaults.
   - The null uses `anti_window_draw`, the screen's formula: 200 draws of each arc's angle inside its `anti`
     window.
   - **Independence:** trials run in batches of 10. Each batch draws its own independent 200-draw null for the
     real arcs. Each trial then subtracts the arcs it replaced and adds fresh draws for its injected arcs.
   - Only the 32 × 32-cell grid blocks that those arcs' windows touch are recomputed. The per-draw maximum outside
     them comes from cached block maxima, so the result is identical to a full recompute. This is tested against
     `line_counts` on the whole grid, and so is the windowed null against the screen's own null with its seed.
   - The 10 trials of a batch share that batch's real-arc draws, so they are not fully independent: a batch with
     an unusually high null maximum lowers all 10 together. With 20 independent batches per field and mass, the
     binomial σ below is approximate (slightly understated).
7. **Recovery.** A lens counts as recovered when a convergence peak with p_random < 0.05 lies within 2″ of the
   injected centre (ASSUMPTION: 4 grid steps, twice the 1″ line tolerance).
   - We used 200 lenses per field and mass, which gives a binomial σ of at most 3.5 %.
8. **Limit.** The real screens found zero surviving detections. The 95 % Poisson upper limit on the surface
   density of W1 lenses of mass |M| is Σ < 2.996 / Σ_fields(ε_f A_f), where ε_f is the field's recovery efficiency
   and A_f its screened area.

Command: `python scripts/inject_radial.py` (all eight fields, the five masses, 200 lenses each). It writes
`trials.ecsv` and `summary.json` per field and `limits.json` to `outputs/inject_radial/`. Abell 370 fetches its
Gaia DR3 stars with `scripts/gaia_stars.py` (7 stars).

### Base screens and screened areas (`derived`)

| Field (model) | z_l | Screened area A_f (arcmin²) | Border excess (upper bound) | `anti` arcs | max lines, p_random | lines needed for p < 0.05 | Wall time (s) |
|---|---|---|---|---|---|---|---|
| SMACS 0723 (ICLv2, DJA photo-z) | 0.39 | 2.99 | 3.3 % | 31 | 4, 0.965 | 7 | 492 |
| El Gordo (Caminha+23, DJA) | 0.87 | 4.83 | 11.1 % | 37 | 4, 0.64 | 6 | 550 |
| Abell 2744 (Bergamini+23, DJA) | 0.31 | 11.51 | 4.9 % | 134 | 5, 0.505 | 6 | 1061 |
| MACS0416 (CATS, CANUCS photo-z) | 0.40 | 6.17 | 11.2 % | 120 | 7, 0.225 | 8 | 672 |
| MACS1149 (CATS, CANUCS) | 0.54 | 6.44 | 10.9 % | 61 | 3, 1.0 | 5 | 405 |
| MACS0717 (CATS, no photo-z) | 0.55 | 7.58 | 16.2 % | 68 | 5, 0.70 | 7 | 595 |
| Abell 370 (CATS, CANUCS, Gaia veto) | 0.38 | 6.12 | 10.7 % | 100 | 5, 0.495 | 7 | 600 |
| Abell S1063 (CATS, no photo-z) | 0.35 | 5.62 | 14.0 % | 46 | 4, 0.435 | 5 | 427 |
| **Total** | | **51.2** (0.0142 deg²) | 46.0 arcmin² at r → 0 | | | | |

**Border excess** is measured, not assumed (`footprint_border`).
- The footprint area A(r) is computed for source distances r = 4, 5, 6 and 8″. Once interior holes are filled it
  grows by about perimeter × dr.
- A linear fit extrapolated to r = 0 gives the area enclosed by the outermost catalogue sources.
- A(4″) exceeds that area by the percentage in the table. This is an upper bound on the overestimate, because the
  true mosaic edge lies beyond the outermost sources.

**Wall time** covers all five masses (1,000 injections, 100 independent nulls), with four fields running in
parallel on 4 cores.

### Recovery efficiency (recovered / 200; `derived` from `simulated` injections)

| Field | 2 × 10¹⁰ M☉ | 2 × 10¹¹ | 2 × 10¹² | 8 × 10¹² | 2 × 10¹³ | θ_E(z_s = 2) at 2 × 10¹² |
|---|---|---|---|---|---|---|
| SMACS 0723 | 0 | 0 | 0 | 1 | 15 | 3.19″ |
| El Gordo | 0 | 0 | 0 | 0 | 6 | 2.05″ |
| Abell 2744 | 0 | 0 | 2 | 12 | 20 | 3.59″ |
| MACS0416 | 0 | 0 | 1 | 7 | 15 | 3.17″ |
| MACS1149 | 0 | 0 | 0 | 4 | 6 | 2.70″ |
| MACS0717 | 0 | 0 | 0 | 8 | 24 | 2.70″ |
| Abell 370 | 0 | 0 | 1 | 16 | 18 | 3.26″ |
| Abell S1063 | 0 | 0 | 4 | 36 | 52 | 3.38″ |
| **All (of 1,600)** | **0** | **0** | **8** | **84** | **156** | |
| Mean lensed sources per lens (range; previous run, same injections) | 0.1–0.4 | 0.7–4.0 | 7–32 | 26–126 | 58–289 | |
| Mean injected arcs entering the screen | ≤ 0.04 | 0.07–0.36 | 0.5–2.7 | 2.3–10.7 | 4.8–23.7 | |
| Lenses with ≥ 3 lines at the centre | 0–3.5 % | 0–4.5 % | 2.5–20.5 % | 9.5–51.5 % | 13–46 % | |

θ_E scales as √|M| (0.1×, 0.32×, 2× and 3.16× the last column for the other masses). El Gordo, at z_l = 0.87, has
the smallest θ_E and the fewest lensed sources.

### Upper limits (95 %, zero detections; `derived`)

**Headline: the six fields with photo-z** (SMACS 0723, El Gordo, Abell 2744, MACS0416, MACS1149, Abell 370;
38.0 arcmin²). In MACS0717 and Abell S1063 every non-star row is lensable at z_s = 2, so cluster members and
foreground galaxies get painted as W1 images and the efficiency is biased high; those two fields enter only the
optimistic all-field set (`combine` in the script).

| \|M\| (M☉) | Σ ε_f A_f, photo-z fields (deg²) | Σ(W1 lenses) < (deg⁻²), headline | all eight fields, optimistic | all eight, border-corrected areas |
|---|---|---|---|---|
| 2 × 10¹⁰ | 0 | no limit (ε = 0 in every field) | no limit | — |
| 2 × 10¹¹ | 0 | no limit | no limit | — |
| 2 × 10¹² | 4.9 × 10⁻⁵ | 6.1 × 10⁴ | 3.7 × 10⁴ | 4.1 × 10⁴ |
| 8 × 10¹² | 4.3 × 10⁻⁴ | 7.0 × 10³ | 3.8 × 10³ | 4.2 × 10³ |
| 2 × 10¹³ | 7.6 × 10⁻⁴ | 4.0 × 10³ | 2.1 × 10³ | 2.4 × 10³ |

θ_E(z_s = 2) across fields (`model_prediction`): 0.21–0.36″, 0.65–1.14″, 2.05–3.59″, 4.11–7.18″, 6.49–11.35″ for the
five masses.

The border-corrected column uses the r → 0 area of each footprint (the conservative version of the all-field set). θ_E(|M| = 10¹², z_l = 0.4,
z_s = 2) = 2.23″ reproduces D-047's 2.2″ (tested).

**Reading.** The radial screen is effectively blind to W1 lenses below about 10¹² M☉: none of 3,200 injections
at 2 × 10¹⁰ and 2 × 10¹¹ M☉ was recovered. At 2 × 10¹² M☉ it recovered 8 of 1,600.
- **Comparison (`derived`, rough).** Takahashi & Asada (2013) limit negative masses above 10¹² M☉ to
  n < 10⁻⁴ h³ Mpc⁻³. Spread over 0 < z < 1 (Planck18 comoving volume, 3.98 × 10⁶ Mpc³ deg⁻²), that is about
  120 deg⁻².
- Our best headline limit, 4.0 × 10³ deg⁻² at 2 × 10¹³ M☉ (2.1 × 10³ optimistic), is about 30× (18×) weaker. It
  also holds only for lenses near the cluster redshift.

**Why the efficiency is low.** The cut that dominates is the screen's elongation cut, not a lack of sources.
- For a round source, the image axis ratio is (x² + 1)/(x² − 1), where x is the outer image's position in θ_E.
  Ellipticity ≥ 0.5 needs β ≲ 2.31 θ_E, about 11 % of the [2, 4] θ_E annulus.
- `anti`, measured against the *cluster's* tangential direction, keeps only images whose direction to the W1 lens
  lies within ±30° of the cluster-radial direction. That is about a third of them.
- The p < 0.05 threshold comes from each field's own null and needs 5–8 converging lines.
- In practice, 0.5–2.7 injected arcs reach the screen at 2 × 10¹² M☉.

### Strong-lensing region

A W1 lens inside a cluster's own radial-arc region is degenerate with the cluster's radial arcs (D-047 "ordinary
mimics"). The screen settles this in favour of the ordinary explanation:
- arcs with a model radial magnification ≥ 3, at any redshift of their class, are vetoed;
- near the critical curves, `anti` arcs point at the cluster centre.

Injected lenses there are still counted, at the efficiency they actually get. The screened area therefore includes
the strong-lensing region, where W1 is unrecoverable by construction. We did not cut it out, so as not to depend on
the model.

### Limits and caveats

- **Thresholds.** Every screen threshold is an ASSUMPTION (D-031): e ≥ 0.5, semimajor σ ≥ 2 px, S/N ≥ 10, 1″ line
  tolerance, 15″ line length, the 60° `anti` window, a global-maximum null over 200 draws, and model-radial
  μ_r < 3. So are the injection choices above. The limits hold only for this screen with these settings.
- **Blending.** Image pairs whose isophotes overlap are painted as one blended row (method, step 5).
  - The isophotal equivalent radius is a proxy for the pipeline's deblending scale. Elongated images, and
    deblending of overlapping segments, are not modelled.
  - Painted images can also overlap real neighbours, which the catalogue would deblend differently.
- **Incompleteness.** Stretched faint images lose surface-brightness-limited area, and a real pipeline may
  fragment or miss them; the S/N ≥ 5 floor approximates detection only.
  - Only images of catalogued sources are painted, so sources below the detection limit that lensing would raise
    into the catalogue are ignored. Outer images mostly have |μ| ≈ 1–1.7, so this effect is small.
  - Shapes come from Gaussian second moments with an isotropic PSF, and no pixel-level rendering was done.
- **Redshifts.** The lens is placed at the cluster redshift only.
  - Rows without a photo-z take z_s = 2; their θ_E is uncertain by the spread of the real redshift distribution.
  - Nearest-match photo-z can be a blend (field docs).
- **Isolated lens.** We assume no coupling to the cluster's shear or magnification (D-047). The observed shapes
  already contain the cluster's shear, and the W1 Jacobian is applied on top of it. Rows beyond 4 θ_E stay
  unlensed, which loses a little weak radial stretching.
- **Area.** The border excess, measured above at 3–16 %, makes the main limits slightly too strong. The
  border-corrected column bounds the effect.
- **Statistics.** Each batch of 10 trials has its own independent null, so the trials are not conditioned on one
  draw set. Efficiency counts of 0/200 give a 95 % upper bound of 1.5 % per field and mass.

### Cases inspected

We plotted catalogue ellipses for two MACS0416 θ_E = 3″ injections (scratch plot, not committed). These came from
the first version of the harness, which used one θ_E for every source and did not merge pairs. Its painting was
otherwise the same.
- **Recovered** (8 lines, p 0.030): one source near the caustic gave a long radial pair, and three more pairs
  pointed at the centre.
- **Missed** (91 sources lensed, 10 injected arcs, 5 lines, p 1.0): most images were small and weakly stretched,
  so the sources' own orientations dominated.
- Both plots look as the simulator predicts: an empty umbra apart from unlensed foreground and member rows, inner
  images bunched near the centre, and outer images stretched radially.

## W1 negative-mass lenses (shear screen)

D-053 (design D-050); `exotic_screens.py shear`, `scripts/inject_shear.py`, tests in `tests/test_exotic_screens.py`
and `tests/test_inject_shear.py`. Same lenses and painting as the radial section (D-049, `paint_lens`), measured
with a catalogue aperture-mass map instead of converging arcs.

### Method (every threshold an ASSUMPTION)
- **Shapes** (`observed` moments → `derived`): F200W pipeline `_cat.ecsv` second moments, isotropic Gaussian PSF
  removed (σ = 1st percentile of the S/N > 50 minor axes), ε = (1 − q)/(1 + q) e^{2i PA}. Rows: lensable (D-049
  rule), not a diffraction-spike segment (`spike_segments`, Gaia stars where the field uses them; D-043), S/N ≥ 10,
  resolved, cluster model |g| < 0.5 and κ < 1 at the row's photo-z.
- **Cluster shear.** The measured ε follows the model's reduced shear with a responsivity R ≈ 0.45, not 1
  (isophotal, unweighted moments). R is fitted per field, Re(ε e^{-i arg g}) against |g| through the origin:
  Abell 2744 0.451 ± 0.018, MACS0416 0.483 ± 0.034, MACS1149 0.476 ± 0.039, Abell 370 0.414 ± 0.023. The screen
  removes R g; the residual along g is then −0.009 ± 0.007 (MACS0416) and +0.006 ± 0.004 (Abell 2744). Removing the
  full g would leave −(1 − R) g, a radial pattern with the W1 sign.
- **Statistic.** S = −Σ Q e_t / √(Σ Q² |e|² / 2) (Schneider 1996) on a 1″ grid, aperture R_ap = 10″, Schirmer et al.
  (2007) Q_TANH filter (eqs. 15–16, x_c = 0.15), sources at 1″ ≤ θ ≤ 10″, ≥ 5 per aperture. S > 0 is radial
  alignment (W1). S_× (cross component) is the B-mode check.
- **Null.** 200 random rotations of every shape (positions and |e| kept); field maximum per draw.
- **Spike veto.** Without it (first run), E and B maps were both wider than the rotation null (std 1.01–1.15 vs
  1.00) with heavy tails (|S| > 3 in 0.06–0.6 % of centres vs 0.02–0.05 %): Abell 370 E max 4.51, B max 4.68 (both
  p < 0.005), Abell 2744 E 3.92 (p 0.01). Spike segments point radially at their star, the W1 sign. Dropping
  77–178 segments per field brought Abell 370 to E 3.48, B 3.46 and Abell 2744 to E 3.39.
- **E/B rule.** B-mode extremes still beat the rotation null in three fields (p 0.005–0.04), so residual shape
  systematics set the floor. A peak counts only if p_random < 0.05 **and** S exceeds the field's largest |S_×|.
  Adopted after the first run's results were seen; it is conservative.
- **Injections.** Painted images keep R of their lens-induced shape change (the catalogue's moments respond to shear
  by R, the painted moments by 1). Recovery: the largest S within 2″ of the lens passes the rule above. 200 lenses
  per field and mass, a fresh rotation null per 10 trials. Spike-segment rows, and rows where the cluster
  shear cannot be removed (κ ≥ 1 or |g| ≥ 1), are never painted (conservative).
- **Choice of filter** (`simulated` only, 50 injections in MACS0416 and Abell 2744, before the E/B rule and the R
  correction): Schirmer 10″ recovered 11/50 and 18/50 at 2 × 10¹² M☉; top-hat 6″ 7 and 17; point-mass 1/x² 10″ 3 and
  9; point-mass 4″ 2 and 8. Schirmer 10″ also led at 8 × 10¹² M☉.

### Real fields (`derived`)

| Field | spike segments dropped | sources used | per arcsec² | S_max (p_rot) | max S_× / max −S_× (p_rot) | E > B? |
|---|---|---|---|---|---|---|
| Abell 2744 | 114 | 4,951 | 0.12 | 3.39 (0.58) | 4.21 (0.005) / 3.52 (0.32) | no |
| MACS0416 | 102 | 1,777 | 0.08 | 3.79 (0.050) | 3.87 (0.04) / 4.04 (0.02) | no |
| MACS1149 | 77 | 1,399 | 0.06 | 3.40 (0.19) | 3.76 (0.015) / 3.45 (0.16) | no |
| Abell 370 | 178 | 2,191 | 0.10 | 3.48 (0.20) | 3.37 (0.31) / 3.46 (0.21) | yes, by 0.02 |

No field has an E-mode peak with p_rot < 0.05 that also beats its B-mode extremes: **null**. Abell 370's maximum
(39.99604, −1.56754, near the screened edge) passes E > B but not the rotation null.

### Recovery efficiency (`derived` from `simulated`; recovered / 200)

| \|M\| (M☉) | θ_E(z_s = 2) | Abell 2744 | MACS0416 | MACS1149 | Abell 370 | radial, all 8 fields (D-049) |
|---|---|---|---|---|---|---|
| 2 × 10¹⁰ | 0.27–0.36″ | 0 | 0 | 0 | 0 | 0 / 1600 |
| 2 × 10¹¹ | 0.85–1.14″ | 0 | 0 | 0 | 1 | 0 / 1600 |
| 2 × 10¹² | 2.7–3.6″ | 19 | 15 | 5 | 10 | 8 / 1600 |
| 8 × 10¹² | 5.4–7.2″ | 86 | 59 | 37 | 61 | 84 / 1600 |
| 2 × 10¹³ | 8.5–11.4″ | 156 | 121 | 81 | 89 | 156 / 1600 |

MACS1149 is lowest (fewest sources per arcsec²). Before the spike veto, with the rotation null alone, the four
fields recovered 66, 309 and 486 of 800 at the top three masses (biased high by the spike systematics).

### Upper limits (95 %, zero detections; `derived`; four photo-z fields, 30.2 arcmin²)

| \|M\| (M☉) | Σ ε_f A_f (deg²) | Σ(W1) < (deg⁻²) | border-corrected | radial headline (D-049) |
|---|---|---|---|---|
| 2 × 10¹⁰ | 0 | no limit | — | no limit |
| 2 × 10¹¹ | 8.5 × 10⁻⁶ | 3.5 × 10⁵ (one recovery) | 4.0 × 10⁵ | no limit |
| 2 × 10¹² | 5.6 × 10⁻⁴ | 5.3 × 10³ | 5.8 × 10³ | 6.1 × 10⁴ |
| 8 × 10¹² | 2.7 × 10⁻³ | 1.1 × 10³ | 1.2 × 10³ | 7.0 × 10³ |
| 2 × 10¹³ | 5.0 × 10⁻³ | 6.0 × 10² | 6.5 × 10² | 4.0 × 10³ |

**Reading.** The shear screen gives limits 6–11× stronger than `radial` from 20 % less area, but is equally blind
at ≤ 2 × 10¹¹ M☉ (θ_E ≲ 1″): there the umbra and images cover too few sources per aperture (S/N ≈ 1–2, as D-050
expected). The best limit is about 5× weaker than Takahashi & Asada's ~120 deg⁻² (radial section). Below 10¹² M☉
the orphan-pair limits (D-051) remain the only ones.

### Caveats
- Isophotal catalogue moments, not a weak-lensing shape pipeline: no PSF anisotropy model, no noise-bias or
  selection calibration beyond the fitted R. The remaining B-mode tails are the price.
- SMACS 0723 and El Gordo are missing: their DJA photo-z tarballs returned 404 on 2026-10-08.
- The E/B rule uses one B-mode map per field; the floor is noisy (3.46–4.21 across fields).
- Injection rows keep the real field's rotation null (only near-lens rows change).

## W2 / dark-deflector pairs (orphan-pair screen)

D-051. The script is `scripts/inject_pairs.py`, with tests in `tests/test_inject_pairs.py`. The screen is the
orphan-pair search of [orphan_pairs.md](orphan_pairs.md), run on six deep fields:
- the five CANUCS NIRCam flanking fields;
- DJA GOODS-North.

The search was null: no pair survived vetting, and the orphan counts are within about 2σ of chance. The injections
turn that into upper limits on the surface density of dark deflectors that make an image pair.

Provenance labels:
- the injected images are **simulated**;
- the efficiencies and limits are **derived** from them;
- masses and throat radii are a **model_prediction** of the lens formulae.

A limit describes the screen's sensitivity to a hypothetical population. It is not evidence for or against
wormholes.

### Lens types (exotic_sim, D-047)

| Type | Model | Images | Lensed sources |
|---|---|---|---|
| `point` (sanity) | dark point mass, n = 1, ε > 0 | two, on opposite sides; separation ≥ 2 θ_E | β ≤ 3 θ_E |
| `W2` | Ellis wormhole, n = 2, ε > 0 | two, on opposite sides, no visible deflector | β ≤ 3 θ_E |
| `W1` | negative mass, n = 1, ε < 0 | two radially stretched images on the source's side; separation 0–3.5 θ_E | β ∈ [2, 4] θ_E; rows at β < 2 θ_E removed (umbra) |

The point-mass and W2 β cut is an ASSUMPTION and a conservative one: beyond 3 θ_E the fainter point-mass image
carries < 1.3 % of the flux.

θ_E ∈ {0.15, 0.3, 0.7, 1.5}″. A point or W2 lens with θ_E = 1.5″ makes pairs ≥ 3″ apart, outside the search
annulus. Its zero efficiency is by construction.

### Method (every threshold an ASSUMPTION)

1. **Sources are the catalogue rows themselves.**
   - Every row with z_best > z_l = 0.4 and some valid photometry is lensable.
   - A lensed row is **removed** and replaced by its images, which take positions and signed μ from
     `exotic_sim.inject_images`.
2. **Painting.**
   - Each image gets |μ| × the row's (already noisy) SED. Fresh Gaussian noise of √max(1 − μ², 0) × the row's
     errors tops a demagnified image up to the row's own noise; a magnified image keeps |μ| × the row's noise
     (slightly conservative). The errors stay the row's (sky-limited).
   - Photo-z and catalogue μ are inherited.
   - Kron radius: √(r_psf² + (r² − r_psf²) s²), with s = max(1/|λ_t|, 1/|λ_r|). r_psf is the 5th percentile of
     the selected sources' radii (0.15–0.17″ in CANUCS, 0.34″ in DJA).
   - Images get no segment box, so the box-overlap rule never applies to them.
3. **Deblending (explicit rule).** Two images closer than d_blend become one row at the flux-weighted centroid
   with the summed flux. d_blend is the 2nd percentile of the catalogue's nearest-neighbour separations, capped at
   0.5″: 0.24–0.27″ (CANUCS) and 0.30″ (DJA). A merged lens cannot be recovered.
4. **Recovery.** The unchanged orphan-pair rules (`select_sources`, `match_table`, `classify_pairs`) run on the
   painted neighbourhood. A lensed source is recovered when both image rows pass the selection, lie 0.3–3″ apart,
   match in SED and photo-z, and are classed `orphan`.
5. **Two efficiencies** per field, type and θ_E:
   - **per source**: 400 rows that the search selects unlensed, each at a random β, uniform in area;
   - **per deflector**: 2000 lenses at uniform random positions in the searched footprint. A lens counts as
     recovered when any of its lensed rows is. Fainter rows lifted into the selection by μ are included.
6. **Limits** use the summed exposure Σ_f ε_f A_f. The footprint A_f is the 1″ grid points that have, within
   4″, a row with valid S/N bands and ≥ 8 valid bands: 108.1 arcmin² in total, 9.8–9.9 arcmin² per flanking
   field and 58.7 arcmin² for GOODS-N.
   - **No-candidate limit:** 2.996 / Σ εA. It assumes every orphan is ordinary.
   - **Background-aware limit:** s₉₅ / Σ εA, where s₉₅ = 72.3 is the CLs 95 % upper limit (Read 2002) on a Poisson
     signal, which never collapses to 0 when fewer events than the background are observed. The inputs are 355 observed orphans over a known background of 315.4 (null (e), summed). This is the
     defensible number: about 40 faint orphans cannot be told apart from lensed pairs.

Reproduce:
- `python scripts/inject_pairs.py` runs all six fields, about 3.5–4.5 min per field on one core.
- Alternatively, run fields in parallel processes with `--fields ...`, then `--combine-only`.

The outputs are `outputs/inject_pairs/limits.json`, plus per field `injection_summary.json` and
`per_source_<type>_<θ_E>.ecsv`.

### Efficiencies (derived, 6 fields)

| Type | θ_E | Per-deflector ε (mean over 108 arcmin²) | Per-source recovered | Both images selected (per source) | Main loss |
|---|---|---|---|---|---|
| point | 0.15″ | 0.02 % | 12 / 2400 (0.5 %) | 11 % | same_galaxy (pairs 0.3–0.5″) |
| point | 0.3″ | 0.54 % | 81 / 2400 (3.4 %) | 12 % | faint counter-image; same_galaxy |
| point | 0.7″ | 0.84 % | 15 / 2400 (0.6 %) | 10 % | visible_lens (a catalogued source inside the pair's circle) |
| point | 1.5″ | 0 | 0 / 2400 | 11 % | separation > 3″ |
| W2 | 0.15″ | 0.01 % | 8 / 2400 (0.3 %) | 10 % | same_galaxy |
| W2 | 0.3″ | 0.64 % | 59 / 2400 (2.5 %) | 9 % | faint counter-image; same_galaxy |
| W2 | 0.7″ | 0.83 % | 6 / 2400 (0.25 %) | 10 % | visible_lens |
| W2 | 1.5″ | 0 | 0 / 2400 | 10 % | separation > 3″ |
| W1 | 0.15″ | 0 | 0 / 2400 | 1 % | merged or faint inner image |
| W1 | 0.3″ | 0.10 % | 7 / 2400 (0.3 %) | 5 % | faint inner image; same_galaxy |
| W1 | 0.7″ | 3.6 % | 62 / 2400 (2.6 %) | 7 % | faint inner image |
| W1 | 1.5″ | 18.6 % | 85 / 2400 (3.5 %) | 7 % | faint inner image; 18–23 sources lensed per deflector |

- **Per field.** The per-deflector ε is similar in every field. For W1 at 1.5″ it ranges from 16.7 to 21.1 %.
- **SED match.** Every injected pair that is selected and 0.3–3″ apart passes the SED and photo-z test (expected:
  χ² p ≥ 0.01 passes 99 % of identical SEDs, and the photo-z intervals are inherited). The
  losses come from geometry and depth, not from the SED test.
- **Magnitude** (per source, all fields; m = AB of the mean F277W/F356W/F444W flux). Recovery is a
  bright-source effect:

  | Type, θ_E | m < 24 | 24–25 | 25–26 | 26–27 | 27–28 |
  |---|---|---|---|---|---|
  | point, 0.3″ | 12.5 % (31/248) | 4.1 % | 4.1 % | 1.8 % | 1.3 % |
  | W1, 0.7″ | 15 % (36/241) | 6.5 % | 1.8 % | 0.3 % | 0.1 % |
  | W1, 1.5″ | 11.5 % (30/260) | 7.3 % | 5.2 % | 1.5 % | 0.7 % |

### Limits (95 %, deg⁻²; z_l = 0.4, z_s = 2, Planck18)

| θ_E | \|M\| (n = 1; model_prediction) | Ellis throat a (pc; model_prediction) | point: no-candidate / background-aware | W2: no-candidate / background-aware | W1: no-candidate / background-aware |
|---|---|---|---|---|---|
| 0.15″ | 4.5 × 10⁹ M☉ | 0.96 | 5.4 × 10⁵ / 1.3 × 10⁷ | 1.1 × 10⁶ / 2.6 × 10⁷ | none (ε = 0) |
| 0.3″ | 1.8 × 10¹⁰ M☉ | 2.7 | 1.9 × 10⁴ / 4.5 × 10⁵ | 1.6 × 10⁴ / 3.8 × 10⁵ | 9.6 × 10⁴ / 2.3 × 10⁶ |
| 0.7″ | 9.8 × 10¹⁰ M☉ | 9.7 | 1.2 × 10⁴ / 2.9 × 10⁵ | 1.2 × 10⁴ / 2.9 × 10⁵ | 2.8 × 10³ / 6.7 × 10⁴ |
| 1.5″ | 4.5 × 10¹¹ M☉ | 30 | none (separation > 3″) | none (separation > 3″) | 5.4 × 10² / 1.3 × 10⁴ |

Reading:
- **What the screen can do.** It is sensitive at θ_E ≈ 0.3–1.5″. That is the regime where the W1 `radial`
  screen is blind (D-049, branch `claude/radial-injection`: no limit at 0.3″ or 1″, < 3.6 × 10⁴ deg⁻² at 3″).
  - The best W1 limits here are < 5.4 × 10² (no candidate) and < 1.3 × 10⁴ deg⁻² (background-aware) at
    θ_E = 1.5″, |M| ≈ 4.5 × 10¹¹ M☉.
  - In the no-candidate case that is two orders of magnitude below the radial screen's 3″ limit.
- **Comparison with Takahashi & Asada (2013).** Their Ellis volume limit, n < 10⁻⁴ h³ Mpc⁻³ for a = 10–10⁴ pc,
  corresponds to about 120 deg⁻² if spread over 0 < z < 1 (Planck18 comoving volume; computed here).
  - The W2 limits at a ≈ 10 pc are 100× (no candidate) to 2400× (background-aware) weaker.
  - At a ≲ 3 pc (θ_E ≤ 0.3″) the screen probes throat radii below their range, but only at > 10⁴ deg⁻².
- **The point-mass column** is the same calculation for ordinary compact dark deflectors (dark subhaloes,
  compact dark objects). It shows that the screen's reach is set by geometry and depth, not by the exotic
  physics.

### Caveats

- **Catalogue-level injection.** No pixels are painted. Blending with neighbours other than the twin image,
  deblending failures and photometric-pipeline effects are not modelled, apart from the explicit merge rule.
- **Photo-z.** The images inherit the source's photo-z intervals, so the z-overlap test is passed by
  construction. This is optimistic for faint counter-images, whose real photo-z could scatter.
- **Noise.** Magnified images carry |μ| × the row's noise, demagnified ones are topped up to it. The PSF-convolved size model is approximate. The DJA
  same_galaxy radius is calibrated on CANUCS (D-051).
- **Single deflector redshift** z_l = 0.4 and single source redshift z_s = 2 for the masses. A real population
  spread in z_l would change the mass scale, not the angular efficiency.
- **No macro-model.** CANUCS flanking-field μ = 1.0–1.4 is ignored for the injected lens.
- **Background-aware limit.** It treats null (e) as exact. Null (e) itself under-predicts the flanking-field
  orphans by ~15 % (orphan_pairs.md), so the true background may be higher and the limit is then conservative.

## W3 inverted microlensing / dimming (multi-epoch)

D-052; `scripts/dimming_screen.py` (`fetch`, `screen`, `inject`, `forced`, `combine`), fields and epochs in
`configs/dimming_screen.yaml`, tests in `tests/test_dimming_screen.py`. The signature is D-047's W3: an n = 1 lens
with ε < 0 crossing a compact source. The lensed flux is 0 inside the umbra (β < 2, for 2 t_E √(4 − u₀²)) between
two caustic spikes (×7.0 for ρ = 0.01, ×2.35 for ρ = 0.1; `exotic_sim` as merged in PR #67). The umbra is the
robust part of the signal; the spike heights depend on ρ.

### Data (observed)

| Field | Epochs | Baseline | Master sources | Monitored compact sources | Overlap area |
|---|---|---|---|---|---|
| NEXUS-Center (5105) | 8 (F200W + F444W): 2024-09-12 … 2026-03-28 | 1.54 yr | 81,136 | 151 | 79.5 arcmin² |
| MACS0416 (1176 o211–o213, 1208 o004; 6882 o054) | 4 F200W + F444W (2022-10-07 … 2023-02-10), plus 6882 o054 in F444W only (2026-01-10) | 0.35 yr (F200W), 3.26 yr (F444W) | 7,708 | 44 | 10.8 arcmin² |
| Abell 2744 (2561 o001/o002/o006) | 3 (F200W + F444W): 2022-11-02 … 2024-07-31 | 1.74 yr | 17,669 | 63 | 14.0 arcmin² |

- Catalogues are level-3 `_cat.ecsv` files, 31 in total (`data/manifests/dimming_*.ecsv`, 339 MB). Every file
  header says jwst 3.0.0 / photutils 3.0.0.
- MAST's product listing still gives `prvversion` 2.0.1, with smaller sizes, for the six PEARLS MACS0416 files
  (obs 211–213). The files were reprocessed and the listing is stale; `acquire` warned about it.
- Caveat: the MACS0416 epochs come from two programmes with different depths and dither patterns (PEARLS 2.9 ks,
  CANUCS 6.4 ks in F200W). Cross-epoch processing differences show up as catalogue flags, which the forced stage
  must remove.
- The master list and the W3 light-curve tests use the detection band (F200W) epochs only. VENUS o054 (F444W only)
  enters MACS0416's F444W light curves (vanish / dim tests in that band). It does not add master sources, and it
  does not lengthen the F200W baseline.
- No `_i2d.fits` was downloaded; forced photometry reads S3 byte ranges. El Gordo has one shared band and was not
  used. JADES is left for later.

### Method (every threshold an ASSUMPTION, `Params` in the script)

1. **Light curves (derived).**
   - Each catalogue is tied to the deepest F200W catalogue by the D-027 global frame shift.
   - Master list: F200W rows with aper50 S/N ≥ 5 from any F200W epoch, matched mutually within 0.3″.
   - Zero point per band and epoch: the median flux ratio of S/N ≥ 30 pairs.
   - Errors are scaled by the robust scatter of epoch-to-epoch differences (1.52–1.73), plus a 3 % floor.
   - A catalogue non-detection counts as flux 0 only where the source would have been ≥ 10σ.
2. **Flags.**
   - `vanish`: ≥ 10σ → < 3σ with a ≥ 5σ drop, in every band testable in that epoch pair.
   - `dim_achromatic`: ≥ 20 % and ≥ 5σ in two bands, with drops equal within 3σ.
   - `rise_dip_rise`: ≥ 3 epochs.
3. **Catalogue-level ordinary tests, cheapest first.**
   - `near_star`: within the D-027 radius of another Gaia DR3 source. The source's own Gaia match (< 0.3″) is
     excluded before the mask is applied, unless it is brighter than G = 18.
   - `bright_neighbour`: a master source ≥ 100× brighter within 1.5″ (new; PSF wings and spikes of saturated
     stars).
   - Edge proxy, blend (< 0.5″), sharper than the PSF (CI_70_30 < 1.9), and `single_epoch` (detected in one F200W
     epoch only; persistence suspect, D-039).
4. **Forced photometry (`forced`).**
   - A 0.15″ aperture at the recentred position in every epoch (`transient_forced.measure`). SCI (MJy/sr) is
     converted to µJy with each epoch's pixel solid angle.
   - Zero points come from controls selected on the reference epoch only (S/N ≥ 10 there), as 1 / median(f_k /
     f_ref), so faded epochs are not dropped. An epoch with < 10 such controls is recorded as uncalibrated.
   - The forced noise scale needs ≥ 50 controls in some epoch pair, else the band is marked uncalibrated.
   - Same flag logic again; a flag is confirmed when the same flag type reappears. Then cutout tests (edge / no
     data, WHT < 0.5, spike statistic ≥ 3), then SIMBAD/NED.
5. **Compact sources monitored** (the population of the limit): 1.9 ≤ CI_70_30 ≤ 2.7 (the F200W stellar locus is
   2.0–2.5), no catalogue veto, and S/N ≥ 10 in ≥ 2 F200W epochs.

### Results (derived)

Counts of the last complete forced run (2026-10-08 06:39–07:03 UTC). Re-runs after the review fixes failed on S3:
s3fs returned "bucket does not exist" through the proxy for multi-target cutout jobs, while single reads worked.
Its forced fluxes were re-calibrated offline with the fixed calibration code; the catalogue screen is the current
code.

| Field | Catalogue flags (vanish / dim / rdr) | After catalogue tests (current code) | Forced-measured | Forced-confirmed | After cutout tests | After the bright-neighbour veto and visual check |
|---|---|---|---|---|---|---|
| NEXUS | 3,793 (3,197 / 264 / 1,276) | 1,051 (12 compact) | 100 (13 compact + 87 random) | 0 | 0 | 0 |
| MACS0416 | 854 (714 / 67 / 310) | 358 (10 compact) | 375 | 6 | 2 | 0 |
| Abell 2744 | 530 (525 / 49 / 6) | 32 (1 compact) | 32 | 0 | 0 | 0 |

- The bright-neighbour veto removes 75 / 828 / 21 catalogue flags (MACS0416 / NEXUS / Abell 2744). It vetoes all
  six forced-confirmed MACS0416 flags automatically. Five of them (`d03232`, `d03844`, `d03845`, `d06533`,
  `d07517`) lie 0.5–1″ from saturated stars; the sixth is the saturated star `d03349` itself (spike statistic ≥ 3).
- Each i2d grid has a different orientation, so spikes and wings cross a fixed aperture differently in each epoch
  (cutouts inspected).
- Before this veto existed, the old Gaia mask failed because it took the source's own faint Gaia match as the
  nearest star. **0 surviving events.**
- Controls flagged by the forced stage: 7/300 (MACS0416), 1/150 (NEXUS), 0/14 (Abell 2744).
- Calibration:
  - MACS0416 is calibrated: noise 2.15 / 1.99 (F200W / F444W). Its forced/catalogue fractional scatter on the same
    300 controls is 1.00 (F200W) and 1.98 (F444W).
  - NEXUS (150 controls) and Abell 2744 (14) are **uncalibrated**: too few controls in any epoch pair.
- Wall time: `screen` 76 s (NEXUS), 5 s (MACS0416), 9 s (Abell 2744); `inject` 140 / 20 / 34 s; `forced` about
  32 min (MACS0416, 675 positions × 9 images), 20 min (NEXUS), 3–9 min (Abell 2744).

### Injection-recovery (simulated)

- Monitored compact sources that are not flagged without an injection (34 / 139 / 62) get 20 copies per model.
- W3 (n = 1, ε < 0) uses `exotic_sim.inject_light_curve` on the real epoch times. Each copy draws u₀ ~ U[0, 2) and
  t₀ ~ U[t_first − 2t_E, t_last + 2t_E].
- Noise model (ASSUMPTION): sky-limited for F ≤ 1 and Poisson-like above.
  f = F f_obs + √max(0, 1 − F²) σ_n z, with quoted error σ_n max(F, 1) plus the floor. The realised scatter
  equals the quoted error, and a vanished epoch keeps sky noise only.
- Recovered = a flag, no static veto, and, per injected copy, detections in ≥ 2 F200W epochs (the `single_epoch`
  veto re-applied). In MACS0416 the event must also stay flagged with errors inflated by the forced/catalogue
  fractional-scatter ratio (1.00 / 1.98).

Efficiency over all magnitudes (ρ = 0.01 / 0.1):

| Model | MACS0416 (calibrated) | NEXUS (uncalibrated) | Abell 2744 (uncalibrated) |
|---|---|---|---|
| W3 t_E = 0.01 yr | 0.02 / 0.01 | 0.04 / 0.03 | 0.00 / 0.01 |
| W3 t_E = 0.03 yr | 0.04 / 0.04 | 0.10 / 0.10 | 0.02 / 0.01 |
| W3 t_E = 0.1 yr | 0.06 / 0.08 | 0.23 / 0.23 | 0.03 / 0.04 |
| W3 t_E = 0.3 yr | 0.16 / 0.12 | 0.23 / 0.24 | 0.04 / 0.06 |
| W3 t_E = 1 yr | 0.20 / 0.22 | 0.12 / 0.16 | 0.04 / 0.06 |
| W3 t_E = 3 yr | 0.13 / 0.12 | 0.07 / 0.10 | 0.04 / 0.05 |
| dimming 20 / 50 / 100 % | 0.09 / 0.59 / 0.35 | 0.13 / 0.70 / 0.41 | 0.02 / 0.41 / 0.02 |

- Dimming by 100 % is now recovered only 2–41 % of the time. Most monitored sources are detected in just 2 F200W
  epochs, so making one of them vanish leaves a single detection, which the `single_epoch` (persistence) veto
  removes. In Abell 2744 almost every source is covered by only two of its three epochs.
- The screen is therefore blind to full vanishes of two-epoch sources until the `single_epoch` veto is replaced by
  the D-039 persistence test for vanish flags (follow-up).

### Upper limits (derived; 0 surviving events, Poisson 95 % = 3.0)

Exposure E = Σ_bins N_bin ε_bin W, in source years, with W = T + 4t_E the t₀ window. N_bin counts the injected
(baseline-unflagged) monitored sources.
- **Rate per compact source per year:** 3 / E.
- **τ:** the fraction of compact sources inside an umbra at a single epoch, equal to that rate × π t_E (the mean
  umbra duration for u₀ ~ U[0, 2) is π t_E).
- **Per deg² per year:** 3 / Σ A ε̄ W, with ε̄ = Σ N ε / N_monitored. This applies at these fields' compact-source
  densities.
- **Per deg² per epoch:** the per-deg²-per-year limit × π t_E, i.e. sources inside an umbra per deg² at one
  epoch.
- |M| = (t_E / 12 yr)² M☉ is a **model_prediction** for z_l = 0.4, z_s = 2 and v⊥ = 1000 km/s (D-047; re-derived
  as 12.09 yr). It does not apply to Galactic stars.

Headline (calibrated field only: MACS0416, 34 injected compact sources; ρ = 0.1):

| t_E (yr) | \|M\| (M☉) | exposure (source yr) | per source per yr | τ | per deg² per yr | per deg² per epoch |
|---|---|---|---|---|---|---|
| 0.01 | 6.9e-07 | 1.2 | 2.6 | 0.082 | 2.9e4 | 927 |
| 0.03 | 6.3e-06 | 4.7 | 0.63 | 0.060 | 7.2e3 | 679 |
| 0.1 | 6.9e-05 | 9.9 | 0.30 | 0.095 | 3.5e3 | 1.1e3 |
| 0.3 | 6.2e-04 | 18.7 | 0.16 | 0.15 | 1.8e3 | 1.7e3 |
| 1 | 6.9e-03 | 53.7 | 0.056 | 0.18 | 635 | 2.0e3 |
| 3 | 6.2e-02 | 62.6 | 0.048 | 0.45 | 545 | 5.1e3 |

All three fields, **including the two uncalibrated ones** (indicative only; `limits_combined_all_fields.ecsv`,
ρ = 0.1):

| t_E (yr) | exposure (source yr) | per source per yr | τ | per deg² per yr | per deg² per epoch |
|---|---|---|---|---|---|
| 0.01 | 7.9 | 0.38 | 0.012 | 2.7e3 | 84 |
| 0.1 | 78.1 | 0.038 | 0.012 | 267 | 84 |
| 0.3 | 120.1 | 0.025 | 0.024 | 179 | 169 |
| 1 | 194.0 | 0.015 | 0.049 | 119 | 375 |
| 3 | 290.6 | 0.010 | 0.097 | 79 | 748 |

**Caveats.**
- The isolated-lens light curve has no cluster macro-magnification or shear (D-047).
- The lensed fraction is 1.
- The injection acts on catalogue fluxes; catalogue re-detection and deblending are not simulated.
- The forced stage is emulated only through MACS0416's measured fractional-scatter inflation.
- The headline rests on 34 compact sources in one field, so it is weak.
- NEXUS extended flags are only audited (100 of 1,051 measured).
- No SN/TNS check was run (no survivors).
