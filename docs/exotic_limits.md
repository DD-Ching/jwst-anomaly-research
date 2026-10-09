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
- **Injections.** A painted image's ε is its source's corrected ε plus R times the change of the *measured* moments
  (image minus source; the catalogue's moments respond to shear by R, the painted moments by 1, and the cluster
  shear cancels in the difference); only resolved sources (finite deconvolved ε) are painted; |ε| is capped at 0.99. Recovery: the largest S within 2″ of the lens passes the rule above. 200 lenses
  per field and mass, a fresh rotation null per 10 trials. Spike-segment rows, and rows where the cluster
  shear cannot be removed (κ ≥ 1, |g| ≥ 1 or |R g| ≥ 1), are never painted (conservative). Painted images that
  land on a spike axis and point at the star are vetoed as on the real field (`ShearInjector.spike_veto`; the
  spike axes are estimated with the field's own vetoed segments). This cost 4 of 409 recoveries at 2 × 10¹³ M☉.
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
| 2 × 10¹² | 2.7–3.6″ | 14 | 7 | 3 | 9 | 8 / 1600 |
| 8 × 10¹² | 5.4–7.2″ | 72 | 50 | 30 | 45 | 84 / 1600 |
| 2 × 10¹³ | 8.5–11.4″ | 146 | 112 | 74 | 73 | 156 / 1600 |

MACS1149 is lowest (fewest sources per arcsec²). Before the spike veto, with the rotation null alone, the four
fields recovered 66, 309 and 486 of 800 at the top three masses (biased high by the spike systematics).

### Upper limits (95 %, zero detections; `derived`; four photo-z fields, 30.2 arcmin²)

| \|M\| (M☉) | Σ ε_f A_f (deg²) | Σ(W1) < (deg⁻²) | border-corrected | radial headline (D-049) |
|---|---|---|---|---|
| 2 × 10¹⁰ | 0 | no limit | — | no limit |
| 2 × 10¹¹ | 8.5 × 10⁻⁶ | 3.5 × 10⁵ (one recovery) | 4.0 × 10⁵ | no limit |
| 2 × 10¹² | 3.9 × 10⁻⁴ | 7.7 × 10³ | 8.4 × 10³ | 6.1 × 10⁴ |
| 8 × 10¹² | 2.2 × 10⁻³ | 1.3 × 10³ | 1.5 × 10³ | 7.0 × 10³ |
| 2 × 10¹³ | 4.6 × 10⁻³ | 6.5 × 10² | 7.1 × 10² | 4.0 × 10³ |

**Reading.** The shear screen gives limits 5–8× stronger than `radial` from 20 % less area, but is equally blind
at ≤ 2 × 10¹¹ M☉ (θ_E ≲ 1″): there the umbra and images cover too few sources per aperture (S/N ≈ 1–2, as D-050
expected). The best limit is about 5× weaker than Takahashi & Asada's ~120 deg⁻² (radial section). Below 10¹² M☉
the orphan-pair limits (D-051) remain the only ones.

### Caveats
- Isophotal catalogue moments, not a weak-lensing shape pipeline: no PSF anisotropy model, no noise-bias or
  selection calibration beyond the fitted R. The remaining B-mode tails are the price.
- SMACS 0723 and El Gordo are missing: their DJA photo-z tarballs returned 404 on 2026-10-08.
- The E/B rule uses one B-mode map per field; the floor is noisy (3.46–4.21 across fields).
- Injection rows keep the real field's rotation null (only near-lens rows change).
- Only resolved, correctable rows are lensed in the shear injections; the radial injections (D-049) also lens
  unresolved rows. The efficiency ratio between the screens mixes screen power with that (painting unresolved rows
  too gave 50 / 238 / 449 instead of 33 / 197 / 405). The linear ε sum is capped at |ε| = 0.99.

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
   - **Background-aware limit:** s₉₅ / Σ εA, where s₉₅ = 54.9 (D-055; 72.3 in D-051) is the CLs 95 % upper limit (Read 2002) on a Poisson
     signal, which never collapses to 0 when fewer events than the background are observed. The inputs are 355 observed orphans over a known background of 334.0 (null (e) as fixed in D-055, summed; D-051 used 315.4). This is the
     defensible number: the ~21 orphans above that background cannot be told apart from lensed pairs.

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
| 0.15″ | 4.5 × 10⁹ M☉ | 0.96 | 5.5 × 10⁵ / 1.0 × 10⁷ | 1.1 × 10⁶ / 2.0 × 10⁷ | none (ε = 0) |
| 0.3″ | 1.8 × 10¹⁰ M☉ | 2.7 | 1.9 × 10⁴ / 3.4 × 10⁵ | 1.6 × 10⁴ / 2.9 × 10⁵ | 9.6 × 10⁴ / 1.8 × 10⁶ |
| 0.7″ | 9.8 × 10¹⁰ M☉ | 9.7 | 1.2 × 10⁴ / 2.2 × 10⁵ | 1.2 × 10⁴ / 2.2 × 10⁵ | 2.8 × 10³ / 5.1 × 10⁴ |
| 1.5″ | 4.5 × 10¹¹ M☉ | 30 | none (separation > 3″) | none (separation > 3″) | 5.4 × 10² / 9.8 × 10³ |

Re-run 2026-10-08 with the D-055 nulls (`inject_pairs.py` per field, then `--combine-only`): the no-candidate
column reproduces D-051 (0.15″ point: 5.46 × 10⁵, printed 5.4 × 10⁵ in D-051); the background-aware column uses
s₉₅ = 54.9 (D-051: 72.3).

Reading:
- **What the screen can do.** It is sensitive at θ_E ≈ 0.3–1.5″. That is the regime where the W1 `radial`
  screen is blind (D-049, branch `claude/radial-injection`: no limit at 0.3″ or 1″, < 3.6 × 10⁴ deg⁻² at 3″).
  - The best W1 limits here are < 5.4 × 10² (no candidate) and < 9.8 × 10³ deg⁻² (background-aware; 1.3 × 10⁴ before D-055) at
    θ_E = 1.5″, |M| ≈ 4.5 × 10¹¹ M☉.
  - In the no-candidate case that is two orders of magnitude below the radial screen's 3″ limit.
- **Comparison with Takahashi & Asada (2013).** Their Ellis volume limit, n < 10⁻⁴ h³ Mpc⁻³ for a = 10–10⁴ pc,
  corresponds to about 120 deg⁻² if spread over 0 < z < 1 (Planck18 comoving volume; computed here).
  - The W2 limits at a ≈ 10 pc are 100× (no candidate) to 1800× (background-aware; D-055) weaker.
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
- **Background-aware limit.** It treats null (e) as exact. After the D-055 fix, null (e) and the companion-aware
  null (f) both match the observed orphans (docs/orphan_pairs.md); (f) predicts ~14 more, so the (e)-based limit
  is the more conservative.

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

## W1/W2 in published lens catalogues

D-056; `src/jwst_anomaly/lenscats.py`, `scripts/w12_lenscats.py`, tests in `tests/test_lenscats.py` and
`tests/test_w12_lenscats.py`. Every number below is in the run's `summary.json`.

Question: do any **published** strong-lens candidates or confirmed lenses have lensed images but **no visible
deflector** (W2, and W1's empty centre), down to the depth of public wide imaging? Catalogue values are **observed**.
Matches, required magnitudes and limits are **derived**. Required magnitudes rest on a **model_prediction** (SIS) and
on **assumptions**. A system without a deflector would be an anomaly to vet, never evidence of exotic physics.

### Data and footprint (observed; SOURCES.md "Published lens catalogues and deep-imaging checks (D-056)")

- **Catalogues:** lenscat 1.1.3 (32,838 entries), Euclid Q1 Discovery Engine (2,584) and SuGOHI (3,961), merged
  within 3″ into 35,862 systems.
- **Galaxy scale:** no group or cluster entry, including the lenscat rows that cite 14 cluster surveys, and θ_E ≤ 3″
  when given. 20,986 systems.
- **Footprint and depth** come from the Legacy Surveys DR10 (DECam) brick summary `ls_dr10.bricks_s` (332,581 bricks
  with survey_primary and z data; 23.5 MB; sha256 in `summary.json`). They do not depend on detected sources.
  - Covered means the brick has nexp_r ≥ 1 and nexp_z ≥ 1: 17,555 systems.
  - Depth is the brick's 5σ galaxy depth, median z = 23.43 (5–95 %: 22.84–24.27; ASSUMPTION: the brick median holds
    at the position).
  - `ls_dr10.tractor` is DECam-only, so the DR9 north bricks are not used.
- **Sources:** Tractor sources in ±5″ boxes, from the Data Lab TAP.
- **W1 geometry is not testable:** no catalogue gives image positions.

### Which systems can test a dark deflector at all

A finder that searches around galaxies cannot select a lens without one. A single-dish sub-mm centroid is uncertain by
more than θ_E. The test is therefore defined per selection class (ASSUMPTIONs in `Params`):

| Selection (galaxy scale / covered) | Test a dark lens could fail | Covered systems by outcome |
|---|---|---|
| galaxy finders or unknown (20,424 / 17,102) | none: insensitive by selection | 17,102 insensitive |
| sub-mm (SPT, Herschel, ACT, Planck; 118 / 110) | none: centroid error > θ_E | 110 insensitive |
| lensed-quasar searches (407 / 325) | quasar pair test (below) | blended 242, too close 65, deflector 7, faint galaxy 1, none 10 |
| radio interferometric (CLASS, JVAS, MG, mJIVE; 37 / 18) | extended source within 1.5″ bright enough for the lens; a fainter one gives "faint galaxy"; else the pair test on optical point sources; else "none" | blended 1, too close 1, deflector 6, faint galaxy 1, none 9 |

**Quasar pair test.**
- The images are the two brightest PSF-typed sources within 3″ of the position. The pair must be ≥ 2″ apart, or the
  system is "too close".
- PSF-typed sources within 3″ of the pair midpoint with g − z within 0.5 mag of the pair's colour (or with no
  colour) are taken as further images (the 3rd and 4th images of a quad), never as the deflector. A colourless pair
  therefore takes a compact lens typed PSF as an image (errs towards "none", weaker limits).
- Deflector candidates lie inside the smallest circle about the images' centroid that holds every image (the circle
  with the pair as diameter for a double; a fold or cusp pair of a quad leaves the lens outside the pair circle), and
  more than 0.5″ from every image. The ±5″ box limits the circle to pairs within 3″ of the position.
- **The LS pair must be the catalogued pair** (D-064 check, `lenscats.pair_match`, D-056 amendment 2026-10-08). The
  catalogues give θ_E, not image positions, so the catalogued separation is 2θ_E (SIS `model_prediction`). A system
  whose status came from the LS pair (`used_pair`) and whose pair separation differs from 2θ_E by > 0.5″ is
  undecided. Without a catalogued θ_E the pair is unchecked and the system stays decided.

**Status rules.**
- A candidate with m_z ≤ typical required magnitude + 2 rms is a "deflector". Where that magnitude is undefined, the
  conservative one is used, and otherwise −∞.
- **"faint galaxy"** (only fainter candidates) is **undecided** and is removed from N. A galaxy sits where the lens
  should be, and whether it is luminous enough depends on the Faber–Jackson scatter beyond the margin. Such a system
  can neither show nor exclude a dark deflector.
- Position and imaging problems are evaluated on every covered system and remove it from N whatever the test said.
  In the 343 covered quasar and radio systems: maskbits 9, rounded positions 5, name/RA–Dec mismatch 1, LS pair
  not the catalogued pair 1.
- A position counts as rounded if either axis is a whole 0.01° or 0.1°, or if both axes show a finer step.
- **One entry per lens** (D-056 amendment 2026-10-08, `dedup_same_lens`). Three radio lenses are listed twice in
  lenscat, about 11″ apart (beyond the 3″ merge): MG0414+0534 (L00042, L04558), B2114+022 (L00194, L10187) and
  B2319+052 (L00211, L04612). Decided systems with the same designation (`lenscats.designation_key`: survey prefix,
  J/B, spaces and suffixes removed, HHMM±DD) within 30″ (ASSUMPTION, `Params.same_lens_radius`) are one lens. A copy
  with a deflector is kept, since a deflector seen at one catalogued position of the lens explains it; otherwise the
  first copy is kept (the catalogues give no lens-galaxy position). Each of the three keeps its "deflector" copy. The
  dropped "none" copies sit on empty sky about 11″ from the lens (L04558 cutout inspected). The merged ids are in
  `summary.json` `vetting.same_lens_merged`.
- **Decided: 25 systems**, of which 12 have a deflector and 13 do not ("none"). Before the pair check and the
  dedup: 29, 13 and 16.
- **Pair check coverage.** 15 decided systems got their status from the LS pair; none of the 15 has a catalogued
  θ_E, so the check removed only a system that has one: 115252+004733 (θ_E 1.67″, so 3.34″ expected; LS pair 4.18″;
  "deflector" before). Its cutout shows an 18.1 mag (z) red galaxy at the centre and two faint PSF sources
  (z ≈ 23.2–23.3) about 2″ either side. The galaxy is plainly there, but the faint pair is not the catalogued one,
  so the system is undecided. Three of the 15 have separations in the SQLS tables pinned for D-064
  (`data/manifests/w12_niq_inputs.json`): J1322+1052 1.88″, J1349+1227 2.99″ and J1515+1511 2.03″, against LS
  2.00″, 3.01″ and 2.01″. All three match (a one-off check, not part of the chain). The other 12 are unchecked.

**Required lens light.**
- SIS σ from θ_E: the catalogue's, else (quasar and radio systems only) half the image separation, else 1″. Then an
  empirical Faber–Jackson calibration, fitted on 605 galaxy-selected lenses that have an extended Tractor source
  within 0.5″ of the catalogued lens, θ_E and z_l. Result (LS z): a = 20.41 at σ = 200 km/s and z_l = 0.5,
  k = 0.67, rms 0.88 mag.
- **Typical:** catalogue z_l (else 0.5) and z_s (else 2). Median required m_z = 19.73.
- **Conservative:** the faintest over z_l, with z_s = 3 when unknown, + 2 rms, and θ_E = 0.5″ when unknown
  (`detectable_floor`). Median 25.05.
- An ordinary lens counts as detectable when required m_z < depth − 0.5. That holds for 99.9 % of covered systems
  in the typical variant and 3.2 % in the conservative variant.

**Completeness is assumed, not measured.** The limits take the test to be complete for decided systems: a dark lens
there would give "none". This was checked only as logic, in the unit tests: a pair with no source between its images
reaches "none". A realistic injection would need Tractor re-run on images, which this offline step cannot do. Not
simulated:
- a faint lens blending with the images;
- Tractor typing a compact lens as PSF (counted as a further image when its colour matches the pair);
- seeing, and depth varying within a brick;
- loss to "faint galaxy": an unrelated faint source in the search circle removes a dark-lens system from N but not
  an ordinary one (a few × 10⁻³ sources arcsec⁻² over 3–10 arcsec² gives an efficiency near 0.98, not 1).

### Vetting of the 13 "none" systems (`candidates_vetted.ecsv`; all cutouts inspected)

| Ordinary explanation | Systems |
|---|---|
| SIMBAD galaxy within 3″ (catalogued lens galaxy: HE1104−1805, 2M1134−2103, J1322+1052, J1349+1227, J1515+1511, MG0751+2716) | 6 |
| lens redshift published (MG1549+305, MG2016+112, MG1131+0456) | 3 |
| literature: HSC J2212−0103, lens galaxy fitted in HSC with i = 22.40 (He et al. 2025, arXiv:2509.03858) | 1 |
| **open in the typical variant:** SuGOHI IX 090434−005328 (A), 091517+040747 (C), 104122−005618 (B), CHITAH pairs of 2.0–2.25″ | 3 |

- 10 of the 13 have a deflector in the literature that the LS test missed. LS DR10 often does not detect real lens
  galaxies next to quasar images, so a "none" is weak evidence.
- All three open pairs need a conservative lens (m_z 23.3–23.5) fainter than the local depth (22.9–23.2). An ordinary
  lens below the LS depth explains them, so none survives every ordinary test, and there is no `docs/candidates/` file.
- Follow-up: their HSC (CHITAH) lens models, and spectra of both images (binary quasar vs lens).

### Limits (95 %; f_dark < s₉₅(k) / N per selection class; completeness assumed 1; `derived`)

| Variant | quasar: k / N / f_dark < | radio: k / N / f_dark < | all |
|---|---|---|---|
| typical | 3 / 15 / 0.52 | 0 / 10 / 0.30 | 3 / 25 / 0.31 |
| conservative | 0 / 5 / 0.60 | 0 / 0 / — | 0 / 5 / 0.60 |

With the pair check and one entry per lens (D-056 amendment 2026-10-08). Before them: typical quasar 3 / 16 / 0.48,
radio 0 / 13 / 0.23 and all 3 / 29 / 0.27; conservative 0 / 6 / 0.50.

f_dark is the fraction of quasar- or radio-selected galaxy-scale lenses whose deflector is dark (fainter than an
ordinary lens of that θ_E).
- These limits are **weak and rest on an assumed complete test.** LS DR10 decides only 25 of the 343 covered quasar
  and radio systems (242 blended, 65 closer than 2″).
- There is **no limit** on galaxy-finder or sub-mm systems, or on W1: there is no image geometry, and W1 pairs are
  unlikely to pass lens finders.
- Earlier versions of this section quoted f_dark < 1.5 × 10⁻⁴, then 0.13. Both are withdrawn. The first counted
  systems that could not fail the test. The second pooled k and the efficiency, counted under-luminous galaxies as
  explained, and inferred coverage from detected sources (which drops exactly the dark configurations).

### Caveats

- Catalogue quality dominates raw "no deflector" counts:
  - cluster-survey rows typed "galaxy";
  - AGEL rows with the declination degrees dropped;
  - SPT rows up to 1.7° off;
  - 468 covered systems with rounded positions (a random 32 inspected);
  - candidates their own papers rejected (MJV16999, Spingola et al. 2019);
  - name-based merges across lists.
- Copies of one lens are merged only when both carry a designation (HHMM±DD) and lie within 30″. Copies with
  different names, or more than 30″ apart, would still count twice.
- HSC-SSP imaging (account required) and HST photometry were not used. They would decide the 307 blended or close
  lensed quasars.
- Relation to D-051: D-051 limits dark deflectors per unit area in JWST deep fields. This is a per-lens fraction in
  published lens lists; the two are not combinable without a lensing cross-section model.
- Prior art: Jackson, Helbig & Browne 1998 found lens galaxies in 12 of 12 JVAS/CLASS lenses (astro-ph/9804136).

## W3 in the OGLE-IV microlensing samples

D-057; `src/jwst_anomaly/ogle.py` (adapter), `scripts/w3_microlensing.py`
(`fit`, `vet`, `revet`, `sheet`, `inject`, `audit`, `limit`, `summary`, `manifest`), tests in
`tests/test_ogle.py` and `tests/test_w3_microlensing.py`. Run outputs live under
`$JWST_ANOMALY_DATA/derived/w3_ogle/` (not in git); data pins in SOURCES.md "OGLE-IV microlensing
samples" and `data/manifests/ogle_mroz.ecsv`.

Question: in the published, homogeneous OGLE-IV samples, does any event fit a negative-mass (n = 1,
ε < 0) or Ellis (n = 2) lens better than ordinary microlensing? Published photometry and fits are
**observed**; our fits, ΔBIC and limits are **derived**; injected events are **simulated**; the
t_E → |M| conversion is a **model_prediction**. A better exotic fit is an anomaly to vet, never
evidence of exotic physics.

### Data (observed)

| Sample | Events fitted | Fields with efficiencies | Monitored sources | ΔT |
|---|---|---|---|---|
| Mróz et al. 2019 bulge (low-cadence) | 5,790 | 112 | 6.13 × 10⁸ (I < 21, Table 7) | 2,741 d |
| Mróz et al. 2020 Galactic plane (GVS) | 460 (Table B1) | 1,982 | 1.86 × 10⁹ (database stars, Table A1; ASSUMPTION ≈ N_s) | 2,650 d |

The nine high-cadence bulge fields of Mróz et al. 2017 are outside this sample (no light curves in
`phot.tar.gz`); the 170 "possible" plane events of Table B2 failed the published selection and are
not used.

### Method (every threshold an ASSUMPTION, `Params` in the script)

1. **One trajectory, one flux solve.** Source positions come from MulensModel's trajectory
   (`Model.get_trajectory`; the annual-parallax basis is affine in π_E and is reproduced to 1e-10,
   test). F = f_s A + f_b from one weighted linear least squares per model, with the published
   bounds f_b ≥ −F_min (I = 20.5) and f_s ≥ 0.
   - Ordinary: PSPL (equal to MulensModel `point_source` to 1e-12); FSPL (uniform disk) when the
     PSPL u₀ < 0.1, cross-checked against `finite_source_uniform_WittMao94`; PSPL + annual parallax
     when t_E ≥ 20 d.
   - Exotic (`exotic_sim`): `N1neg` (n = 1, ε < 0), `E2pos` (Ellis, n = 2), `E2neg` (n = 2, ε < 0),
     each with a uniform-disk source of free ρ, so the caustic spikes are capped. The exact disk
     integral runs within 10ρ of a singular radius, a tabulated point-source magnification
     (2 × 10⁻⁴ accurate against `exotic_sim`, test) elsewhere.
   - Nelder–Mead from grids (published fit; PSPL bump; absolute t_E; both caustic spikes placed on
     pairs of light-curve maxima), best starts plus a restart. ΔBIC = BIC(exotic) − min BIC(ordinary),
     flag at ΔBIC < −10.
2. **Vetting, cheapest first** (`vet`, then `revet` for the last three, which need the survivors' fits):
   refit with FSPL and parallax for every flag; isolated 4σ outliers removed and errors rescaled to
   χ²/dof = 1; baseline variability; free blend per observing season; free blend **and** linear drift
   per season; binary source (xallarap proxy); binary lens (VBMicrolensing through MulensModel);
   VSX and Gaia DR3 variability matches (CDS XMatch, 1″); arXiv mentions of both names; then
   **feature coverage** (are there ≥ 3 epochs where the exotic and the best ordinary model differ by
   > 3σ, and does the Δχ² come from them?); a **jackknife** (drop the most influential epochs and
   refit; at most 3, and never fewer than 3 left inside the feature, so a short, well-sampled W3 event
   survives it; skipped with exactly 3 feature epochs); and **two unrelated events** (two independent PSPL bumps, the second started at the
   epoch that favours the exotic fit most outside ±2 t_E). The `revet` tests are checked on synthetic
   W3 events in the unit tests but are not run in the injection loop.
3. **Selection emulation** (`published_selection`): Mróz et al. 2019 Table 2. Not emulated: n_DIA ≥ 3
   (difference-image centroids), the s < 0.4 artifact statistic, neighbours brightened together, the
   moving search window (one window on the brightest 3-point mean here), the human inspection, and the
   unpublished bump-counting rule. **Audit on the 5,790 published bulge events, which all passed the
   real selection: 63.9 % pass the emulation** (per cut: F_b 78 %, χ²_fit,tE 87 %, χ²_out 98 %, every
   other cut ≥ 98 %), so the emulation is somewhat stricter than the original, mostly on blending.
4. **Injection-recovery** (`inject`): 400 bulge events with χ²/dof ≤ 1.5 give real cadences and noise
   (PSPL subtracted, residuals rescaled to the baseline). 600 W3 events (n = 1, ε < 0; u₀ ~ U[0, 2),
   t₀ ~ U over the efficiency window, t_E ∈ {3, 10, 30, 100, 300} d, ρ ∈ {0.01, 0.1}, 60 per cell)
   and 300 PSPL controls (u₀ ~ U[0, 1)) on the same light curves.

### Results (derived)

| Sample | Events | Flags (ΔBIC < −10) | After vetting |
|---|---|---|---|
| bulge 2019 | 5,790 (0 fit failures) | 127 (E2pos 96, E2neg 22, N1neg 9; ΔBIC −10.2 … −2788; 45 fields) | **0** |
| plane 2020 | 460 | 6 (all E2pos; −12.6 … −33.8) | **0** |

- Best ordinary model, bulge: PSPL 5,377, parallax 401, FSPL 12. ΔBIC (best exotic) quantiles
  5/25/50/75/95 % = −3.6 / 3.8 / 5.7 / 6.6 / 12.0; 554 events (9.6 %) have any ΔBIC < 0.
- Bulge vetting, cumulative: 127 → 113 (refit with all ordinary models) → 80 (robust errors,
  isolated outliers) → 73 (variable baseline) → 14 (**free blend per season**) → 9 (season drifts)
  → 7 (binary source; binary lens removed none further) → 7 (9 flags matched VSX or Gaia DR3
  variables, none of them still alive) → 1 (feature coverage) → 1 (jackknife) → **0** (two unrelated
  events).
- The last two tests matter: of the 7 flags that passed everything else, five differ from the best
  ordinary model by < 3σ anywhere, or put their caustic spike in an observing gap (0–2 epochs inside
  the feature; BLG603.25.29679 has 2, and loses the preference when its most influential epoch is
  dropped, ΔBIC −15.2 → +0.4). BLG519.21.110304 (OGLE-2011-BLG-0589; 4 epochs in the feature, ΔBIC −14.2 after
  one dropped) is an N1neg fit with f_s ≈ 0.01 that puts one caustic spike on the published 2011 event
  and the other on three points of a 1-day brightening in September 2015 (JD 2457277.5–2457278.7, ~6σ),
  with an umbra too shallow to see. Two unrelated PSPL bumps (the 2011 event plus a t_E ≈ 2.7 d bump)
  fit better by ΔBIC 24.4: an ordinary second brightening (a flare or a second lens), not W3.
  The survivors' light curves with every model curve were inspected (scratch contact sheet).
- Plane sample: every flag loses its preference once each season gets a free baseline offset (D-057).
- Wall time: bulge `fit` 10,413 s (5,790 events, 4 cores, 1.8 s per event), `vet` 6,592 s (127 flags,
  2 cores, including the binary-lens fits), plane `fit` 1,423 s, `inject` 1,974 s (900 injections,
  2 cores). ~5 CPU hours in total.

### What the published samples can contain (simulated)

| t_E (d) | \|M\| (M☉, model_prediction) | flagged and cheaply vetted (ρ = 0.01 / 0.1) | flagged and vetted (both ρ) | **passes the published selection** | PSPL control passes |
|---|---|---|---|---|---|
| 3 | 1.7 × 10⁻³ | 0.35 / 0.52 | 0.43 | **0.00** | 0.15 |
| 10 | 1.8 × 10⁻² | 0.62 / 0.38 | 0.50 | **0.00** | 0.27 |
| 30 | 0.17 | 0.88 / 0.75 | 0.82 | **0.00** | 0.35 |
| 100 | 1.8 | 0.93 / 0.78 | 0.86 | **0.00** | 0.43 |
| 300 | 17 | 0.93 / 0.78 | 0.86 | **0.00** | 0.15 |

- **0 of 600** injected W3 events pass the emulated published selection, in every t_E, ρ and u₀ bin
  (u₀ < 1, 1–1.8, 1.8–2). The fitter finds them (42–97 % flagged; 35–93 % also survive the cheap
  vetting), so the loss is the selection, not the search. The cuts they fail, in order of how often:
  one bump (70 %), χ²_fit/dof ≤ 2 of the PSPL fit (66 %), χ₃₊ ≥ 32 (60 %), three consecutive 3σ
  points (51 %), F_b > −F_min (49 %); the median injected event fails four cuts at once.
- This is the D-054 selection caveat, measured: a PSPL-shaped finder cannot select a light curve
  whose signal is an umbra between two spikes. It also means the 5,790 + 460 published events are not
  a W3-complete sample, and the null above is a statement about *PSPL-selected* events only.

### Limit (derived)

**No W3 event-rate limit follows from these samples.** With zero recovered injections, the measured
W3 efficiency is 0 and the 95 % limit is formally infinite. What can be stated:

- The 95 % Poisson upper bound on the recovery fraction is 3/60 per cell, i.e. ε_W3 / ε_PSPL < 0.12–0.33
  depending on t_E. Even at that bound the rate limit would be no stronger than 0.8–3.9 × 10⁻⁸ per
  monitored star per year (`rate95_floor_per_star_yr` in `limits_bulge2019.ecsv`), against a measured
  ordinary event rate of 5–25 × 10⁻⁶ per star per year (Mróz et al. 2019, Table 7). A real W3 limit
  therefore needs a search run on the OGLE light curves **before** the PSPL selection, not on the
  published event lists.
- Scale: for n = 1, θ_E = √(κ |M| π_rel), D_L = 4 kpc, D_S = 8 kpc, μ_rel = 5 mas/yr (ASSUMPTION),
  t_E = 73.7 d (|M|/M☉)^½ — the grid above covers 1.7 × 10⁻³ to 17 M☉ of |negative| mass.
- D-052's JWST limits (above) are per compact source per year in three deep fields; they are not
  combinable with a Galactic per-star rate without a lens-population model.

### Caveats

- The emulated selection is stricter than the published one (64 % pass rate on real events), so the
  "0 of 600" is an upper bound on W3 selectability only up to that factor; the injected events fail
  four cuts on average, so a factor ~1.6 cannot change the conclusion.
- Detection efficiencies are the published PSPL ones; they are not W3 efficiencies, which is the point
  of this section.
- The injection noise model (residual rescaling by √F with a 0.3 floor) is an ASSUMPTION; blending uses
  each host event's own f_s.
- This bulge run predates D-058, so its parallax fits are unbounded and a few reach |π_E| ~ 10³. That is
  conservative for the null (a more flexible ordinary model can only remove exotic flags), but the flag list under
  the bounded fitter may be longer; the D-059 chunk tables are the re-fit under the current `Params`.
- Binary-lens fits are a 54-start grid with a 600-evaluation local search, not a global search; they
  are a vetting test, not a characterisation.
- Only the published samples were used; OGLE EWS seasons wait for the owner's decision on their terms
  (D-054).

## W3 in the Gaia DR3 microlensing candidates

D-061; `src/jwst_anomaly/gaia_mulens.py` (adapter), `scripts/w3_gaia.py` (`fetch`, `fit`, `inject`, `summary`,
`manifest`), tests in `tests/test_w3_gaia.py`. Tracked tables: `results/w3_gaia/`; inputs:
`data/manifests/gaia_dr3_mulens.ecsv`.

- **No limit.** The emulated Sample A selection (Wyrzykowski et al. 2023, Table C.1; passes 143 / 163 real Sample A
  events) passes 2 / 240 injected W3 events and none that the fitter also flags (PSPL controls 17 / 120). The
  catalogue cannot constrain W3 at any t_E (`derived`; ASSUMPTIONs in `SelParams`).
- Fits of all 363 candidates: one flag (4053892503992268288, ΔBIC −40.3); truncated event on a variable baseline,
  not a candidate.

## W3 in MOA-II (pilot: gb22)

D-062; `src/jwst_anomaly/moa.py` (adapter), `scripts/w3_moa.py` (`prescreen`, `fit [--chunk]`, `merge-chunks`, `vet`, `sheet`,
`inject`, `limit`, `manifest`), tests in `tests/test_moa.py` and `tests/test_w3_moa.py`. Outputs under
`$JWST_ANOMALY_DATA/derived/w3_moa/` (not in git); data pins in SOURCES.md "MOA-II 9-year bulge release" and
`data/manifests/moa_ii.ecsv`. The gb22 fit table is tracked as `results/w3_moa/fits_gb22_chunk1of1.ecsv.gz`
(`merge-chunks --n 1`); `fits_gb22_chunk1of8.ecsv.gz` belongs to the superseded 1,058-pass pre-screen.

Why MOA: the OGLE samples above cannot limit W3 because their PSPL selection rejects every injected W3
event. The MOA-II 9-year release publishes a light curve for **every Cut-0 object** — 2,409,061 variable
objects found on difference images with positive *or negative* PSF profiles (Koshimoto et al. 2023, Table 2:
S/N > 2.7, N_continue,8 ≥ 3) — before any bump cut. Light curves and metadata are **observed**; scan
statistics, fits, ΔBIC and limits are **derived**; injected events are **simulated**; the t_E → |M| scale is a
**model_prediction**; every threshold is an **ASSUMPTION** (`Params` and module constants in the script). A
flag is an anomaly to vet; it is not evidence of exotic physics.

### Data (observed)

gb22 (l, b) ≈ (10.0°, −7.5°), the smallest field: 18,599 Cut-0 objects (248 have too few good epochs to
scan), ~3,100 epochs each over HJD 2453824–2456970 (8.61 yr). Difference flux in counts; the de-trended
`cor_flux` is used where it is complete (6,801 objects), else `flux`; epochs with `included = False` are
dropped. gb22 is not in Nunota et al. 2024 (no clear red clump), so it has no published N_s.

### Method

1. **Pre-screen** (`deficit_scan`, every light curve, 207 s on 4 cores): for boxes of width 1, 3, 10, 30, 100,
   300 d (≥ 3 epochs on ≥ 2 nights), the significance of the box mean below its flanks *and* below the median
   flux (the weaker of the two), with errors scaled to the point-to-point scatter and each statistic divided by
   its own robust spread over the light curve (red noise). Pass: z < −10, white-noise S < −5, no second
   deficit (z < −8) more than 4 box widths away, and the deficit epoch not shared with improbably many other
   objects — same width in the whole field, or any width on the same chip (Poisson p < 10⁻³ against deficits
   spread uniformly in time). Thresholds were set from the real distribution and the injections (below).
2. **Fit** (`fit`): PSPL, FSPL and the exotic N1neg, E2pos, E2neg models of `w3_microlensing.py` on one
   trajectory, blend flux unbounded (difference flux), extra repulsive-lens starts on the deficit, t_E ≤ 1,000 d
   and ρ ≤ 0.3. Parallax is fitted in vetting only (it can only remove flags). Flag: ΔBIC < −10.
3. **Vetting, cheapest first, stopping at the first failure** (`vet_one`): repeated deficit outside the exotic
   feature (z < −6); variable baseline (χ²/dof > 2 beyond 60 d of the feature); a neighbouring Cut-0 object
   (< 12 px) with |S| > 5 over the same window; **an eclipse or occultation** (baseline − depth × trapezoid,
   no caustic spikes); isolated outliers and rescaled errors; free baseline per season; plus drift per season;
   baseline linear in seeing, airmass and sky; refit with FSPL and parallax; binary source; binary lens
   (MulensModel + VBMicrolensing grid); the exotic feature sampled on **≥ 3 nights** with Δχ² from inside it;
   epoch jackknife; **night jackknife** (drop the 1–2 most influential nights); two unrelated events; VSX and
   Gaia DR3 variability matches (2″).
4. **Injection-recovery** (`inject`): 600 W3 events (n = 1, ε < 0; u₀ ~ U[0, 2); t₀ ~ U over the data span;
   t_E ∈ {3, 10, 30, 100, 300} d × ρ ∈ {0.01, 0.1}, 60 per cell; source MOA-Red magnitude I_s ~ U[14.2, 21.4],
   flux from the archive zero point) added as F_s(A − 1) to the real light curves of 300 quiet gb22 objects
   (|z| < 4), plus 100 PSPL controls (u₀ < 1). Each goes through a light-curve-level Cut-0 emulation (≥ 3
   epochs with |signal|/σ > 2.7, each ≤ 8 d after the previous), the pre-screen (against the real field's
   shared-epoch population), the fit and the vetting (binary lens and VSX/Gaia excepted).

### Results (derived)

| Stage | gb22 |
|---|---|
| Light curves scanned | 18,599 (18,351 with enough epochs) |
| Pass the shape cuts | 1,058 |
| … and not at a shared epoch (field or chip) | **30** (1,022 removed by the field test, 6 more by the chip test) |
| Flags (ΔBIC < −10) | 30 (N1neg 22, E2neg 7, E2pos 1; ΔBIC −95 … −3,889) |
| After vetting | **0** |

- Vetting funnel: 30 → 17 (repeated deficit) → 15 (variable baseline) → 13 (neighbour) → 5 (eclipse dip) → 3
  (isolated outliers) → 0 (feature on ≥ 3 nights). Every flag fails because the dip is ordinary: an
  eclipse-like flat-bottomed dip (gb22-R-6-0-4889, -10-3-16229, -10-0-24345: 5–8 d deep dips that a trapezoid
  fits better than any lens by ΔBIC 36–142, with the predicted caustic spikes absent where sampled), a
  periodic or repeated dimming, a noisy baseline, or a dip confined to one or two nights.
- The first vetting pass (without the eclipse model, the chip-level shared-epoch test and the ≥ 3-night rule)
  left 9 survivors. All nine were inspected on a contact sheet: four were one- or two-night drops (the
  predicted spikes fell in gaps), five were flat or V-shaped dips without spikes. The three tests were added
  because these are the ordinary explanations the chain lacked (eclipsing and dipping stars are a listed W3
  mimic; MOA takes several exposures a night, so one bad night is several bad epochs); their cost to real W3
  events is in the injection numbers below, which were run after the change.
- The pre-screen is W3-specific: **0 of 100** PSPL controls pass it (they have no deficit below the baseline),
  so no control can become a W3 survivor.
- Wall time: prescreen 203 s, fit 335 s (30 light curves), vet 22 s, inject 888 s (700 injections), all on
  4 cores; the abandoned first fits of 1,058 shape passes ran ~26–35 s each before the shared-epoch cut existed.

### Injection-recovery and limit (simulated / derived)

| t_E (d) | \|M\| (M☉, model_prediction) | Cut-0 (ρ = 0.01 / 0.1) | pre-screen | recovered | ε per star (LF-weighted) | **Γ₉₅ per star per yr** |
|---|---|---|---|---|---|---|
| 3 | 1.7 × 10⁻³ | 0.43 / 0.37 | 0.22 / 0.20 | 0.07 / 0.13 | 0.009 / 0.016 | 11 / 6.2 × 10⁻⁶ |
| 10 | 1.8 × 10⁻² | 0.57 / 0.62 | 0.32 / 0.33 | 0.18 / 0.25 | 0.032 / 0.068 | 3.1 / 1.5 × 10⁻⁶ |
| 30 | 0.17 | 0.83 / 0.72 | 0.48 / 0.33 | 0.28 / 0.23 | 0.051 / 0.070 | 2.0 / 1.4 × 10⁻⁶ |
| 100 | 1.8 | 0.83 / 0.90 | 0.30 / 0.38 | 0.17 / 0.20 | 0.028 / 0.081 | 3.5 / 1.2 × 10⁻⁶ |
| 300 | 17 | 0.90 / 0.90 | 0.22 / 0.35 | 0.07 / 0.10 | 0.027 / 0.023 | 3.6 / 4.4 × 10⁻⁶ |

- **First W3 rate limit**: with zero survivors, the 95 % upper limit on the rate of umbra crossings (u₀ < 2 in
  |ε| Einstein radii) by an n = 1, ε < 0 lens is **Γ₉₅ ≈ 1.2–4.4 × 10⁻⁶ per monitored star per year** in gb22 for
  t_E = 10–300 d (|M| ≈ 0.02–17 M☉ in the stated geometry) and 0.6–1.1 × 10⁻⁵ at t_E = 3 d. Γ₉₅ = 3 / (N_s T ε);
  N_s = 3.5 × 10⁶ stars with 10 ≤ I ≤ 21.4, T = 8.61 yr. With the smallest N_s of the star-count model
  (2.4 × 10⁶) the limits are 1.43× weaker (`rate95_conservative`). 101 of 600 W3 injections are recovered, 0 of
  100 PSPL controls end as W3 survivors.
- Fractions are of all injections; ε per star weights them by a luminosity function ∝ 10^(0.319 I) and counts
  stars brighter than I = 14.2 (0.5 % of N_s) as undetectable. The weights favour faint sources, so the
  effective number of injections per cell is ~20–25 (`n_eff_lf`) and ε is uncertain by ~±25–75 % per cell
  (binomial); the ρ = 0.01 / 0.1 differences at fixed t_E are within that noise.
- Where efficiency is lost: Cut-0 (faint sources and short events); the shape cuts, mostly long events whose
  umbra spans seasons (second-deficit rule); then the **variable-baseline test (73 of 188 vetted injections)**,
  the eclipse model (6), season offsets (3) and the repeated-deficit test (2). The baseline test is the largest
  loss because 35 % of the 300 carrier light curves have χ²/dof > 2 about a constant on their own (red noise of
  difference photometry), so it is conservative for the limit. It is also sensitive to the carrier definition:
  a first run whose quiet-carrier cut required *both* white-noise significances above −5 (before the review fix
  that made `s_min` the weaker one) lost 34 of 231 to it, recovered 171 / 600, and gave Γ₉₅ ≈ 0.5–5 × 10⁻⁶.

### Assumptions and caveats

- **N_s** (no published value for gb22): N_s per Cut-0 object in Nunota et al.'s 20 fields (median 188, range
  131–240) × 18,599 = 3.5 × 10⁶ (2.4–4.5 × 10⁶). A Gaia DR3 RP power-law extrapolation over the same magnitude
  range gives ~6 × 10⁶ for a 2.18 deg² field, so the adopted N_s is not optimistic.
- **Luminosity function**: slope 0.319 dex/mag from Gaia DR3 RP counts at RP 13–17.5 within 0.3° of the field
  centre, extended to I = 21.4; MOA-Red ≈ RP ≈ I to ~0.3 mag. A steeper faint end would lower ε.
- **Cut-0 fidelity**: emulated on the light curve, not the images — the light-curve S/N stands in for the
  image S/N, and the spurious-detection filters (σ_x,y, moving objects, PSF χ²) are not emulated. A source
  that passes the emulation but not the real filters would be missing from the release, so ε is an upper
  bound on that step, and so the limit is optimistic by that unknown factor.
- **Carriers**: quiet Cut-0 light curves stand in for the difference light curves of constant stars, which are
  not in the release; the source's own photon noise is not added (sky- and blend-dominated noise assumed).
- **Event rate definition**: per star per year for umbra crossings (u₀ < 2); events with t₀ outside the data
  span are not injected. Mass scale as D-057 (n = 1, D_L = 4 kpc, D_S = 8 kpc, μ_rel = 5 mas/yr).
- Three vetting tests were added after the first pass left nine survivors (see Results); the injections
  measure their cost, and their removal would raise ε but leave nine ordinary-looking dips unexplained.
- One field (0.8 % of the release's light curves). The other 21 fields need the same chain (bulk tars of
  7–474 GB each; per-object files are also served).
- Not a statement about OGLE or the Mróz samples, and not combinable with the D-052 JWST limits without a lens
  population model.

## W5 count deficits

D-063; `src/jwst_anomaly/countmap.py` (count-map adapter, matched filter, injector), `scripts/w5_counts.py`
(`fetch`, `numcounts`, `predict`, `screen`, `vet`, `sheet`, `inject`, `limit`), tests in `tests/test_countmap.py`
and `tests/test_w5_counts.py`. Small result tables in `results/w5_counts/`; maps and injections under
`$JWST_ANOMALY_DATA/derived/w5_counts/` (not in git). Galaxy counts are **observed**; maps, filter amplitudes,
nulls, vetting and limits are **derived**; injected deficits are **simulated**; the deficit profile and θ_E(|M|) are
**model_prediction**s; every threshold and the lens geometry are **ASSUMPTION**s (module constants). A flag is an
anomaly to vet; it is not evidence of exotic physics.

### Prediction (model_prediction)

An n = 1, ε < 0 lens puts both images on the source's side, so every image radius x = θ/θ_E is reached: inside the
radial critical curve (x < 1) by the inner images of distant sources, demagnified as |μ| ≈ x⁴. Magnification bias
with the measured counts of the galaxy sample (`numcounts.ecsv`, d log N/dm = 0.37 at r = 23.5, 0.6 brighter
than r = 16 assumed) gives N_obs/N̄ (`profile.ecsv`, `countmap.deficit_profile`, |μ| capped at 30 near x = 1):

| x | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 | 0.9 | 1.05 | 1.2 | 1.5 | 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N_obs/N̄ | 0.02 | 0.14 | 0.56 | 0.86 | 1.00 | 1.06 | 0.94 | 0.86 | 0.94 | 0.98 | 0.99 |

An almost empty core (x ≲ 0.25), a 6 % excess at x ≈ 0.6–0.7 and a shallow dip at the critical curve. Net, about 10 %
of the galaxies inside θ_E are missing. Counts with α = 2.5 d log N/dm = 1 would give no signal at all; α < 1 at the
limit makes the dip at x ≈ 1 a deficit, and the steep bright end empties the core. With 37,400 galaxies deg⁻²
(r < 23.5), θ_E = 10″ leaves 0.09 missing galaxies (untestable per lens), θ_E = 2′ 13, 8′ 214, 32′ 3,400.

θ_E vs |M| (`theta_mass.ecsv`; Planck18, sources at z_s = 1, geometries ASSUMPTIONs; θ_E ∝ |M|^½):

| θ_E | 10″ | 2′ | 8′ | 32′ | 1° |
|---|---|---|---|---|---|
| D_L = 1 kpc | 1.2 × 10⁷ | 1.8 × 10⁹ | 2.8 × 10¹⁰ | 4.5 × 10¹¹ | 1.6 × 10¹² M☉ |
| D_L = 1 Mpc | 1.2 × 10¹⁰ | 1.8 × 10¹² | 2.8 × 10¹³ | 4.5 × 10¹⁴ | 1.6 × 10¹⁵ M☉ |
| z_L = 0.3 | 1.8 × 10¹³ | 2.6 × 10¹⁵ | 4.2 × 10¹⁶ | 6.7 × 10¹⁷ | 2.4 × 10¹⁸ M☉ |

θ_E between 10″ and 1° therefore means |M| ≈ 10⁷–10¹² M☉ at 1 kpc and 10¹⁰–10¹⁵ M☉ at 1 Mpc; a cosmological lens
needs ≳ 10¹³ M☉ (≳ 10¹⁶ M☉ for the θ_E this screen can test). For the 1 kpc geometry the "background" r < 23.5
population still is galaxies; Milky Way stars are excluded by the extended-model selection.

### Data and method

- **Count maps** (observed): Legacy Surveys DR10 Tractor via Data Lab TAP, aggregated server-side per HEALPix
  nest4096 pixel (0.74 arcmin²): galaxies (`type` ≠ PSF/DUP, `maskbits` = 0, dereddened r < 23.5), all primary
  sources and sources with any maskbit (unmasked fraction w = 1 − n_bad/n_all, ASSUMPTION), mean galaxy depth,
  E(B−V). Two DES-wide regions, desA (RA 20–40°) and desB (RA 50–70°), Dec −30° to −20°, b < −45°: 100 chunks,
  ~170 MB of FITS, fetched in ~2 h because Data Lab returned HTTP 502 for ~30 min and serves ~1 chunk/min.
  Pixels used: w ≥ 0.5, 5σ galaxy depth ≥ 24.3 (median 24.84), E(B−V) ≤ 0.1, not cut by the region border:
  **171.8 + 168.7 = 340.5 deg²**, 6.30 + 6.15 million galaxies.
- **Matched filter**: the pixel counts are painted on a 0.25′ equal-area (sinusoidal) raster; around every cell,
  a local least-squares fit D = w a (1 + A T) within 2.5 θ_E, T = profile − 1, so A = 1 is the predicted deficit
  and the local mean a absorbs large-scale gradients. Filter scales θ_E = 2′, 4′, 8′, 16′, 32′ (raster binned to
  θ_E/8). Z = A/σ_Poisson; positions need ≥ 80 % unmasked area in the aperture and in the inner disc.
- **Null** (no Gaussian assumption): local maxima of Z in each region. Every peak of both regions goes through the
  ordinary vetting tests below; each region's surviving peaks are the null of the other region at the same scale.
  A flag is a peak above the null region's highest peak; an exponential fit to the top 1 % of null peaks gives
  N_false, the expected number of null peaks at least as high among the search region's peaks. Detection needs
  N_false < 0.01 per region and scale (`cosmic_variance` otherwise). The Poisson σ underestimates the scatter by
  ×1.1 (2′) to ×5.8 (32′): galaxy clustering, which the null carries.
- **Vetting** (cheapest first; the same code vetoes injections): `mask` (> 20 % masked or missing area inside θ_E),
  `depth` (> 0.3 mag shallower than the region median), `depth_edge` (p90 − p10 of pixel depth > 0.5 mag within
  2 θ_E: a coverage boundary), `dust` (E(B−V) > median + 0.03), `bright_star` (Gaia DR3 G < 9 whose halo,
  2′ × 10^0.15(9−G), reaches the core 0.5 θ_E and could empty ≥ 10 % of it), `large_galaxy` (HyperLEDA D25 ≥ 1′
  within 0.5 θ_E + 2 D25: sky over-subtraction), `cluster` (Wen & Han 2024 M500 ≥ 3 × 10¹⁴ M☉, z ≤ 0.6, within
  0.5 θ_E + 3′: magnification-bias depletion), `cosmic_variance` (N_false ≥ 0.01).
- **Injection-recovery**: the predicted deficit (binomial thinning where the ratio < 1, Poisson additions where
  > 1, ratio averaged over each pixel) painted into the real pixel counts at random footprint positions, ≥ 2 × 3.5
  max(θ_E, 8′) apart, 3 realisations per region; θ_E = 2–32′, 77–884 injections each (5,232 in all). Recovered = a flag within
  max(θ_E/2, 1′) at any filter scale that survives every vetting test.

### Results (derived)

- 247,361 peaks, **40 flags** (7 desA, 33 desB), **0 survivors** (`vetting.ecsv`; `contact_sheet.png` shows the
  eight with the smallest N_false). The six with N_false < 0.01 are all explained:
  - three 2′/4′ flags 7–8′ from NGC 1398 (D25 = 7.2′): sharp-edged, CCD-sized deficits and excesses of sky
    over-subtraction (`large_galaxy`; first seen by a parallel cloud run, `ngc1398_dr10_cutout.jpg`);
  - the 8′ flag at (64.07°, −23.43°): masked holes next to a G < 9 star (`bright_star`);
  - the 16′ flag at (53.37°, −28.05°): the edge of the DES-SN C3 deep field (`depth_edge`);
  - the 4′ flag at (63.75°, −24.95°), N_false = 0.004, which survived the first vetting pass: inside 2′ both the
    galaxies (1.0–1.6 per pixel against 7.6) and *all* sources (14–26 against 57) drop along a straight edge of a
    deep-coverage tile (depth 24.9 → 25.9). A lens would not remove foreground stars or follow a tile edge; the
    `depth_edge` test was added for it, and the injections were run after the change.
- The other 34 flags are `cosmic_variance` (N_false ≥ 0.01), 12 of them around NGC 1398 and 3 in the C3 field.

### Efficiency and limits (95 %, zero survivors, Poisson 3.0; derived)

| θ_E | 2′ | 3′ | 4′ | 6′ | 8′ | 12′ | 16′ | 24′ | 32′ |
|---|---|---|---|---|---|---|---|---|---|
| injected | 884 | 871 | 866 | 879 | 873 | 421 | 241 | 120 | 77 |
| flagged | 0 | 7 | 155 | 441 | 782 | 395 | 235 | 117 | 76 |
| after vetting | 0 | 0 | 0 | 4 | 433 | 284 | 172 | 88 | 53 |
| ε | 0 | 0 | 0 | 0.005 | 0.50 | 0.67 | 0.71 | 0.73 | 0.69 |
| n₉₅ (deg⁻²) | — | — | — | 1.9 | 0.018 | 0.013 | 0.012 | 0.012 | 0.013 |

- **Sky density of n = 1, ε < 0 lenses with θ_E = 8–32′: n₉₅ ≈ 0.012–0.018 deg⁻²** (fewer than one per ~55–85 deg²)
  over 340.5 deg² of DES-depth DR10 sky. In the stated geometries this is |M| ≈ 3 × 10¹⁰–5 × 10¹¹ M☉ at 1 kpc and
  3 × 10¹³–5 × 10¹⁴ M☉ at 1 Mpc (`limits.ecsv`).
- **Blind below θ_E ≈ 6′** (|M| ≲ 1.6 × 10¹⁰ M☉ at 1 kpc): the empty core (x ≲ 0.25) holds < 15 galaxies and the
  required N_false < 0.01 is deep in the clustered tail; 18 % (4′) and 50 % (6′) of injections are flagged but end as
  `cosmic_variance`.
- Where efficiency goes at θ_E ≥ 8′: `cosmic_variance` (8′: 264 of 873), then `bright_star`, `depth_edge`,
  `large_galaxy`, `cluster` (each 1–6 %). Large injections are found mostly by the 8′ filter, which matches their
  empty core.

### Assumptions and caveats

- Point lens, n = 1, isolated, all r < 23.5 galaxies behind the lens (for z_L = 0.3 a foreground fraction would
  dilute the deficit; not modelled); demagnified galaxies are assumed to keep their extended Tractor type.
- The profile uses the counts of a pilot area of the same selection; the bright end (r < 16) slope 0.6 is
  assumed. The core is what the screen sees, and it depends on the counts at r ≲ 18.
- w = 1 − n_bad/n_all is a proxy for the unmasked area (DR10 randoms are not on Data Lab; D-063).
- The `depth_edge` test was added after inspecting a survivor; its cost (2–9 % of injections) is in ε.
- Each region calibrates the other; desB holds more artefacts (NGC 1398, DES-SN C fields, deep tiles). Peaks that
  the ordinary tests remove are left out of both the null and the search sample.
- No limit below θ_E ≈ 6′ or above 32′ (the regions are 10° high; 1° needs a larger contiguous area).
- Not combinable with the W1 radial/shear limits (different mass and θ_E ranges and selection).

### Euclid Q1: deeper counts do not open θ_E < 6′; shapes would (2026-10-09, D-065)

`scripts/w5_euclid_feasibility.py` → `results/w5_counts/euclid_q1_feasibility.json`. Galaxy densities are
**observed** (IRSA TAP counts in 0.1° discs inside EDF-F, EDF-S and EDF-N); the gains are **model_prediction**.
- Extended VIS detections (`vis_det`, clean flags, `point_like_prob` < 0.1; ASSUMPTION): 2.8 × 10⁴ (VIS < 23.5),
  6.8 × 10⁴ (< 24.5) and 1.07 × 10⁵ deg⁻² (< 25.0); field-to-field spread ±10 %. VIS < 24.5 is 1.81× the DR10
  r < 23.5 density; the count slope d log N/dm = 0.378 (DR10 0.365), so the predicted deficit profile is the same.
- The screen's scatter is mostly galaxy clustering, which does not shrink with depth: DR10 Z scatter is 1.1× (2′),
  1.6× (4′) and 2.4× (8′) the Poisson value. With 1.81× the galaxies and the same (half the) clustering variance
  (ASSUMPTION range), the count S/N gains only ×1.24 (1.36) at 2′, ×1.10 (1.39) at 4′, ×1.04 (1.40) at 8′.
  DR10 reached ε ≈ 0.5 at 8′; its 4′ S/N is 0.72 and its 2′ S/N 0.52 of that, so Euclid Q1 counts reach
  0.79–1.0 at 4′ and 0.65–0.71 at 2′: at best a floor near 4–6′ instead of 6–8′, on 63 deg² instead of 340 deg²
  (n₉₅ ≥ 5× weaker above 8′). Not worth a ~3.5 h row fetch.
- A radial-shear test (point-mass γ = (θ_E/θ)², mean 0.46 over θ_E–2θ_E; σ_γ = 0.3 per component and S/N ≥ 6,
  ASSUMPTIONs) reaches θ_E ≈ 0.3′ ≈ 18″ with VIS < 24.5 shapes (≈ 4 × 10⁷ M☉ at 1 kpc, 4 × 10¹⁰ M☉ at 1 Mpc),
  an order of magnitude below the count floor in θ_E. The sign of the tangential shear separates a negative-mass
  lens (radial) from every ordinary foreground mass (tangential). MER `ellipticity` / `position_angle` are
  SExtractor image moments without PSF correction (`position_angle` is CCW from the image x axis), so the test
  needs a PSF-anisotropy check on stars first.

## W1/W2 in rejected lensed-quasar pairs

D-064; `scripts/w12_niq.py screen --sheet` (`--repin` only after inspecting changed inputs); outputs `results/w12_niq/`
(`systems.ecsv`, `summary.json`, contact sheets); inputs pinned in `data/manifests/w12_niq_inputs.json`.

**Question.** Lens searches reject same-redshift quasar pairs when no lens galaxy is seen: Lemon et al. 2023's
"unclassified quasar pairs" (UQP, "akin to NIQs") and "QSO pair" classes, and the SQLS "no lens(ing) object",
"QSO pair" and "binary QSO" rejections. A dark deflector (W2, W1's empty centre) would hide in exactly these lists.
The ordinary explanations are binary quasars, unrelated pairs and lens galaxies below the depth.

**Method.**
- Inputs:
  - rejected: 123 pairs (Lemon "UQP (?)" is undecided and left out);
  - control: 106 real quasar lenses (Lemon lens/quad, SQLS "SDSS lens"/"known lens"; Lemon "lens (?)" is
    undecided and left out);
  - entries within 3″ are merged transitively. A group containing a catalogued lens is a control. A rejection that
    Lemon classifies as a non-pair (QSO+star, projected, …) is dropped: this applied to J0947+0247, which Lemon
    classifies as "QSO + star". Undecided Lemon classes ("?") and SQLS non-pair companion rows do not veto.
- Exclusions:
  - 181 entries with catalogued separations > 3″: a catalogue position is one image, so the second image must
    fall inside the 3″ image search;
  - lenses the pair test cannot decide: the SQLS "component" rows of one cluster lens and Lemon's lensed galaxies.
    These still promote a merged rejection to control.
- The D-056 chain unchanged: DR10 brick coverage and depth, Tractor boxes, the quasar pair test (two PSF images
  ≥ 2″ apart, a deflector between them), and the required lens magnitude from the D-056 Faber–Jackson calibration.
- A system counts as decided only if the LS image pair is the catalogued pair: separations agree within 0.5″.
- Vetting columns (thresholds are ASSUMPTIONs, `NiqParams`, recorded in `summary.json` with every count below):
  - image colours: |Δ(g − z)| ≤ 0.5;
  - two quoted redshifts: |Δz| / (1 + z) > 0.01 means two quasars;
  - a Hennawi et al. 2006 binary-quasar match within 3″.

**Result (derived).**

| Sample | covered | blended | too close (< 2″) | decided | deflector | none | none, colours match |
|---|---|---|---|---|---|---|---|
| rejected | 88 | 17 | 46 | 24 | 0 | 24 | 14 |
| control (real lenses) | 79 | 45 | 29 | 5 | **0** | 5 | 4 |

- **The test does not find real lens galaxies at these separations.** All 5 decided control lenses (1.9–2.6″;
  J0628−7448, J1550+0221, J2308+3201, SDSS J1322+1052, SDSS J1515+1511) give "none". Their required typical
  m_z ≈ 19.2–19.8 is ~3 mag brighter than the depth, yet the contact sheet shows the lens light blended into the
  images, and Tractor fits the blend as two point sources. Measured efficiency for an ordinary lens: **0 / 5 (95 %
  upper bound 0.45)**. A further 45 control lenses are "blended" (fewer than two PSF images). A "none" at ≤ 3″ in LS
  DR10 therefore carries no information about a dark deflector.
- Ordinary explanations among the 24 rejected "none" pairs:
  - 10 have mismatched colours (blue+orange on the contact sheet: unrelated objects, not lens images);
  - J0740+2926 and J1035+0752 are catalogued binaries (Hennawi et al. 2006);
  - J1212+0912 has two redshifts (1.686, 1.600).
- **No limit and no candidate.** The 11 remaining colour-matched pairs are untestable here: J0130+0725, J0728+2607,
  J0941−2443, J1428+0500, J2355−4553, J0927+2113, J1242+2543, J0904+1134, J0942+2310, J1324+2823 and J1711+2929.
  Their discovery papers' deeper follow-up already found no lens galaxy, and SQLS calls J1242 and J0942 quasar pairs
  from their spectra. Telling a binary from a dark lens needs spectral comparison or HST/Euclid/HSC imaging.
- Most rejected pairs are closer than LS can resolve (46 "too close"; the Lemon UQP median separation is 1.22″).
- One pair (J0041−5350) is not decided: its LS pair (3.3″) is not the catalogued one (1.1″).

**Consequence for D-056.** D-056's quasar-class limits assume that a dark lens gives "none" (true). Its "none"
systems were explained by literature lens galaxies, never by LS. This measurement confirms the D-056 vetting finding
(13 of 16 "none" systems had a literature lens galaxy): at ≤ 3″, LS DR10 cannot show a lens galaxy, so only
literature or deeper imaging decides.
