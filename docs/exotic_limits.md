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

## W3 inverted microlensing / dimming (multi-epoch)

D-052; `scripts/dimming_screen.py` (`fetch`, `screen`, `inject`, `forced`, `combine`), fields and epochs in
`configs/dimming_screen.yaml`, tests in `tests/test_dimming_screen.py`. The signature is D-047's W3: an n = 1 lens
with ε < 0 crossing a compact source. The lensed flux is 0 inside the umbra (β < 2, for 2 t_E √(4 − u₀²)) between
two caustic spikes (×7.0 for ρ = 0.01, ×2.35 for ρ = 0.1; `exotic_sim` as merged in PR #67). The umbra is the
robust part of the signal; the spike heights depend on ρ.

### Data (observed)

| Field | Epochs (F200W + F444W) | Baseline | Master sources | Monitored compact sources | Overlap area |
|---|---|---|---|---|---|
| NEXUS-Center (5105) | 8: 2024-09-12 … 2026-03-28 | 1.54 yr | 81,136 | 154 | 79.5 arcmin² |
| MACS0416 (1176 ×3, 1208, 6882 F444W only) | 5: 2022-10-07 … 2026-01-10 | 3.26 yr | 7,708 | 44 | 10.8 arcmin² |
| Abell 2744 (2561 o001/o002/o006) | 3: 2022-11-02 … 2024-07-31 | 1.74 yr | 17,669 | 63 | 14.0 arcmin² |

All 31 catalogues are level-3 `_cat.ecsv` with jwst 3.0.0 (`data/manifests/dimming_*.ecsv`, 339 MB). No `_i2d.fits`
was downloaded; forced photometry uses S3 byte-range cutouts. El Gordo has one shared band across epochs and was
not used. JADES is left for later.

### Method (every threshold an ASSUMPTION, `Params` in the script)

1. **Light curves (derived).** Each catalogue is tied to the deepest detection-band catalogue by the D-027 global
   frame shift. Master list: rows with aper50 S/N ≥ 5 from any epoch. Mutual matches within 0.3″. Per band and
   epoch, the zero point is the median flux ratio of S/N ≥ 30 pairs. Errors are scaled by the robust scatter of
   epoch-to-epoch differences (1.52–1.73 in all bands and fields; D-027 found 1.2–1.5), plus a 3 % floor.
2. **Catalogue non-detections.** A covered epoch without a match counts as flux 0 only where the source's peak
   flux would have been ≥ 10σ there. Catalogues are assumed complete at S/N ≥ 10.
3. **Flags.**
   - `vanish`: an epoch pair with S/N ≥ 10 in one epoch and < 3σ in the other, with a ≥ 5σ drop, in every band
     testable in that pair (achromatic).
   - `dim_achromatic`: one epoch ≥ 20 % below the median of the others at ≥ 5σ in both bands, with drops equal
     within 3σ.
   - `rise_dip_rise`: with ≥ 3 epochs, an epoch ≥ 20 % and ≥ 5σ fainter than an earlier and a later one.
4. **Catalogue-level ordinary tests, cheapest first.** Gaia DR3 star mask (D-027 radii; a faint Gaia match to the
   source itself is kept as `gaia_star`), edge proxy (neighbour count < 0.5 × median), blend (neighbour < 0.5″),
   single-epoch detection (persistence suspect, D-039), sharper than the PSF (CI_70_30 < 1.9).
5. **Forced photometry (`forced`).** For every flag that passes the catalogue tests: a 0.15″ aperture at the
   recentred position in every epoch (`transient_forced.measure`), zero points and noise scale from random
   controls, then the same flag logic again. A flag is confirmed when the same flag type reappears. Cutout tests
   come next: edge / no data, WHT < 0.5 × typical, hexagonal spike statistic ≥ 3 (D-018/D-021). The last step is
   a SIMBAD/NED match.
6. **Compact sources monitored** (the population of the limit): 1.9 ≤ CI_70_30 ≤ 2.7 in the reference epoch (the
   F200W stellar locus is 2.0–2.5), no catalogue veto, and S/N ≥ 10 in ≥ 2 detection-band epochs. These are mostly
   Galactic stars in NEXUS and a mix of stars and compact sources in the clusters.

### Results (derived)

| Field | Catalogue flags (vanish / dim / rdr) | After catalogue tests | Forced-measured | Forced-confirmed | After cutout tests | Survivors after visual check |
|---|---|---|---|---|---|---|
| NEXUS | 3,793 (3,197 / 264 / 1,276) | 1,095 (13 compact) | 100 (all 13 compact + 87 random) | 0 | 0 | 0 |
| MACS0416 | 854 (714 / 67 / 310) | 375 (10 compact) | 375 | 6 | 2 | 0 |
| Abell 2744 | 530 (525 / 49 / 6) | 32 (1 compact) | 32 | 0 | 0 | 0 |

- Most catalogue flags are catalogue effects. Sources are missing from one epoch's catalogue (deblending,
  segmentation) but present in its image; others sit at mosaic edges (cutouts inspected).
- The two MACS0416 flags that passed every automatic test, `d06533` and `d07517` (RA 64.0206–64.0207), lie
  0.5–1″ from the saturated star `d03349`. Each i2d grid has a different orientation, so the star's PSF wings and
  spikes cross the fixed aperture differently in each epoch, and the F200W flux jumps up and down by ×10.
  The cutout spike statistic is centred on the target's own peak, so it does not catch this. **Ordinary
  (instrumental); no candidate.**
- The controls flagged 7 of 300 (MACS0416), 1 of 150 (NEXUS) and 0 of 14 (Abell 2744). This is the forced stage's
  false-positive rate on ordinary sources.
- Wall time: `screen` 76 s (NEXUS), 4 s (MACS0416), 9 s (Abell 2744); `inject` 136 / 20 / 34 s; `forced` about
  32 min (MACS0416, 675 positions × 9 images), 20 min (NEXUS, 250 × 16), 3–9 min (Abell 2744).
- The forced noise scale is calibrated for MACS0416 (2.15 / 1.99, from 300 controls). NEXUS (150 controls) and
  the final Abell 2744 run (14 controls) had too few controls covering ≥ 2 epochs, so their forced errors stay
  ERR-based (scale 1). An earlier Abell 2744 run with 300 controls gave 2.45 / 3.27 and also confirmed 0 flags.

### Injection-recovery (simulated)

- Every monitored compact source that is not already flagged gets 20 copies per model.
- W3 (n = 1, ε < 0) uses `exotic_sim.inject_light_curve` over the real epoch times. Each copy draws
  u₀ ~ U[0, 2) and t₀ ~ U[t_first − 2t_E, t_last + 2t_E].
- Injected flux: F·f_obs + (1 − F)·N(0, σ_noise), the same factor F in both bands. A vanished epoch therefore
  keeps sky noise only.
- Recovered: any flag and no catalogue veto. For MACS0416, the event must also still be flagged with errors
  inflated by the forced/catalogue noise ratio (1.30 / 1.15).
- Plain achromatic dimming of 20/50/100 % in one random epoch is the sensitivity floor.

Efficiency, all magnitudes (source-weighted):

| Model | NEXUS | MACS0416 | Abell 2744 |
|---|---|---|---|
| W3 t_E = 0.01 yr, ρ = 0.01 / 0.1 | 0.06 / 0.05 | 0.03 / 0.02 | 0.04 / 0.05 |
| W3 t_E = 0.03 yr | 0.17 / 0.16 | 0.07 / 0.07 | 0.10 / 0.09 |
| W3 t_E = 0.1 yr | 0.44 / 0.43 | 0.14 / 0.15 | 0.30 / 0.30 |
| W3 t_E = 0.3 yr | 0.60 / 0.59 | 0.30 / 0.26 | 0.55 / 0.59 |
| W3 t_E = 1 yr | 0.36 / 0.36 | 0.48 / 0.45 | 0.36 / 0.36 |
| W3 t_E = 3 yr | 0.17 / 0.18 | 0.26 / 0.22 | 0.18 / 0.16 |
| dimming 20 % / 50 % / 100 % | 0.13 / 0.71 / 1.00 | 0.07 / 0.62 / 1.00 | 0.02 / 0.40 / 1.00 |

Pooled over fields, by detection-band AB magnitude (ρ = 0.1):

| Model | 15–22 | 22–24 | 24–25 | 25–26 | 26–27 | 27–29 |
|---|---|---|---|---|---|---|
| W3 t_E = 0.1 yr | 0.39 | 0.42 | 0.40 | 0.34 | 0.28 | 0.24 |
| W3 t_E = 0.3 yr | 0.56 | 0.61 | 0.57 | 0.53 | 0.48 | 0.48 |
| W3 t_E = 1 yr | 0.45 | 0.42 | 0.40 | 0.34 | 0.29 | 0.27 |
| dimming 20 % | 0.20 | 0.14 | 0.06 | 0.01 | 0.00 | 0.00 |
| dimming 50 % | 0.88 | 0.79 | 0.91 | 0.56 | 0.16 | 0.09 |
| dimming 100 % | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

Short events fall between epochs, and events longer than the baseline cover every epoch, so no epoch pair
differs. The screen is most sensitive at t_E ≈ 0.3 yr. Efficiency depends weakly on ρ, because the umbra, not
the spike, drives recovery.

### Upper limits (derived; 0 surviving events, Poisson 95 % = 3.0)

- Rate per source per year: 3 / Σ N ε (T + 4t_E), with t₀ drawn over T + 4t_E.
- Per source per epoch pair: 3 / Σ N ε n_pairs, with n_pairs = 28 / 10 / 3.
- Per deg² per epoch pair: 3 / Σ A n_pairs ε̄.
- |M| = (t_E / 12 yr)² M☉ is a **model_prediction** for z_l = 0.4, z_s = 2 and v⊥ = 1000 km/s (D-047; re-derived
  here as 12.09 yr). It does not apply to the Galactic stars that dominate the NEXUS compact sample.

| t_E (yr) | \|M\| (M☉) | ρ | exposure (source yr) | per source per yr | per source per epoch pair | per deg² per epoch pair |
|---|---|---|---|---|---|---|
| 0.01 | 6.9e-07 | 0.01 | 23.8 | 0.13 | 0.011 | 77 |
| 0.03 | 6.3e-06 | 0.01 | 65.9 | 0.045 | 0.0038 | 27 |
| 0.1 | 6.9e-05 | 0.01 | 193.5 | 0.016 | 0.0015 | 10.8 |
| 0.3 | 6.2e-04 | 0.01 | 413.2 | 0.0073 | 0.0011 | 7.8 |
| 1 | 6.9e-03 | 0.01 | 582.0 | 0.0052 | 0.0017 | 12.6 |
| 3 | 6.2e-02 | 0.01 | 682.9 | 0.0044 | 0.0034 | 26 |

ρ = 0.1 gives the same limits within 10 % (`limits_combined.ecsv`).

**Assumptions and caveats.**
- The isolated-lens light curve has no cluster macro-magnification or shear (D-047).
- The lensed fraction is 1 (blend = 1).
- The injection is applied to catalogue fluxes. Catalogue re-detection and deblending of the injected epoch are
  not simulated. The forced stage is emulated only through the MACS0416 error inflation.
- The monitored population is small (261 sources), so the limits are per source and weak. A rate of ≲ 0.005
  umbra crossings per compact source per year for t_E ≈ 0.3–3 yr is the main result.
- The forced stage measured only 100 of NEXUS's 1,095 catalogue-passing flags (all 13 compact ones). The limit
  concerns compact sources only, so it is unaffected, but extended-source flags in NEXUS are not fully vetted.
- No SN/TNS check was run (no survivors). SIMBAD/NED is wired into `forced` but was not needed.
