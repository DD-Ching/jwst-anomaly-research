# Exotic-lens limits from injection-recovery

What the null exotic screens (docs/exotic_lensing.md "Screens and results") rule out. We inject simulated signals
into the real catalogues and run the unchanged screens on them. Every injected object is **simulated**. Recovery
efficiencies and limits are **derived** from those simulations. Masses are a **model_prediction** of the point-lens
formula. A limit describes how sensitive a screen is to a hypothetical lens population. It is not evidence that any
such lens exists or does not exist beyond what the stated efficiency covers.

## W1 negative-mass lenses (radial screen)

D-049; `scripts/inject_radial.py`, tests in `tests/test_inject_radial.py`. The signature is D-047's W1: an n = 1 lens
with ε < 0. There are no images of sources within β < 2 θ_E (the umbra), and two radially stretched images on the
source's side of the lens for every source beyond that.

### Method

1. **Fields and inputs.** The eight fields where `exotic_screens.py radial` came back null. Each uses the
   catalogue, photo-z, `--max-radius` and spike veto from its field doc (`FIELDS` in the script).
   - The models are those of the field docs, in the JWST frame.
   - The base screen reproduces each documented result: number of `anti` arcs, maximum lines and p_random (table
     below). Abell 2744 is the exception: its doc predates the D-034 spike veto. The current code drops 42 spike
     segments and gets 134 arcs, max 5 lines, p 0.505 (null).
   - Abell S1063 runs without photo-z, the configuration of D-043 that needs no 75 MB tarball. The documented
     result with photo-z is also null.
2. **Lens placement.** Each lens sits at a uniform random position in the **screened footprint**, at the cluster
   redshift (ASSUMPTION). The footprint is the set of 1″ grid points within `--max-radius` of the model centre
   that have a catalogue source within 4″ (ASSUMPTION).
3. **Lensed sources.** The lenses act on the field's own catalogue rows behind the lens position, so the number of
   lensed sources follows the real local surface density of background sources.
   - A row counts as background under the screen's own rule: either it has no photo-z, or its photo-z z160 is
     above z_l + 0.1. Point sources brighter than AB 20 (stars) are never lensed.
   - Rows with β < 2 θ_E are removed (umbra).
   - Rows with 2 ≤ β ≤ 4 θ_E (the D-047 range, ASSUMPTION) are replaced by their two images from
     `exotic_sim.inject_images(n = 1, sign = −1)`.
   - Rows beyond 4 θ_E are left unchanged.
   - Expected number lensed: Σ_bg × 12π θ_E². Σ_bg is the lensable-row density in the screened area,
     0.046–0.094 arcsec⁻². The realised means agree (table).
4. **Painting each image** (ASSUMPTIONs):
   - **Shape.** Start from the source's catalogue second moments (`semimajor_sigma`, `ellipticity`,
     `sky_orientation`). Remove an isotropic PSF, map the result by A⁻¹, then add the PSF back. A⁻¹ has the signed
     eigenvalues 1/λ_r along the lens direction and 1/λ_t across it. The PSF σ is the 1st percentile of
     `semiminor_sigma` among S/N > 50 sources: 1.0–1.44 px (**derived**).
   - **Brightness and noise.** Magnitude − 2.5 log₁₀|μ|; isophotal area × |μ|, because lensing conserves
     surface brightness; S/N × √|μ|, assuming background-limited noise.
   - **Detection.** Images below the catalogue's own S/N floor are dropped as undetected.
   - **Orientation.** A round source gives the simulator's radial position angle exactly (tested).
5. **Screen.** The injected catalogue goes through `exotic_screens.radial_candidates`, the code `cmd_radial` itself
   uses, and then `line_counts` and `convergence_peaks`, with every default threshold.
   - The null is the screen's: 200 draws of each arc's position angle inside its `anti` window.
   - The real arcs' draws are cached once per field, with the screen's seed and draw order. Each injection then
     subtracts the removed arcs and adds fresh draws for the injected ones. This gives the same distribution as a
     full re-run, not bit-identical numbers.
6. **Recovery.** A lens counts as recovered when a convergence peak with p_random < 0.05 lies within 2″ of the
   injected centre (ASSUMPTION: 4 grid steps, twice the 1″ line tolerance).
   - We used 200 lenses per field and θ_E, which gives a binomial σ of at most 3.5 %.
   - θ_E ∈ {0.3″, 1″, 3″} is the D-047 recommendation. 6″ and 10″ were added to find where the screen starts to
     respond.
7. **Limit.** The real screens found zero surviving detections. The 95 % Poisson upper limit on the surface
   density of W1 lenses is therefore Σ < 2.996 / Σ_fields(ε_f A_f), where ε_f is the field's recovery efficiency
   and A_f its screened area. The sum runs over the fields with photo-z only (see "Upper limits").

Command: `python scripts/inject_radial.py --fields smacs0723 elgordo abell2744 macs0416 macs1149 macs0717
abell370 abells1063 --n-inject 200 --theta-e 0.3 1 3 6 10`. It writes `trials.ecsv` and `summary.json` per field and
`limits.json` to `outputs/inject_radial/`. Abell 370 fetches its Gaia DR3 stars with `scripts/gaia_stars.py`
(7 stars).

### Base screens and screened areas (`derived`)

| Field (model) | z_l | Screened area (arcmin²) | Σ_bg (arcsec⁻²) | `anti` arcs | max lines, p_random | lines needed for p < 0.05 | Wall time (s) |
|---|---|---|---|---|---|---|---|
| SMACS 0723 (ICLv2, DJA photo-z) | 0.39 | 2.99 | 0.055 | 31 | 4, 0.965 | 7 | 447 |
| El Gordo (Caminha+23, DJA) | 0.87 | 4.83 | 0.046 | 37 | 4, 0.64 | 6 | 572 |
| Abell 2744 (Bergamini+23, DJA) | 0.31 | 11.51 | 0.067 | 134 | 5, 0.505 | 6 | 877 |
| MACS0416 (CATS, CANUCS photo-z) | 0.40 | 6.17 | 0.094 | 120 | 7, 0.225 | 8 | 506 |
| MACS1149 (CATS, CANUCS) | 0.54 | 6.44 | 0.056 | 61 | 3, 1.0 | 5 | 344 |
| MACS0717 (CATS, no photo-z) | 0.55 | 7.58 | 0.075 | 68 | 5, 0.70 | 7 | 531 |
| Abell 370 (CATS, CANUCS, Gaia veto) | 0.38 | 6.12 | 0.091 | 100 | 5, 0.495 | 7 | 411 |
| Abell S1063 (CATS, no photo-z) | 0.35 | 5.62 | 0.088 | 46 | 4, 0.435 | 5 | 322 |
| **Total** | | **51.2** (0.0142 deg²) | | | | | |

Wall time covers all five θ_E (1,000 injections), with four fields running in parallel on 4 cores. The base screen
took 14–145 s of it. A sequential SMACS run of the three D-047 θ_E took 128 s.

### Recovery efficiency (recovered / 200; `derived` from `simulated` injections)

| Field | θ_E = 0.3″ | 1″ | 3″ | 6″ | 10″ |
|---|---|---|---|---|---|
| SMACS 0723 | 0 | 0 | 0 | 4 | 13 |
| El Gordo | 0 | 0 | 0 | 3 | 7 |
| Abell 2744 | 0 | 0 | 0 | 13 | 26 |
| MACS0416 | 0 | 0 | 2 | 8 | 13 |
| MACS1149 | 0 | 0 | 1 | 8 | 21 |
| MACS0717 | 0 | 0 | 0 | 11 | 30 |
| Abell 370 | 0 | 0 | 2 | 7 | 14 |
| Abell S1063 | 0 | 0 | 5 | 26 | 57 |
| Mean lensed sources per lens (range over fields) | 0.2–0.3 | 1.7–3.5 | 15–30 | 55–103 | 133–268 |
| Mean injected arcs entering the screen | ≤ 0.03 | 0.14–0.32 | 1.2–2.8 | 4.7–8.6 | 11.8–23.9 |
| Lenses with ≥ 3 lines at the centre | 0–3 % | 0–4.5 % | 7–22 % | 18–43 % | 28–52 % |

### Upper limits (95 %, zero detections; `derived`)

| θ_E | \|M\| (`model_prediction`; z_l = cluster, z_s = 2, Planck18) | Σ ε_f A_f (deg²) | Σ(W1 lenses) < (deg⁻²) |
|---|---|---|---|
| 0.3″ | 1.4–4.3 × 10¹⁰ M☉ | 0 | no limit (ε = 0 in every field) |
| 1″ | 1.6–4.7 × 10¹¹ M☉ | 0 | no limit |
| 3″ | 1.4–4.3 × 10¹² M☉ | 4.3 × 10⁻⁵ | 7.0 × 10⁴ |
| 6″ | 5.6–17 × 10¹² M☉ | 4.4 × 10⁻⁴ | 6.8 × 10³ |
| 10″ | 1.6–4.7 × 10¹³ M☉ | 9.3 × 10⁻⁴ | 3.2 × 10³ |

The limits use only the six fields with photo-z (38.1 arcmin²). In MACS0717 and Abell S1063 every non-star row
counts as lensable, so cluster members and foreground galaxies near an injected lens are painted as W1 images, and
the efficiency is biased high (Abell S1063 has the highest efficiency of all eight fields). `limits.json` reports
the all-field values (8.2 × 10⁻⁵, 7.6 × 10⁻⁴ and 1.7 × 10⁻³ deg²; < 3.6 × 10⁴, 3.9 × 10³ and 1.8 × 10³ deg⁻²)
separately as `optimistic_all_fields`.

Mass conversion: |M| = θ_E² c² D_L D_S / (4 G D_LS). The low end of each range is Abell 2744 (z_l = 0.31), the
high end El Gordo (0.87). `theta_e_to_mass(2.2″, z_l = 0.4)` reproduces D-047's 10¹² M☉ (tested).

**Reading.** The radial screen is effectively blind to W1 lenses at the D-047 amplitudes. At 0.3″ and 1″ nothing
was recovered, and at 3″ 10 of 1,600 lenses were. The limits that do exist, at 3–10″, are weak.
- **Comparison (`derived`, rough).** Takahashi & Asada (2013) limit negative masses above 10¹² M☉ to
  n < 10⁻⁴ h³ Mpc⁻³. Spread over 0 < z < 1 (Planck18 comoving volume, 3.98 × 10⁶ Mpc³ deg⁻²), that is about
  120 deg⁻².
- Our best limit, 3.2 × 10³ deg⁻² at about 2 × 10¹³ M☉, is about 27× weaker. It also holds only for lenses near
  the cluster redshift.

**Why the efficiency is low.** The cut that dominates is the screen's elongation cut, not a lack of sources.
- For a round source, the image axis ratio is (x² + 1)/(x² − 1), where x is the outer image's position in θ_E.
  Ellipticity ≥ 0.5 needs β ≲ 2.31 θ_E, about 11 % of the [2, 4] θ_E annulus.
- `anti`, measured against the *cluster's* tangential direction, keeps only images whose direction to the W1 lens
  lies within ±30° of the cluster-radial direction. That is about a third of them.
- The p < 0.05 threshold comes from each field's own null and needs 5–8 converging lines. The ≥ 3 lines that
  define a peak are not enough.
- In practice, 1–3 injected arcs per 3″ lens reach the screen.

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
- **Blending.** The two images of a source are √(β² − 4) θ_E apart, i.e. 0–3.5 θ_E. At θ_E = 0.3″ they would
  blend, and the faint inner image (|μ| below 0.1 for β ≳ 2.4) would merge with the outer one. Painting them as separate rows is optimistic. This does
  not matter, because 0.3″ has zero efficiency anyway. Painted images can also overlap real neighbours, which the
  catalogue would deblend differently.
- **Incompleteness.** Stretched faint images lose surface-brightness-limited area, and a real pipeline may
  fragment or miss them. Only images of catalogued sources are painted, so sources below the detection limit that
  lensing would raise into the catalogue are ignored. Outer images mostly have |μ| ≈ 1–1.7 (it diverges only
  at β → 2), so this effect is small.
  Shapes come from Gaussian second moments with an isotropic PSF, and no pixel-level rendering was done.
- **Isolated lens.** We assume no coupling to the cluster's shear or magnification (D-047). The observed shapes
  already contain the cluster's shear, and the W1 Jacobian is applied on top of it.
- **Redshift.** The lens is placed at the cluster redshift only; other redshifts are not covered. Sources without
  a photo-z count as background, which is why the two fields with no photo-z are left out of the limits. Rows beyond 4 θ_E stay unlensed, which loses a little weak radial stretching.
- **Area.** The footprint extends up to 4″ past the catalogue edge, so the area is overestimated by up to about
  10 %, which makes the limits slightly too strong.
- **Statistics.** The null is computed by an incremental update, statistically equivalent to a full re-run but not
  bit-identical. Efficiency counts of 0/200 give a 95 % upper bound of 1.5 % per field and θ_E.

### Cases inspected

We plotted catalogue ellipses for two MACS0416 θ_E = 3″ injections (scratch plot, not committed).
- **Recovered** (8 lines, p 0.030): one source near the caustic gave a long radial pair, and three more pairs
  pointed at the centre.
- **Missed** (91 sources lensed, 10 injected arcs, 5 lines, p 1.0): most images were small and weakly stretched,
  so the sources' own orientations dominated. Only 5 lines met within 2″, short of the 8 that MACS0416's null
  requires.
- Both plots look as the simulator predicts: an empty umbra apart from unlensed foreground and member rows, inner
  images bunched near the centre, and outer images stretched radially.

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
