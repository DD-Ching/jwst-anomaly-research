# Exotic-lens limits from injection-recovery

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
   - Each image gets |μ| × the row's SED, plus fresh Gaussian noise drawn from the row's own errors in every band.
     The errors stay the row's (sky-limited).
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
   - **Background-aware limit:** s₉₅ / Σ εA, where s₉₅ = 72.2 is the classical 95 % upper limit on a Poisson
     signal. The inputs are 355 observed orphans over a known background of 315.4 (null (e), summed). This is the
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
- **Noise.** Fresh noise uses the source's errors. The PSF-convolved size model is approximate. The DJA
  same_galaxy radius is calibrated on CANUCS (D-051).
- **Single deflector redshift** z_l = 0.4 and single source redshift z_s = 2 for the masses. A real population
  spread in z_l would change the mass scale, not the angular efficiency.
- **No macro-model.** CANUCS flanking-field μ = 1.0–1.4 is ignored for the injected lens.
- **Background-aware limit.** It treats null (e) as exact. Null (e) itself under-predicts the flanking-field
  orphans by ~15 % (orphan_pairs.md), so the true background may be higher and the limit is then conservative.
