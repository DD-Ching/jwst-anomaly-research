# W3 / Kepler: physical mapping of the tested grid and the value of a larger search — one expected detection needs a negative-mass density ≥ 0.87 M☉ pc⁻³ (about 10× the local baryonic density); recommendation C, pause

- **Plan:** `PLAN.md` was committed at 10:43:25 UTC (a7c8afe), before the numbers. Downloads: 11,118,549 B of the
  150 MB cap (`sources.json`).
- **Inputs:**
  - the DR25 radius and distance of the 2,000 searched stars (official `keplerstellar` table, model values;
    `inputs/stellar_dr25.csv.gz`). 1,971 have both, with errors; 29 lack a value and are left out;
  - the measured efficiencies of `s2000/combined.json` (`round_plus_546`).
- **No:** new light curve, image, injection or threshold change.
- **Labels:** catalogue radii and distances are model values. The mapping and the population requirement are
  derived. Scenario speeds are assumptions. Nothing here is an abundance, a physical event-rate limit or a
  discovery probability.
- **Deviations:** none. Notes:
  - The QMASK comment in `scripts/w3_kepler.py` is corrected (section 4). The mask itself is unchanged.

## 0. What already existed (checked first)
- **Model:** `exotic_sim` (D-047):
  - θ_E = √(4G|M|/c² · D_LS/(D_L D_S)) (`einstein_radius`);
  - β(t) = √(u0² + ((t − t0)/t_E)²) (`impact_track`);
  - ρ = θ_*/θ_E (`finite_source_magnification`);
  - for n = 1, sign −1, the umbra is β < 2 (`caustic_beta(1)` = 2).
  - The test `test_cell_map_is_consistent_with_the_model_einstein_radius` feeds this round's |M| and D_L back
    into `einstein_radius` and reproduces θ_E to 0.2 %.
- **Literature already read in full (SOURCES "Exotic-lensing predictions"):** Safonova, Torres & Romero 2002 (PRD
  65, 023001) and Cramer et al. 1995 for the light-curve shape.
  - Their abstracts (re-read now through the arXiv API, ledger) give no population or rate model.
- **Takahashi & Asada 2013 (ApJL 768, L16):** SDSS quasar-lens bounds for negative masses at |M| > 10¹² M☉
  (galaxy scale) and for Ellis wormholes with throat radius 10–10⁴ pc. This is not the mass scale tested here, and
  it is not transferred.
- **Torres, Romero & Anchordoqui 1998 (PRD 58, 123001):** a cosmological bound of O(10⁻³⁶) g cm⁻³ on the
  negative-mass density, from GRB data. That is about 1.5 × 10⁻¹⁴ M☉ pc⁻³ (conversion not stored).
  - Its mass range is not in the abstract, and applying it to the local Galaxy would need a clustering model, so it
    is not used.
- **What this tests:** the W3 grid tests only the n = 1 negative-mass point lens. Ellis wormholes (n = 2) and
  general wormholes are other models and are not tested.
- **Prior work:** no earlier physical mapping of the W3 light-curve searches exists in the repository. The OGLE and
  MOA limits (docs/exotic_limits.md) stop at events per star-year.

## 1. What each tested cell means physically (`mapping.json`)
- **Identity:** ρ · t_E = θ_*/μ_rel = t_*, the source-radius crossing time.
  - For a lens at D_L = x D_S moving at transverse speed v: t_* = x R_*/v, so x = v t_*/R_*.
  - Then R_E = v t_E and |M| = R_E² c² / (4G D_L (1 − x)).
  - So ρ and t_E are not free once R_*, D_S and v are given: each cell fixes x for a given speed.
- **Stars:** median R_* 1.08 R☉ (5–95 %: 0.7265–1.64); median D_S 587.9 pc (265.6–918.9).
- **Table** (blend 1.0; the blend 0.7 cells map identically). Values are medians over stars, with the 5–95 % range
  in brackets. "Reachable" is the share of stars for which 0 < x < 1.

| cell | ρ | t_E (d) | t_* (d) | θ_E (µas) | μ_rel (mas/yr) | 220 km/s: reachable | D_L (pc) | \|M\| (M☉) | 30 km/s: reachable | D_L (pc) | \|M\| (M☉) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.01 | 0.1 | 0.001 | 876 | 3,200 | 1.0 | 14 (8.7–20) | 0.0014 (0.001–0.0024) | 1.0 | 2 (1.2–2.7) | 0.00019 (0.00014–0.00031) |
| 1 | 0.01 | 0.3 | 0.003 | 876 | 1,070 | 1.0 | 44 (26–59) | 0.0044 (0.0033–0.0076) | 1.0 | 5.9 (3.5–8.1) | 0.00057 (0.00042–0.00095) |
| 2 | 0.01 | 1 | 0.01 | 876 | 320 | 1.0 | 140 (87–200) | 0.018 (0.013–0.033) | 1.0 | 20 (12–27) | 0.0019 (0.0014–0.0033) |
| 3 | 0.01 | 3 | 0.03 | 876 | 107 | 0.863 | 450 (270–590) | 0.15 (0.078–1.2) | 1.0 | 59 (35–80) | 0.0062 (0.0046–0.011) |
| 4 | 0.01 | 10 | 0.1 | 876 | 32 | 0.0 | — | — | 1.0 | 200 (120–270) | 0.029 (0.02–0.056) |
| 5 | 0.1 | 0.1 | 0.01 | 87.6 | 320 | 1.0 | 140 (87–200) | 0.00018 (0.00013–0.00033) | 1.0 | 20 (12–27) | 1.9 × 10⁻⁵ (1.4–3.3 × 10⁻⁵) |
| 6 | 0.1 | 0.3 | 0.03 | 87.6 | 107 | 0.863 | 450 (270–590) | 0.0015 (0.00078–0.012) | 1.0 | 59 (35–80) | 6.2 × 10⁻⁵ (4.6–11 × 10⁻⁵) |
| 7 | 0.1 | 1 | 0.1 | 87.6 | 32 | 0.0 | — | — | 1.0 | 200 (120–270) | 0.00028 (0.0002–0.00056) |
| 8 | 0.1 | 3 | 0.3 | 87.6 | 10.7 | 0.0 | — | — | 0.4444 | 600 (370–810) | 0.0037 (0.0016–0.035) |
| 9 | 0.1 | 10 | 1 | 87.6 | 3.2 | 0.0 | — | — | 0.0 | — | — |

- **Reading:**
  - The tested grid corresponds to negative masses of about 2 × 10⁻⁵ to 0.15 M☉ (medians; about 6 Earth masses to
    160 Jupiter masses) at distances of a few pc to about 600 pc.
  - With halo-like speeds (220 km/s), the cells with t_* ≥ 0.1 d (cells 4, 7, 8, 9) are not reachable.
  - With disk-like speeds (30 km/s), only cell 9 is not reachable, and cell 8 only for 44 % of the stars.
  - Small t_* (ρ 0.01 with t_E ≤ 0.3 d) needs very nearby lenses: 2–59 pc.
  - The catalogue errors move θ_* and therefore |M| (`M_msun_median_theta_star_high/low`). The scenario speeds are
    assumptions, and v moves x linearly.
- **No events between grid points:** the 20 cells are measured points. No efficiency is claimed between or
  outside them.

## 2. Would any population give an expected signal? (`population.json`)
- **Assumptions:**
  - a uniform lens mass density |ρ_neg| along each sight line;
  - umbra crossings (u0 < 2);
  - all events have one tested t_E and a ρ in 0.01–0.1;
  - efficiency between the lowest and highest measured cell at that t_E.
- **Expected detections:** N = Σ_i ε (4/π) τ_i T_i / t_E, with τ_i = (4πG/c²) |ρ_neg| D_S,i²/6 the Einstein-radius
  optical depth.
  - This does not depend on |M| or v once t_E is fixed. Σ D_S² T = 3.089 × 10⁹ pc² yr for the 1,971 stars.
- **Density needed, in M☉ pc⁻³:**

| t_E (d) | ε (lowest–highest cell) | N = 1, now (2,000) | N = 3, now | N = 3, +500 | N = 1, +5,000 | N = 3, +5,000 |
|---|---|---|---|---|---|---|
| 0.1 | 0.7485–0.7985 | 0.8696–0.9277 | 2.605–2.779 | 2.084–2.223 | 0.2485–0.2651 | 0.7444–0.7941 |
| 0.3 | 0.7573–0.796 | 2.617–2.751 | 7.84–8.241 | 6.272–6.593 | 0.7477–0.7859 | 2.24–2.355 |
| 1 | 0.734–0.8103 | 8.569–9.46 | 25.67–28.34 | 20.54–22.67 | 2.448–2.703 | 7.335–8.098 |
| 3 | 0.7128–0.794 | 26.24–29.22 | 78.6–87.56 | 62.88–70.04 | 7.496–8.35 | 22.46–25.02 |
| 10 | 0.6205–0.6885 | 100.9–111.9 | 302.2–335.3 | 241.7–268.2 | 28.81–31.97 | 86.33–95.79 |

- **Scale references, not population models (Read 2014):**
  - local dark-matter density extrapolated from the rotation curve: 0.005–0.015 M☉ pc⁻³;
  - local baryonic density: 0.0914 ± 0.009 M☉ pc⁻³ (Flynn et al. 2006).
- **How many more stars would be needed** (`star_factor_needed_vs_now`), even at the best t_E (0.1 d):
  - N = 1 at the baryonic density needs 9.514× the present 2,000 stars;
  - N = 1 at 0.015 M☉ pc⁻³ needs 57.97×;
  - N = 3 needs 28.54× and 173.9×.
  - The remaining unused FGK pool is 17,740 stars, so exhausting it would give about 9.9× (sum not stored).
- **Reading:**
  - For one expected detection, this search needs a negative-mass density at least about 9.5× the local baryonic
    density and at least about 58× the highest dark-matter reference.
  - A negative-mass component that large would need an equally large additional positive density to keep the
    measured local dynamical density. Read 2014's dark-matter estimates (0.005–0.015) show no such component.
  - This is reasoning from a review's numbers, not a dynamical model built here.
- **Not modelled:**
  - spatial clustering;
  - mass or t_E distributions (events spread over t_E are bounded by the per-t_E rows);
  - t_E < 0.1 d: a higher rate per density (∝ 1/t_E), but efficiency is not measured there, and Kepler's 29.4-min
    cadence resolves only a few cadences in the umbra.

## 3. Sky to SAP: what is not verified (Kepler Archive Manual, KDMC-10008-006, Table 2-3; sha256 in `sources.json`)
1. **Cosmic-ray correction (bit 128).** "Cosmic Ray was found and corrected in optimal aperture pixel" is set
   from the data, and `QMASK` excludes these cadences.
   - A real caustic spike is a 1–2 cadence brightening, and it might be corrected and flagged. That would remove
     shoulder cadences that the local shoulder test needs.
   - An injection into SAP cannot reproduce this. The pipeline algorithm (Kepler Data Processing Handbook) was not
     read here.
2. **SPSD (bit 1024).** It is set by PDC "on the last non-gapped cadence before the maximum positive change due
   to the detected SPSD". A sharp umbra ingress might trigger it, which would mask one cadence. Not verified.
3. **Unmasked bits.** Manual exclude (256), impulsive outlier (2048), mod/out Argabrightening (4096), collateral
   cosmic ray (8192) and not-in-fine-point (32768) cadences are not masked.
   - The #500 plan and the `QMASK` comment said manual exclude and impulsive outliers were masked. Those labels
     followed a shifted bit table; the comment is corrected.
   - The mask itself, and therefore every search and efficiency result since #500, is unchanged.
   - The labels matter for ordinary false-positive control, not for W3 recovery: the injections ran through the
     same mask.
4. **Aperture photometry.** Crowding and aperture changes between quarters are modelled only by a constant blend.
   Limb darkening is not modelled (uniform disk).
- **Effect on the conclusion:** every listed step can remove or dilute a real signal, and none adds one. If they
  matter, real efficiencies are lower, and the densities in section 2 are lower bounds on what is needed.
  - This is reasoning, not measured.
  - The product-level conclusion stands. No abundance or physical rate limit is quoted.
- **Three separate questions:** (i) the model produces this light-curve shape (Safonova et al. 2002); (ii) data
  support the model (no: zero events); (iii) a controllable signal or communication channel (not addressed by
  any of this). None implies the next.

## 4. Recommendation: C, pause the W3 / Kepler branch
- **What binds is physics, not the pipeline.** With 2,000 stars, one expected detection needs |ρ_neg| ≥ 0.87 M☉
  pc⁻³ at the best tested t_E (section 2).
  - +500 stars lowers this by a factor 0.8.
  - +5,000 lowers it by a factor 0.286, to 0.25 M☉ pc⁻³. That is still about 2.7× the local baryonic density (not
    stored).
  - Even the whole remaining FGK pool would only reach the baryonic density for N = 1.
- **Why not B (fill a gap first):** the open steps (section 3) can only lower efficiency. Settling them would not
  move the requirement towards plausible densities, so it would not change this judgement.
  - The one documentation error found (QMASK labels) is corrected here and changes no result.
- **What is kept:**
  - the zero result in 2,000 FGK stars;
  - the measured product-level sensitivity (#546, #554);
  - this mapping.
  - They remain valid as an exploratory, model-light check for 1-mag W3-shaped dimming.
  - "No unique discovery probability" does not mean "no exploration value". It means the value per added star is
    linear and small against the physical gap.
- **Not done:** no other branch is opened or switched to by this recommendation.

## 5. Plain answers
1. **What we are looking for:** a compact object with negative mass (the n = 1 negative-mass point lens) passing
   in front of a Kepler star. Its "umbra" hides the star completely for about 4 t_E, with bright spikes at both
   edges.
   - The tested grid means lenses of about 6 Earth masses to 0.15 solar masses at a few to about 600 pc, moving at
     disk-like or halo-like speeds.
   - It is not an Ellis wormhole or wormholes in general.
2. **What the data can reach:** the efficiency is measured at 20 points of that space. With halo-like speeds,
   four cells (t_* ≥ 0.1 d) cannot occur at all.
   - The physically interesting quantity is the density of such objects.
   - For even one expected event, they would have to outweigh all the ordinary matter near the Sun about tenfold
     (0.87–0.93 M☉ pc⁻³ against 0.0914). That is far beyond any dark-matter estimate.
3. **What another batch adds:**
   - 500 stars cut the required density by 20 %;
   - 5,000 stars cut it to about 29 %.
   - Neither reaches densities that the local mass budget allows, so the main gain would be a smaller
     product-level number.
4. **Recommendation:** pause this branch (C). The limit comes from how rarely such lenses can line up with stars
   at plausible densities, not from the pipeline. The open pipeline questions could only make it worse. Keep the
   results; no other branch is opened automatically.
