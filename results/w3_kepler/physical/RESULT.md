# W3 / Kepler: physical mapping of the tested grid and the value of a larger search — for a uniform negative-mass population with t_E ≥ 0.1 d, one expected detection needs |ρ_neg| ≥ 0.69 M☉ pc⁻³ even if every umbra crossing were detected (about 7.6× the local baryonic density); recommendation C, pause

- **Plan:** `PLAN.md` was committed at 10:43:25 UTC (a7c8afe), before the numbers.
- **Downloads:** 11,118,549 B of the 150 MB cap (`sources.json`).
- **Inputs:**
  - the DR25 radius and distance of the 2,000 searched stars (official `keplerstellar` table, model values;
    `inputs/stellar_dr25.csv.gz`). 1,971 have both, with errors; 29 lack a value and are left out;
  - the measured efficiencies in `s2000/combined.json` (`round_plus_546`).
- **No:** new light curve, image, injection or threshold change.
- **Labels:** catalogue radii and distances are model values. The mapping and the population requirement are
  derived. Scenario speeds are assumptions. Nothing here is an abundance, a physical event-rate limit or a
  discovery probability.
- **Deviation:** scenario S2 (30 km/s, disk-like) is an unsourced assumption. The owner asked for sourced
  scenarios, and S1's 220 km/s is from Read 2014. S2 only moves the mass and distance mapping and the
  restricted-ρ rows, not the ε = 1 floor.
- **Notes:**
  - The `QMASK` comment in `scripts/w3_kepler.py` is corrected (section 3); the mask is unchanged.
  - After the independent review, `population.json` gained the ε = 1 floor, the tested-ρ restriction, the
    unread-pool count and the short-t_E numbers. One test was made independent of the code's constants (astropy and
    a numerical integral).
  - After `/code-review`, the θ_* error bracket in `mapping.json` was fixed: DR25 stores `dist_err2` positive, so
    the old "D low" corner moved D up. Only the two bracket keys changed; the medians and ranges in the table did
    not.

## 0. What already existed (checked first)
- **Model:** `exotic_sim` (D-047):
  - θ_E = √(4G|M|/c² · D_LS/(D_L D_S)) (`einstein_radius`);
  - β(t) = √(u0² + ((t − t0)/t_E)²) (`impact_track`);
  - ρ = θ_*/θ_E (`finite_source_magnification`);
  - for n = 1, sign −1, the umbra is β < 2 (`caustic_beta(1)` = 2).
  - The test `test_cell_map_is_consistent_with_the_model_einstein_radius` feeds this round's |M| and D_L back
    into `einstein_radius` and reproduces θ_E to 0.2 %.
- **Literature, and how much of it was read:**
  - **Safonova, Torres & Romero 2002** (PRD 65, 023001): read in full on 2026-10-08 (SOURCES "full text") for
    the umbra, the caustic at 2θ_E and A(u). Its abstract, re-read now, gives no population or rate model.
  - **Cramer et al. 1995** (PRD 51, 3117): abstract only (SOURCES).
  - **Takahashi & Asada 2013** (ApJL 768, L16): abstract only, through the arXiv API in this round. It bounds
    negative masses with |M| > 10¹² M☉ (galaxy scale) and Ellis wormholes with throat radii 10–10⁴ pc, from SDSS
    quasar lenses. That is not the mass scale tested here, and it is not transferred (as in `docs/exotic_limits.md`).
  - **Torres, Romero & Anchordoqui 1998** (PRD 58, 123001): abstract only. It gives a cosmological bound of
    O(10⁻³⁶) g cm⁻³ on the negative-mass density from GRB data, about 1.5 × 10⁻¹⁴ M☉ pc⁻³ (conversion not stored).
    Its mass range is not in the abstract, and applying it locally needs a clustering model, so it is not used.
- **What this tests:** the W3 grid tests only the n = 1 negative-mass point lens. Ellis wormholes (n = 2) and
  general wormholes are other models and are not tested.
- **Prior work:** no earlier physical mapping of the W3 light-curve searches exists. The OGLE and MOA limits stop
  at events per star-year.

## 1. What each tested cell means physically (`mapping.json`)
- **Identity:** ρ · t_E = θ_*/μ_rel = t_*, the source-radius crossing time.
  - For a lens at D_L = x D_S moving at transverse speed v: t_* = x R_*/v, so x = v t_*/R_*.
  - Then R_E = v t_E and |M| = R_E² c² / (4G D_L (1 − x)).
  - So ρ and t_E are not free once R_*, D_S and v are given: each cell fixes x for a given speed.
- **Stars:** median R_* 1.08 R☉ (5–95 %: 0.7265–1.64); median D_S 587.9 pc (265.6–918.9).
- **Table** (blend 1.0; the blend 0.7 cells map identically). Values are medians over stars, with the 5–95 % range
  in brackets. "Reachable" is the share of stars with 0 < x < 1.
  - In partly reachable cells (3 and 6 at 220 km/s, 8 at 30 km/s), D_L and |M| cover only the reachable stars.

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
    160 Jupiter masses, conversion not stored), at a few pc to about 600 pc.
  - At exactly 220 km/s, cells 4, 7, 8 and 9 (t_* ≥ 0.1 d) are not reachable for any star. With both blends that
    is 8 of the 20 measured cells.
  - At 30 km/s, only cell 9 is not reachable, and cell 8 only for 44 % of the stars.
  - Small t_* (ρ 0.01, t_E ≤ 0.3 d) needs nearby lenses: medians of 2–44 pc (5–95 %: 1.2–59 pc).
  - The catalogue errors move θ_* and therefore |M| (`M_msun_p5_50_95_theta_star_high/low`). For cell 2 at
    220 km/s, the median |M| of 0.01828 M☉ becomes 0.01287 (θ_* low) to 0.02641 (θ_* high). The bracket ignores the
    correlation between the DR25 R and D errors, so it is conservative.
- **No efficiency between grid points:** the 20 cells are measured points.

## 2. Would any population give an expected signal? (`population.json`)
- **Expected detections:** N = Σ_i ε (4/π) τ_i T_i / t_E, for umbra crossings (u0 < 2, the injected
  distribution). τ_i = (4πG/c²) |ρ_neg| D_S,i²/6 is the Einstein-radius optical depth of a uniform lens mass density
  along each sight line. Σ D_S² T = 3.089 × 10⁹ pc² yr for the 1,971 stars.
- **The catch:** for one speed, ρ = x R_*/(v t_E) varies along the sight line. Most lenses have ρ outside the
  tested 0.01–0.1, where ε is not measured. At t_E = 0.1 d only 13.91 % of τ (220 km/s) and 0.3052 % (30 km/s)
  comes from lenses inside the tested ρ range (`tested_rho_only.share_of_tau`).
- **Three bounds are reported** (density in M☉ pc⁻³ for N = 1). The sums run over the 1,971 stars with a distance,
  unscaled; scaling to all 2,000 would lower every number by about 1.5 % (not stored).
  - **(a)** the measured ε applied at every distance (an extrapolation);
  - **(b)** ε = 1 at every distance: a floor that holds whatever the efficiency is outside the tested range;
  - **(c)** only lenses inside the tested ρ range, with the measured ε.

| t_E (d) | (a) measured ε, all lenses | (b) ε = 1 floor | (c) tested ρ only, 220 km/s | (c) tested ρ only, 30 km/s |
|---|---|---|---|---|
| 0.1 | 0.8696–0.9277 | 0.6944 | 6.252–6.67 | 284.9–303.9 |
| 0.3 | 2.617–2.751 | 2.083 | 3.576–3.759 | 100.1–105.2 |
| 1 | 8.569–9.46 | 6.944 | 9.973–11.01 | 35.85–39.58 |
| 3 | 26.24–29.22 | 20.83 | 103.4–115.2 | 29.06–32.38 |
| 10 | 100.9–111.9 | 69.44 | none (no lens in range) | 133.1–147.6 |

- **Scale references, not population models (Read 2014):**
  - local dark-matter density extrapolated from the rotation curve: 0.005–0.015 M☉ pc⁻³;
  - local baryonic density: 0.0914 ± 0.009 M☉ pc⁻³ (Flynn et al. 2006).
- **More stars:**
  - **+500 / +5,000:** under row (a), assuming a zero result and similar efficiency, the requirement scales by 0.8 and
    by 0.286 (`plus_500`, `plus_5000` in `population.json`).
  - **The whole unread FGK pool:** 16,240 stars after #554 (`pool_unread_after_554`), 9.12× the 2,000 (not stored).
    One expected detection at t_E 0.1 d would still need 0.07614 M☉ pc⁻³ under the floor (b). Under (a) it would
    need 0.09535, which is above the baryonic density.
  - Either way it is at least 5× the highest dark-matter reference.
- **Shorter t_E, the one direction that helps:** the rate per density grows as 1/t_E.
  - With ε = 1 and these stars, one expected detection at the baryonic density needs t_E ≈ 0.01316 d (19 min). Its
    umbra lasts at most 2.58 long cadences.
  - At 0.015 M☉ pc⁻³ it needs t_E ≈ 0.00216 d, with an umbra shorter than one cadence (0.423)
    (`eps1_tE_needed_for_N1_at_reference`).
  - The efficiency below t_E = 0.1 d is not measured. An umbra of 0.423 cadence cannot give three consecutive
    1-mag cadences. At 2.58 cadences it could only at a favourable phase and small u0, and the local shoulder test
    also needs 2 cadences on each side (`LOCAL_NCAD`).
- **Reading, as reasoning and not a dynamical model built here:**
  - A negative-mass density of 0.69 M☉ pc⁻³ or more would have to be cancelled by an equally large unseen positive
    density to keep the local dynamical density.
  - Measurements of the net density (Read 2014) cannot by themselves exclude such a cancelling pair. But it would be
    a fine-tuned positive component at least 46× the highest dark-matter reference (0.6944 / 0.015; not stored).
- **Not modelled:** spatial clustering; mass and speed distributions (a single speed per row (c)).

## 3. Sky to SAP: what is not verified (Kepler Archive Manual, KDMC-10008-006, Table 2-3; copied into `sources.json`)
1. **Cosmic-ray correction (bit 128).** "Cosmic Ray was found and corrected in optimal aperture pixel" is set
   from the data, and `QMASK` excludes these cadences.
   - A real caustic spike (a 1–2 cadence brightening) might be corrected and flagged, which would remove shoulder
     cadences that the local shoulder test needs.
   - An injection into SAP cannot reproduce this. The pipeline algorithm (Kepler Data Processing Handbook) was not
     read.
2. **SPSD (bit 1024).** "SPSD detected (last non-gapped cadence before the maximum positive change)", set by
   PDC. A sharp umbra edge might trigger it and mask one cadence. Not verified.
3. **Unmasked bits:** manual exclude (256), impulsive outlier (2048), mod/out Argabrightening (4096), collateral
   cosmic ray (8192) and not-in-fine-point (32768) cadences are not masked.
   - The #500 plan listed the masked bits by number but mislabelled kept bits: 2048 and 4096 as "straylight", 8192
     as "planet-search exclude".
   - The old code comment labelled masked bits as aperture cosmic ray (64), manual exclude (128), impulsive outlier
     (512) and so on. Those labels are corrected now; the old ones resemble the TESS flag table.
   - The mask, and therefore every result since #500, is unchanged. The labels matter for ordinary false-positive
     control. W3 recovery is not affected, because the injections ran through the same mask on the archive's real
     flags.
4. **Aperture photometry and source model:**
   - crowding and per-quarter aperture changes are modelled only by a constant blend;
   - the source is a uniform disk without limb darkening.
- **Effect on the conclusion:**
  - For bits 128 and 1024 and for crowding (blend 1 is already the top of the bracket), a real signal is plausibly
    removed or diluted, not added. This is reasoning, not measured.
  - For limb darkening, the sign of its effect on the local shoulder test is not known.
  - Hence the robust statement uses floor (b), ε = 1, which no efficiency question can undercut. The conclusion
    stays at product level.
- **Three separate questions:** (i) the model produces this light-curve shape (Safonova et al. 2002); (ii) data
  support the model (no: zero events); (iii) a controllable signal or communication channel (not addressed here).
  None implies the next.

## 4. Recommendation: C, pause the W3 / Kepler branch
- **What binds is the event rate per density, not the pipeline.** Even if every umbra crossing along the sight
  lines were detected (ε = 1), one expected detection in these 2,000 stars needs |ρ_neg| ≥ 0.6944 M☉ pc⁻³ at the
  best tested t_E (0.1 d).
  - That is about 7.6× the local baryonic density (not stored).
  - Restricted to the tested ρ range with the measured ε, it needs 3.576 or more.
- **Why not A (expand):**
  - +500 stars multiplies the requirement by 0.8, and +5,000 by 0.286.
  - Even the whole unread pool (9.12×, not stored) stays at the baryonic scale (0.07614–0.09535, against 0.0914) and
    at least 5× the dark-matter references.
- **Why not B (fill a gap first):**
  - The efficiency questions (section 3) cannot move floor (b).
  - The only direction that lowers the requirement is shorter t_E. Reaching dark-matter-scale densities needs
    t_E ≈ 3 min, an umbra shorter than one long cadence, so no existing Kepler long-cadence step resolves it.
  - The documentation error found (QMASK labels) is corrected here and changes no result.
- **What is kept:**
  - the zero result in 2,000 FGK stars;
  - the measured product-level sensitivity (#546, #554);
  - this mapping.
  - They stay valid as an exploratory, model-light check for 1-mag W3-shaped dimming.
  - "No unique discovery probability" is not "no exploration value". Here, the value of each added star is linear
    and small against the physical gap.
- **Not done:** no other branch is opened or switched to by this recommendation.

## 5. Plain answers
1. **What we are looking for:** a compact object with negative mass (the n = 1 negative-mass point lens) passing
   in front of a Kepler star. Its "umbra" hides the star completely for up to 4 t_E (at u0 = 0), with bright
   spikes at both edges.
   - The tested grid means lenses of about 6 Earth masses to 0.15 solar masses at a few to about 600 pc, at
     disk-like or halo-like speeds.
   - It is not an Ellis wormhole or wormholes in general.
2. **What the data can reach:**
   - Efficiency is measured at 20 points. At exactly 220 km/s, 8 of them cannot correspond to any real lens
     (slower lenses can reach some of them).
   - The physically decisive quantity is how many such objects there are.
   - Even if every passage were caught, one expected event in these 2,000 stars needs the negative-mass objects to
     outweigh all ordinary matter near the Sun about 7.6-fold, assuming a uniform population and t_E ≥ 0.1 d.
   - That is far beyond any dark-matter estimate.
3. **What another batch adds:**
   - 500 stars lower the required density by 20 %, and 5,000 to about 29 %.
   - Even the whole unread pool stays at the scale of all ordinary matter.
   - So the main gain would be a smaller product-level number.
4. **Recommendation:** pause this branch (C). The limit comes from how rarely such lenses can line up with stars
   at plausible densities. It does not come from the pipeline, and the open pipeline questions cannot remove it.
   Keep the results; no other branch is opened automatically.
