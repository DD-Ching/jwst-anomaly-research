# B on A1: destruction report for the "Relational Transition World" (RTW)

System B ("Physics Destruction"). Target: `A1-relational.md`. Labels: **[lit]** = from a reference I opened
(list at the end); **[B-est]** = my own order-of-magnitude estimate; **[std]** = textbook physics.

## 1. Internal consistency

1. **Ω cannot be countable.** A6(i) needs infinitely many loci (polynomial growth "at large n"). But a product
   Π_{x∈V} Σ_x over infinite V, with |Σ_x| ≥ 2, is uncountable. So A1 contradicts A4 unless Ω is a restricted
   (finite-support) product. A1 never says this, and the restriction would bring back a preferred "vacuum" state.
2. **The finest factorization is not defined.** Product decompositions of a set do not form a lattice that is closed
   under refinement, so a "finest" one satisfying the k-locality condition need not exist or be unique. A8 then changes
   the factorization over time, while A4 defines it once for all of Ω. Loci are therefore not well defined.
3. **The redshift mechanism gives z = 0.** A1 says a carrier's wavelength "counted in loci, is fixed in transit while
   loci split". Emitter atoms, observer atoms and the carrier would all keep a fixed number of loci, so
   λ_obs/λ_emit = 1. Cosmological redshift needs the carrier to span *more* loci over time, i.e. splits inside the
   carrier. Splits are a-transitions (A8 ⊂ T_a), so they conserve h (A7). Then the carrier's h is unchanged while its
   frequency falls, so h is not photon energy, and A7's "energy-like conservation" constrains nothing observable.
4. **Γ_b has inconsistent units and meaning.** Γ_b ≡ εκ_b f_b is a rate per unit *time*, but A6(iii) ties leakage to
   "crossing" a b-endpoint, which is per unit *path length* (≈ εf_b/ℓ_a, with κ_b dropping out). Under the first
   reading, stationary matter leaks h at the same rate as light (§2). Under the second, D1's parameter is a different
   one from the parameter in D3/P2.
5. **The small-α branches have infinite degree.** The expected b-degree is Σ_n n²·p₀n^{−α}, which diverges for α ≤ 3.
   D2's "α < 3" branch and D1's "spread over the whole volume when α is small" therefore contradict "b-edges are
   sparse" (A6ii) in an infinite universe. In a finite one the degree grows as N^{(3−α)/3} as loci split, so Γ_b depends
   on time and the "linear in χ" form of η(z) fails. A1 also leaves undefined what happens to b-edges when a locus splits.
6. **Circularity (partly admitted).** 3-D geometry, isotropy, diffusive gravity and the factor 2 are imposed. Duality
   is "derived iff manifold-like", which is the assumption itself. Minor slips: the Lieb–Robinson v is missing the
   degree factor and is a bound, not the front speed 1/μ; §5b omits the (1+z) dilation of Δt_emit.
7. **Checked and correct:** Δ = log2/log(2d/α) [math/0304418]; 1 µs at 0.8 Gpc gives ℓ_a ≲ 1.0 nm [B-est].

## 2. Conservation, stability, causality

- **Momentum is not conserved.** A7 conserves only a scalar h. A b-transfer moves h to a distant locus with no rule for
  its momentum. If the direction is kept, every source gets displaced ghost images. If it is randomized, local momentum
  conservation is violated. There is no stress-energy tensor, so nothing guarantees local Lorentz covariance.
- **Preferred frame.** Global Markov time plus a bounded-degree G_a pick out a rest frame; Lorentz-invariant
  discreteness needs Poisson-type randomness with no equivariant direction [gr-qc/0605006]. The global foliation forbids
  causal loops, so b-hops signal superluminally without paradoxes, at the price of explicit LIV.
- **Generic dispersion.** Excitations hopping on a lattice of step ℓ_a generically get ω(k) corrections of order
  (kℓ_a)². A1 asserts that every effect is "achromatic" but derives no dispersion relation. From LHAASO GRB 221009A,
  E_QG,2 > 6×10⁻⁸ E_Pl and E_QG,1 > 10 E_Pl [lit 2402.06009]. With E_QG,2 ≈ ħc/ℓ_a this gives
  **ℓ_a ≲ 2.7×10⁻²⁸ m** [B-est]. Fermi GRB 090510 gave the earlier linear bound of about E_Pl [0908.1832].
- **Matter leakage (first reading of §1.4).** If bound matter (mostly relativistic constituent energy) also
  leaks, LLR's Ġ/G = (7.1 ± 7.6)×10⁻¹⁴ yr⁻¹ (quoted in [1812.11181]) caps it near 10⁻¹³ yr⁻¹, forcing
  η − 1 ≲ 8×10⁻⁴ at z = 1 [B-est, conditional].
- **Gravity.** Without the factor 2, light bending is half the GR value. Cassini: γ = 1 + (2.1 ± 2.3)×10⁻⁵ (Bertotti+ 2003).

## 3. Already excluded? (per prediction)

**P1 (achromatic leakage plus diffuse excess).** Cosmic-opacity bounds: Δτ < 0.012 (95%) between z = 0.2 and 0.35,
and ε = −0.04 (+0.08/−0.07) in d_L = d_A(1+z)^{2+ε} [1004.2053]. In RTW variables, τ = Γ_bχ/c gives
**Γ_b ≲ 7×10⁻¹² yr⁻¹, so τ(z=1) ≲ 0.07** [B-est]. Recent DESI+SN analyses disagree with each other: one reports a 2σ
violation (6σ with SH0ES+BBN calibration), but no deviation in its model-independent reconstruction [2504.01750], and
others find consistency. These tensions depend on calibration, so the SN Ia absolute-magnitude or evolution systematics
are an ordinary mimic for any η ≠ 1. The diffuse side: New Horizons measures a COB of 11.16 ± 1.65 against an IGL of
8.17 ± 1.18 nW m⁻² sr⁻¹, an excess of 2.99 ± 2.03 that is not significant [2407.06273]. Galaxy-count EBL is in
[1605.01523]. The allowed excess, ≲ 0.07 × 8 ≈ 0.6 nW m⁻² sr⁻¹, is about 0.3σ of the current excess error. **The
consistency relation cannot fail with present data.** The idea is known: photon-number violation as a distance-duality
probe [astro-ph/0312443] and cosmic-opacity tests [1004.2053].

**P2 (pre-leading-image hybrid arrivals).** I found no published search at the lags s_i − s_j. The geometry is right:
the lag relative to image j is (s_i + o_j) − (s_j + o_j), and a discontinuous path can beat the global Fermat minimum.
**But the amplitude cannot be reached under A1's own scaling** [B-est]. The total leaked fraction along a path is
≤ 0.07, and the part leaked within ~100 kpc of a lens at ~Gpc is ≲ 10⁻⁵. Its partner must land inside tube j, whose
cross-section at the lens plane is about (source size)² ~ (10¹⁵ m)². The routing probability is ~10⁻¹⁰ for large α
and ~10⁻²⁸ for α < 3, so a_ij ≪ 10⁻¹⁰. For comparison, COSMOGRAIL light curves (15 yr, 23 systems)
[2002.05736] reach about a_ij ~ 10⁻² per system [B-est]. SN Refsdal pre-explosion imaging [1411.6009, 1512.04654]
reaches a few ×10⁻². The only route to a detectable signal is structured, non-random b-edges, which is a new
unconstrained assumption. A1's claim that "microlensing is chromatic" is only weakly true; microlensing is
achromatic for compact sources.

**P3 (ghost gravity).** Large α: the Bullet Cluster has its total-mass centroid offset from the baryonic (mostly gas)
peak at 8σ [astro-ph/0608407]. Gas dominates the baryons, so a depletion-weighted centroid sits near the gas, which
this excludes. The 72-collision sample also finds dark matter tracking the galaxies [1503.07675]. Small α: the
component is smooth and does not cluster, so it cannot form halos. That branch is essentially the
"disordered-locality dark energy" proposal [0903.5303], and both branches resemble nonlocal gravity
[0706.2151, 0812.1059]. **Excluded or known.**

**P4 (achromatic broadening ∝ L^{1/3}).** RTW's δl = ℓ_a^{2/3}L^{1/3} is exactly the "α = 2/3" spacetime-foam
scaling with ℓ_a in place of l_P. Perlman et al. claim γ-ray image data exclude α = 2/3 even at the Planck length
[1411.7262], though that method is debated. My independent estimate [B-est]: with KPZ 1+1 exponents, the difference
in arrival time across a baseline b saturates at ~√(ℓ_a b)/c, independent of L. Keeping interferometric coherence
needs ℓ_a ≲ λ²/(4π²b). That is 2.6×10⁻¹⁵ m for the diffraction-limited HST (0.5 µm, 2.4 m), and 4×10⁻¹⁵ m for EHT
fringes at 230 GHz on ~10⁷ m baselines, including on 3C 279 at z = 0.54 (Kim et al. 2020, A&A 640, A69,
doi:10.1051/0004-6361/202037493). Then σ_T(1 Gpc) ≲ 2×10⁻¹⁰ s, and only 4×10⁻¹⁹ s under the LIV bound. FRB widths
go down to ~60 ns [2105.11446] and µs bursts [2307.02303, 2002.12539], so the FRB channel is 10²–10¹¹ times too
coarse. Intrinsic emission timescales and instrument resolution mimic any "floor". **Excluded or untestable.**

**P5 (direction-dependent speed).** If light and gravity both propagate on G_a, they share μ(n̂). GW170817's bound,
−3×10⁻¹⁵ < Δc/c < 7×10⁻¹⁶ [1710.05834], then constrains only a *differential* anisotropy. **A1's claim that it
implies δ_μ ≲ 10⁻¹⁵ is wrong inside its own model.** If the effect is universal and ellipsoidal it is a coordinate
choice. If it differs between light and matter, rotating optical resonators already bound it at 10⁻¹⁸
[1412.6954], about 10³ below anything astrophysical timing can reach. Excluded, or unobservable by construction.

**D3 (out-of-cone event dependence).** No direct limit exists, but it is unfalsifiable as stated (see §7).

## 4. Renaming check

Emergent-locality graph plus sparse links "identifying far-away points" = Markopoulou–Smolin **disordered locality**
[gr-qc/0702044] in **quantum graphity** [hep-th/0611197]. Its dark-energy use [0903.5303] = D2 at small α, and its
dispersion effects are modelled in [1201.3206]. Diamond-property causal order = causal invariance in the Wolfram model
[2004.14810]; see also causal sets [gr-qc/0605006]. D1 = cosmic opacity / achromatic photon loss
[astro-ph/0312443, 1004.2053]. D4 = α = 2/3 foam [1411.7262]. D2 at large α = nonlocal gravity [0812.1059]. New:
routing to a specific partner (P2) and the admissibility-margin bookkeeping.

## 5–6. Novelty and survivor triage

| Pred. | Novelty | Verdict | Reason |
|---|---|---|---|
| P1 | Known (opacity / photon-number violation) | **FAIL** | Γ_b ≲ 7×10⁻¹² yr⁻¹. The required diffuse excess (≲0.6 nW m⁻² sr⁻¹) is below the COB error (2.0). Untestable consistency relation; SN calibration mimics it. |
| P2 | Partially novel (precursor at lag s_i − s_j) | **Conditional PASS (limit only)** | Self-consistent and geometry-clean. The model-natural amplitude ≪10⁻¹⁰ vs ~10⁻² sensitivity, so it can only set a first limit. |
| P3 | Known (disordered-locality dark energy, nonlocal gravity) | **FAIL** | Bullet Cluster 8σ plus 72 collisions (large α). Small α has no halos and is not novel. |
| P4 | Known (α = 2/3 foam) | **FAIL** | ℓ_a ≲ 10⁻¹⁵ m from interferometric coherence (B-est) and ≲3×10⁻²⁸ m from LIV, so broadening ≪ ns. Intrinsic widths mimic it. |
| P5 | Partially known (SME anisotropy) | **FAIL** | Internal error (GW and light share μ); lab bound 10⁻¹⁸. |
| D3 | Known in spirit (non-local links) | **FAIL** | Not a prediction until τ_h and the edge statistics are fixed. A1 itself says ε can shrink without limit. |

**P2 mimics to exclude:** quasar red-noise self-correlation at any lag (null from off-model lags and
phase-randomized curves); achromatic microlensing; PSF crosstalk (at t_i − t_j, not s_i − s_j); lens-model systematics
(the source/observer split is unobservable and mass-sheet degeneracy rescales delays, so test both conventions with a
trials correction); host transients for SNe. Light echoes arrive later and cannot fake a precursor. **Data:**
COSMOGRAIL public light curves [2002.05736]; HST/JWST archival imaging of SN Refsdal [1411.6009, 1512.04654].
**Sensitivity:** about 10⁻² stacked, a few ×10⁻² for SN precursors. That gives a first limit on structured routing,
not A1's natural amplitude. Low priority, cheap null test.

## 7. Is the "blindness theorem" correct?

**No, as stated.** A1 requires Δχ to ~10 ly. But out-of-cone dependence only needs
|Δt_emit| < |Δr|/c, and |Δr| ≥ D_⊥ (angle × distance), which is known to about 1%. For a pair 1° apart at z = 0.5,
D_⊥ ≈ 33 Mpc ≈ 108 Mly [B-est]. Any dependence seen with |Δt_obs| < 10 yr has
|Δt_emit| ≤ 10 yr + σ_χ/c. Here σ_χ ≈ 1 Mly from σ_z = 10⁻⁴, or ≈ 14 Mly once 300 km/s peculiar velocities are
included. The pair is therefore spacelike with a margin of ~7–100σ. The admissible region *is* computable for every
pair with D_⊥ ≫ σ_χ. Only pairs closer than ~0.1° at Gpc distances are genuinely undetermined. A1's S ~ 10⁻⁵ is
correct, but it is an *efficiency*: only pairs with |Δχ| < c × baseline produce lags inside the survey window. It is
not an in-principle limit. The real obstacles are statistical: chance coincidences, trials factors, and common-mode
causes on the observer side (Earth, Sun, instrument), which are timelike to both detections and mimic correlation.

**Clean-geometry loopholes:** (a) same-source messengers: GW170817/GRB 170817A bounds Δv/c at ~10⁻¹⁵
[1710.05834], and a b-shortcut would show as a precursor; (b) line-of-sight precursors from one unlensed FRB, GRB or
SN, arriving early by (skipped length)/c, which needs no lens model; (c) lensed images (P2); (d) A1's binary-pulsar
and solar-flare tests. The "blindness theorem" holds only for close pairs and does not apply to events from one source.

## References (each opened and checked)

- astro-ph/0312443: Bassett & Kunz, "Cosmic distance-duality as probe of exotic physics and acceleration" (2003)
- 1004.2053: Avgoustidis et al., "Constraints on cosmic opacity and beyond the standard model physics from cosmological distance measurements" (2010)
- 2504.01750: Keil et al., "Probing the Distance Duality Relation with Machine Learning and Recent Data" (2025)
- 2407.06273: Postman et al., "New Synoptic Observations of the Cosmic Optical Background with New Horizons" (2024)
- 1605.01523: Driver et al., "Extra-galactic background light measurements from the far-UV to the far-IR from deep ground and space-based galaxy counts" (2016)
- 2002.05736: Millon et al., "COSMOGRAIL XIX: Time delays in 18 strongly lensed quasars from 15 years of optical monitoring" (2020)
- 1411.6009: Kelly et al., "Multiple Images of a Highly Magnified Supernova Formed by an Early-Type Cluster Galaxy Lens" (2014)
- 1512.04654: Kelly et al., "Deja Vu All Over Again: The Reappearance of Supernova Refsdal" (2015)
- 2105.11446: Nimmo et al., "Burst timescales and luminosities link young pulsars and fast radio bursts" (2021)
- 2307.02303: Snelders et al., "Detection of ultra-fast radio bursts from FRB 20121102A" (2023)
- 2002.12539: Cho et al., "Spectropolarimetric analysis of FRB 181112 at microsecond resolution" (2020)
- 1411.7262: Perlman et al., "New Constraints on Quantum Gravity from X-ray and Gamma-Ray Observations" (2014)
- 2402.06009: LHAASO Collab., "Stringent Tests of Lorentz Invariance Violation from LHAASO Observations of GRB 221009A" (2024)
- 0908.1832: Fermi GBM/LAT Collabs., "Testing Einstein's special relativity with Fermi's short hard gamma-ray burst GRB090510" (2009)
- 1710.05834: LIGO/Virgo/Fermi/INTEGRAL, "Gravitational Waves and Gamma-rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A" (2017)
- 1412.6954: Nagel et al., "Direct Terrestrial Test of Lorentz Symmetry in Electrodynamics to 10⁻¹⁸" (2014)
- 1812.11181: Belgacem et al., "Testing nonlocal gravity with Lunar Laser Ranging" (2018). The Ġ/G value is quoted in its text, not its abstract.
- astro-ph/0608407: Clowe et al., "A direct empirical proof of the existence of dark matter" (2006)
- 1503.07675: Harvey et al., "The non-gravitational interactions of dark matter in colliding galaxy clusters" (2015)
- gr-qc/0605006: Bombelli, Henson & Sorkin, "Discreteness without symmetry breaking: a theorem" (2006)
- gr-qc/0702044: Markopoulou & Smolin, "Disordered locality in loop quantum gravity states" (2007)
- hep-th/0611197: Konopka, Markopoulou & Smolin, "Quantum Graphity" (2006)
- 0903.5303: Prescod-Weinstein & Smolin, "Disordered Locality as an Explanation for the Dark Energy" (2009)
- 1201.3206: Caravelli & Markopoulou, "Disordered locality and Lorentz dispersion relations: an explicit model of quantum foam" (2012)
- 2004.14810: Gorard, "Some Relativistic and Gravitational Properties of the Wolfram Model" (2020)
- math/0304418: Biskup, "On the scaling of the chemical distance in long-range percolation models" (2003)
- 0706.2151: Deser & Woodard, "Nonlocal Cosmology" (2007)
- 0812.1059: Hehl & Mashhoon, "Nonlocal Gravity Simulates Dark Matter" (2008)
- Not on arXiv (verified via publisher record): Bertotti, Iess & Tortora, Nature 425, 374 (2003); Kim et al., A&A 640, A69 (2020).
