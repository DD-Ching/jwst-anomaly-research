# B-on-A2: destruction review of "Carry" (no identity, no locality)

Reviewer: System B. Target: `A2-no-identity.md`. All references were checked against arXiv
metadata or abstracts on 2026-10-09 (list at the end). Numbers marked [calc] are this reviewer's own
order-of-magnitude estimates.

## 1. Internal consistency

**1a. Carry is not well-defined (circular).** AX1 says cells are unnamed and "the same cell at t and t+1"
is meaningless. AX3 then defines c_uv = I(s_u(t); s_v(t+1)). A mutual information between a cell at t
and a cell at t+1 needs a joint distribution over a *specific pairing* (u at t, v at t+1), which means a
cross-time identification of cells. That is exactly the identity AX5 claims to *derive* from carry.
Relabeling-equivariance only lets you define the carry graph up to isomorphism *after* a pairing exists.
Identity is presupposed, not derived.

**1b. Carry is not conserved, so the "Kirchhoff law" (AX4) is false.** Pairwise mutual information is not
an additive, conserved flow under reversible maps. Two counterexamples [calc]:
- *XOR (synergy):* (a,b) uniform bits, reversible map (a,b)→(a⊕b, b). Then I(a;a′)=I(a;b′)=0, so the
  carry out of cell a is 0, yet a's bit survives (jointly, in a′ and b′). Carry out < information content.
- *Fan-out (copy):* (a,0,0)→(a, a⊕0, a⊕0) is a bijection on full configurations (CNOTs). Then
  c_{a→v}=1 for three targets, Σ_v c_av = 3 > b = 1. This violates the capacity axiom.
Reversibility conserves *joint* entropy, not the sum of pairwise MIs. Since §2 derives energy
conservation from "total carry is conserved", **energy conservation is not derived**.

**1c. The Margolus–Levitin step is misused.** ML (quant-ph/9710043) is an *inequality*: the time to reach an
orthogonal state satisfies t ≥ πħ/(2E), with E the mean energy *above the ground state*. It counts
orthogonal-state transitions, not bits of mutual information. One orthogonal transition can carry
anywhere from ~0 to log₂d bits. Setting E = (πħ/2)Ṅ with Ṅ in bits/s therefore (i) assumes the bound is
saturated, (ii) equates bits with orthogonal transitions, and (iii) uses rest energy where ML needs energy
above the ground state. This is Lloyd's operations-per-second identification (quant-ph/9908043;
quant-ph/0110141) under a new name. The electron figure, 2m_ec²/(πħ) = 4.94×10²⁰ s⁻¹, is arithmetically
right [calc], but it is a rate of orthogonal transitions, not bits.

**1d. Lorentz invariance is not derived.** The hop metric on a finite-valence graph is not isotropic in
general (on a cubic lattice it is the L¹ metric). Ignatowski's argument needs a continuum relativity
principle, which is the thing to be shown. A theorem forbids it outright: there is no way to associate a
finite-valency graph to a Poisson sprinkling consistently with Lorentz invariance (Bombelli, Henson &
Sorkin, gr-qc/0605006). AX4 (finite capacity, hence finite effective valence) sits squarely inside that
no-go. The model therefore has a preferred frame, and a hop length near the Planck scale generically gives
linear-in-E dispersion, bounded at E_QG,1 > 1.2 E_Pl (GRB 090510, 0908.1832) and > 10 E_Pl (LHAASO, 2402.06009). E = mc² also does not follow: the 4-vector is
*declared* ("p ∝ hop-carry rate has invariant norm mc"), not derived.

**1e. Gravity sector.** D2's algebra checks out [calc]. n = 1/(1−2m/r) gives a photon sphere at r = 4m and
b_c = 8m, versus 3√3 m, a factor 1.54. That is excluded: M87* gives a crescent diameter of 42±3 µas,
consistent with GR (1906.11243), and Sgr A* gives an image size within ~10% of Kerr (2311.09484).

## 2. Conservation, stability, causality

- **Energy.** E_out = E_in/s with s = (1−f_in)/(1−f_out) is a ratio of lapse-like factors. It is the
  gravitational redshift between two static potentials, which *conserves* Killing energy. So the claim
  "photon number conserved, energy not" is not new physics, and it fixes the size of s. Using the model's
  own f ≈ 2GM/(rc²): galaxy halo f ~ (200 km/s / c)² ≈ 4×10⁻⁷, cluster ≈ 10⁻⁵ [calc]. **So |s−1| ≲ 10⁻⁵**
  unless a chord mouth sits within ~10³ Schwarzschild radii of a compact object. The model has no cosmic
  expansion, so it has no other source of s ≠ 1. If expansion is added, s_t·s_ν = 1 holds for *any*
  redshift (Doppler, cosmological or gravitational) and is not chord-specific.
- **Momentum.** A photon leaving at a remote mouth in a new direction changes its linear and angular
  momentum by ~D·E/c, and nothing absorbs the recoil (chords are massless edges). Chords break the very
  homogeneity that §2 uses to derive momentum conservation, so momentum is violated at every transfer.

- **Causality.** A chord crossing a bulk distance D in one update is a spacelike signal. With exact Lorentz
  invariance, spacelike signalling plus boosts gives closed causal loops (the tachyonic antitelephone).
  The only escape is a preferred frame (the graph frame), which contradicts §2's Lorentz claim. The model
  must pick one: Lorentz violation or causal paradoxes.
- **Stability of geometry (quantitative kill of the chord amplitude).** In a d-dimensional small-world
  network, random shortcuts set a crossover length ξ such that there is about one shortcut per ξ^d
  volume. Above ξ, distances collapse to ~log N (Newman & Watts, cond-mat/9904419). Keeping Euclidean
  geometry out to the Hubble scale L ≈ 1.3×10²⁶ m therefore allows **≲ 1 long chord per Hubble volume**,
  out of (L/ℓ_Pl)³ ≈ 5×10¹⁸² cells [calc]. Then τ_ch ≈ n_ch σ_ch L ≲ σ_ch/L², which is fatal for every rate-based prediction (P2, P3, T). Reaching τ_ch = 10⁻³ needs
  σ_ch^{1/2} ≳ 0.03 L ≈ 140 Mpc [calc]. That is a mouth subtending several degrees at z ~ 1, which would
  appear as a degree-scale hole or duplicated patch in galaxy surveys and the CMB. So: either τ_ch is
  unobservably small, or the mouths are huge and should already have been seen.

## 3. Already excluded or already tested?

- **P1 (twins).** Lensed-GRB searches look only for *same-position, same-spectrum* repeats.
  Ahlgren & Larsson (2006.07095; ~2,700 GBM bursts, 11 yr) cut on position first, then on spectra and
  duration, then on light-curve cross-correlation. They found no candidates, and the most similar pairs
  were single-pulse bursts explained by population similarity. A 2025 re-analysis (2508.21413) rejects six
  millilensing candidates with χ², count-hardness and pulse-shape tests. CHIME (2204.06014; 172 bursts)
  searched for coherent same-burst copies with ns–100 ms delays. **None of these tested wide separation,
  arbitrary delay, or time-stretch s ≠ 1.** So P1 is untested rather than excluded.
- **P2 (sub-floor DM).** The Macquart relation (2005.13161) and the larger host-localized samples (e.g.
  2409.16952) show the DM_cosmic distribution with a sharp lower envelope. No host-localized FRB is known
  to sit significantly below it. With O(10²) localized FRBs and zero outliers, τ_ch ≲ 0.03 (95%). P2 is
  also **internally inconsistent**: a chord-transported burst arrives from the exit mouth, not from the
  host, so it would fail host association. The host-localized sample is blind to chord events except for
  chords whose two mouths lie on the same sightline. The predicted incidence is therefore τ_ch × (alignment
  probability ≪ 1), not τ_ch.
- **P3 (achromatic η).** This is exactly the cosmic-opacity / distance-duality test. Avgoustidis et al.
  (1004.2053) give ε = −0.04 (+0.08/−0.07) at 2σ and Δτ < 0.012 (95%) for 0.2 < z < 0.35. With DESI
  data, Dhawan & Mörtsell (2506.22599) find no DDR violation, a maximum Δm ~ 0.05 mag (τ ≲ 0.05) at the
  highest z. DESI DR2 + Pantheon+/Union3 are consistent to ≲ 1σ (2509.19899). Grey opacity at τ_ch ≳ 0.02–0.05
  is excluded; below that, P3 is identical to existing photon-number-violation opacity tests.
- **P4 (shadow/ring).** The linear law is excluded (above). See §5 for why the nonlinear version is empty.
- **T (transfer covariance).** No wide-angle signed photon-number covariance study is known to me. But it
  is ill-posed (§5).

## 4. Renaming check

- **Chords = "disordered locality".** Graph models in which space emerges from a network, with leftover
  non-local links as defects, are quantum graphity (Konopka, Markopoulou & Smolin, hep-th/0611197) and
  disordered locality (Markopoulou & Smolin, gr-qc/0702044). Cosmological consequences
  (Prescod-Weinstein & Smolin, 0903.5303) and dispersion consequences (Caravelli & Markopoulou, 1201.3206)
  have been worked out. A2's D1 is this idea; only the photon-transfer phenomenology is a new framing.
- **Energy = (πħ/2)·rate** is Lloyd (quant-ph/9908043, quant-ph/0110141).
- **Capacity → refractive-index gravity** is in the entropic/information-gravity family (Verlinde,
  1001.0785) and the optical-metric analogue. Chords as microscopic wormholes connect to ER=EPR-type
  ideas (1306.0533).
- **Genuinely new:** the wide-angle, photon-number-matched twin criterion and the signed T statistic.

## 5. Per-prediction analysis and triage

**P1 — partially novel. CONDITIONAL PASS (degraded).** As written, it fails the internal check: s ≈ 1
within 10⁻⁵ (§2), so the s_t·s_ν test is vacuous. "Photon fluence agrees after correcting each mouth's
beam-loss fraction" has a free correction for every pair, so it tests nothing. It survives only as a
**model-agnostic, cheap search**: an all-sky, any-delay search for spectro-temporal twins at separations
well beyond the localization error, with s fixed at 1 (primary) and free s as a secondary scan. This is
the one combination prior lensing searches cut away.
- Mimics to exclude:
  - single-FRED pulse degeneracy (most of Ahlgren & Larsson's false positives);
  - pulse width ∝ E^-0.4 (astro-ph/9504075). A true redshifted twin seen in a fixed band stretches as
    ~s^0.6, not s, so the comparison must map energy channels by s;
  - Amati E_pk–E_iso clustering (astro-ph/0205230);
  - duplicate catalogue entries and re-triggers of the same burst; GBM/Swift/Konus cross-instrument
    duplicates;
  - FRB repeaters (same source, similar morphology);
  - ordinary strong lensing (should sit at < localization error; use it as a positive control);
  - trials factor: N²/2 ≈ 4×10⁶ GBM pairs, 10⁷ CHIME pairs, times delay and s grids.
- Data: Fermi GBM burst catalogue with TTE and spectral fits (2002.11460 and later, ~4,000 bursts); CHIME/FRB
  Catalog 2 (2601.09399; 4,539 bursts, 0.983 ms dynamic spectra); Swift-BAT as a cross-check.
- Sensitivity [calc]: restrict to the ~500 brightest multi-pulse GRBs (~1.2×10⁵ pairs) and calibrate
  a false-match rate ≲ 10⁻⁶ per pair on time-scrambled and pulse-shuffled catalogues. Zero twins then gives
  a twin fraction ≲ 6×10⁻³ per burst (95%). Reaching A2's 10⁻³ target needs the full catalogue and
  CHIME. A twin must also fall within the instrument's dynamic range (flux ratio ≲ 30), and the model
  supplies no reason to expect that.

**P2 — partially known. FAIL.** Internally inconsistent with host localization (§3); already bounded at
τ_ch ≲ 0.03; mimicked by underdense sightlines and host-DM misestimates (the lower tail is model-dependent).

**P3 — known. FAIL.** It is the existing cosmic-opacity test, already at τ ≲ 0.012–0.05. It is the photon-number-violation opacity class of 1004.2053. A2 adds no
discriminating feature, because "no P1/P2 counterpart" depends on free ρ_ch and σ_ch.

**P4 — known/empty. FAIL.** Any static, spherically symmetric null-geodesic structure can be reproduced by
an isotropic n(r). Choosing f = 1 − 1/n_GR, with n_GR = (1+m/2r)³/(1−m/2r) (monotone, f → 2m/r, f → 1
at the horizon), reproduces GR's shadow, photon-ring spacing and α₂ = 15π/4 (m/b)² *exactly*. The kill
criterion can therefore never fire: there is always a fitting f. P4 is unfalsifiable as stated.

**T — novel form. FAIL.** The sign is wrong by the model's own physics. Photons a mouth captures were headed
toward the mouth, not toward the observer, so the direct image loses a *constant* fraction and still
tracks the intrinsic variability. The copy then correlates **positively** with the source at lag Δt. T < 0
needs a time-variable chord throughput, which A2 never introduces. It is also ill-posed: the copy appears
at the exit mouth, not on another catalogued quasar, so pair-wise covariance between catalogued sources
has no expected signal.

| Item | Novelty | Verdict | Main reason |
|---|---|---|---|
| P1 twins | partial | CONDITIONAL PASS (s=1 wide-angle twin search only) | Untested combination; s≈1 per model; free beam-loss factor |
| P2 sub-floor DM | partial | FAIL | Chord bursts lose host association; τ_ch ≲ 0.03 already |
| P3 achromatic η | known | FAIL | Existing opacity tests; τ ≲ 0.01–0.05 |
| P4 shadow/ring | known | FAIL | Free f(r) reproduces GR exactly; unfalsifiable |
| T covariance | novel form | FAIL | Wrong sign by own physics; ill-posed pairing |
| Foundations | — | FAIL | Carry undefined (circular), not conserved; Lorentz no-go; small-world ⇒ τ_ch tiny |

## References (verified via arXiv on 2026-10-09)

- quant-ph/9710043 — Margolus & Levitin, *The maximum speed of dynamical evolution*, 1997.
- quant-ph/9908043 — Lloyd, *Ultimate physical limits to computation*, 1999.
- quant-ph/0110141 — Lloyd, *Computational capacity of the universe*, 2001.
- gr-qc/0605006 — Bombelli, Henson & Sorkin, *Discreteness without symmetry breaking: a theorem*, 2006.
- cond-mat/9904419 — Newman & Watts, *Scaling and percolation in the small-world network model*, 1999.
- 0908.1832 — Abdo et al. (Fermi), *Testing Einstein's special relativity with Fermi's short hard GRB 090510*, 2009.
- 2402.06009 — LHAASO, *Stringent Tests of Lorentz Invariance Violation from LHAASO Observations of GRB 221009A*, 2024.
- 1906.11243 — EHT, *First M87 EHT Results. VI. The Shadow and Mass of the Central Black Hole*, 2019.
- 2311.09484 — EHT, *First Sagittarius A\* EHT Results. VI: Testing the Black Hole Metric*, 2023.
- 2006.07095 — Ahlgren & Larsson, *A search for lensed gamma-ray bursts in 11 years of observations by Fermi GBM*, 2020.
- 2508.21413 — *A re-identification of six Candidate Gravitationally Lensed Gamma-Ray Bursts*, 2025.
- 2204.06014 — Kader et al. (CHIME/FRB), *A High-Time Resolution Search for Compact Objects using FRB Gravitational Lens Interferometry with CHIME/FRB*, 2022.
- 2005.13161 — Macquart et al., *A census of baryons in the Universe from localized fast radio bursts*, 2020.
- 2409.16952 — Connor et al., *A gas-rich cosmic web revealed by the partitioning of the missing baryons*, 2024.
- 1004.2053 — Avgoustidis et al., *Constraints on cosmic opacity and beyond the standard model physics from cosmological distance measurements*, 2010.
- 2506.22599 — Dhawan & Mörtsell, *Implications for dark energy of cosmic transparency in light of DESI data*, 2025.
- 2509.19899 — *Calibration-independent consistency test of DESI DR2 BAO and SNIa*, 2025.
- hep-th/0611197 — Konopka, Markopoulou & Smolin, *Quantum Graphity*, 2006.
- gr-qc/0702044 — Markopoulou & Smolin, *Disordered locality in loop quantum gravity states*, 2007.
- 0903.5303 — Prescod-Weinstein & Smolin, *Disordered Locality as an Explanation for the Dark Energy*, 2009.
- 1201.3206 — Caravelli & Markopoulou, *Disordered locality and Lorentz dispersion relations: an explicit model of quantum foam*, 2012.
- 1001.0785 — Verlinde, *On the Origin of Gravity and the Laws of Newton*, 2010.
- 1306.0533 — Maldacena & Susskind, *Cool horizons for entangled black holes*, 2013.
- astro-ph/9504075 — Fenimore et al., *Gamma-Ray Burst Peak Duration as a Function of Energy*, 1995.
- astro-ph/0205230 — Amati et al., *Intrinsic spectra and energetics of BeppoSAX GRBs with known redshifts*, 2002.
- 2002.11460 — von Kienlin et al., *The Fourth Fermi-GBM Gamma-Ray Burst Catalog: A Decade of Data*, 2020.
- 2601.09399 — CHIME/FRB, *The Second CHIME/FRB Catalog of Fast Radio Bursts*, 2026.
