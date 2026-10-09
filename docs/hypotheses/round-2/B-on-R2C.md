# B on R2-C "Two-Slot World": destruction report (System B, round 2)

Labels: **[V]** verified against the cited source, **[C]** my own calculation, **[J]** my judgement. Every arXiv ID was
opened.

## Verdict

R2-C's correct kinematic core is textbook material:
- spinor → null vector;
- SL(2,C) → SO⁺(1,3);
- s_ij = |ε_ij|².

It is also nearly word for word von Weizsäcker's ur-theory, which already had k^μ = σ^μ ū u **and** N ≈ 10¹²⁰ urs tied
to the Planck length [10].

The new parts fail in four ways:
- The axioms produce no fermions, no positions, no Λ > 0 background and no redshift.
- The new conservation law contradicts neutrino mixing.
- The one sharp number (Σm_ν = 58.8 meV) is the minimal-seesaw prediction [12–14].
- Together with the model's own w = −1, that number sits at 2.2–3.4σ tension with DESI+CMB [3, 21].

This project cannot test any of it.

## 1. Internal consistency

**1a. 3+1 from C² + ε: yes for momenta, no for positions.**
- [V] Herm(2) with det is Minkowski space, and uu† is null. This is Penrose [8, 9], Lyre's eq. (1) [10], and the
  R/C/H/O ladder [11].
- It yields null *momenta*, not events. A2 makes the translation-invariant ε_ij the only observables, so no
  observable depends on position. The incidence ω = i x u is introduced and then made unobservable.
- [J] As a result:
  - "affine distance r" (§3) is undefined, so the inverse-square law is assumed, not derived;
  - there is no locality;
  - the "exact zero dispersion" holds only because nothing propagates.

**1b. A4 fails.** [C]
- The lowest-degree, non-zero, symmetric holomorphic function of the ε_ij is the **constant**: an empty world.
- Excluding constants, the first admissible choice is Ψ = Σε_ij², because Σε_ij is antisymmetric. It fails three ways:
  - It is not homogeneous per element (degrees 0 or 2), so the claimed per-element helicity conservation is false.
  - It gives helicity 0 or −1, never −½.
  - It is exchange-**symmetric**, so the elements are bosons. Calling them Weyl fermions violates spin-statistics.
    Weizsäcker's urs are explicitly "Bose urs" [10].
- The antisymmetric option, Pf(ε_ij), vanishes for N ≥ 4 because rank ε = 2.
- SL(2,C) orbits have infinite volume, so any invariant Ψ has ∫|Ψ|² = ∞ and the Born frequencies of A3 are undefined.

**1c. Mass.**
- [V/C] det Σ u_i u_i† = Σ_{i<j} |ε_ij|² (Cauchy–Binet) is the spinor-helicity invariant mass s [15].
- That is the invariant mass of *free* null momenta: a continuum, with no binding, spectrum or stability.
- **Weinberg–Witten** [16] forbids massless |h| > 1 states when a Lorentz-covariant, conserved stress tensor exists.
  R2-C asserts both (§2, P = Σuu†), so its "graviton-like composites", and with them "c_gw = c_γ identically", cannot
  exist.

**1d. Λ from N.**
- [C] The arithmetic is right: N = 4.69×10¹²². Two small slips: log₂N = 407.5, not 409; ½ln N = 141.2 ✓.
- It is one parameter fitted to one datum, and borrowed:
  - "finite N ↔ Λ" is Banks and the Bousso N-bound [17, 18];
  - "N = 10¹²⁰ ↔ ℓ_P = 10⁻⁶⁰R" is Weizsäcker and Görnitz [10].
- Ω_Λ(z) = S_H/S_dS = (R_H/R_dS)² is an identity of flat ΛCDM, not a result of this model.
- **Contradiction:**
  - §2 asserts exact Poincaré invariance, but Λ > 0 has the de Sitter group, with no global translations.
  - Variant 4c.1 itself shows the only SL(2,C)-commuting flow gives no observable redshift.
  - So 1+z, the §3 distances and "exact Etherington" have nothing to apply to inside A1–A4. Etherington is inherited
    from metric-theory photon conservation [19].

**1e. m_lightest = 0 is a choice.**
- [J] A1–A4 say nothing about confinement, lepton flavour, weak coupling or mixing.
- "Massless + left-handed + neutral ⇒ ν₁" is a pattern match, and the model has no charge at all.
- ν₃ (IO) or a dark Weyl fermion would fit equally well. NO is chosen from data.
- The parity "explanation" is contradicted by observation. Charged leptons and quarks are composites in the model,
  yet only their left-handed parts couple to the W.

## 2. Conservation, stability, causality

- **Causality** ✓: a sum of future-null vectors is future-causal. This is a standard fact, not a new mechanism.
- **Stability**: cannot be assessed, because there is no Hamiltonian.
- **Element-number conservation refutes itself** [J].
  - ν₁ has 1 element, while ν₂ and ν₃ have ≥ 3.
  - A conserved N creates superselection sectors, so 1+3 superpositions are forbidden.
  - Yet flavour states are such superpositions: |U_e1|² ≈ 0.68, sin²θ₁₂ = 0.308 ± 0.012 [1]. All oscillation data are
    interference between those sectors.
  - The law survives only mod 2, which is ordinary fermion parity.
- If anti-elements count −1, element number is a lepton-number-like U(1). That forbids 0νββ and contradicts the model's
  own m_ββ row.
- **CPT and antiparticles** are asserted, not derived. A3 gives every element negative helicity, so the right-helicity
  reactor ν̄ used in the θ₁₂ and Δm²₂₁ fit [1] has no bare-element counterpart.

## 3. Already excluded?

**Oscillation floor** [V 1, C]. NuFIT 6.0 (IC24 + SK-atm): Δm²₂₁ = 7.49 ± 0.19×10⁻⁵, Δm²₃₁ = +2.513 (NO) and
Δm²₃₂ = −2.484 (IO) ×10⁻³ eV², sin²θ₁₂ = 0.308, sin²θ₁₃ = 0.02215.

| m_lightest = 0 | Σm_ν | m_β | m_ββ |
|---|---|---|---|
| NO | **58.78 meV** (3σ: 57.8–59.8) | 8.84 meV | 1.50–3.72 meV |
| IO | **98.9 meV** (IC19: 99.5) | 48.8 meV | 18.2–48.2 meV |

- DESI quotes 58.78 ± 0.23 and 98.92 ± 0.41 meV [3]. R2-C's NO values are right. Its IO value of 99.9 is 1 meV high,
  and its ±0.5 is about twice the true 1σ.
- "Data disfavour IO" depends on the dataset: with SK-atm, IO has Δχ² = 6.1; without it (IC19), IO is the best fit,
  with NO at Δχ² = 0.6 [1]. The rejection of IO comes from cosmology.

**Cosmology in ΛCDM, which w = −1 commits R2-C to** [V]:
- DESI DR2 + CMB [3]:
  - Σm_ν < 64.2 meV (95%), σ = 20 meV;
  - Feldman–Cousins < 53 meV;
  - Σm_ν,eff = −101 (+47, −56) meV, in **3.0σ** tension with the NO floor (2.7σ GoF, 4.0σ parameter-shift);
  - profile-likelihood parabola μ₀ = −36, σ = 43 meV, which puts 58.8 meV **2.2σ** high and the IO floor 3.1σ high
    [C].
- DR1 gave −160 ± 90 meV [4]. Excess CMB lensing drives the preference [5].
- CMB-SPA + DESI DR2 + DES: < 52 meV [20].
- 2026, with DES-Dovekie: −75 (+39, −53) meV [21], which puts 58.8 meV at ≈ 3.4σ (Gaussian extrapolation) [C].
- Escape routes:
  - w₀w_a: Σm_ν < 163 meV [3, 6];
  - a sign-switching Λ [21];
  - neutrino decay [22];
  - a high optical depth, τ = 0.11, which gives Σ = 100 (+40, −50) meV [23].

  R2-C forbids the first three.
- **Result:** in tension at 2.2–3.4σ, but **not excluded at ≥ 3σ under the author's rule**. The rule needs a ≤ 10 meV
  systematics budget, and τ alone moves Σ by ~40 meV.

**Other experiments** [V]:
- KATRIN: m_ν < 0.45 eV (90%) [7]. The final 1000-day reach is ≈ 0.3 eV, so the "m_β > 20 meV" kill rule cannot fire.
- KamLAND-Zen: m_ββ < 28–122 meV [25].
- β:
  - Planck + WMAP: 0.342° (+0.094, −0.091), 3.6σ [26];
  - ACT DR6: 0.215° ± 0.074°, 2.9σ, with "systematics … not understood" [27];
  - see also [28].

  None is at 5σ.
- w₀w_a:
  - DESI + CMB: 3.1σ; with SNe, 2.8–4.2σ [2];
  - DES-Dovekie: 3.2σ [29];
  - full DES + CMB: 3.0σ [30];
  - DESI DR2 + Lyα AP: 2.7σ (3.2σ with SNe) [31].
- c_gw: |Δc/c| ≲ 10⁻¹⁵ [32], which GR also satisfies.

## 4. Renaming check

| R2-C element | Prior art |
|---|---|
| spinor ↦ null vector; SL(2,C) → Lorentz | ur-theory [10]; Penrose [8, 9] |
| N ≈ 10¹²² ↔ ℓ_P/R_H | Weizsäcker–Görnitz N = 10¹²⁰ [10] |
| Λ ↔ finite N | Banks [17]; Bousso [18] |
| M² = Σ\|ε_ij\|²; ω = i x u | spinor-helicity [15]; twistor incidence [8, 9] |
| R/C/H/O → 3/4/6/10 | [11] |
| Etherington | [19] |
| m_lightest = 0 | minimal seesaw, rank 2 [12–14]; ~10⁻¹³ eV at two loops [33] |

The only differences from ur-theory are a fixed N (Weizsäcker's N grows and serves as cosmic time [10]) and the
neutrino identification.

## 5. Novelty

- **P1** is not new. It is numerically identical to the minimal seesaw [12–14]; only the *reason* is new, and no
  observation can separate the two.
  - m₁ = 5 meV gives Σ = 65.4 meV, and m₁ = 10 meV gives 74.3 meV [C].
  - Separating m₁ = 0 from m₁ ≤ 5 meV at 3σ needs σ(Σ) ≈ 2 meV, beyond any planned survey.
  - "Exactly zero" versus 10⁻¹³ eV [33] is observationally empty.
- **P2** is not new. Every zero is the GR + ΛCDM default.
- **Λ ↔ N** is not new and is not a prediction [10, 17, 18].
- **P3a** is new but empty. Δt = (d/2c)(m/E)² at 10 kpc and 10 MeV is 5 ns for m = 1 meV, 0.4 µs for m₂ and 13 µs for
  m₃ [C]. Burst structure lasts milliseconds, and detectors see flavour, not mass eigenstates.
- **P3b** has the same content as m₁ = 0 and no separate observable.

## 6. Triage

| # | Prediction | Verdict | Data, sensitivity | Mimics, degeneracies |
|---|---|---|---|---|
| P1 | Σm_ν = 58.8 meV (NO) | **Conditional PASS** (falsifiable, live; not novel) | Planck PR4, ACT DR6, SPT-3G, DESI DR2 likelihoods (Cobaya/CLASS); σ ≈ 20 meV; 2.2–3.4σ tension [3, 21] | τ [23], lensing excess [5], dark energy [6, 21], ν decay [22], prior edge [20] |
| P1′ | Σ ≈ 99 meV (IO) | **FAIL** (~3σ disfavoured; author has abandoned it) | as P1 | as P1 |
| P1-β | m_β = 8.8 meV | **FAIL as a test** (reach ≈ 300 meV) | [7] | — |
| P1-ββ | m_ββ = 1.5–3.7 meV | **FAIL** (beyond reach; contradicts element conservation) | [25] | nuclear matrix elements |
| P2a | zero dispersion, c_gw, tensor-only, 1/d_L, Etherington | **FAIL** (= GR, non-discriminating) | [32] | — |
| P2b | β = 0 | **Conditional PASS** (= ΛCDM default) | Planck, ACT, SPT EB/TB; σ ≈ 0.07° [27] | angle miscalibration, dust EB, I→P leakage |
| P2c | w = −1 | **Conditional PASS** (= Λ) | DESI + CMB + SNe; 2.7–3.2σ against [29–31] | SN calibration (4.2 → 3.2σ) [29] |
| Λ–N | Λ = 3π/(N ln2 ℓ_P²) | **FAIL** (calibration; contradicts Poincaré) | — | — |
| P3a | ν₁ zero delay | **FAIL** (ns versus ms) | — | emission physics |
| P3b | relic ν₁ w = 1/3 | **FAIL** (no separate observable) | — | — |

**Sharpest test (joint):** "Σ = 58.8 ∧ w = −1 ∧ no new ν physics" is the ΛCDM-at-the-floor configuration, now at
~3σ. The model has no escape parameter, so a ≥ 5σ result from any of these kills it:
- the DESI DR3 + CMB Σm_ν,eff tension;
- w₀w_a evidence;
- β.

**Can this project add anything?** No [J]. P1, P2b and P2c run on CMB, BAO and SN likelihoods with nuisance models
that this JWST imaging repository does not have. Re-running public chains would only reproduce [3, 20, 21, 29]. The
right action is to log R2-C as "tracked externally: DESI DR3, ACT/SPT β calibration, final KATRIN result".

## References (all opened 2026-10-09)

1. Esteban et al., NuFIT 6.0, arXiv:2410.05380 (Table 1; IO Δχ² = 6.1 with SK-atm, NO Δχ² = 0.6 without).
2. DESI Collaboration, DESI DR2 Results II, arXiv:2503.14738, PRD 112, 083515 (2025).
3. Elbers et al. (DESI), Constraints on neutrino physics from DESI DR2 BAO and DR1 full shape, arXiv:2503.14744
   (v2 full text: Σm_ν,eff = −0.101 (+0.047, −0.056) eV; FC < 0.053 eV; floors 58.78 ± 0.23 and 98.92 ± 0.41 meV).
4. Craig, Green, Meyers, Rajendran, "No νs is Good News", arXiv:2405.00836.
5. Green & Meyers, "The Cosmological Preference for Negative Neutrino Mass", arXiv:2407.07878.
6. Elbers et al., "Negative neutrino masses as a mirage of dark energy", arXiv:2407.10965.
7. KATRIN (Aker et al.), arXiv:2406.13516, Science 388, 180 (2025).
8. Penrose, "Twistor Algebra", J. Math. Phys. 8, 345 (1967), doi:10.1063/1.1705200.
9. Atiyah, Dunajski, Mason, "Twistor theory at fifty", arXiv:1704.07464.
10. Lyre, "C. F. von Weizsäcker's Reconstruction of Physics", arXiv:quant-ph/0309183 (eq. 1; N = 10¹²⁰; Bose urs).
11. Baez & Huerta, "Division Algebras and Supersymmetry I", arXiv:0909.0551.
12. Frampton, Glashow, Yanagida, "Cosmological Sign of Neutrino CP Violation", arXiv:hep-ph/0208157.
13. Guo, Xing, Zhou, minimal seesaw review, arXiv:hep-ph/0612033, IJMPE 16, 1 (2007).
14. Xing & Zhao, "The minimal seesaw and leptogenesis models", arXiv:2008.12090, Rep. Prog. Phys. 84, 066201 (2021)
    (§2: "rank two, leading us to a massless Majorana neutrino at the tree level").
15. Elvang & Huang, "Scattering Amplitudes", arXiv:1308.1697.
16. Weinberg & Witten, "Limits on massless particles", Phys. Lett. B 96, 59 (1980), doi:10.1016/0370-2693(80)90212-9.
17. Banks, "Cosmological Breaking of Supersymmetry?", arXiv:hep-th/0007146.
18. Bousso, "Positive vacuum energy and the N-bound", arXiv:hep-th/0010252.
19. Ellis, "On the definition of distance in general relativity: I. M. H. Etherington (1933)", GRG 39, 1047 (2007),
    doi:10.1007/s10714-006-0355-5.
20. Hou et al., "Constraints on the Sum of Neutrino Masses from ACT DR6 and DESI DR2 …", arXiv:2606.17994.
21. Kıbrıs, Elbers, Akarsu et al., "Negative neutrino mass or negative dark energy?", arXiv:2605.21456.
22. Franco Abellán, "Neutrino decays as a natural explanation of the neutrino mass tension", arXiv:2601.04312,
    PRD 113, 123527 (2026).
23. Sullivan, de Belsunce, Ivanov, "Cosmological Concordance in an Especially Opaque Universe", arXiv:2606.30903.
24. Calabrese et al. (ACT), DR6 extended models, arXiv:2503.14454.
25. KamLAND-Zen (Abe et al.), arXiv:2406.11438, PRL 135, 262501 (2025).
26. Eskilt & Komatsu, arXiv:2205.13962.
27. Diego-Palazuelos & Komatsu, "Cosmic Birefringence from ACT DR6", arXiv:2509.13654.
28. Minami & Komatsu, arXiv:2011.11254, PRL 125, 221301 (2020).
29. Popovic et al. (DES), DES-Dovekie, arXiv:2511.07517, MNRAS (2026).
30. DES Collaboration, full-DES dynamical dark energy, arXiv:2605.27221.
31. DESI Collaboration, DESI DR2 Results IV (Lyα AP), arXiv:2607.27410.
32. LIGO/Virgo/Fermi/INTEGRAL, GW170817 / GRB 170817A, arXiv:1710.05834.
33. Davidson, Isidori, Strumia, "The smallest neutrino mass", arXiv:hep-ph/0611389.

Also consulted, but not used as evidence: the SPT-3G D1 abstract (arXiv:2506.20707), ACT DR6 ΛCDM (arXiv:2503.14452),
King (arXiv:2502.07877). KATRIN's 1000-day completion is reported in a search-result summary of KIT/MPIK news, not a
primary source. It is labelled as such and does not affect any verdict.
