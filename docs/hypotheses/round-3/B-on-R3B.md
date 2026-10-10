# System B on R3-B "Phase-Locked Dilaton" (PLD), round 3

Labels: **[obs]** published measurement, **[der]** my derivation from published numbers, **[A]** System A's claim,
**[hyp]** hypothesis. Arithmetic was done in a scratch script (not committed).
References verified via the arXiv API, 2026-10-09.

**Verdict in one line.** PLD is a non-universally coupled ω = 0 Brans–Dicke scalar whose charge sits only in neutron
stars. Published NS–WD timing already excludes its predicted band: the upper half (s ≥ 8×10⁻³) at more than 6σ, and
the self-consistent minimum (α² = 5.3×10⁻⁶) at about 2σ. Only the corner that A itself calls non-self-consistent
survives. The one sharp surviving prediction is a laboratory measurement, outside this project's scope.

## 1. Internal consistency

**Checks that pass [der].**
- *Dipole formula.* `Ṗ_b^D = −(4π²/P_b) T_⊙ (m_p m_c/M) Δα²` is exactly Freire et al. 2012 eq. 21 and Zhao et al. 2022
  eq. 2, in DEF normalisation. A's D values (1.0×10⁻⁹ for J1738, 3.5×10⁻⁹ for J0348) and the dipole/quadrupole ratio
  for J1738 (3.7×10⁴ α²) reproduce to 2%.
- *BCS condensate fraction.* Integrating Σ_k u_k²v_k² gives n₀/n = 3πΔ/(8ε_F). That is the weak-coupling limit of
  Salasnich, Manini & Parola 2005, and A's table entries follow from it.
- *Feedback algebra.* For a uniform sphere with a canonical σ, the self-energy is −(6/5)β²G(f_cM)²/R. Comparing it with
  ½N(0)Δ² gives K′/N(0) = 1.11 β² C x_i (mc²/E_F), which is A's expression. Its values (0.17–0.35 for protons, 0.78
  for neutrons) are also right.

**Errors and unforced assumptions.**
1. **The β² normalisation is off by a factor of 2.** With J_σ = −β√(8πG) f_c T and a canonical σ, the force between
   two bodies is G(1 + 2β²q_Aq_B). The DEF α that enters the dipole formula is therefore α = √2 β s, not βs. ω = 0
   Brans–Dicke (metric f(R)) has α_DEF² = 1/(2ω+3) = 1/3, which corresponds to a *canonical* β² of 1/6 (Sotiriou &
   Faraoni 2010). A's "alternative counting, 1/6" is in fact the ω = 0 value. A mixes the two conventions. The
   dipole calculation in §4.3 is consistent only if β_canonical² = 1/6. In that case:
   - G_AB = G(1 + α_Aα_B), not G(1 + 2α_Aα_B);
   - the He-II excess is (1/3)f_c² ≈ 1.6–2.7×10⁻³, not 3–5×10⁻³.

   If A really means β² = 1/3 canonical, every dipole α² doubles and the minimum becomes 1.06×10⁻⁵, which is excluded
   at about 4σ (§2).
2. **A5 does not fix β.** A Cauchy-type equality of stiffnesses says nothing about the scalar–matter coupling, so ω
   stays a free parameter. "β² = 1/3" is a choice, not a derivation.
3. **The trace is not −ρc².** In NS cores T^μ_μ = −ρ + 3p, with p/ρ ≈ 0.1–0.3 at 1.4–2 M_⊙. This lowers s by about
   20–50%, and more for heavy stars. A ignores it. It pushes PLD toward the surviving corner, but it also weakens the
   claimed s(M) shape.
4. **The crust lies outside weak coupling.** Salasnich 2014 finds that neutron-matter n₀/n peaks at 0.42 in the inner
   crust (BCS–unitarity crossover). The (3π/8)Δ/E_F formula is not valid there. The crust term is uncertain by about
   ×3, in an unknown direction.
5. **The ODLRO-only charge is load-bearing at about 10¹²** (A's weak point 1, made quantitative). Finite nuclei carry
   the same local pair amplitude as the NS crust (Δ ≈ 12/√A MeV). The same formula gives f ≈ 0.055 for Ti, 0.027 for Pt
   and about 0.07 for an Earth-mix. MICROSCOPE would then see η ≈ 2β²f_E Δf ≈ 7×10⁻⁴ to 1.3×10⁻³. The measured value
   is η = (−1.5 ± 2.3 ± 1.5)×10⁻¹⁵ (Touboul et al. 2022), so the overshoot is 10¹¹–10¹². PLD survives only if the
   σ-source is the *thermodynamic-limit* ODLRO, which is non-local. No local Lagrangian has that property, so A4 is a
   rule, not a coupling.
6. **The feedback is O(1).** The σ self-energy is 1.6× the BCS condensation energy for neutrons (2K′/N(0)), so the
   claimed "enhancement, not runaway" relies on a mean-field, uniform-sphere functional far beyond the regime where it
   holds. Since Δ ∝ e^{K′/N(0)} with K′ ∝ C, the exponent itself carries the ×2 uncertainty A admits. The feedback
   also makes s increase with compactness, which works against the claimed non-monotonic s(M).
7. **"α_WD = α_BH = 0 exactly" is not distinctive.** In DEF theory, WDs already have α ≈ α₀ (≲ 3×10⁻³ by Cassini) and
   black holes have zero scalar charge (no-hair). See Freire et al. 2012 §5 and Takeda et al. 2023.

## 2. Already excluded?

Published excess decay Ṗ_b^xs = Ṗ_b^int − Ṗ_b^GR **[obs]** is converted to α² = −Ṗ_b^xs/D **[der]**. The 95% upper
limit is computed with a truncated-Gaussian posterior on α² ≥ 0. z is the exclusion significance of each PLD value.

| System (m_p) | Ṗ_b^xs (fs s⁻¹) | α² from data | 95% UL α² (α) | z at 5.3e-6 / 2.1e-5 / 8.5e-5 |
|---|---|---|---|---|
| J1738+0333 (1.47) [Freire 2012, as tabulated in Zhao 2022] | 3.15 ± 3.69 | (−3.1 ± 3.6)e-6 | 5.3e-6 (2.3e-3) | 2.3 / 6.7 / 24 |
| J1738+0333, 2026 update [der from Vaglio 2026] | ≈ 0.5 ± 3.1 | (−0.5 ± 3.0)e-6 | 5.6e-6 (2.4e-3) | 1.9 / 7.1 / 28 |
| J1012+5307 (1.72) | 4.8 ± 5.1 | (−8.6 ± 9.1)e-6 | 1.3e-5 (3.6e-3) | 1.5 / 3.2 / 10 |
| J2222−0137 (1.83) | −6.3 ± 7.6 | (+8.9 ± 10.8)e-6 | 2.8e-5 (5.3e-3) | — / 1.1 / 7.1 |
| J0348+0432 (2.01) | −15 ± 46 | (+4 ± 13)e-6 | 2.9e-5 (5.4e-3) | — / 1.3 / 6.1 |
| J1909−3744 (1.49) | −1.7 ± 7.8 | (+6 ± 29)e-6 | 6.1e-5 (7.8e-3) | — / 0.5 / 2.7 |
| **Combined, mass-independent α** | | (−2.2 ± 3.1)e-6 | **4.8e-6 (2.2e-3)** | **2.4 / 7.5 / 28** |

The Vaglio et al. 2026 row is my derivation: observed −18.2 ± 2.5 fs s⁻¹, minus their Galactic (−0.30) and
Shklovskii (+9.3) terms, minus GR −27.7 (+1.5 −1.9) from Freire 2012. It is a preprint.

- **K1 is met for most of PLD now.**
  - s ≥ 8×10⁻³ is excluded at 6.7σ by J1738 alone.
  - The "self-consistent minimum" α² = 5.3×10⁻⁶ is disfavoured at 1.9–2.4σ. That is **not yet 95%-decisive** in
    a two-sided sense, but it sits exactly at the 95% upper limit.
  - With the β² = 1/3-canonical reading (§1.1) the minimum is excluded at 4.0–4.1σ.
  - The no-feedback corner α² = 1.3×10⁻⁶ is allowed (z = 0.6–1.1).
  - Including the trace correction (§1.3) moves the minimum to about 2.6×10⁻⁶ (z ≈ 1.1–1.6).
- **A's recalled J1738 ±15% is roughly right.** The ratio Ṗ_b^int/Ṗ_b^GR is 0.94 ± 0.13 (Freire 2012). A's inference
  that "most of the range is excluded" holds.
- **Combined DEF analyses.** Zhao et al. 2022 combine 7 pulsars and get |α_A| ≲ 6×10⁻³ (90%) for all masses and EOSs
  within DEF. Batrakov et al. 2023 (DDSTG, J2222) and Miao et al. 2026 (J1913+1102, a double NS, so only Δα) tighten
  this further. These are model-specific, but in DEF the per-system Δα bounds above carry over to any α(M) when the
  companion is a WD. Shao et al. 2017 covers the earlier 5-system version.
- **NS–BH gravitational waves (K4).** Santos et al. 2025, using GW230529, GW200105, GW200115 and GW190814, find
  |α̂| ≤ 0.028 (90%, κ₀ = ½), so α̂² ≲ 8×10⁻⁴. That is about 150× weaker than the pulsar bound. Takeda et al. 2023 is
  comparable for GW200115. **K4 cannot bite** before third-generation detectors.
- **Laboratory.**
  - Tajmar et al. (2006–2009) reported "gravitomagnetic"/gyroscope anomalies near rotating Nb and He. These are
    vector, rotation-dependent and unreplicated, and they later shifted to a non-superconducting "helium" origin.
    Tajmar's own 0707.3806 reports a Canterbury ring-laser null (Graham et al.). I did not open the Graham paper
    itself, so treat that null as second-hand.
  - Tajmar et al. 2004 found no weight change of YBCO/BSCCO through T_c.
  - PLD is **blind to all of these**: Earth's charge is zero, and PLD predicts a static scalar force, not a
    gravitomagnetic one. I found no published Cavendish-type measurement between two He-II masses. The PLD lab
    channel is untested, not refuted.

## 3. Renaming check

- **The skeleton is DEF/Brans–Dicke.** One massless conformally coupled scalar with ω = 0 (α₀² = 1/3, f(R)), per-body
  effective charges α_A, and dipole radiation ∝ (α_p − α_c)². The dipole and G_AB formulas are DEF's verbatim.
- **The "feedback" is a mean-field cousin of DEF spontaneous scalarisation.** The NS's own scalar field raises its
  charge by e^{K′/N(0)}. A says this itself.
- **The selective charge belongs to the "dense-object-only scalar charge" class.** Hook & Huang 2018 (QCD axion sourced
  only by NSs because of finite-density effects, "evades fifth-force constraints" because Earth and Sun carry none) is
  the same phenomenological design, as are its binary-pulsar and GW bounds (Poddar et al. 2020; Zhang et al. 2021;
  Huang et al. 2018). The composition-selective coupling is the Damour–Polyakov non-universal-dilaton idea, and
  condensate-mediated forces appear in superfluid-DM models (Berezhiani & Khoury 2015); phase-as-time is Volovik's
  superfluid-vacuum programme.

## 4. Novelty verdict

**Low.** The only new element is the *source term*: charge ∝ condensate fraction f_c ≈ (3π/8)Δ/E_F, which gives an
s(M) template. That source has no local-operator realisation (§1.5), it rests on a normalisation slip (§1.1), and its
amplitude is uncertain by ×2–3 from gap models, the trace, the crust regime and the feedback exponent. The template is
therefore not sharp enough to be distinguished from a generic α(M), which is K2's discriminant. Phenomenologically
PLD is "DEF with α₀ = 0 for all non-NS matter", a case the existing per-mass α_A bounds already cover.

## 5. Triage

| # | Prediction | Verdict | Reason |
|---|---|---|---|
| P1 | NS–WD dipole decay, α_NS = 2.3–9×10⁻³ | **FAIL** | Excluded: >6σ for s ≥ 8e-3, about 2σ at the minimum (§2) |
| P2 | Stacked Ṗ_b residual regressed on s(M) (K1/K2) | **Conditional PASS** (bound-setting only) | Cheap and data are public, but the template is soft (×2–3), so it can only bound the surviving corner |
| P3 | α_WD = α_BH = 0 exactly | **FAIL** (as a discriminant) | Same as DEF with α₀ → 0 plus BH no-hair |
| P4 | NS–BH GW −1PN dipole (K4) | **FAIL** | Current limit α̂² ≲ 8e-4, about 150× above the target |
| P5 | He-II–He-II attraction excess ≈ 1.6–2.7×10⁻³ (A says 3–5e-3) tracking n₀(T)², off at T_λ | **Conditional PASS** | Sharp, never done as far as found, and a single experiment decides it; conditional on the non-local A4. Not this project's domain (no archive data) |
| P6 | NS gaps enhanced ×1.4–2 (cooling, Cas A) | **FAIL** | Degenerate with the ×3 nuclear-theory gap uncertainty |
| P7 | J0337+1715 triple and EP tests null | **FAIL** (no test) | Consistency statement with zero discriminating power |
| P8 | DM must be incoherent | **FAIL** (liability) | A condensate DM would get +1/3 self-gravity. Sgr streams already limit DM–baryon acceleration differences to about 10% (Kesden & Kamionkowski 2006), so PLD forbids BEC/fuzzy DM ad hoc |

**P2 in detail (what this project could add).**
- **Data.** Published Ṗ_b^xs and masses for J1738, J1012, J1909, J2222 and J0348 (Zhao 2022, tables 1–2) and the J1738
  update (Vaglio 2026). Raw TOAs and timing models are in NANOGrav 15-yr and EPTA DR2.
  ATNF psrcat supplies P_b, PBDOT, PX and PM, but **not intrinsic Ṗ_b**, so you must use the published kinematic
  corrections.
- **Mimics.**
  - Shklovskii term: 9.3 ± 0.6 fs s⁻¹ in J1738, about 3× its σ; it depends on parallax.
  - Galactic-potential model.
  - WD mass from atmosphere models, which sets Ṗ_b^GR (the 1.5–1.9 fs s⁻¹ error in J1738).
  - Ġ/G, red noise, an unseen third body.
  - Any universal α(M).
- **Sensitivity today.** σ(α²) ≈ 2.7×10⁻⁶ combined, about 90% from J1738. Reaching 5σ on the minimum needs
  σ ≈ 1×10⁻⁶, i.e. about ×3. That is plausible in about 5–10 years only if the J1738 Shklovskii and GR-mass errors
  shrink along with timing. A's T^(−5/2) scaling ignores those floors.
- **Novel?** Strictly, nobody has fitted a condensate-fraction template. But mass-resolved combined α(M) fits exist
  (Shao 2017, Zhao 2022, Batrakov 2023), and the per-mass table in §2 *is* that regression to the precision the
  template allows. **Recommendation: do not open project work on it.** Record the §2 table as the PLD bound
  (α_NS < 2.2×10⁻³ at 95% for a mass-independent α, 1.4–2.0 M_⊙) and close R3-B unless someone runs P5 in a lab.

## References (all opened and verified)

- Freire, Wex, Esposito-Farèse, Verbiest et al. 2012, MNRAS 423, 3328, "PSR J1738+0333 II", arXiv:1205.1450
- Antoniadis et al. 2013, "A massive pulsar in a compact relativistic binary", arXiv:1304.6875
- Zhao, Freire, Kramer et al. 2022, CQG 39, 11LT01, "Closing a spontaneous-scalarization window with binary pulsars", arXiv:2201.03771
- Shao, Sennett, Buonanno et al. 2017, PRX 7, 041025, arXiv:1704.07561
- Batrakov, Hu, Wex et al. 2023, "A new pulsar timing model for scalar-tensor gravity… PSR J2222−0137", arXiv:2303.03824
- Guo et al. 2021, "PSR J2222−0137 I. Improved physical parameters", arXiv:2107.09474
- Vaglio, Carleo, Susobhanan et al. 2026, "Constraints on Einstein-aether gravity from … PSR J1738+0333" (preprint; Ṗ_b, Ṗ_b^Shk, Ṗ_b^Gal used), arXiv:2605.01436
- Miao, Freire, Wex et al. 2026, A&A 713, A1, "Improved proper motion and gravity tests with PSR J1913+1102", arXiv:2606.19276
- Kramer et al. 2021, "Strong-field gravity tests with the Double Pulsar", arXiv:2112.06795
- Bhat, Bailes & Verbiest 2008, "Gravitational-radiation losses from … PSR J1141−6545", arXiv:0804.0956
- Bussieres, Caldarola & Nesseris 2025, "Updated constraints on modified gravity from binary pulsars" (preprint; its κ_D results look degeneracy-limited and are not used), arXiv:2507.18188
- Santos, Nunes & de Araujo 2025, PRD 111, 084087, "Constraining scalar charge … NSBH mergers", arXiv:2504.00782
- Takeda, Tsujikawa & Nishizawa 2024, PRD 109, 104072, "GW constraints on scalar-tensor gravity from … GW200115", arXiv:2311.09281
- Touboul et al. 2022, "MICROSCOPE mission: final results", arXiv:2209.15487
- Damour & Esposito-Farèse 1996, PRD 54, 1474, "Tensor-scalar gravity and binary-pulsar experiments", arXiv:gr-qc/9602056
- Damour & Polyakov 1994, "The string dilaton and a least coupling principle", arXiv:hep-th/9401069
- Sotiriou & Faraoni 2010, RMP 82, 451, "f(R) theories of gravity", arXiv:0805.1726
- Hook & Huang 2018, "Probing axions with neutron star inspirals and other stellar processes", arXiv:1708.08464
- Huang, Johnson, Sagunski et al. 2018, "Prospects for axion searches with Advanced LIGO through binary mergers", arXiv:1807.02133
- Poddar, Mohanty & Jana 2020, PRD 101, 083007, arXiv:1906.00666
- Zhang et al. 2021, "First constraints on nuclear coupling of axionlike particles from … GW170817", arXiv:2105.13963
- Berezhiani & Khoury 2015, "Theory of dark matter superfluidity", arXiv:1507.01019
- Volovik 2001, "Superfluid analogies of cosmological phenomena", arXiv:gr-qc/0005091
- Salasnich, Manini & Parola 2005, PRA 72, 023621, "Condensate fraction of a Fermi gas in the BCS-BEC crossover", arXiv:cond-mat/0506074
- Salasnich 2014, J. Phys. Conf. Ser. 497, 012026, "Fermionic condensation in … neutron stars", arXiv:1308.0922
- Diallo et al. 2012, PRB 85, 140505, "BEC in liquid ⁴He near the liquid-solid line" (gives n₀ = 7.25 ± 0.75% at SVP), arXiv:1111.2284
- Kesden & Kamionkowski 2006, "Galilean equivalence for galactic dark matter", arXiv:astro-ph/0606566
- Tajmar et al. 2006, "Experimental detection of the gravitomagnetic London moment", arXiv:gr-qc/0603033; Tajmar et al. 2007, "Search for frame-dragging-like signals…", arXiv:0707.3806; Tajmar et al. 2008, arXiv:0806.2271; Tajmar & Plesescu 2009, arXiv:0911.1033
- Tajmar et al. 2004, "Weight measurements of high-temperature superconductors during phase transition", arXiv:gr-qc/0404005
- Manchester et al. 2005, "The ATNF Pulsar Catalogue", arXiv:astro-ph/0412641; NANOGrav 15-yr data set, arXiv:2306.16217; EPTA DR2, arXiv:2306.16224
- Not opened, cited second-hand only: Graham et al. 2008 (Canterbury ring-laser null on a spinning Pb disc), via arXiv:0707.3806.
