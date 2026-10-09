# R3-B, System A: "Phase-Locked Dilaton" (PLD)

Everything below is a **hypothesis**. Numerical inputs marked *(recalled)* come from memory with no
network access and must be checked before anyone uses them. The derivations are order-of-magnitude
unless stated otherwise.

## 0. What the seed forces: a derived no-go

The seed asks for a collective state that changes how geometry responds, while energy conservation
and the Bianchi identity still hold. The most direct version has no new degrees of freedom: the
substrate's curvature stiffness depends on the local order parameter, `S = ∫√-g [F(ψ) R / 16πG + L_m]`
with `F = 1 + φ(order)`. Working through it gives two results.

1. **At O(φ) the effect is equivalent to GR.** `R` is slaved to the local trace (`R = 8πG(ρ−3p)`). So
   `φR/16πG` is a local energy density `φ(ρ−3p)/2`, and that is not suppressed by G. In the Einstein
   frame it is only a shift in the matter equation of state. If φ is chosen to "store the condensation
   energy in geometry" (`φ = 2ε_ord/(ρ−3p)`), the result cannot be told apart from standard BCS
   thermodynamics plus GR. This is round-2 failure mode 6.
2. **At O(φ²) laboratory superfluids exclude it.** In the Einstein frame there is a gradient term
   `(3/16πG)(∇φ)²`. For He-II, `φ ≈ ε_ord/ρc² ≈ 5×10⁻¹⁴`. A wall where φ jumps over the healing length
   (about 1 Å) would cost about 10²⁵ J m⁻², while the condensation energy available per unit area is about
   10⁻⁴ J m⁻². Balancing the two forces the wall to spread over w ~ M_P φ/√ε_ord ~ 10⁵ m. Superfluid He
   exists in centimetre-sized containers, so this is false. Any local, non-Planck-suppressed F(order) is
   excluded the same way.

So a "response material" without a new mode is either invisible or already excluded. The surviving
structure is **a different coupling of an existing source through a second substrate response
channel**. It is built below.

## 1. Axioms

- **A1 (substrate).** The primitive is a substrate of relational events with no background space, time
  or metric. Observables are the substrate's linear-response (correlation) kernels to small
  perturbations.
- **A2 (single quadratic form).** At long wavelength every response kernel depends on separations only
  through one quadratic form η. That quadratic form *is* the emergent metric, and the substrate supplies
  one stiffness-to-inertia ratio, which is c.
- **A3 (two massless channels).** At long wavelength the substrate's response to stress has exactly two
  gapless channels:
  - **shear (spin-2):** it couples to the full energy-momentum T_μν of every excitation. This is the
    usual metric response.
  - **dilatation (spin-0, σ):** it is the substrate's conformal or phase mode. The time-like direction is
    the advance of the substrate's collective phase, so σ is the local rate of that phase.
- **A4 (phase locking).** An excitation couples to σ only if it carries a *macroscopic phase of its own*
  that can lock to the substrate phase. In the language of matter physics that means off-diagonal
  long-range order (ODLRO). The source is the trace of the condensed component only:
  `J_σ = −β √(8πG) f_c T^μ_μ`. Here f_c is the condensate fraction, the weight of the long-range part
  of the one- or two-body density matrix. Incoherent matter, normal Fermi liquids, crystals and
  radiation carry no σ-charge. Coherent radiation such as laser or maser light has a traceless T, so it
  carries none either.
- **A5 (isotropic stiffness).** The substrate's dilatational stiffness equals its shear stiffness, a
  Cauchy-type relation for a medium with a single elastic constant. With canonical normalisation of the
  conformal mode this fixes **β² = 1/3**. That is the value a conformal mode gets with no extra kinetic
  weight (ω = 0). The alternative counting in which only the trace couples gives 1/6, so the theory
  carries a factor-2 uncertainty here.

**Resemblances, listed honestly.** A3 and A5 have the structure of Brans–Dicke with ω = 0, which is
also f(R)'s α² = 1/3. The self-consistent feedback in §4 resembles Damour–Esposito-Farèse
scalarisation. The composition-selective charge resembles fifth forces coupled to B or L. "Superfluid
vacuum" programmes share the phase-as-time picture. The difference from all of these: the scalar
couples **only to ODLRO matter**, so the Sun, the Earth, white dwarfs and black holes have zero
charge.

## 2. Lorentz invariance and the bounds

- By A2 the two channels and all matter propagate on one light cone, so c_gw = c_σ = c_γ at long
  wavelength. That satisfies GW170817 (10⁻¹⁵). Violations come only from higher-gradient terms
  suppressed by (k ℓ_sub)². With ℓ_sub at or below the Planck length that is about 10⁻⁴⁰ at GeV, well
  below lab bounds (10⁻¹⁸).
- The action is generally covariant: the tensor term, the canonical σ, and J_σ built from a scalar
  density. There is no preferred frame, and both total T_μν conservation and the Bianchi identity hold
  identically.
- **Honest caveat.** "Derived" means derived from A2. A2 is a symmetric primitive, which is exactly
  lesson 5's warning. The risk sector it creates (§5) is *not* already excluded at O(1). The predicted
  charge is about 10⁻³, at the edge of current bounds.

## 3. Recovering GR and the equivalence-principle tests

- **Ordinary matter** (f_c = 0) couples only to the shear channel, which is exactly GR. Solar-system
  PPN parameters are γ = β_PPN = 1 because the Sun has no charge. That covers Cassini, perihelia and
  lunar laser ranging.
- **Eöt-Wash and MICROSCOPE.** A σ force needs a charge on *both* bodies. The Earth has no ODLRO matter,
  so the anomalous force on any test mass is zero, superconducting masses included.
- **The pulsar–WD triple (J0337+1715).** The anomalous differential acceleration is proportional to
  α_NS·α_outerWD = α_NS·0 = 0. The strongest strong-EP test is **exactly blind** to this hypothesis.
  That is a derived statement, not a tuned one.
- **Double neutron stars.** G_AB = G(1 + 2α_Aα_B), which is about 10⁻⁵ (see §4). It is absorbed into the
  masses except for a Shapiro-versus-Keplerian mismatch of order 10⁻⁵, below current precision.
- **Gravitational-wave speed and polarisation.** The tensor sector is unchanged. Scalar emission in
  binary neutron stars goes as (α_A − α_B)², which is negligible.

Exterior predictions differ from GR only where **two phase-locked bodies** interact, or where one such
body radiates dipole σ-waves.

## 4. The derived effect

### 4.1 Neutron-star charge

The σ-charge of a neutron star relative to its mass is `α_NS = β s`, with
`s = (1/M) Σ_i ∫ f_c,i ρ_i dV`.

For weak-coupling BCS, the fraction of particles in the pair condensate is
**f_c = (3π/8) Δ/E_F** *(recalled standard result)*. That is a ratio of two microscopic scales.
Using x_p ≈ 0.05–0.1, E_F,p ≈ 20 MeV and E_F,n ≈ 80 MeV at about 2n₀:

| Component | Gap Δ (MeV) | f_c | Mass weight | Contribution to s |
|---|---|---|---|---|
| Core protons, ¹S₀ | 0.5–1 | 0.03–0.06 | × x_p × region fraction | (0.7–6)×10⁻³ |
| Core neutrons, ³P₂ | ≈ 0.05–0.1 | ≈ 10⁻³ | × ~0.9 | ≈ 10⁻³ |
| Inner-crust neutrons, ¹S₀ | 1.5–2 | ≈ 0.1 | × ~0.01 | ≈ 10⁻³ |

This gives **s₀ ≈ 2–8×10⁻³** before feedback.

### 4.2 Self-consistent feedback (derived, finite)

The σ self-energy, roughly −(6/5)β²Gρ₀²V/R, is quadratic in Δ. It renormalises the BCS functional,
`F = N(0)Δ²[ln(Δ/Δ₀) − ½] − K′Δ²`. Minimising gives `Δ = Δ₀ exp(K′/N(0))`, which is enhancement and
not runaway. Here `K′/N(0) ≈ 1.1 β² C x_i (m c²/E_F,i)`, where C = GM/Rc² ≈ 0.2. That evaluates to
about 0.35 for protons and about 0.75 for neutrons. Gaps grow by a factor of 1.4–2, so
**s ≈ 4×10⁻³ to 1.6×10⁻²**, which gives

**α_NS = s/√3 ≈ 2.3×10⁻³ to 9×10⁻³**, and the white-dwarf charge α_WD is exactly 0.

Side prediction: neutron-star pairing gaps exceed nuclear-many-body values by
e^{K′/N(0)}. This bears on cooling, for example Cas A.

### 4.3 Observable: dipole σ-radiation in NS–WD binaries

At leading order *(standard dipole form; coefficient to be re-verified)*:

`Ṗ_b^σ = −(4π²/P_b) T_⊙ (m_p m_c / M) α_NS²`, with T_⊙ = GM_⊙/c³.

Relative to the GR quadrupole term this is `0.0153 (T_⊙M/P_b)^(−2/3) α²`. For J1738+0333 that ratio is
about 3.7×10⁴ α², so a predicted α² of 5×10⁻⁶ to 8×10⁻⁵ means a **20%–300% excess in orbital decay**.

The regression variable is new. The template is not a universal α(M) as in standard scalar-tensor
fits. It is **α_NS ∝ s(M)** from gap models, which is non-monotonic in mass: the proton ¹S₀ gap closes
at high density while the ³P₂ region grows. It is also exactly zero for the companion.

## 5. Dataset, N, 5σ N and the kill criterion

**Quantity.** Intrinsic Ṗ_b residuals of NS–WD binaries: observed minus Shklovskii minus Galactic
acceleration minus the GR quadrupole. Sources are the IPTA, NANOGrav, EPTA and PPTA data releases,
the ATNF catalogue, and MeerKAT and FAST timing. Usable N today is roughly 10–30 systems with
Ṗ_b significance above 3σ *(estimate; count to be made)*. As far as I know, nobody has stacked these on
a condensate-fraction template.

**Per-system sensitivity.** `σ(α²) = σ(Ṗ_b)/D` with `D = 4π² T_⊙ (m_p m_c/M)/P_b`.

| System | D | σ(Ṗ_b) *(recalled)* | σ(α²) |
|---|---|---|---|
| J1738+0333 (m_p ≈ 1.46, m_c ≈ 0.18, P_b ≈ 0.355 d) | 1.0×10⁻⁹ | ≈ 3.7×10⁻¹⁵ | ≈ 3.6×10⁻⁶ |
| J0348+0432 (2.01 M_⊙) | 3.5×10⁻⁹ | — | ≈ 1.3×10⁻⁵ |

**5σ N.** The per-system significance in J1738-equivalents is z = α²/3.6×10⁻⁶, and
N₅σ = (5/z)².

| s | α² | z | N₅σ |
|---|---|---|---|
| 4×10⁻³ (lowest gaps with feedback) | 5.3×10⁻⁶ | 1.5 | ≈ 11 |
| 8×10⁻³ | 2.1×10⁻⁵ | 5.9 | < 1 (J1738 alone decides) |
| 1.6×10⁻² | 8.5×10⁻⁵ | 24 | already decisive |
| 2×10⁻³ (no feedback, minimum gaps; not self-consistent under PLD) | 1.3×10⁻⁶ | 0.37 | ≈ 180 |

Ṗ_b precision improves roughly as T_span^(−5/2), so ten more years of timing on the existing systems
gives about a 5× gain.

**Laboratory cross-check (N = 1, decisive).** Two He-II masses should attract with an excess
`2β² f_c² ≈ (2/3)(0.07–0.09)² ≈ 3–5×10⁻³` of the Newtonian value. The excess should switch off at T_λ
and track f_c(T)². Earth's own field shows no change.

**Kill criteria**, any one of which ends the hypothesis:

- **K1.** The stacked NS–WD Ṗ_b residual gives a template amplitude α² below 5×10⁻⁶ (the
  self-consistent minimum) at 95% CL. If J1738's published Ṗ_b agrees with GR to about ±15%
  *(recalled; verify)*, then α² is already below about 6×10⁻⁶, and **most of the PLD range is excluded
  now**. Only the lowest-gap corner survives.
- **K2.** A non-zero excess appears but does not follow s(M). For example, J0348 (2.0 M_⊙) and J1738
  (1.46 M_⊙) show the ratio of a universal α(M) instead of the gap-model ratio. Or a white-dwarf-only
  system shows any σ signature.
- **K3.** A He-II–He-II attraction experiment finds less than 10⁻³ excess, with a T_λ on/off
  modulation.
- **K4.** A neutron-star–black-hole gravitational-wave stack bounds the −1PN dipole term below
  α_NS² ≈ 5×10⁻⁶. Black holes carry no charge, so this is a clean test.

## 6. Systematics that could mimic the signal

- **Kinematics.** The Shklovskii term μ²d/c and Galactic-acceleration corrections, which depend on the
  potential model and on parallax versus DM distance. Their P_b scaling differs from the dipole's 1/P_b,
  so a regression on P_b helps separate them.
- **Mass errors.** Errors in the white-dwarf mass from atmosphere models propagate into the GR
  prediction.
- **Other orbital effects.** Unmodelled companions, mass loss or tides (exclude black-widow and redback
  systems), red timing noise, and solar-system ephemeris errors.
- **Theory systematic.** The gap models: Δ(n) is uncertain by about a factor of 3 and E_F,p by the
  equation of state. This shifts s(M), and it dominates.
- **Look-alike physics.** Any other dipole mechanism, such as a universal scalar with α(M), or a
  massive-graviton or extra-field model. These are distinguished only by the s(M) shape and by
  α_WD = α_BH = 0 exactly.

## 7. My three weakest points

1. **A4 couples to a non-local quantity.** ODLRO, meaning the long-range part of the density matrix, is
   not a local operator. Normal matter has short-range pair fluctuations. If those leak any σ-charge q_n,
   torsion-balance and MICROSCOPE tests require roughly 2β²Δq_n² below 10⁻¹⁵ in the relevant
   combination, so q_n must be about 10⁻⁸ or smaller with no mechanism given. I asserted the
   long-wavelength filter. I did not derive it.
2. **β² = 1/3 depends on A5, and A2 assumes Lorentz invariance.** The Cauchy-type equality of stiffnesses
   is posited, and the trace-only counting would halve α². A2 puts Lorentz invariance in through one
   quadratic form rather than deriving it from anything deeper (lesson 5). The feedback estimate
   e^{K′/N(0)} uses a uniform-sphere, mean-field approximation, and an O(1) error in the exponent moves
   s by a factor of 2.
3. **It is nearly dead already, and it has a cosmological hostage.** Typical gap values are already
   excluded by J1738 if my recalled Ṗ_b number holds. Only the minimum-gap corner survives. Separately,
   if dark matter is a coherent condensate (axion-like or ultralight), PLD makes dark-matter
   self-gravity stronger by 1 + 2β²f_c, up to 5/3. Structure-formation and tidal-stream bounds would
   reject that. So PLD requires dark matter to be incoherent, a requirement the hypothesis does not
   explain.

**Bottom line.** The seed's literal form, local F(order)·R, is ruled out by He-II walls or is
equivalent to GR. The surviving version is a σ-charge carried only by ODLRO matter, with an amplitude
derived from Δ/E_F and β² = 1/3. It predicts NS–WD dipole orbital decay at α_NS ≈ 2–9×10⁻³, and a
J1738-anchored stack decides it with about 11 equivalent systems or fewer.
