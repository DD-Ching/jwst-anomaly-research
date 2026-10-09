# R2-C — "Two-Slot World": a minimal-axiom hypothesis (System A, round 2)

Status labels: **[A]** axiom, **[D]** derived inside the model, **[I]** interpretive identification (not derived),
**[O]** observed input, **[H]** hypothesis/prediction. No web or repository material was used; numerical inputs
quoted from memory are marked "(approx., from memory)" and must be checked before use.

## 0. Stance

Ask what is the least structure a world can have and still contain relations and interference. My answer is this.
There is a finite collection of things, and each thing is a pair of complex numbers. The only structure is the
complex "area" between two such pairs. The rest is meant to be derived.

## 1. Axioms (four, plus one integer)

- **A1 (Carrier).** There is a finite set S with |S| = N. Each element i carries a non-zero vector u_i ∈ C².
- **A2 (Relation).** C² carries exactly one structure: the antisymmetric bilinear form
  ε(u,v) = u¹v² − u²v¹. Physical content is invariant under every map that preserves ε, so observables are functions
  of the relations ε_ij = ε(u_i,u_j) alone.
- **A3 (Amplitude).** The state of the world is one function Ψ of the ε_ij, and Ψ is **holomorphic**, meaning no
  complex conjugates appear in it. Relative frequencies are given by |Ψ|² (normalized). Complex conjugation enters
  only through this Born square.
- **A4 (Indifference).** Elements carry no labels. Ψ is the lowest-degree non-zero function that is symmetric under
  permutations of S and allowed by A2–A3. No further numbers are specified.

**Information content.** A1–A4 can be written in roughly 200 bits of formal notation. The one free datum is the integer
N, about log₂N ≈ 409 bits if N ≈ 5×10¹²². The dimensionful unit is a convention. N is the only dimensionless input.

Why C² over C: R² with ε gives SL(2,R) ≅ SO(2,1), a 2+1 world. Quaternions give 5+1 and octonions give 9+1. C is the
smallest algebraically closed field, and A3 needs a continuous interference phase. I treat this as a minimality
argument, not a theorem.

## 2. Lorentz invariance, derived

**[D] The symmetry group is exact.** The maps of C² that preserve ε are the 2×2 complex matrices with det M = 1, the
group SL(2,C), because ε(Mu,Mv) = det M · ε(u,v). This holds by construction, not approximately. A1–A4 contain no
other structure that could break it.

**[D] Spacetime emerges as a Born bilinear.** The Born square makes conjugates available, so each element defines
p_i = u_i u_i†, a 2×2 Hermitian matrix. Write X = t·1 + x·σ. The space Herm(2) is R⁴, det X = t² − |x|² is a
Minkowski form, and M acts as X ↦ M X M†, which preserves det. So SL(2,C)/{±1} ≅ SO⁺(1,3). Dimension 3+1 and Lorentzian
signature both follow from "pairwise complex relation on C²".

**[D] Every elementary object is exactly null.** det(u u†) = 0 is an algebraic identity: a rank-1 matrix has zero
determinant. It is not a dispersion relation with corrections. No power series in E/E_Planck exists that could deform
it, because A1–A4 contain no length scale.

**[D] Mass is relational.** By the Lagrange identity, det(u u† + v v†) = |ε(u,v)|². Therefore

  M² (any composite) = Σ_{i<j} |ε_ij|².

A composite is massless only if every pair has ε_ij = 0, which means all its constituents are collinear. Massless
composites of any kind (photon-like, graviton-like) therefore move exactly on the null cone of their constituents.
This gives c_gw = c_γ = c_ν(massless) **identically**, a zero rather than "below 10⁻¹⁵".

**[D] A speed limit exists.** Positive-semidefinite Hermitian matrices form a convex cone. A sum of p_i stays inside
it, so no composite can be spacelike or superluminal. The speed limit is convexity.

**[D] Time orientation.** u u† is never negative-definite, so no past-pointing vector can be built from elements. This
gives a microscopic time orientation. It does not give the thermodynamic arrow (see §6).

**[D] Gravity is universal.** By Schur's lemma, the only SL(2,C)-equivariant real quadratic map C² → Herm(2) is
u ↦ c·uu†. A4 makes c the same for every element. Any universal, conserved, vector-valued "charge" must therefore be
P = c Σ u_i u_i†. Every system couples to the same P with the same coefficient, which is the equivalence principle
without a separate postulate.

**How the finite-valency no-go is evaded.** The no-go applies to discrete point sets embedded in Minkowski space with a
neighbour relation, because those pick a frame. Here the elements are not spacetime points. They are null momenta: a
point on the celestial sphere CP¹ times a scale. They have no neighbour graph. Every pair has a relation ε_ij, but it
is a Lorentz scalar (an invariant mass) and not a link in position space. Discreteness lives in the **count of
momenta**, not in positions, so a boost has nothing to distort.

**Why this is not "disordered locality" (lesson 2).** No metric is defined from ε-links. There is no graph distance,
so no distance can collapse. Positions are a conjugate variable. In the holomorphic representation of A3, ω_i ≡ ∂/∂u_i
is conjugate to u_i. A spacetime point x ∈ Herm(2) is the Lorentz-covariant incidence ω = i x u. Ψ depends only on the
ε_ij, which are invariant under ω → ω + i a u, so translations are an exact symmetry.

**Bounds.** The 10⁻¹⁵ GW170817 bound and the 10⁻¹⁸ laboratory bounds are satisfied with margin "∞". The model predicts
exactly zero vacuum dispersion, zero anisotropy of c, zero speed difference between massless messengers, and zero extra
GW polarizations. Finite N introduces no frame into the *laws*. The cosmic frame P_total is a property of the *state*,
just as in FRW cosmology.

## 3. Inverse square, distances, conservation

**[D] Conservation laws.**
- Six Lorentz charges come from SL(2,C) invariance (Noether).
- Four-momentum comes from translation invariance (ω-shift).
- Each element's helicity, −½·(degree of Ψ in u_i), is conserved because Ψ is homogeneous in each u_i.
- **Element number N is conserved.** This is new: it is a count conservation underlying all the others.

**[D] Inverse square.** Element number is conserved along a ray. Directions of null elements fill CP¹, which has real
dimension 2. The emergent incidence geometry gives a sphere at affine distance r an area of 4πr². Flux is therefore
F = L/4πr². The exponent 2 is dim_R CP¹ = 2·dim_C(C²) − 2, so it is exact and tied to the same fact that fixes 3+1. A
"leakage" exponent is zero, so in the GW-damping test h ∝ 1/d_L with γ_extra = 0.

**[D, given an expansion history] Four distances.** Define 1+z as the ratio of a free element's energy to a bound
composite's standard energy. Bound composites are fixed sets of ε_ij; free elements are not. With a homogeneous state
(A4) and spatial flatness:
- D_C = c∫dz/H
- D_M = D_C
- D_A = D_M/(1+z)
- D_L = (1+z)D_M
- light-travel distance c∫dz/[(1+z)H]

The **Etherington relation D_L = (1+z)² D_A holds exactly**. It follows from element-number conservation together with
Liouville's theorem, because SL(2,C) preserves the holomorphic area form, so phase-space density along rays is conserved.
There is no opacity term at all, which satisfies lesson 3. H(z) itself is **not** derived (§6).

## 4. Where this world differs from GR + ΛCDM

The axioms reproduce relativistic kinematics exactly and have no length scale. The model therefore cannot generate
Planck-suppressed deviations. Any difference must come from **structural selections**: things ΛCDM leaves free that the
axioms fix. I found two of these, and one of them carries an amplitude that is fully derived.

### 4a. The lightest neutrino is exactly massless [D + I]

A single element is a massless helicity −½ object (§2). It is a left-handed Weyl fermion with no right-handed partner,
because A3's holomorphy gives every element the same chirality. Every massive particle is a composite of at least two
elements. A spin-½ massive particle needs an odd count, so at least three.

**[I]** The bare element, if it exists unconfined, is a massless, purely left-handed, electrically neutral fermion.
That is exactly the description of a neutrino mass eigenstate with m = 0. The identification also "explains" the least
explained regularity I chose, **maximal parity violation**: only bare elements are chiral. Photon- and graviton-like
objects are Born bilinears (u_A ū_A′, u_A u_B ū_A′ ū_B′) and are parity-symmetric. Parity violation is therefore
confined to fermion couplings, which is what is observed.

**Derived amplitude.** m_lightest = 0 exactly. Combined with the measured oscillation splittings [O]
(Δm²₂₁ ≈ 7.4×10⁻⁵ eV², |Δm²₃₁| ≈ 2.5×10⁻³ eV², sin²θ₁₂ ≈ 0.30, sin²θ₁₃ ≈ 0.022; approx., from memory):

| Quantity | Normal ordering (m₁ = 0) | Inverted ordering (m₃ = 0) |
|---|---|---|
| Σm_ν | **58.8 meV** | **≈ 99.9 meV** |
| m_β (β-decay endpoint) | 8.8 meV | ≈ 49 meV |
| m_ββ (if Majorana) | 1.5–3.7 meV | ≈ 19–49 meV |

There are no adjustable parameters. Standard 3ν + ΛCDM leaves m_lightest free in [0, ~20 meV]. The Planck "baseline"
0.06 eV is a convention, not a prediction.

### 4b. Λ fixed by the count, with w = −1 exactly [D given one calibration]

The only dimensionless datum is N. If the de Sitter horizon's capacity is identified with the element count (one
two-state slot each, S_dS = N ln 2), then Λ = 3π/(N ln2 · ℓ_P²). Observed Λℓ_P² ≈ 2.9×10⁻¹²² gives N ≈ 4.7×10¹²². This
is **calibration, not prediction.** The consequences are:
- N is conserved, so **w = −1 exactly**.
- Ω_Λ(z) = S_H(z)/S_dS, which is an identity in this model. The "coincidence problem" becomes the statement that the
  Hubble horizon now spans ~70 % of all slots.
- G and Λ are not independent constants.

### 4c. Variants I tried and killed

1. A global C* "expansion flow" acting on every element (the only operation that commutes with SL(2,C)) leaves every
   dimensionless ratio unchanged, so redshift would be unobservable. Rejected.
2. Adding a phase to that flow with "one topological winding over the ½lnN ≈ 141 e-folds of history" gives a
   wavelength-independent CMB polarization rotation β ≈ 9°. Observed |β| ≲ 0.5°, so rejected. Minimality then forces
   the phase to zero, which gives **β = 0 exactly**.
3. Giving all elements equal energy in the cosmic frame (a "less information" option) forces redshifted photon energy
   into a bath that scales as ρ ∝ a⁻³/a_initial. That diverges. Rejected.

## 5. Predictions and kill criteria

**P1 (main; existing public data).** Σm_ν = 58.8 ± 0.5 meV (NO) or ≈ 100 meV (IO). The ± comes only from the
oscillation-parameter uncertainty.
- (a) No ordinary astrophysical process makes a mass eigenvalue *exactly* zero. Systematics can shift a cosmological
  Σm_ν estimate, but they cannot mimic a sharp parameter-free value that is cross-checked by β-decay and 0νββ.
- (b) Testable now with public data: Planck PR4 / ACT DR6 CMB lensing and DESI DR2 BAO likelihoods (the current
  ΛCDM 95 % upper limit is ~64 meV; approx., from memory), KATRIN, KamLAND-Zen and LEGEND.
- (c) **Kill if any of these happen:**
  - A cosmological analysis with a systematics budget ≤ 10 meV excludes 58.8 meV at ≥ 3σ.
  - IO is established and Σm_ν ≠ 99.9 meV at ≥ 3σ.
  - A β-decay experiment measures m_β > 20 meV (NO) or > 60 meV (IO).
  - m_ββ > 4 meV with NO established.
  - A light right-handed (sterile) neutrino is shown to be elementary.

  Existing data already disfavour IO, so the model effectively commits to **NO, 58.8 meV**. The persistent cosmological
  preference for Σm_ν → 0 or "negative" is the live threat.

**P2 (exact zeros; existing data).** All of the following are predicted to be zero:
- vacuum dispersion at every order;
- c_gw − c_γ;
- non-tensor GW polarizations;
- GW-damping leakage (h ∝ 1/d_L exactly);
- Etherington violation η(z) − 1;
- isotropic cosmic birefringence β;
- dark-energy evolution (w₀ = −1, w_a = 0).

ΛCDM is silent or extendable on several of these; this model is not. **Kill if any one of them is established at ≥ 5σ
with independent calibration.** Two candidates are already live: the reported β ≈ 0.3° (approx., from memory) and the
DESI DR2 w₀w_a preference. The model says both must be systematics.

**P3 (never measured).** Two quantities:
- **The mass-eigenstate-resolved neutrino time of flight.** The ν₁ component of a Galactic supernova burst should show
  exactly zero delay relative to the gravitational-wave or optical onset. Flavour-averaged SN1987A timing never
  isolated it.
- **The late-time equation of state of the relic ν₁ background.** It should stay exactly radiation-like (w = 1/3) today,
  not become non-relativistic as a ~1 meV neutrino would. No survey measures this component separately.

## 6. Three weakest points

1. **Dynamics is not derived; only kinematics is.** I derive the Lorentz group, 3+1, the null cone, the equivalence
   principle, the inverse-square law and Etherington. I do **not** derive:
   - Einstein's equations or any H(z);
   - spatial flatness;
   - G separately from Λ;
   - α or any mass ratio;
   - how localized wave packets and local interactions emerge from the holomorphic incidence ω = i x u.

   The claim "everything follows" is far from earned. The model explains *why the stage is Minkowskian*, not what plays
   on it.
2. **The neutrino identification [I] is the load-bearing step and is interpretive.** If bare elements are confined, or
   if antiparticles do not arise from the Born square as claimed, P1 loses its derivation. CPT is asserted via the
   Born square's conjugation symmetry, not proved. P1 is also *compatible* with standard physics (m₁ = 0 is allowed),
   so a confirmation would be weak evidence, while a refutation would be decisive.
3. **The Λ–N relation borrows horizon-entropy counting.** "One slot per element" is not derived from A1–A4, so Λ is
   calibrated rather than predicted. The thermodynamic arrow and the low-entropy beginning are **not explained**.
   Positivity of every p_i gives a time *orientation*, not an entropy gradient. Strict w = −1 may already be under
   pressure.

## Resemblances I noticed

- von Weizsäcker's "ur-alternatives" (binary alternatives as C² → 3+1): the closest, and nearly identical in A1.
- Penrose's two-spinor calculus and twistor incidence (ω = i x u).
- The spinor-helicity formalism of scattering amplitudes (|ε_ij|² = 2p_i·p_j).
- The division-algebra ↔ spacetime-dimension correspondence (R, C, H, O → 3, 4, 6, 10).
- Celestial-sphere Möbius structure (celestial-holography-adjacent).
- Self-dual / holomorphic gravity variables (A3's holomorphy).
- de Sitter horizon entropy counting (§4b).
- Relational/Page–Wootters time.
- Dirac's large numbers (avoided: N is conserved, so no Ġ).
- Sorkin-type Λ ~ 1/√N (avoided: I use 1/N with N conserved, which gives no fluctuation).
