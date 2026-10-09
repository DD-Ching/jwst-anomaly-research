# R2-B: "Rate-Cone World": a speculative world where geometry is inferred from a spectrum of rates

Status: **speculative hypothesis** (System A, round 2, group R2-B). Nothing here is an observation or
an established result. I did no web or repository lookups, and every "as far as I know" is flagged.

---

## 0. One-paragraph summary

The primitive is a continuous spectrum of rates. It is the positive cone of a rank-2 Euclidean Jordan
algebra, which by classification is a Lorentz cone, with a spectral measure and complex amplitudes on
it. Nothing has a position. An observer *infers* positions as stationary points of phase, x = ∂φ/∂k.
Lorentz symmetry is the automorphism group of the primitive cone, so it is exact and is not an
approximation. In the scalar sector this world matches special relativity plus FRW kinematics:
inverse-square flux, Etherington reciprocity, and **zero** phase-versus-amplitude distance mismatch.
The difference is in the *orientation* channel. Without a background connection, the only map that
carries a polarisation frame from source to observer is the orthogonal (polar) factor of the
stationary-phase Jacobian. For a negative-parity (saddle) gravitational image that factor is a
reflection. So **saddle images have reversed circular-polarisation handedness and a linear
polarisation angle mirrored about the lens's principal axis.** The amplitude is O(1) and exact, with
no free parameter. GR forbids both effects, because helicity is conserved along null geodesics.

---

## 1. Axioms

**A1 (spectrum).** There is a finite-dimensional Euclidean (formally real) Jordan algebra J. Its
cone of squares Ω = {a∘a} is the set of possible *rates*. Physical content is a positive Radon
measure μ on the closure Ω̄ (the "spectral weight") together with a complex half-density
ψ with |ψ|² dν = dμ, where ν is the Aut(Ω)-invariant measure class. Nothing else is primitive:
there is no manifold, metric, time or position.

**A2 (rank 2).** Every element of Ω̄ is a sum of at most two elements of the boundary ∂Ω (rank J = 2).
Gloss: every rate splits into at most two "pure" rates, so a massive label is two null rates.

**A3 (phase carriers).** The pure rays (rank-1 idempotents, the extreme rays of ∂Ω) form a manifold
P. Amplitudes on rays carry a U(1) phase that can interfere, so P must carry a complex structure with
a non-trivial U(1) bundle. The bundle's Chern number defines the *helicity* h ∈ ½ℤ.

**A4 (additivity).** Composite signals combine by adding rates (phases multiply), so the rate
labels form the additive group V = J. An observer is a pair (u, Φ). Here u ∈ int Ω is the observer's
own clock rate (an element of the spectrum). Φ is the phase that the observer records on a band of
rates.

**A5 (inference, no connection).** The observer assigns the source an inferred event
x^μ = ∂φ/∂k_μ (stationary phase), its flux from |ψ|² against ν, and its transverse frame from the
Jacobian of that inference map. No other structure exists to transport frames.

**A6 (gravity is inference, not scattering).** A mass label M ∈ int Ω changes *only the phase
function* φ(k) that observers reconstruct (a Fermat-type functional). It has no local vertex.
Material media act through local, Aut(Ω)-covariant vertices with mass labels (absorption and
re-emission). They are not the inference map.

## 2. Deriving Lorentz invariance (and why the bounds hold)

- Rank-2 Euclidean Jordan algebras are exactly the spin factors J_n = ℝ ⊕ ℝ^{n−1} (Jordan–von
  Neumann–Wigner classification). Their cone is Ω_n = {(t, **v**) : t > |**v**|}, the forward cone
  of ℝ^{1,n−1}. Its automorphism group is ℝ₊ × SO⁺(1, n−1) (Koecher–Vinberg). **The Lorentz group is
  the symmetry group of the primitive itself.**
- Choice of dimension (A3). The pure rays of Ω_n form the sphere S^{n−2}. Of the spheres, only S²
  carries an integrable complex structure (CP¹) with non-trivial U(1) bundles, since S⁶ is only
  almost complex. So n = 4, which gives 3+1 inferred dimensions, J ≅ Herm(2, ℂ), and
  Aut ≅ ℝ₊ × SL(2,ℂ)/ℤ₂. The rays are rank-1 matrices k = ξξ† with a spinor ξ ∈ ℂ². Helicity h is
  the homogeneity in ξ versus ξ̄.
- Evading the finite-valency no-go: there is no graph, lattice or cutoff in the axioms. The spectrum
  is a continuum, and the symmetry is an exact automorphism, not an emergent average. So no
  Lorentz-violating operator arises at any order.
- **Speed of gravity = speed of light, exactly.** Photons (|h| = 1) and gravitons (|h| = 2) are both
  supported on the *same* boundary ∂Ω, because there is only one cone. Their group velocity is
  ∂ω/∂|**k**| = 1 identically, so Δc/c = 0. This meets the 10⁻¹⁵ bound (GW170817) and the 10⁻¹⁸
  laboratory bounds trivially, not by tuning.
- Preferred frames (for example the CMB frame) are properties of the *state* μ, not of the laws. They
  break no dispersion relation.

## 3. Recovering standard kinematics

**Mass and energy are labels.** m² = det k (the Jordan determinant, which is the Minkowski norm),
and E_u = ⟨u, k⟩ = tr(u∘k). Massless ⇔ k ∈ ∂Ω.

**Translations and conservation.** The Pontryagin dual of the additive group V is the space of
characters e^{i⟨k,x⟩}. This is where the inferred position x *lives*: positions are a derived dual
of rates. Translation invariance means only phase differences are observable. 4-momentum conservation
is A4 (rates add). Angular momentum and boost charges are the Noether charges of SL(2,ℂ). Dilations
are automorphisms of Ω but not of the mass shells (det k ↦ λ² det k), so mass labels break them.
The only remnant is conformal invariance of free radiation.

**Inverse square.** For a source of spectral weight concentrated near |**k**| = ω, the amplitude at
inferred distance r is an integral over the ray sphere CP¹:
A(r) = ∫_{S²} ψ(**n̂**) e^{iωr **n̂**·**ê**} d²n̂. Stationary phase on S² gives two poles, each
contributing 2π/(iωr). So |A|² ∝ r⁻², the inverse-square law. It holds exactly in the far field
because the invariant measure d³k/|k| is exact.

**Cosmology as dilation.** An expanding universe is a state whose null-sector weight is transported
by the dilation k ↦ k/(1+z), while mass shells are mapped to themselves:
(E, **p**) ↦ (√(m² + p²/(1+z)²), **p**/(1+z)). This is the same momentum redshift as FRW. Inferred
time intervals scale as ∂φ/∂ω ∝ (1+z), which reproduces supernova time dilation.

**Four distances** (inferred flat-FRW, with comoving χ = ∫dz/H):
- d_M = χ. This is also the wavefront-curvature ("parallax") distance d_φ = (∂²φ/∂k_⊥²)/…, so d_φ = d_M.
- d_A = d_M/(1+z).
- d_L = (1+z) d_M. This follows from photon-number conservation (μ is conserved because the
  dilation preserves ν up to the Jacobian (1+z)²) plus energy and arrival-rate dilation.
- Hence d_L = (1+z)² d_A **exactly** (Etherington).

**Derived phase-versus-amplitude mismatch.** Time delay (phase) and magnification (amplitude) are
the value and the inverse Hessian determinant of one function φ, evaluated at the same stationary
point. Phase-derived distances (time-delay distance D_Δt, wavefront curvature) and amplitude-derived
distances (d_L, magnification) therefore agree in modulus to all orders: **Δd/d = 0**. This world
does not predict wavelength-independent dimming (round-1 lesson 3). What differs is the *sign*.
Amplitude sees |det A| and phase sees sgn det A. Standard physics uses the sign only for the scalar
Morse phase (−π/2 per negative eigenvalue). Here it also acts on polarisation (Section 4).

## 4. Where this world departs from GR + ΛCDM

In GR, polarisation is parallel-transported along each null geodesic. Helicity is conserved, and
weak-field lensing rotates polarisation by zero at first order. In this world (A5) there is no
connection. The observer relates the source's transverse frame to their own using the
stationary-phase Jacobian A = ∂β/∂θ (lens map). The unique orthogonal map canonically defined by A is
its polar factor R = A (AᵀA)^{−1/2}.

- **Weak lensing, minimum and maximum images** (A symmetric, det A > 0): R = 1, with no
  polarisation change at first order. The CMB and cosmic shear are therefore safe at first order. At
  second order the antisymmetric (rotation) part ω of multi-plane A rotates polarisation by ω.
  Typical rms values are about 10⁻⁴ to 10⁻³ (my order-of-magnitude estimate, not computed), which is
  below current limits.
- **Saddle images** (det A < 0, symmetric A with eigenvalues of opposite sign): R is the
  **reflection** across the eigen-axis ê₊ with the positive eigenvalue. Derived consequences, exact,
  with no free parameter:
  1. **Circular polarisation flips:** V_sad/I_sad = −V_min/I_min.
  2. **Linear polarisation angle mirrors:** χ_sad = 2θ₊ − χ_min (mod π), where θ₊ is the position
     angle of ê₊ from the lens model.
  3. **GW helicities swap:** h_R ↔ h_L. For a compact binary, inclination ι maps to π − ι between the
     type-I and type-II images, on top of the standard Morse phase π/2.
  4. The lens absorbs 2ħ|h| of angular momentum per quantum, like a half-wave plate. Angular momentum
     is conserved.

**Point-lens corollary (time domain, parameter-free).** A point lens has μ₊ − |μ₋| = 1 identically.
For an unresolved source, the net Stokes values during an event are I = (μ₊ + |μ₋|) I₀ = μ I₀ and
V = (μ₊ − |μ₋|) V₀ = V₀. So **the circular-polarisation fraction falls as V/I = (V/I)₀ / μ(t)**,
following the photometric light curve with no free parameter. Linear polarisation behaves the same
way: the component perpendicular to the source–lens axis keeps its unmagnified value, while the
component along the axis is magnified by μ. The polarisation angle therefore drifts as the axis turns
during the event. In GR, V/I and the polarisation angle stay constant for a point source.

## 5. Predictions, tests and kill criteria

**P1 (primary): polarimetry of strongly lensed radio sources resolved by image parity.**
- *Data that already exist (public archives):* VLBI, VLA and ALMA full-Stokes multi-frequency
  observations of lensed radio-loud quasars with known saddle images. Examples: B0218+357,
  PKS 1830−211, MG J0414+0534, B1422+231, CLASS and JVAS quads. Public lens models give θ₊ for each
  image.
- *Procedure:* fit RM per image across frequencies to remove Faraday rotation (∝ λ²) and get the
  intrinsic χ₀. For variable sources, compare at time-delay-shifted epochs. Test
  χ₀,sad − θ₊ = −(χ₀,min − θ₊) against GR's χ₀,sad = χ₀,min. Where circular polarisation is detected
  (≈0.1% in quasar cores), test its sign by parity.
- *Not mimicked by ordinary astrophysics:* Faraday rotation is ∝ λ² and Faraday conversion has a
  steep λ dependence. This effect is achromatic. Differential source substructure and lens-galaxy
  Faraday screens produce χ offsets that are random relative to θ₊. They cannot produce a
  *systematic* mirror about the lens eigen-axis that holds only for negative-parity images across a
  sample.
- **Kill criterion:** take ≥3 saddle images with |χ₀,min − θ₊| > 15° and RM-corrected uncertainties
  < 5°. If they satisfy χ₀,sad = χ₀,min (GR) and reject the mirror at > 3σ, this world is dead. A
  single coherent source (FRB or pulsar) seen in both images with the *same* circular handedness at
  > 5σ also kills it.
- *Honest caveat:* published VLBI position angles for some of these systems may already decide this.
  I could not check. That is the first thing System B should look up.

**P2: millilensing echoes of FRBs (CHIME/FRB baseband with full Stokes).** A compact-object lens
produces a primary and a delayed, demagnified saddle echo. Prediction: the echo has opposite Stokes V
sign and a mirrored polarisation angle, with flux ratio |μ₋|/μ₊ fixed by the delay fit. Existing
echo searches that correlate total intensity or same-handedness voltages would miss the polarised
coherent part. Re-run them with an R↔L cross-correlation template. Kill: one confirmed lensed echo
with same-sign V.

**P3: lensed GW pairs (public LVK strain).** For candidate type-I/type-II pairs, compare the
evidence for ι₂ = ι₁ (GR) against ι₂ = π − ι₁ (here). The difference can be separated only through
higher modes or precession, so the test is weak now and stronger in O5.

**Quantities current searches never measure:**
1. **Cross-image off-diagonal coherence.** Take the 2×2 cross-coherence matrix
   C = ⟨E_img1 E_img2†⟩ between time-delay-aligned images of a coherent source, in the circular
   basis. GR predicts a diagonal C (R with R, L with L, times the Morse phase). This world predicts
   an anti-diagonal C (R with L).
2. **Stokes V fraction versus magnification within a microlensing event.** The prediction is
   V/I ∝ 1/μ(t).
3. **Polarisation angle measured in the lens-eigenframe and split by image parity.** Image- and
   flux-based searches drop the sign of det A. This world puts its whole signal there.

## 6. The three weakest points

1. **A6 is put in by hand.** Laboratory optics already shows that a negative-parity focus (for
   example a cylindrical-lens mode converter) reverses *orbital* angular momentum but *not* spin
   handedness. To survive, this world has to say that refraction is local scattering while gravity is
   pure inference. That distinction is principled within the axioms but not derived. Also, Newton's
   constant G is not derived. I only *assume* that the phase response to mass labels reproduces GR's
   weak field (Shapiro delay, bending), so the scalar sector is fitted, not predicted. If gravity is
   just another refractive phase, P1–P3 vanish.
2. **The transport rule.** The observer's inferred geometry has an emergent metric and so a
   Levi-Civita connection. I have not proved that the observer *cannot* use it, and A5 forbids it by
   fiat. The polar-factor rule is canonical only for symmetric A. Strong multi-plane lenses with large
   antisymmetric parts need an extension.
3. **This world may already be falsified, and its selection arguments are thin.** VLBI polarimetry of
   saddle images, or microlensing polarimetry, may already show χ_sad = χ_min. I could not check
   either. The choice of rank 2 and n = 4 rests on a "complex phase on the ray sphere" argument that
   is heuristic. Nothing here produces Λ, dark matter, or the particle mass spectrum.

## Resemblances I noticed (honest list)

- Jordan–von Neumann–Wigner Jordan-algebraic quantum mechanics and the Koecher–Vinberg classification
  of symmetric cones.
- Penrose spinor and twistor ideas: the celestial sphere as CP¹ and helicity as homogeneity in ξ
  versus ξ̄. Wigner's little group.
- Connes-style spectral or noncommutative geometry ("geometry from a spectrum"). Momentum-space
  approaches and relative locality (positions inferred from momenta).
- Standard wave-optics lensing (Fermat potential, Morse and Maslov phases, type-II GW images).
  Pancharatnam–Berry and Rytov–Vladimirskii geometric phases. The gravitational spin Hall effect.
  Half-wave plates, and mode converters that reverse orbital angular momentum.
- Etherington reciprocity and conformal invariance of Maxwell fields (recovered, not new).
