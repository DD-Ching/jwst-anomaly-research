# A1 — Relational Transition World (RTW)

System A ("Physics Invention"), group A1. Status: **speculative hypothesis**. Nothing here is established physics.
I did no web searches and read no repository files. Every number is either a derivation inside the model
or an order-of-magnitude estimate, and is labelled as such.

---

## 1. Axioms

**A1 (States).** A countable set Ω of global states.

**A2 (Rules).** A finite set of rule families R = {a, b}. Each family r ∈ R is a relation T_r ⊂ Ω × Ω
with a rate function κ_r : T_r → (0, ∞).

**A3 (Dynamics).** The world is a continuous-time Markov jump process on Ω with generator
Q(ω, ω′) = Σ_r κ_r(ω, ω′) · 1[(ω, ω′) ∈ T_r]. A *history* is a sample path. This is the only source
of probability.

**A4 (Derived loci, i.e. factorization).** Among all product decompositions Ω ≅ Π_{x∈V} Σ_x
(finite alphabets Σ_x) choose the **finest one in which every a-transition changes at most k coordinates**,
with k fixed. Its index set V is the set of *loci*. Loci are not primitive. They are the coordinates
in which rule a is local. Write σ_x for the x-component of a state.

**A5 (Adjacency per rule).** For each r, the graph G_r on V has x ~_r y iff some transition in T_r
changes σ_x and the rate of that transition depends on σ_y, or the transition changes both.
The *r-distance* d_r(x, y) is the graph distance in G_r. The distance for a set S ⊂ R is
d_S(x, y), the graph distance in ∪_{r∈S} G_r.

**A6 (Two regimes.)** These conditions are imposed, not derived.
(i) G_a has bounded degree and polynomial growth |B_a(x, n)| ≍ n³ (quasi-isometric to R³ at large n).
(ii) b-edges are sparse. They link pairs at *arbitrary* a-distance, with probability
P(x ~_b y) = p₀ · d_a(x, y)^(−α) (α > 0), and with total rate κ_b ≪ κ_a.
(iii) A b-transition copies the local "excitation" content only partially. A fraction ε ≪ 1 of the
a-conserved quantity (A7) crossing a b-endpoint is transferred to the partner locus. The rest stays.

**A7 (Conserved label).** There is a function h : Σ → ℝ≥0 such that every a-transition preserves
Σ_{x∈supp} h(σ_x) on its support, and every b-transition preserves h(σ_x) + h(σ_y) over its endpoint pair.

**A8 (Growth rewrite).** A subfamily a_split ⊂ T_a replaces one locus by two adjacent loci
(it refines the factorization) at rate H per locus. No other new primitive is introduced.

That is the whole world: a set, two transition relations with rates, a Markov process, and a factorization
criterion. Space, time, causality, energy and matter are all derived below.

## 2. What an internal observer recovers

**Events and causality.** Two transitions are *independent* if applying them in either order reaches the
same state (diamond property). An *event* is a transition occurrence modulo swaps of adjacent independent
transitions. Event e *precedes* f (e ≺ f) if e occurs before f in every representative history.
Causal order is therefore derived: it is the non-commutation order of rewrites.

**Metric distance.** Physical "time" for an observer is the count of a-transitions in its own
neighbourhood (a local clock). Let T(x, y) be the first-passage time for a perturbation at x to change
the law of σ_y using only a-transitions. With state-dependent (effectively random) rates, first-passage
theory gives a shape theorem: T(x, y) / d_a(x, y) → μ(direction) almost surely. The emergent metric is
d(x, y) = lim T(x, y) · v, where v is defined below.

**Light cone.** With bounded rates and bounded degree on G_a, the probability that a-influence travels
graph distance n within time t is ≤ (e κ_max t / n)^n. So influence is confined to n ≤ v t with
v ≈ e κ_max ℓ_a, where ℓ_a is the coarse-grained length of one a-step. This is the maximum signal speed:
a *derived*, statistical light cone. With b present the cone **leaks**. The causal future of x is the
a-cone plus thin a-cones rooted at every b-partner of x, delayed by one hop latency τ_h ≈ 1/κ_b.

**Why two places 1 Mly apart for light can be adjacent.** Light consists of excitations that propagate
only by a-transitions (they are the h-carrying patterns that a moves). Their distance is d_a. A b-edge
gives d_{a∪b}(x, y) = 1 regardless of d_a(x, y). Adjacency is rule-relative. "1 Mly apart" is a statement
about G_a only. If α lies between 3 and 6, the long-range-edge result says d_{a∪b} grows only
polylogarithmically in d_a, ~ (log d_a)^Δ with Δ = log 2 / log(6/α). Light does not use this geometry
because each b crossing carries only the fraction ε of the h-content.

**Energy-like conservation.** By A7 and the boundedness of supports, for any region A,
dH_A/dt = −Φ_a(∂A) − Φ_b(A), where H_A = Σ_{x∈A} h(σ_x). The first term is an ordinary boundary flux,
which gives a local continuity equation. The second is a *non-local* transfer to b-partners outside A.
Global H is conserved. Local H is conserved only up to O(ε κ_b n_b), where n_b is the b-edge density.

**Matter and gravity (weakest derivation).** Persistent local patterns ("matter") occupy their loci and
lower the rate of available a-transitions nearby. Let φ(x) = 1 − (local a-rate)/(vacuum a-rate). Clocks
run slow by (1 − φ), and first-passage times grow by ∫φ dl/v. This gives Shapiro-like delays and,
through Fermat's principle on T, deflection. If the depletion relaxes diffusively in a 3-D graph with
pattern sources, then quasi-statically ∇²φ ∝ ρ_pattern. That reproduces a Newtonian-like potential.
I have **not** derived the GR factor of 2 in light bending. It must be imposed by requiring that
spatial step-lengths dilate together with clock rates.

**Agreement of the four distances.** On a quasi-isometric R³ graph with expansion from A8
(a(t) = e^{Ht} at constant H; more generally a(t) follows the split history), the observer defines:
- **time-delay distance**: from first-passage times T;
- **angular-diameter distance d_A**: from transverse extent divided by the fraction of the observer's
  incoming a-directions, so d_A ∝ |∂B_a(n)|^{1/2} / (1+z);
- **luminosity distance d_L**: from conservation of h-carrier number across ∂B_a(n), so F = L/|∂B| with
  two factors of (1+z) from count-dilution and clock-rate;
- **dynamical distance**: from the depletion field φ, which is sourced and relaxed on the same graph.

Redshift arises because a carrier's wavelength, counted in loci, is fixed in transit while loci split.
The three conditions together give d_L = (1+z)² d_A, and all four agree, **iff**:
(1) G_a is manifold-like;
(2) carrier number is conserved along a-paths;
(3) φ propagates on the same graph as light.
The b-family breaks (2) at order ε, and breaks (3) if φ also crosses b-edges. Every prediction below
comes from these breaks.

## 3. Differences from GR + ΛCDM

**D1 — Achromatic carrier leakage.** Each b-endpoint met along a path diverts a fraction ε of the
carriers. The surviving fraction is exp(−Γ_b χ/c), where Γ_b ≡ ε κ_b f_b and f_b is the fraction of loci
with a b-edge. Hence
η(z) ≡ d_L / [(1+z)² d_A] = exp(+½ Γ_b χ(z)/c) ≈ 1 + ½ Γ_b χ(z)/c.
The leakage is frequency-independent and linear in comoving distance at low z.
The leaked carriers are **not destroyed**. They reappear at the partners, which are spread over the
whole volume when α is small, and form a diffuse, isotropic background. **Consistency relation**: the
fractional excess diffuse background satisfies I_excess / I_resolved ≈ ⟨1 − e^{−Γ_b χ/c}⟩, averaged
over the luminosity-weighted source population. The excess has the integrated spectral shape of
galaxies, redshifted by the extra path. One parameter fixes both effects.

**D2 — Ghost gravity (conditional branch).** If depletion also crosses b-edges, then
φ = φ_local + K_b * φ_local, where K_b is the b-kernel (∝ d^{−α}). For α < 3 the convolution is dominated
by distant matter. It is nearly uniform, so it acts as an extra smooth component that drives expansion
but does not cluster. For large α it is short-range and halo-like, tracking baryons smoothed over a
length ℓ_b.

**D3 — Out-of-cone dependence (the owner's causal event network).** Event pairs (e₁, e₂) can show
statistical dependence with emission-time separation ≈ τ_h, which is much less than |r₁ − r₂|/c.

**D4 — Arrival-time roughness.** First-passage times on a random-rate graph fluctuate. Standard growth
universality gives σ_T ∝ L^χ ℓ_a^{1−χ}/c with χ ≈ 1/3 (or somewhat smaller in 3-D). This is an
**achromatic**, distance-growing intrinsic broadening, unlike plasma scattering (∝ ν^{−4}) or
dispersion (∝ ν^{−2}). An estimate: FRB microstructure of ~1 μs at L ~ 1 Gpc requires
ℓ_a ≲ (c σ_T)^{3/2} / L^{1/2} ≈ 1 nm for χ = 1/3. That bound is not yet a kill, but it is a hard constraint.

**D5 — Anisotropic limit shape.** The shape theorem's μ(direction) is generically **not** a sphere.
Isotropy of v must come from a symmetry of the a-rule. Any residual effect is a direction-dependent but
energy-independent signal speed, of fractional size δ_μ.

## 4. Predictions with kill criteria

**P1 (D1): leakage–background consistency.**
- *Data*: SN Ia luminosity distances together with BAO angular-diameter distances, or strongly lensed SNe
  that have both time-delay distance and calibrated flux; plus direct extragalactic-background-light
  measurements and deep galaxy counts.
- *Why ordinary astrophysics struggles*: dust is chromatic and re-emits in the far-IR. Calibration drifts
  do not produce a matching optical diffuse excess. Photon–particle mixing is chromatic or depends on
  magnetic fields.
- *Kill*: fit η(z) − 1 = ½ Γ_b χ/c across 0 < z < 2. Independently, fit Γ_b from the diffuse excess.
  (a) If η is consistent with 1 to 10^{−2} at z ≈ 1, then Γ_b < ~10^{−2} c/χ(1). If the diffuse excess
  then exceeds what that Γ_b allows, the leakage mechanism cannot be its origin.
  (b) If η − 1 is found to be nonzero but chromatic, or not linear at low z, D1 is falsified.

**P2 (D3, geometry-clean version): hybrid images in strongly lensed transients.**
For a lens, split each image's light-travel time into a source-side leg s_i (source to lens plane) and an
observer-side leg o_i (lens plane to observer). A b-edge joining tube i to tube j near the lens plane
produces a faint **hybrid arrival** at image j's position, at time s_i + o_j. Relative to image j's main
signal, it appears at lag δ_ij = s_i − s_j, which is computable from the lens model. Its flux ratio is
≈ ε_ℓ × (magnification factor), where ε_ℓ is the b-crossing probability within the bundle. When
s_i < s_j and o_j < o_i, the hybrid arrives **before the first Fermat image**. No a-path can do this,
because the leading image is the global minimum of arrival time.
- *Data*: decade-long quasar-lens light curves from lens-monitoring campaigns; HST/JWST light curves and
  pre-explosion imaging of multiply imaged supernovae, especially cluster-lensed ones with image
  separations of 10⁵–10⁶ ly.
- *Test*: model each image as main + Σ_i a_ij · (intrinsic curve shifted by δ_ij), with a_ij ≥ 0, then
  stack over systems. Microlensing is chromatic and uncorrelated with the δ_ij. Lens-model errors move
  all delays coherently.
- *Kill*: stacked a_ij consistent with 0 at < 10^{−3} bounds ε_ℓ. Excess power at lags that do **not**
  match the model's δ_ij falsifies this mechanism as its cause. A precursor before the leading image,
  achromatic and at a model-predicted δ_ij, is the positive signal; vetting would still be required.

**P3 (D2): ghost-gravity offsets.**
- *Data*: merging clusters with weak/strong-lensing maps, X-ray gas maps and galaxy distributions;
  SZ, X-ray and lensing mass comparisons.
- *Prediction*: the lensing excess equals K_b convolved with the baryon depletion field. For large α,
  the lensing peaks sit at the depletion-weighted baryon centroid smoothed by ℓ_b, *between* the galaxies
  and the gas.
- *Kill*: lensing peaks that coincide with collisionless galaxies to within the smoothing length, with
  significance > 3σ across several merging systems, kill large-α D2. A dark/baryon ratio that varies by
  > 10× at a fixed smoothed environment kills small-α D2. Honest note: dark-matter-dominated dwarfs
  already put D2 under heavy pressure.

**P4 (D4): achromatic broadening floor.**
- *Data*: FRB catalogues (burst widths, scattering times, DM-inferred and host-redshift distances); GRB
  minimum variability times.
- *Prediction*: after removing the ν^{−4} scattering component, the minimum achieved intrinsic width
  w_min(L) follows the envelope ∝ L^χ.
- *Kill*: if w_min shows no upward trend with L across two decades of distance, while sub-μs structure
  is seen at Gpc, then ℓ_a < 1 nm. Any chromatic residual floor is not D4.

**P5 (D5): direction-dependent arrival offsets.**
- *Data*: multi-messenger events (GW standard sirens with EM counterparts, GRB–GW lag) and pulsar timing
  arrays.
- *Kill*: D5 predicts an achromatic offset of fractional size δ_μ(n̂) with a fixed angular pattern
  (cubic-harmonic for a lattice-like G_a). Random-sign residuals with no ℓ = 4 pattern bound δ_μ. Note
  that a GW–GRB coincidence at ~40 Mpc already implies δ_μ ≲ 10^{−15}.

## 5. Quantities current searches never look at

**(a) Pre-leading-image flux at hybrid lags δ_ij.** Lens searches look for images and for delays at
Fermat stationary points. Nobody fits light-curve power at the cross-terms s_i − s_j, or *before* the
leading image.
*Measure*: from the lens model, decompose the geometric and potential delays into source-side and
observer-side parts. In thin-lens geometry the potential term sits on the lens plane, so assign it with
both conventions and treat the result as two candidate lag sets. Cross-correlate each image's light curve
against the leading image's curve at those lags only, so that the trials factor is small.

**(b) A blindness theorem for the causal event network, and the resulting admissibility statistic.**
For events e₁ and e₂ at comoving positions r₁ and r₂, the observed lag is
Δt_obs = Δt_emit + (χ₁ − χ₂)/c. A-admissibility (light-cone causality) requires
|Δt_obs − (χ₁ − χ₂)/c| ≥ |r₁ − r₂|/c.
A b-linked pair has Δt_emit ≈ τ_h, so Δt_obs ≈ (χ₁ − χ₂)/c. To test admissibility you must know
χ₁ − χ₂ to better than c × (survey baseline) ≈ 10 ly. With spectroscopic σ_z ~ 10^{−4} the radial
uncertainty is ~1 Mly. The geometric acceptance is therefore S ≈ 10 ly / 1 Mly ~ 10^{−5}.
**Consequence: timing-coincidence searches between distinct extragalactic sources cannot, even in
principle, demonstrate out-of-cone dependence. Their admissible region is not computable.** Only
event sets whose geometry is known to better than the cadence are testable:
- images of one lensed source (P2);
- binary pulsar systems (orbital geometry known to light-microseconds): look for dependence between the
  companion's orbital phase and pulsar timing residuals with a lag shorter than the projected light time;
- Solar System pairs (solar flares versus terrestrial detector channels in the window −8.3 min < lag < 0
  before the flare light arrives).

The quantity to record for every candidate pair is the **admissibility margin**
m = (|Δt_obs − Δχ/c| − |Δr|/c) / σ_m. Claims are allowed only when |m| > 5. Current anomaly searches
never compute this number.

## 6. Self-critique: three weakest points

1. **Geometry is put in by hand.** A6(i) assumes G_a is quasi-isometric to R³. Isotropy of the derived
   speed needs a symmetry I have not justified (D5 shows the problem). The model explains neither why the
   dimension is 3 nor why the limit shape is round to 10^{−15}.
2. **Gravity is not derived.** "Depletion" gives a Newtonian-like Poisson field only under an assumed
   diffusive relaxation. The GR light-bending factor, gravitational waves and the cosmic expansion law
   (beyond a constant split rate H) are imposed or missing. P3 tests a branch I cannot derive from the
   axioms.
3. **Retreat to small ε.** The b-family has free parameters (p₀, α, ε, κ_b). The blindness theorem shows
   that most b-edges are untestable by timing. The model can therefore always shrink ε below every bound,
   which makes it nearly unfalsifiable. Only P1's two-observable consistency relation and P2's
   model-fixed lags tie the parameters to something that could fail. If both come out null, Occam's
   razor should remove the b-family. (Also missing: any account of quantum interference. The process is
   classical-stochastic.)

## Resemblances I noticed (not built from)

- The diamond-property event order resembles trace theory and event structures, rewriting systems with
  "causal invariance", and the partial orders of causal-set approaches.
- Emergent geometry from graph growth resembles graph-based emergent-space programmes.
- Using rule-relative adjacency to make far-apart places adjacent resembles "non-local link"/wormhole-like
  ideas and entanglement-geometry proposals.
- The mathematics used is standard: long-range percolation (polylog distances for 3 < α < 6),
  first-passage percolation and its shape theorem, KPZ-type fluctuation exponents, and Lieb–Robinson-type
  velocity bounds.
- D1 is phenomenologically similar to standard "distance-duality / photon-number violation" tests, and its
  diffuse excess is reminiscent of extragalactic-background-light excess debates.
- D2 resembles non-local or smoothed-source modified-gravity proposals.
