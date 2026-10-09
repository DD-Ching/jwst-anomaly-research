# A2 — "Carry": a world with no fundamental identity and no fundamental locality

Status: **speculative hypothesis for adversarial review (System B).** Nothing here is established physics.
Provenance labels: [AX] axiom, [DER] derived inside the model, [ASM] added assumption/conversion,
[PRED] prediction, [OBS] existing observation used as a check.

## 0. Seed and one-line summary

No object has identity; nothing is fundamentally "near" anything. What exists is a relational state that
updates. The single basic quantity is **carry**, X ≡ C: the mutual information, in bits per update, that
one part of the relational state passes to another part across one update. Persistence, motion, mass,
energy, distance, time and gravity are all read off the carry matrix.

## 1. Axioms

**[AX1] Relational state.** At update index t ∈ ℤ the world is an equivalence class [S_t] of a finite
structure S_t = (V, {s_v}), where V is a finite set of *cells* and s_v ∈ Σ (a finite alphabet) is a cell's
local state. Only quantities invariant under all relabelings π: V → V are physical. Cells have no names,
so "the same cell at t and t+1" is meaningless a priori. There is no adjacency in this axiom: no locality
is assumed.

**[AX2] Update.** [S_{t+1}] is produced from [S_t] by a Markov kernel K that is (a) relabeling-equivariant
(it commutes with every π), and (b) *reversible*: K is a bijection on configurations of the full
state (so total Shannon information is conserved exactly; ensembles enter only through observers'
ignorance).

**[AX3] Carry (the quantity X).** For two sub-configurations A, B (sets of cells defined by their content,
not their names), define

  C(A→B; t) := I( s_A(t) ; s_B(t+1) )   [bits / update]

where I is mutual information over the observer's ensemble of microstates consistent with coarse data.
The carry matrix c_uv := C({u}→{v}) is non-negative and in general non-symmetric. **Carry is X**: units
bits·update⁻¹, structure = a non-negative weighted directed graph defined only up to isomorphism.

**[AX4] Capacity.** Each cell has finite outgoing capacity: Σ_v c_uv ≤ b, with b = log₂|Σ| bits/update.
Reversibility (AX2b) turns the data-processing inequality into a balance: summed over the whole state,
carry in = carry out each update (a Kirchhoff law for information).

**[AX5] Pattern and lineage.** A *pattern* P is an isomorphism class of sub-configuration. Its
*retention* is Π(P) := C(P_t → P_{t+1}) bits/update, where P_{t+1} is the occurrence at t+1 maximising
that quantity. Two occurrences are *the same pattern* iff they are linked by a chain of positive carry
along which the conditional mutual information I(P_t ; P_{t+k} | rest_{t+k}) stays positive. Identity is
lineage of carry, not a substance.

**Dynamics** is fully given by K. The empirical content I commit to is a coarse statement about K: its
carry graph G_C = {(u,v): c_uv > 0} is, at coarse scales, a statistically homogeneous graph whose ball
volumes grow as N(r) ∝ r³ (r = hop count), **plus a sparse set of "chords"**: edges with c_uv > 0
between cells whose bulk hop distance is ≫ 1. Chord density ρ_ch and chord capacity are free parameters.

## 2. Recovering familiar physics

**Time.** [DER] An observer is itself a pattern; its clock is its own retention events. Proper time
elapsed = number of retention updates the observer's pattern completes. The update index t is unobservable;
only ratios of retention counts are.

**Distance.** [DER] Define the signal time τ(A,B) = min{k : I(s_A(t); s_B(t+k)) ≥ ε}. Because carry
moves only along G_C edges, τ is a hop distance on G_C. Symmetrise d(A,B) := ½ c[τ(A,B)+τ(B,A)]
with c := ℓ/t₀ (hop length per update, both unobservable separately). On a coarse-grained homogeneous
3D graph, d satisfies the triangle inequality and approximates a Riemannian metric (geodesics = carry
paths). Locality is thus *derived*: "near" means "high carry, few hops".

**Maximum signal speed.** [DER] Information cannot cross more than one carry edge per update, so
d(A,B)/τ ≤ c everywhere in the bulk. With homogeneity and relabeling-equivariance, coarse dynamics has
no preferred frame on the graph scale; the Ignatowski argument (relativity principle + isotropy +
homogeneity) then gives either Galilean or Lorentz kinematics; a finite invariant c selects Lorentz.
(This step needs *emergent* isotropy, see Self-critique 1.)

**Energy and mass.** [ASM] Fix one conversion constant by identifying carry rate with the orthogonal-
transition rate of quantum dynamics: E := (πħ/2)·Ṅ, where Ṅ = total carry a pattern commands per unit
proper time, in bits/s. A pattern at rest spends all of it on re-creating itself in place: Ṅ₀ = Π/t₀.
Define m := πħṄ₀/(2c²). Under the Lorentz kinematics just derived, the 4-vector (E/c, p) with
p ∝ hop-carry rate has invariant norm mc, giving **E = mc² at rest** and E² = p²c² + m²c⁴ in general.
A massless pattern (light) is one with Π_internal = 0: all its carry goes into displacement, so it
moves one hop per update, at c. Numerically an electron has Ṅ₀ = 2m_ec²/(πħ) ≈ 4.9×10²⁰ bits/s.

**Inverse-square dilution.** [DER] A source emits carry at rate Ṅ_s. Reversibility conserves it; on a
graph with N(r) ∝ r³, the shell at hop radius r contains dN/dr ∝ r² cells, so carry per cell (flux)
∝ Ṅ_s / (4πr²). The exponent 2 is the ball-growth dimension minus one; it is a measurement of G_C.

**Conservation laws.** [DER] Total carry is conserved (AX2b, AX4) ⇒ energy conservation. Statistical
homogeneity of G_C under coarse translations ⇒ (Noether-like, approximate) momentum conservation; isotropy
⇒ angular momentum. Pattern-number conservation for "quanta" follows when patterns carry fixed, indivisible
information content.

**Gravity (sketch).** [ASM] A massive pattern consumes capacity b locally (its retention occupies the cells
it lives on), leaving less carry for transit. Transit slows: an effective refractive index
n(r) = 1/(1 − f(r)), f = locally occupied capacity fraction. Matching weak-field light bending requires
f ≈ 2GM/(rc²) at large r (G enters as a second conversion constant; not derived).

**Why must continuity pass through intermediate positions?** Because "position" *is* placement in G_C and
continuity *is* a carry lineage; a lineage can only advance along edges with c_uv > 0. In the bulk those
edges join hop-neighbours, so continuity sweeps through intermediate positions tautologically. **But no
fundamental locality** means nothing forbids chords. A lineage that crosses a chord is *genuine continuity*
(positive carry chain) with no intermediate bulk positions. This is the second continuity channel.

**Continuity vs copy.** Three operational criteria, all following from reversibility:
1. *Closure*: continuity conserves carry, so the source region shows a deficit equal (in bits, i.e. in
   quanta) to the excess at the destination. A copy requires an independent carry supply: no deficit.
2. *Screening*: for continuity, I(A_past ; B_future | A_future) > 0. Knowing what remains at A does not
   explain B. For a common-cause copy this conditional information is zero.
3. *No persistence*: after a transfer the original occurrence's retention drops; after a copy it does not.

**Paired deficit/excess ("missing light").** [DER] A photon is a pure-transport pattern. If its lineage
meets a chord mouth, it exits at a remote sky region. Bits (quanta) are conserved, **energy is not
one-to-one**: the energy per bit/s is fixed per unit *local* proper time, and the local update-to-proper-
time ratio differs at the two mouths (differing capacity occupancy f). Writing s := (1−f_in)/(1−f_out):

  N_out = N_in (photon number);  ν_out = ν_in / s;  Δt_out = s·Δt_in;  E_out = E_in / s.

So the excess is the deficit's spectrum and light curve, rescaled by one common factor s, with
**s_t · s_ν = 1** and photon fluence matching.

## 3. Where this world differs from GR + ΛCDM

**D1. Non-local chords.** GR has no transport between distant regions without a traversable geometry and
exotic stress-energy; here chords are generic defects with density ρ_ch.

**D2. Capacity saturation in strong fields.** With the simplest law n = 1/(1−2m/r) (m = GM/c²), the
expansion n = 1 + 2m/r + 4(m/r)² + … differs from GR's isotropic-coordinate value 15/4 at second order:
second-order deflection α₂ ≈ 4π(m/b)² vs 15π/4(m/b)² (assuming PPN-style cross-terms carry over).
More decisive: the photon sphere solves d/dr[r·n(r)] = 0 ⇒ r = 4m, critical impact parameter
b_c = 8m, vs GR's 3√3 m ≈ 5.20m. That 54% larger shadow is **[OBS] already excluded** by horizon-scale
imaging of nearby supermassive black holes, which match GR to roughly ten to twenty percent. So the linear
saturation law is dead; f(r) must be non-linear, and that function is the model's main tunable freedom (see
Self-critique 2).

**D3. Achromatic beam loss ⇒ distance-duality violation.** Chords take photons out of the bulk beam:
D_L = D_A(1+z)²·η(z), with η(z) = exp[τ_ch(z)/2], τ_ch(z) = ∫ρ_ch σ_ch dl. Wavelength-independent,
unlike dust.

**D4. Bulk-path-only dispersion.** Plasma dispersion accrues only on bulk segments. A signal that used a
chord carries less dispersion than its distance implies.

**D5. Diffuse re-emission.** Light lost at chord mouths reappears elsewhere, rescaled by s. That
redistributes, but does not change, the integrated background photon count; the angular power of the
background gains a non-local cross-term.

## 4. Predictions (each: not easily mimicked, existing data, kill criterion)

**P1 (D1, top prediction): Time-stretched spectro-temporal twins.** Pairs of transients at
widely separated sky positions (≫ lensing scales, no intervening mass concentration) whose light curves
match after a time stretch s, whose spectral peak energies satisfy E_pk,1/E_pk,2 = s (the *same* s), and
whose photon fluences (not energy fluences) agree after correcting each mouth's beam-loss fraction.
*Why hard to mimic:* gravitational lensing gives images at one redshift (s = 1) separated by arcseconds;
unrelated bursts have no reason to obey s_t·s_ν = 1 with matching photon count.
*Data:* GRB catalogues with time-tagged light curves and peak-energy fits (gamma-ray monitors), FRB
catalogues (burst morphology + spectral occupancy), optical transient surveys.
*Kill:* do a blind all-pairs search with a pre-registered match statistic, calibrated on time-scrambled
catalogues. If the number of pairs passing both the shape-after-stretch and the s_t·s_ν = 1 tests is
consistent with the scrambled background, with an upper limit on the pair fraction below 10⁻³ per burst,
then chord transport of transients at τ_ch ≳ 10⁻³ is falsified. If the s_t·s_ν relation fails for every
shape-matched pair, the mechanism in §2 is falsified outright.

**P2 (D4): Sub-floor dispersion.** Localized FRBs whose dispersion measure lies *below* the Galaxy +
minimal-IGM floor for their host redshift, at an incidence equal to τ_ch(z) and growing with z as the
integral of ρ_ch. *Why hard to mimic:* host and halo contributions only add DM; the IGM floor can be
lowered only by rare underdense sightlines, whose scatter has a known, z-dependent width.
*Data:* catalogues of FRBs with host-galaxy redshifts.
*Kill:* if N localized FRBs show zero sub-floor events (beyond 3σ of the modelled IGM scatter), then
τ_ch < 3/N (95%). If that bound is inconsistent with a τ_ch measured in P3 or P1, the model is falsified.

**P3 (D3): Achromatic distance-duality offset.** η(z) − 1 ≈ τ_ch(z)/2 > 0, identical in every band.
*Data:* supernova Ia luminosity distances vs baryon-acoustic angular distances; strong-lens time-delay
distances; multi-band SN photometry to check achromaticity.
*Kill:* η(z) = 1 within 1% up to z ≈ 1 ⇒ τ_ch(1) < 0.02. Any detected η ≠ 1 that is chromatic is dust,
not chords. A non-zero achromatic η with no P1/P2 counterpart at the implied rate falsifies the
"same chords for all signals" claim.

**P4 (D2): Shadow/photon-ring ratio.** Whatever f(r) is fitted to the shadow diameter, the same f fixes
the photon-ring sub-image spacing and the weak-field α₂ coefficient. *Data:* published horizon-scale images,
plus binary-pulsar/solar-system bending constraints. *Kill:* no single monotone f(r) with f → 2m/r that
fits both the shadow size and the second-order solar-system bending limits.

## 5. A quantity current searches never measure: the transfer covariance

Lens morphology, microlensing curves and source-count deficits all analyse **one sky location at a time**.
The diagnostic here is non-local and signed:

  T(θ, Δt, s) := ⟨ δN_i(t, ν) · δN_j(t + Δt, ν/s) ⟩ over pairs (i, j) at angular separation θ,

where δN is the *photon-number* residual (observed − predicted count, not flux) of a variable source.
Chord transfer predicts **T < 0** (one source dims as the other brightens), at a fixed lag Δt unrelated to θ,
and only at one s per pair. Calibration systematics give T > 0 (common-mode) at Δt = 0, s = 1. Lensing
gives T > 0 with s = 1 at arcsecond θ.
*How to measure:* take multi-epoch quasar variability surveys (repeated all-sky photometry), convert to
photon-number residuals per band, remove per-epoch and per-field zero points, then compute T on a
(θ, Δt, s) grid using the multi-band light curves to test the s-rescaling (band i at mouth A ↔ band shifted by
s at mouth B). Null tests: time-reversed and band-shuffled pairs. *Kill:* T consistent with zero at
negative sign for all θ > 1°, Δt within the survey baseline and s in [0.5, 2], with per-pair amplitude limits
below 1% of quasar variance.

A second unmeasured quantity: **photon-number closure** between a dimmed region and its partner, i.e.
comparing counts N rather than energies E. Every existing deficit search is energy- or magnitude-based.

## 6. Self-critique: the three weakest points

1. **Emergent Lorentz invariance is asserted, not shown.** Discrete graphs usually pick out preferred frames,
   and chords make it worse: random long edges turn a 3D lattice into a small-world graph, so distances
   collapse unless ρ_ch is tiny. I give no reason ρ_ch is tiny other than observation. Lorentz-violation limits
   from high-energy photon timing may already bound the graph scale and chord density severely.
2. **Constants and f(r) are inserted, not derived.** ħ (via the Margolus–Levitin identification), G and the
   saturation function f(r) are all fitted. The simplest f is already excluded by black-hole shadow sizes
   (§3 D2). Once f is free, the gravity sector predicts little. Only the chord sector (P1–P3, §5) makes sharp
   structural claims, and even there amplitudes depend on unconstrained ρ_ch and σ_ch.
3. **Energy rescaling at chord mouths is postulated.** Rule s = (1−f_in)/(1−f_out) assumes energy per bit
   is set by local capacity occupancy, while quanta are indivisible information packets. Neither is derived
   from K. Also, ordinary lensing *does* produce paired magnification/demagnification with global flux
   conservation; my separation from lensing relies on large θ, s ≠ 1 and absent lens mass, any of which
   instrumental or selection effects could imitate in small samples.

## Resemblances I noticed (not built from)

Information-first ("it from bit") programmes; graph-based emergent-space models and rewriting-rule
universes; the Margolus–Levitin speed limit (used explicitly as a conversion); entropic/capacity pictures of
gravity and analogue-gravity refractive-index treatments; Einstein–Rosen-bridge-like shortcuts (the chords
look like microscopic wormholes, phenomenologically); no-cloning (the continuity/copy criterion); opacity-
based distance-duality tests. The chord idea also resembles "tired light" in producing achromatic loss, but here
photon number is conserved globally, not destroyed.
