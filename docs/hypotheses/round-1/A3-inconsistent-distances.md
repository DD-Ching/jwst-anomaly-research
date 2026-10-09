# A3: Bridged Markov substrate, where distances measured different ways disagree in a fixed pattern

Status: **speculative hypothesis** written for adversarial review by System B. Nothing here is established physics. I did no web lookups and read no repository files.

---

## 1. Minimal axioms

**A1 (Substrate).** There is a countable set of cells `V`, a finite alphabet `Σ`, and a strictly positive probability measure `P` on configurations `s ∈ Σ^V`. That is the only primitive: it says which joint configurations are consistent with each other, and with what weight. Space, time, matter and fields are not assumed.

**A2 (Locality is derived, not assumed).** Let `G = (V, E)` be the minimal graph for which `P` is Markov, meaning that every cell is conditionally independent of all other cells given its neighbours. By Hammersley–Clifford, `P(s) ∝ ∏_C ψ_C(s_C)` over the cliques `C` of `G`. Give each edge a weight `w_ij`, the coupling strength (for example, the maximum conditional mutual information of the pair).

**A3 (Bulk plus bridges).** `G = G₀ ∪ B`.
- `G₀` has bounded degree. At coarse scales it is quasi-isometric to `Z^4`, and its weights depend smoothly on the local configuration.
- `B` is a sparse set of *bridges*. Each bridge is an edge whose endpoints are at least `ℓ_b` apart in `G₀`, and it carries a small complex coupling amplitude `ε e^{iφ_b}` with `ε ≪ 1`. Bridges have comoving number density `n_b(x) = n̄_b (1+δ(x))^α`, where `δ` is the local record-density contrast defined in A4. In other words, bridges attach preferentially to structure.
- `B = B_f ∪ B_c`. Bridges in `B_f` are frozen (always present). A bridge in `B_c` is *conditional*: its potential is `ψ_b = 1 + ε·χ_R(s_{N(x)})·χ_R(s_{N(y)})`, so it exists only when the neighbourhoods of both endpoints lie in a resonance set `R ⊂ Σ^{N}`.

**A4 (Records).** A *record* of region `A` in region `B` is a function `f(s_B)` with `I(f; s_A) ≥ (1−δ)H(s_A)`. The redundancy `Red(A)` is the number of disjoint regions that each hold a record of `A`. Coarse-grained "matter" means configurations with high redundancy, and `δ(x)` is their density contrast.

**A5 (Boundary asymmetry).** `P` is conditioned on one low-entropy boundary of `G₀`. This is a contingent fact about the one realised `P`, not a law, and I treat it as an assumption.

What makes disagreement between distance measures natural: each distance measurement reads a *different functional* of the same coupling operator. Local edges and bridges enter those functionals at different orders in `ε` and with different weightings. With `B = ∅`, every functional collapses onto a single emergent metric. With `B ≠ ∅`, they come apart in a fixed pattern.

---

## 2. What an observer recovers

**Ordering and arrow.** Define `A ⊲ B` ("A precedes B") when `B` holds records of `A` but `A` holds none of `B`. By data processing, information about any region can only decrease across nested separating sets `S_1 ⊃ S_2 ⊃ …`. Redundancy is not symmetric, though: under A5, regions near the low-entropy boundary end up copied into many far regions, and the reverse does not happen. The separating sets therefore form a foliation, and `⊲` is a strict partial order that agrees with the redundancy gradient. That foliation is the derived "time" direction.

The ordering **fails** in two places:
- (i) In equilibrium patches with no redundancy gradient, `⊲` is empty: there is no local time.
- (ii) Across a bridge, each endpoint can hold records of the other, so `⊲` becomes symmetric. The order is *ambiguous*, not reversed. Conditional bridges produce brief, local episodes of this ambiguity.

**Dynamics and maximum signal speed.** Cut `G₀` along the foliation. `P` then factorises into transfer operators between neighbouring layers, `P = ∏_k T_k`, and the generator `H ≡ −log T` acts as an emergent Hamiltonian. With bounded degree and couplings in the Dobrushin-uniqueness regime, the effect of conditioning at `x` on cell `y` falls off as `λ^{d_{G₀}(x,y)}`. Records advance at most one `G₀` step per layer, so a light cone with maximum speed `c = (graph steps)/(layer step)|_max` follows. Bridges break this only locally, with amplitude `ε` and between fixed endpoints.

**Conservation laws.** Any symmetry of the potentials `ψ` that preserves the layering gives a charge whose layer sum `Σ_{x∈S_k} q(s_x)` is the same on every layer. If the `T_k` do not depend on `k`, `H` is conserved; that is "energy". Bridges keep *global* conservation but break *local* conservation: a flux can leave at one endpoint and reappear at the far one.

**Metric and consistent distances (the generic case).** Matter (`δ`) modulates the `G₀` weights. The coarse graph metric then converges, in the Gromov–Hausdorff sense, to a curved Lorentzian metric `g`. When `B = ∅`, linearising `T` gives a wave operator `□_g`, and all four probes are functionals of that one operator:
- angular size reads transverse ray geometry;
- flux reads the conserved norm spread over a ray bundle;
- time delay reads the stationary phase (Fermat);
- dynamics reads the static (ω→0) Green's function.

Geodesic optics plus norm conservation gives Etherington reciprocity, `d_L = (1+z)² d_A`. Fermat's principle on the same `g` gives the standard time-delay distance. The static Green's function of the same operator gives the potential that both dynamics and lensing infer. **All four measures agree.**

---

## 3. Where the world differs from GR + ΛCDM

Treat bridges as sparse scatterers in the coarse operator, `L = L₀ + K_B`. A coherent disturbance crossing bridge depth

`τ(z) = ∫₀^{χ(z)} n_b σ_b dχ = n̄_b σ_b ∫ (1+δ)^α dχ`

feels the forward amplitude of every bridge it passes. The analogue of the optical theorem then fixes how each probe responds:

| Probe | Functional | Bridge response | Sign |
|---|---|---|---|
| Angular diameter `d_A` | transverse ray geometry | no bias at O(ε). At O(ε²) a faint diffuse halo appears, but the image centroid and size do not change | 0 |
| Time-delay path `d_T` | stationary phase | `δ_T ≡ Δln d_T = −ε sinφ_b·τ` (an effective index offset) | set by `φ_b` |
| Luminosity `d_L` | conserved norm in the beam | `δ_L ≡ Δln d_L = +½ ε cosφ_b·τ` (power leaks to remote endpoints and returns as an isotropic background) | dimmer |
| Dynamics `M_dyn` | static Green's function | the static potential leaks through bridges: `M_dyn/M_lens ≈ 1 − κ_s (r/ℓ_b)²` for `r ≲ ℓ_b`, with `κ_s ∝ ε n_b ℓ_b^4`. It saturates at `r ~ ℓ_b` | deficit |

(Here `ε cosφ_b` plays the role of the "imaginary part" of the forward amplitude, the part that removes power, and `ε sinφ_b` the "real part", which shifts phase.)

**The specific pattern:**
1. `d_A` is the clean reference.
2. `δ_L` and `δ_T` are both proportional to the same `τ`. Their ratio `Q ≡ δ_L/δ_T = −½ cot φ_b` is therefore **a universal constant**: it is the same on every sightline, at every redshift, at every wavelength, and for every messenger (photons, gravitational waves, neutrinos), because bridges belong to the substrate and not to any one field.
3. Static probes (dynamics, hydrostatics) show a deficit relative to propagating probes (lensing). The deficit grows as aperture² and saturates at `ℓ_b`.
4. `τ` follows matter through a **flat comoving kernel** weighted by `(1+δ)^α`. Lensing convergence instead uses the geometric kernel `D_l D_ls/D_s`.

**Fiducial sizes.** A CMB constraint forces `α ≳ 2`. Bridges then sit mostly in collapsed structure at z ≲ 2, otherwise the damping of anisotropies, `e^{−ε cosφ_b τ(1100)}`, would show up as a mismatch with the reionisation optical depth. With that in place:
- `ε cosφ_b τ(z=1) ~ 0.01–0.04`, so SNe are about 1–4% dimmer than `(1+z)² d_A`;
- `|δ_T| ~ 1–5%`;
- `κ_s(r_500) ~ 0.1–0.2` in clusters.

These are knobs, not results. **Resemblance (honest):** if `φ_b` makes `δ_T < 0`, time-delay H₀ comes out higher than inverse-ladder (BAO+CMB) H₀. That points the same way as the existing tension. It is not evidence for this model.

---

## 4. Predictions and kill criteria

### P1. Lensed transients: dimming that tracks convergence, independent of parity and wavelength
Images that lie at larger convergence `κ_i` pass through more structure, so `Δτ_ij ∝ Σ_crit(κ_i − κ_j)`.

**Prediction.** After the macro model, the magnitude residual obeys `Δm_ij = β(κ_i − κ_j)` with:
- `β > 0`;
- the **same** `β` for minima and saddle images;
- the **same** `β` in radio, mid-IR, narrow-line and optical flux ratios;
- constant in time.

The time-delay residual follows `Δt_ij^{res}/Δt_ij ≈ (2δ_T/δ_L)·…`, i.e. it is proportional to the flux residual with the universal ratio `1/Q`.

**Why it is hard to mimic:**
- Subhalos and microlensing preferentially demagnify saddle images.
- Microlensing is chromatic and changes with time.
- Dust is chromatic.
- A mass-sheet transform with `λ<1` *shortens* delays and *brightens* images (`μ→μ/λ²`). Bridges with `δ_T<0` shorten delays and *dim* images. **The sign of the flux–delay correlation is opposite to the mass-sheet case.**

**Data:**
- public time-delay quasar lenses with published delays and lens models;
- radio, mid-IR (including JWST MIRI) and narrow-line flux ratios of quads;
- the handful of multiply imaged SNe that have both delays and standardisable fluxes.

**Kill:** any one of the following.
- `|β| < 0.02 mag` per unit Δκ (95%) across ≥30 quads in microlensing-immune bands.
- `β` differs significantly between saddles and minima once subhalo models are marginalised.
- `β` differs between bands.
- In multiply imaged SNe Ia, the flux–delay residual correlation has the mass-sheet sign.

### P2. Messenger-universal Etherington violation on a flat kernel
**Prediction.** `η_E(z) ≡ d_L/((1+z)² d_A) − 1 > 0`, and it is the same for SNe Ia and for gravitational-wave standard sirens at matched z. Per object, SN Hubble residuals correlate with a *flat-kernel* foreground column `T_flat = ∫(1+δ_g)^α dχ`, built from foreground galaxy catalogues, after controlling for the lensing-kernel convergence `κ_lens`. Lensing predicts `Δm ≈ −2.17κ_lens` (brighter). Bridges predict `+γ T_flat` (fainter), with no colour excess.

**Why it is hard to mimic:**
- Dust and axion-like photon loss affect only photons and are generally chromatic.
- Grey intergalactic dust would not dim gravitational waves.
- Selection effects do not prefer a flat-kernel column over a lensing-kernel one.

**Data:** public SN Ia compilations with host and sightline coordinates; photometric-redshift galaxy catalogues for foreground density; BAO `d_A(z)`; public bright and dark siren distance posteriors.

**Kill:** any one of the following.
- Sirens give `η_E^{GW}` consistent with 0 while SNe give `η_E^{EM} > 0` at the same precision (this rules out universality).
- The SN residual's flat-kernel coefficient `γ` is consistent with 0 at <0.005 mag per unit normalised column.
- The dimming is chromatic.

### P3. Static vs propagating mass, scaling as aperture squared
**Prediction.** X-ray hydrostatic, SZ-pressure-based and galaxy-velocity-dispersion masses are *all* below weak-lensing masses by the *same* fraction `κ_s(r/ℓ_b)²` at the same aperture. The deficit:
- is independent of the cluster's relaxation state;
- increases with environment density (through α);
- vanishes inside strong-lens Einstein radii (`r ≪ ℓ_b`), so lens-galaxy kinematics and lensing agree there.

**Why it is hard to mimic:** non-thermal pressure gives a larger bias in disturbed clusters and affects gas-based and galaxy-based dynamics differently.

**Data:** public SZ cluster catalogues, X-ray hydrostatic masses, weak-lensing masses, spectroscopic member velocities.

**Kill:** any one of the following.
- In relaxed clusters, the hydrostatic/lensing ratio is consistent with 1 at r_500 (±5%).
- The deficit tracks dynamical state rather than aperture.
- Galaxy-dynamics and gas-dynamics deficits differ significantly.

### P4. Transient connectivity: resonance-gated events
Conditional bridges open when both endpoints' neighbourhoods are in `R`. Coarse-grained, the local mode scale `Ω(x)` (the spectral gap of `T`, which I take to be tied to local record density) must lie near a universal value `ω_b`. Taking the lowest-order matching form,

`p_open(x) = Γ² / [(Ω(x) − ω_b)² + Γ²]`, `Ω ∝ ln ρ_m(x)`.

**Observable consequences:**
- (a) An image whose sightline crosses a patch that has just become resonant shows an **achromatic dip**. The depth is `ε cosφ_b × (number of open bridges)`, and the duration `t_open ~ Γ^{-1}` is fixed in the rest frame of the patch.
- (b) The power removed in dips comes back as achromatic **flashes** at other resonant sites. Global conservation means the summed fluence of all flashes equals the summed fluence of all dips.

**Statistical prediction:**
- The event rate per unit volume is `Rate ∝ n_b(x) p_open(x)`. It peaks as a Lorentzian in local density (in strong lenses, in surface density `κΣ_crit` at the image).
- Once local density is conditioned on, the rate is flat on the sky (no Galactic or ecliptic dependence beyond exposure).
- Observed durations scale exactly as `(1+z_patch)`, with a scatter that does not depend on source size or transverse velocity. Microlensing durations do depend on these, through `√M/v_⊥` and source size.

**Data:** archival multi-band light curves of lensed quasars and SNe; repeated JWST/HST deep-field imaging for flashes.

**Kill:** any one of the following.
- Achromatic short dips in lensed-quasar light curves show no dependence on image `κ`.
- After rescaling by `(1+z_l)`, their durations do not cluster more tightly than a microlensing model predicts (compare by KS test).
- No excess of achromatic dips over a chromatic-microlensing model across ≥20 well-sampled systems.

---

## 5. Quantities current searches never examine

1. **The per-object distance-closure invariant `Q = δ_L/δ_T`.** Searches based on lens morphology, microlensing curves or source counts treat flux, delay and size separately. Bridges predict that `Q` is the same in every system. To measure it, take each lens with ≥2 well-measured delays and microlensing-immune flux ratios, fit the macro model, and form the flux residual and the delay residual. Plotting one against the other should give a **single line through the origin, with a slope of universal sign across systems**. Under ordinary astrophysics the residuals are scattered and parity-dependent.
2. **The kernel-shape coefficient.** Regress Hubble residuals on two foreground columns at once: the geometric lensing kernel and a flat comoving kernel. No standard pipeline builds the flat-kernel column, since no standard mechanism predicts it.
3. **The rest-frame duration distribution of achromatic dips, rescaled by `(1+z_l)`.** Microlensing analyses fit durations to masses and velocities. They do not test whether durations collapse onto a single universal value.

---

## 6. Self-critique: the three weakest points

1. **Gravity is asserted, not derived.** I claim that matter modulating the `G₀` weights gives a curved emergent metric whose static and propagating limits agree, but I did not derive anything like Einstein's equations. The static–propagating split in P3 rests on an assumed difference between ω→0 and eikonal bridge response. A critic can fairly say that P3 is a Yukawa-like modification with extra words.
2. **The arrow, and causality, are fragile.** The ordering depends on a contingent low-entropy boundary (A5); it is not derived from the axioms. Conditional bridges are worse. If an agent could make a patch resonant on purpose, it could send signals across `ℓ_b` outside the cone, and causal paradoxes could follow. Ruling this out needs an extra no-signalling condition: bridges carry norm but no controllable phase. That condition is currently ad hoc.
3. **Too many parameters, and some degeneracies.** The parameters are `ε, φ_b, n̄_b, σ_b, ℓ_b, α, ω_b, Γ`. The CMB already forces α large and the effects small, and further tuning could shrink every signal below detectability. κ-proportional dimming is partly degenerate with lens-model mass normalisation. Only the sign and parity discriminants in P1 break that degeneracy, and the samples are small (few multiply imaged SNe).

---

## Resemblances I noticed (none were used as a starting point)

- **Statistical mechanics and Euclidean field theory:** Markov random fields and transfer matrices.
- **Network science:** small-world networks (shortcut edges).
- **Quantum Darwinism:** records defined through redundancy.
- **Photon-number-violating tests of Etherington reciprocity:** e.g. axion–photon mixing and grey dust. Mine differs because it claims messenger universality.
- **Massive-gravity or Yukawa modifications:** these resemble the static deficit.
- **Wormhole-like or ER=EPR-like nonlocal shortcuts:** bridges look superficially like these.
- **Emergent and entropic gravity programmes:** the idea that geometry comes from information.
- **The Hubble tension:** P2/P1 share the same direction, by coincidence of sign.
