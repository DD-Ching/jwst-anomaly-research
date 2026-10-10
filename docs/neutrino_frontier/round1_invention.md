# Neutrino Frontier — Round 1: Invention (team N-A)

Literature-blind round (no web or literature search; no citations). A separate adversarial team checks
prior art afterwards; any framework here may turn out to be known. Neutrinos are probes only. Nothing
here is a detection or a claim of new physics.

Status labels used on statements: **[A]** assumption, **[D]** derived (algebra/arithmetic from stated
inputs), **[P]** model_prediction, **[H]** hypothesis. Numbers are order-of-magnitude unless stated.

Common inputs [A] (standard values from memory, not re-checked this round):
Δm²₂₁ ≈ 7.4×10⁻⁵ eV², |Δm²₃₁| ≈ 2.5×10⁻³ eV²; oscillation phase φ = 1.27 Δm²[eV²] L[km] / E[GeV].
SN1987A: D ≈ 51 kpc (D/c ≈ 5.3×10¹² s), ~24 events, E ≈ 7–40 MeV, spread ≈ 13 s, neutrinos arrived
~hours before the optical rise (allowed extra ν–γ offset |Δt| ≲ 10⁴ s). TXS 0506+056: z ≈ 0.34, take
D ≈ 1.4 Gpc (D/c ≈ 1.4×10¹⁷ s); 2017 alert (~290 TeV) during a months-long γ flare (window ≲ 10⁷ s);
2014–15 neutrino flare ≈ 13 events over ≈ 110 d at ~10–100 TeV. ICECAT-1: ~350 tracks / ~12 yr
(alert rate ≈ 9.2×10⁻⁷ s⁻¹ all-sky). 10-yr track sample: ~1.1×10⁶ events. Overlap GBM bursts ≈ 3,500.
GWTC events overlapping ICECAT-1 ≈ 150.

Mechanisms are chosen to be structurally different: (F1) stochastic distance, (F2) hidden partner
states, (F3) a second causal cone, (F4) nonlocal topology. None is a variant of energy-dependent speed.

---

## F1 — Granular locality (stochastic distance)

**Core idea [H].** If spacetime and distance emerge from a discrete or fluctuating substrate, "distance"
is not sharp: a propagating quantum picks up a random walk in direction and path length whose size
grows with how far it travels. The distinctive mark is a **time spread ∝ D²** (random-walk path
excess), not the ∝ D of any speed change. Neutrinos are good probes because nothing else scatters them.

**Assumptions.**
1. [A] Each propagation step of length ℓ adds an independent random direction kick; the kicks are
   isotropic about the instantaneous direction (no preferred frame on average).
2. [A] Kick variance per length depends on energy as a power law, index n (sign free).
3. [A] Photons are affected far less (e.g. the substrate couples to weak-isospin or to the neutrino's
   long coherence); otherwise gamma-ray image sharpness already rules out the region of interest.
4. [A] No absorption, no flavor effect beyond the path-length jitter.

**Minimal math [D from A1–A2].**
- Angular diffusion: dθ²/dL = κ (E/E₀)^{2n} / L₀, so θ²_rms(E,D) = κ (D/L₀)(E/E₀)^{2n}.
  Free parameters: κ (dimensionless), n (dimensionless); fixed reference L₀ = 1 Gpc, E₀ = 100 TeV.
- Arrival delay (small-angle multiple scattering): Δt ≈ c_g D θ²_rms / c, with c_g = O(¼) geometric
  factor; hence **Δt ∝ κ D² E^{2n}** and the burst is also *spread* by ~Δt.
- Phase jitter on oscillations: δL ≈ D θ²/4, δφ = 1.27 Δm² δL/E (same units as above).

**Recovery of known physics.**
- Three-flavor oscillations and Δm² values: lab baselines (≤ 10⁴ km) give δL ≲ 10⁻²⁰ km for any κ
  allowed below; δφ ≈ 0 [D]. Standard oscillations are untouched.
- Special relativity: mean speed is c; only the *path* is longer; local Lorentz invariance holds on
  average (A1). Lab time-of-flight tests see Δt ≈ 0 because Δt ∝ D² [D].
- SN1987A: spread ≲ 10 s ⇒ θ²_rms(15 MeV, 51 kpc) ≲ 4×10/5.3×10¹² ≈ 7.5×10⁻¹² [D].
- TXS 0506+056 (stronger for all n ≥ 0): coincidence within ≲ 10⁷ s ⇒ θ²_rms(TXS) ≲ 4×10⁷/1.4×10¹⁷
  ≈ 2.8×10⁻¹⁰, θ_rms ≲ 1.7×10⁻⁵ rad ≈ 3.4″ [D]. For n < 0, SN1987A dominates instead.
- Allowed region [D]: κ (E/100 TeV)^{2n} ≲ 2×10⁻¹⁰ at ~30–300 TeV; and the SN1987A bound
  scaled by (51 kpc/1 Gpc)(15 MeV/100 TeV)^{2n} for n < 0.

**Distinguishable prediction [P].**
- Angle: maximal allowed broadening ≲ 4″ ≪ IceCube track PSF (~0.5–1°): **undecidable in angle** with
  any listed dataset [D]. Stated honestly: F1 cannot be seen as source "fuzz".
- Time: minimum flare width τ_min(D,E) = τ_TXS (D/1.4 Gpc)² (E/E_TXS)^{2n}. At the edge of the allowed
  region τ_TXS ≈ 10⁷ s; then a BNS merger at 40 Mpc would have τ_min ≈ 10⁷ × (40/1400)² ≈ 8×10³ s at
  30 TeV, and a GRB at z ≈ 2 (D ≈ 5 Gpc) τ_min ≈ 1.3×10⁸ s ≈ 4 yr [D].
- Consequence [P]: standard ±500 s GRB–neutrino stacking (GBM × ICECAT) loses its signal for distant
  bursts while a GW–neutrino search with ±10⁴ s windows does not. The test is the **D²** scaling of
  the best-fitting coincidence window across source classes with known distances (GWTC distances, GBM
  bursts with redshifts, TXS). Decidability: requires ≥ 2 neutrino-detected transients at different D;
  current public data hold one (TXS). Background for a GWTC × ICECAT spatial (combined 5°) match in a
  10⁴ s window: 350 × 150 × 1.9×10⁻³ × 10⁴ / 3.8×10⁸ ≈ 2.6×10⁻³ pairs — clean, but signal unknown.

**Failure conditions.** A neutrino transient at D ≳ 1 Gpc with width ≪ 10⁷ s at ≳ 100 TeV kills the
n ≥ 0 region at the TXS edge; a sub-second Galactic-SN burst at ~10 kpc kills n < 0 at ~10⁻¹³ level
(not in the public list). Any gamma-ray image broadening that tracks the neutrino one kills A3.

**Ladder: C** (consistent in the stated region). D only in principle (needs more neutrino transients).

---

## F2 — Twin-sheet flavor (hidden partner states)

**Core idea [H].** Flavor identity is incomplete: each neutrino mass state νᵢ has a hidden partner νᵢ′
living on a "second sheet" of the emergent spacetime (a hidden degree of freedom with no gauge
charges). The pair is split by a tiny δm²ᵢ, so over cosmological baselines a neutrino oscillates
half-way into the sheet and becomes invisible. Locally, nothing changes.

**Assumptions.**
1. [A] Each νᵢ is maximally mixed with its partner νᵢ′ (mixing angle π/4); the partner is sterile.
2. [A] Splittings δm²ᵢ ≪ |Δm²₂₁|; simplest case equal δm² for i = 1,2,3 (then flavor ratios unchanged).
3. [A] Partners do not decay and do not interact with matter at lab energies.

**Minimal math [D].**
- Survival of the visible component: Pᵢ(E,L) = 1 − sin²(1.27 δm²ᵢ L/E) (units as above); averaged
  over source size/energy bins → ½.
- Free parameters: δm²ᵢ (eV²), three numbers (one in the minimal case).
- First (highest-energy) flux dip of a source at distance L: E₁ = 1.27 δm² L/(π/2) ≈ 0.81 δm² L [GeV,
  eV², km]. Redshift integration shifts this by O(1+z) [A, not worked here].

**Recovery of known physics [D].**
- Three-flavor oscillations: Δm²₂₁, Δm²₃₁ are untouched (A2). Solar (L ≈ 1.5×10⁸ km, E ≈ 1 MeV):
  phase = 1.9×10¹¹ δm² ≪ 1 for δm² ≲ 10⁻¹³ eV². Reactor/atmospheric/accelerator: even smaller.
- Special relativity: standard massive dispersion; no Lorentz violation.
- SN1987A: phase = 1.27 δm² × 1.6×10¹⁸ / 0.015 ≈ 1.4×10²⁰ δm² ⇒ averaged to ½ for δm² ≳ 10⁻¹⁹ eV².
  A factor-2 deficit in total emitted energy is within the ~factor-2–3 uncertainty of the ~24-event
  energetics; timing unchanged. Compatible window: **10⁻²⁰ ≲ δm² ≲ 10⁻¹³ eV²** (plus δm² → 0).

**Distinguishable prediction [P].**
- TXS 0506+056 (L ≈ 4.3×10²² km) at 100 TeV: phase = 5.5×10¹⁷ δm² [D]; NGC 1068-like source
  (14 Mpc, 4.3×10²⁰ km) at 1 TeV: also 5.5×10¹⁷ δm² [D]. Both probe **δm² ≈ 10⁻¹⁹–10⁻¹⁷ eV²**.
- Per-source spectrum: a dip to 0 at E₁ ≈ 0.81 δm² L (single δm²), averaging to ½ below ~E₁/3. For
  δm² = 2×10⁻¹⁸ eV², TXS E₁ ≈ 70 TeV, i.e. inside the 2014–15 flare band.
- Size vs resolution [D]: IceCube track energy resolution ≈ factor 2 (muon-energy proxy), so the dip
  (width ≈ factor 3 in E) is smeared to a ≈ 30% depression; with ≈ 13 TXS flare events the Poisson
  error per half-decade bin is ≳ 50%. **Undecidable with current public single-source statistics.**
- Diffuse: flux suppressed by ½ below E_c ~ 0.81 δm² L_H, unsuppressed far above ⇒ apparent index
  hardening Δγ ≈ log₁₀2 / 1.5 ≈ 0.2 across ~1.5 decades [D]; comparable to current index uncertainty
  and to astrophysical model freedom ⇒ not decidable alone.
- With unequal δm²ᵢ the visible flavor ratio becomes energy-dependent [P]; ICECAT tracks carry no
  flavor, so not testable with the listed data.

**Failure conditions.** A resolved single-source spectrum with no dip between ~0.3 and ~3 × E₁ (at
≥ 5× current statistics) excludes the corresponding δm²; a full-strength flux from a future Galactic
SN matching model energetics to ≲ 20% excludes δm² ≳ 10⁻¹⁹ eV² (not public now).

**Ladder: C.** Testability D requires either a stacked distance-binned spectral analysis of many
sources or flavor data; with the listed public data it is not decidable.

---

## F3 — Second causal cone (neutrino-specific geometry)

**Core idea [H].** Neutrinos propagate on their own effective metric g̃, built from the same matter
but with a slightly different coupling to the gravitational potential. Far from mass, g̃ = g (special
relativity holds); near mass the neutrino cone is slightly wider or narrower than the light cone. This
is nonstandard causal structure: two cones; for ε < 0 neutrinos are locally "warp-like" (outrun light
through potential wells) without any energy dependence.

**Assumptions.**
1. [A] g̃_μν = g_μν − ε (2Φ/c²) u_μ u_ν (weak field, Φ Newtonian potential, u the matter rest frame).
2. [A] ε is a universal dimensionless constant, energy-independent and flavor-universal.
3. [A] Photons and gravitational waves follow g.
4. [A] No closed causal curves: |ε| ≪ 1 keeps the union of cones globally hyperbolic.

**Minimal math [D from A1].**
- Neutrino Shapiro delay relative to photons: Δt_νγ = ε T_S(n̂, D), with T_S = (2/c³)∫|Φ| dl the
  photon Shapiro delay along the line of sight n̂. Free parameter: ε (dimensionless; sign free).
- Milky Way T_S ≈ 10⁷ s for extragalactic directions (M ~ 5×10¹¹ M_⊙, GM/c³ ≈ 2.5×10⁶ s, × 2 × log
  factor) and **varies with direction** by O(1) (largest towards the Galactic centre) [D, rough].
- Extragalactic structure along 1.4 Gpc: rms T_S ~ 2 × 3×10⁻⁶ × √(1400 × 30 Mpc)/c ≈ 10¹¹ s (random
  sign, coherence 30 Mpc assumed) [D, rough].

**Recovery of known physics [D].**
- Oscillations and Δm²: g̃ changes all mass states' paths identically ⇒ no phase difference
  (A2); oscillations standard.
- Special relativity: in a lab on Earth the local Φ/c² ≲ 10⁻⁶ (Galactic) and its gradient is tiny;
  speed anomaly |v−c|/c ≲ ε × 10⁻⁶, below any time-of-flight test for |ε| ≲ 1.
- SN1987A: |ε| × 10⁷ s ≲ 10⁴ s ⇒ |ε| ≲ 10⁻³.
- TXS 0506+056: |ε| × (10⁷ + 10¹¹) s ≲ 10⁷ s ⇒ |ε| ≲ 10⁻⁴ (if the extragalactic rms holds; otherwise
  10⁻³ from the Milky Way term alone).

**Distinguishable prediction [P].**
- Transient neutrino–photon lags are **direction-dependent and energy-independent**: Δt(n̂) = ε T_S(n̂).
  For |ε| = 10⁻⁴ and a Galactic-only term, lags ≈ 10³ s (range ~3×10² to 3×10³ s across the sky);
  standard ±500 s searches would partially miss them.
- Test: scan ICECAT × GBM pairs within a 5° combined radius over lags ±10⁴ s; fit one parameter ε with
  each pair's expected lag ε T_S(n̂). Background [D]: 350 × 3,500 × 1.9×10⁻³ / 3.8×10⁸ ≈ 6×10⁻⁶ pairs
  per second of lag ⇒ ≈ 0.12 chance pairs over the whole ±10⁴ s scan; ≈ 6×10⁻³ in a 1,000 s band.
  So **2 pairs lying on one ε line would be significant (p ≈ 10⁻²)**; 1 pair is suggestive only.
- Same with GWTC (distance known, extragalactic term computable from large-scale-structure maps):
  background over ±10⁴ s ≈ 350 × 150 × 1.9×10⁻³ × 2×10⁴ / 3.8×10⁸ ≈ 5×10⁻³ pairs (GW sky maps are
  often ≫ 5°, so this is optimistic by up to ~10×).
- Signal size is unknown: it equals the number of GRB/GW neutrinos IceCube detects, which current null
  results cap at a few. The test is decidable in background terms, not guaranteed in signal.

**Failure conditions.** A neutrino–photon coincidence from a source behind a large T_S (e.g. a
Galactic-centre-direction transient) with lag ≪ ε T_S for the claimed ε; two coincidences whose lags
need ε of opposite sign; any energy dependence of the lag at fixed direction (kills A2).

**Ladder: C, with a D-level test** (the lag scan is runnable on ICECAT-1 + GBM now; outcome likely a
bound on ε in the ~10⁻⁴–10⁻³ range, conditional on a few true coincidences existing).

---

## F4 — Sparse nonlocal links (wormhole-like topology)

**Core idea [H].** If locality is emergent from a network of connections, the network may contain rare
"long links": two distant regions joined by a short route, like a wormhole with two mouths. Suppose
the links admit only neutral, colorless, very weakly coupled quanta (neutrinos; possibly gravitons).
A neutrino from a source in direction n̂₁ can then reach us from a different direction n̂₂ — a
**ghost image** — with a delay set by the mouth geometry rather than by energy.

**Assumptions.**
1. [A] Links are static in the cosmic rest frame (prevents closed causal curves — see ladder).
2. [A] A neutrino crossing a mouth of cross-section σ_M passes through with probability ≈ 1; photons,
   charged particles and hadrons do not (no gamma-ray ghost images).
3. [A] Mouth pairs are randomly placed; the ghost direction n̂₂ is uncorrelated with n̂₁.
4. [A] The relevant links have one mouth near the observer side at distance r_M; ghost delay
   Δt_g ≈ r_M/c (O(1) geometric factor), energy-independent, flavor-independent.

**Minimal math [D].**
- Ghost ratio per source R_g = P_h/(1 − P_h), with P_h = n_M σ_M D_eff the hop probability (n_M mouth
  number density [m⁻³], σ_M [m²], D_eff effective path [m]). Free parameters: P_h (or n_M σ_M) and r_M [m].
- Ghost flux = R_g × direct flux, same spectrum, same light curve shifted by Δt_g, at random n̂₂.

**Recovery of known physics [D].**
- Oscillations: links do not act on flavor; the hopped neutrino still averages to the standard
  astrophysical flavor mix; lab baselines have negligible P_h (n_M σ_M L ≪ 1).
- Special relativity: ordinary propagation is unchanged; only the topology is altered.
- SN1987A: ghosts arrive from other directions after ~r_M/c; the ~24-event burst is the direct image,
  so P_h(51 kpc) ≲ 0.5 (no visible deficit beyond energetics uncertainty).
- Point-like associations (TXS, steady hotspots) and diffuse isotropy: ghosts dilute sources and add
  an isotropic component, so any P_h ≲ 0.5 at Gpc paths is compatible.

**Distinguishable prediction [P].**
- Ghost of the TXS 2014–15 flare: if r_M/c ≪ 110 d, a second cluster of ≈ 13 R_g events with the same
  spectrum (index ≈ 2.2) appears somewhere in the 10-yr track sample in the same 110 d window.
  Background [D]: 1.1×10⁶ × (110/3,650)/41,253 deg² ≈ 0.8 events deg⁻², i.e. ≈ 2.5 in a 1°-radius
  disc; ≈ 1.3×10⁴ independent discs ⇒ need p ≲ 8×10⁻⁷ for a ~1% global false alarm ⇒ ≥ 14 events
  (excess ≈ 11.5 ⇒ **R_g ≳ 0.9**). With energy weighting cutting background ~10× (to ~0.25): ≥ 6
  events ⇒ **R_g ≳ 0.45** decidable. Smaller R_g is not decidable with this sample.
- Transient ghosts (r_M ≲ 1.5×10¹¹ m, Δt_g ≲ 500 s): any-sky ICECAT × GBM coincidences in ±500 s:
  background 350 × 3,500 × 10³ / 3.8×10⁸ ≈ 3.2 pairs; ICECAT × GWTC: ≈ 0.14 pairs [D]. Signal ≈ R_g ×
  (direct coincidences), and direct ones are near zero ⇒ **undecidable** unless some sources are seen
  only as ghosts.
- If gravitons also pass (A2 variant): GWTC would contain repeat events with matching intrinsic
  parameters at unrelated sky positions with lag ≈ r_M/c [P].

**Failure conditions.** No ghost of the TXS 2014–15 flare at the energy-weighted threshold excludes
R_g ≳ 0.5 for r_M/c ≲ 0.3 yr; a gamma-ray ghost of any bright flare kills A2; a neutrino ghost
arriving *before* its direct image kills A1/A4 (it would signal causal-order violation, a stronger and
separately vetted claim).

**Ladder: B (fragile).** Self-consistency rests on A1: moving mouths generically allow closed causal
curves, and A2's selective transmission has no derivation. Its TXS-ghost test is runnable (D-type), but
the framework does not pass B cleanly, so a null result would constrain parameters, not the idea.

---

## Summary table

| ID | Mechanism | Key free parameters | Sharpest public-data handle | Effect vs resolution | Ladder |
|----|-----------|--------------------|-----------------------------|----------------------|--------|
| F1 | stochastic distance (Δt ∝ D²) | κ, n | flare widths vs distance (TXS; future GW/GRB ν) | angle ≲ 4″ vs ~1° PSF: undecidable; time: needs ≥ 2 transients | C |
| F2 | hidden partner states | δm² (10⁻¹⁹–10⁻¹⁷ eV²) | per-source spectral dip at E₁ ≈ 0.81 δm² L | ~30% smeared dip vs ≳ 50% stat error: undecidable now | C |
| F3 | second causal cone | ε (abs ε ≲ 10⁻⁴–10⁻³) | direction-dependent ν–γ lag scan, ICECAT × GBM | 0.12 chance pairs over ±10⁴ s: decidable if ≥ 2 true pairs | C (D test) |
| F4 | nonlocal links / ghosts | P_h, r_M | ghost of TXS 2014–15 flare in 10-yr sample | needs R_g ≳ 0.45–0.9 | B |

Next step for the program (not done here): the adversarial team checks prior art for each mechanism;
F3's lag scan and F4's TXS-ghost search are the two that are cheap to run on public data.
