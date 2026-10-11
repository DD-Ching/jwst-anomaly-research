# Goal round 8: can Einstein–æther's superluminal spin-0 / spin-1 modes carry a controllable bit?

Owner brief (2026-10-11). The Krasnikov line stays paused and G1 is not started. The new bounded unit is
Einstein–æther (EA) spin-0 and spin-1 modes, approached from the action, the matter coupling and a measurable
response:
1. applicability screen;
2. sender operation → propagation → receiver readout;
3. one minimal computation, after reproducing published controls;
4. separate value (i), a testable signal beyond known physics, from value (ii), seconds-scale interstellar
   latency;
5. a definite judgement.

Code `scripts/goal_ea_response.py`; results `results/goal_round8/ea_response.json`; tests
`tests/test_goal_ea_response.py`.

Labels: [lit] = stated in the cited paper (with its own status: proved, data-constrained or conjectured);
[derived] = computed in this round; [ASSUMPTION] = a choice made here.

**Repository check.** EA appears only in round 5's "not covered" list and as a data citation (2605.01436) in
docs/hypotheses/round-3/B-on-R3B.md. No earlier communication analysis exists. This differs from the Bancal
hidden-influence model (rounds 2–3): here a definite Lagrangian fixes both the source coupling and the readout.

## 1. Applicability screen [lit]

**Action** (signature −+++, matter coupled to g only):

S = (1/16πG)∫√−g [R − c1(∇u)² − c2(∇·u)² − c3 ∇_a u_b ∇^b u^a + c4 a²] + λ(u² + 1)

Conventions differ across papers: signature, the sign of the c4 term, and c_a = c14, c_σ = c13,
c_θ = c13 + 3c2, c_ω = c1 − c3 in 2304 and 2605. The mode speeds below were checked against our own
derivation (§3).

| Item | Statement | Status |
|---|---|---|
| Speeds | c_T² = 1/(1 − c13); c_V² = (2c1 − c1² + c3²)/(2c14(1 − c13)); c_S² = c123(2 − c14)/(c14(1 − c13)(2 + c13 + 3c2)) (Jacobson–Mattingly gr-qc/0402005; 1802.04303 Eq. 3.4) | derived in the papers; reproduced here |
| Tensor speed | \|c13\| ≲ 10⁻¹⁵ (GW170817) | data |
| PPN | \|α1\| ≤ 10⁻⁴, \|α2\| ≤ 10⁻⁷, with α1 = −4c14 for c13 = 0 (1802.04303 Eq. 1.2) | data |
| Allowed region | (i) 0 < c14 ≤ 2 × 10⁻⁷ with c14 ≲ c2 ≲ 0.095; (ii) 2 × 10⁻⁶ ≲ c14 ≲ 2.5 × 10⁻⁵ with c2 − c14 ≲ 2 × 10⁻⁷ (1802.04303 Eq. 3.23). Region (ii) has c_S ≈ 1. **Region (i) lets c_S ≈ √(c2/c14) → ∞.** c_ω, and hence c_V, is unconstrained (2304.06801; 2605.01436). | data plus stability and Cherenkov assumptions |
| Cherenkov | Constrains only subluminal modes; "c_N > 1 is allowed. This does not violate causality" (2304.06801 fn. 9) | theory |
| GW data | GW170817 and GW190425 posteriors on c_a and c_θ follow the priors (2304.06801) | data, uninformative |
| Pulsar | PSR J1738+0333 gives α1 > −7.2 × 10⁻⁵ at 90 % (2605.01436v1) | data |
| Superluminal exclusion | 2505.09544 (v3, Oct 2026; PRD 112, 024052) conjectures that gravitational memory with superluminal modes excludes the superluminal parameter space. In its own words this is a "conjecture". It has no quantitative amplitude and neglects finite-size and nonlinear effects. No rebuttal among 13 INSPIRE citations. Zosso 2610.06528 (an EM analogue) finds that finite superluminal sources give finite bursts. | conjecture |
| EFT range | No paper read here gives a strong-coupling or cutoff scale as c14 → 0. At c14 = 0 the æther equation loses its time derivatives (instantaneous limit; 1802.04303 fn. 4). | open |
| Causality | Signals are causal with respect to the widest cone, i.e. the æther frame. A global time is "only guaranteed whenever the æther field is hypersurface orthogonal … in EÆ theory it is not true in general" (2505.09544 App. C). | stated in the paper |

**Screen result.** Data do **not** bound c_S or c_V from above. A superluminal spin-0 mode (region i) is not
excluded by data; the exclusion that exists is a conjecture. "Not excluded" is not the same as viable, so the
communication question is examined next.

## 2. The communication chain

- **Sender.** The sender can only move ordinary matter, which couples to the metric alone. The source is a
  conserved T^{μν}. Its mass monopole and dipole are fixed by energy and momentum conservation, so the lowest
  radiating scalar moment is Q = I_kk = ∫ρr² (and T^{ij} via the virial identity ∫T^{ij} = ½Ï_ij).
- **Propagation.** Linearised EA about the æther rest frame has three cones: c_T = 1 (with c13 = 0), c_V and c_S.
- **Receiver.** A freely falling detector at rest in the æther frame reads the tidal tensor E_ij = R_{0i0j}
  through geodesic deviation. At linear order this is a coordinate- and gauge-invariant local quantity.
- **Timing.** A component on the c_S cone arrives at D/c_S, before light (D/c), if c_S > 1. Three different times
  must be kept apart:
  - the mode speed (here the characteristic speed, so also the front speed, because the linearised equations are
    hyperbolic);
  - the arrival of the front;
  - the time to accumulate enough signal-to-noise ratio to read the bit.
- **Frequency range.** All results are low-frequency EFT statements (frequencies ≪ any cutoff). Nothing is
  extrapolated to high frequency.

## 3. Calculation [derived; reproduces published controls]

Quadratic EA action in Fourier space, plane wave along z, gauge v^i = 0 = h_0z, general conserved source
(ρ, j_⊥, transverse stresses; longitudinal parts fixed by conservation). The response is the gauge-invariant
E_ij(ω, k).

**Positive controls** (all reproduced exactly; `controls()`):
- The poles of the response sit at the published c_S², c_V² and c_T².
- The static limit gives E_zz = 4πG_N ρ with G_N = G/(1 − c14/2).
- The tensor far field is E_xy = −Gω⁴I_xy/r, as in GR, for c13 = 0.
- **Yagi, Blas, Barausse & Yunes (1311.7144v4) Eq. 106**, breathing polarisation, weakly gravitating source, c13 = 0:
  - our I_kk coefficient equals theirs exactly, including their Z = 2(c14 + 3c2)/(3c2(2 − c14));
  - our n̂n̂Q̈ coefficient is **½ of theirs**. That term is O(c14²) and unresolved. It does not affect any
    conclusion below.
- **Vector sector.** With c13 = 0, E_xz from a current has no c_V pole. The spin-1 mode does not reach a metric
  detector from a matter source. This agrees with 2304.06801 (β1 = −2c_σ/c_V = 0).

**New results:**
1. **Exact exterior field of a pulsating spherical source** (Q = ∫ρr², c13 = 0). The exterior tidal field
   contains only the c_S pole; the remainders are local contact terms, so there is no instantaneous piece.
   Spherical sources also have no exterior field in GR (Birkhoff), so the whole exterior change is spin-0:

   E_ij(t, r) = −κ G [δ_ij Q⁽⁴⁾(t_S)/r + ∂_i∂_j(Q̈(t_S)/r)], with t_S = t − r/c_S and
   κ = c14(c14 + 3c2)/(3c2(2 − c14)²) = c14/4 + O(c14²).

2. **General slow source, far field:** E^S_ij = −(c14/4) Gω⁴ I_kk δ_ij / r + O(c14²).
3. **Trade-off.** c_S² = c2(2 − c14)/(c14(2 + 3c2)) gives c14 ≈ 2c2/((2 + 3c2)c_S²) ≤ 0.083/c_S² for
   c2 ≤ 0.095. **The faster the spin-0 mode, the weaker it couples to matter, as 1/c_S².** This holds in both the
   near zone and the far zone, since both terms carry κ. Nothing diverges as c14 → 0 in the linear theory
   (consistent with the fluxes in Foster gr-qc/0602004 and Yagi 1311.7144).

**Link budget** (ASSUMPTION, deliberately generous):
- Sender: a 10¹² kg shell (a large mountain) at R = 1 km, pulsating by ±10 m at 1 Hz, so ΔQ = 2 × 10¹⁶ kg m².
- Receiver: tidal noise 10⁻²⁰ s⁻²/√Hz.
- SNR 5.
- c_S = D/latency, with c14 at its largest allowed value (region i).

| D | Latency | c_S needed | κ | Tidal signal | Integration time for SNR 5 |
|---|---|---|---|---|---|
| 20 ly | 1 s | 6.3 × 10⁸ | 5 × 10⁻²⁰ | 7 × 10⁻⁶² s⁻² | 5 × 10⁸³ s |
| 20 ly | 1 day | 7.3 × 10³ | 4 × 10⁻¹⁰ | 5 × 10⁻⁵² s⁻² | 9 × 10⁶³ s |
| 20 ly | 1 yr | 20 | 5 × 10⁻⁸ | 7 × 10⁻⁵⁰ s⁻² | 5 × 10⁵⁹ s |
| 100 ly | 1 s | 3.2 × 10⁹ | 2 × 10⁻²¹ | 6 × 10⁻⁶⁴ s⁻² | 8 × 10⁸⁷ s |

The age of the universe is 4 × 10¹⁷ s.
- Even GR's own quadrupole field from this source at 20 ly (1.4 × 10⁻⁴² s⁻², the κ = 1 scale) would need about
  10⁴⁴ s.
- The barrier is the weakness of gravitational coupling to any controllable source. That is shared with GR, and
  EA adds a further factor κ ≤ 0.083/c_S².

## 4. The two kinds of value, kept separate

**(i) A testable signal beyond known physics.** Yes, conditionally. It is already a published prediction, not new
here.
- For c14 > 0, EA predicts a breathing (scalar) tidal component of amplitude ≈ (c14/4)GÏ_kk-type/r from sources
  whose I_kk changes (e.g. core collapse, eccentric or asymmetric systems). Circular binaries have constant I_kk at
  leading order.
- It arrives D(1 − 1/c_S) before the tensor wave.
- Existing tests (2304.06801; pulsar timing 2605.01436) bound c14 through α1. The amplitude is ≤ 2.5 × 10⁻⁵ of the
  GR scale, and ≤ 5 × 10⁻⁸ in the superluminal region (i).
- **Observable prediction:** a breathing-polarisation precursor to a GW event, with lead time D(1 − 1/c_S), at a
  fraction ~c14/4 of the GR amplitude.
- **Main false signals:** non-Gaussian detector transients (glitches); polarisation leakage of tensor noise into
  the breathing null-stream; the unknown lead time (for c_S ≫ 1 the precursor arrives ~D/c before the event: years
  to millennia for Galactic sources), which makes a blind search with an enormous trials factor.
- **Smallest verification step:** a breathing null-stream search in archival multi-detector data at the scalar
  frequency band of one well-localised event, scanning lead times. Its sensitivity is capped at ~c14/4 ≲ 5 × 10⁻⁸
  of the event's GR strain, so it is not worth running for this question.

**(ii) Seconds-scale interstellar latency.** **Excluded under the stated assumptions:**
- matter coupled minimally to the metric;
- linearised EA in the allowed region with c13 = 0;
- a detector reading the tidal tensor;
- any controllable source up to mountain scale.

The required c_S forces κ ≤ 0.083/c_S², so the signal sits 40 or more orders of magnitude below detectability even
at a latency of one year. The trade-off is quantitative: raising c_S by 10 lowers the signal by 100 and multiplies
the integration time by 10⁴.

## 5. Judgement

**Excluded under the stated assumptions** as a low-latency channel. The arrival time of the front is not the
obstacle: the c_S front does reach the receiver before light. The obstacle is that the amplitude a matter source
can put into that front falls as 1/c_S².

What could change this (inputs not supplied here, and not pursued):
- a non-minimal matter–æther coupling. That would be a new model, constrained by Lorentz-violation tests.
- a nonlinear or strong-coupling regime as c14 → 0. No controlled theory of it exists in the sources read.
- a resolution of the Taherasghari–Will (2506.03843) dipole dispute. It concerns only bodies with sensitivities,
  i.e. strongly self-gravitating ones, which a sender cannot build.

No candidate passes, so no further theory chain or large computation follows.

## 6. Bookkeeping

- **Read** (full text via a literature subagent; key formulas re-checked here):
  - 1802.04303v3, 2304.06801v2, 2505.09544v3, 2605.01436v1;
  - Foster gr-qc/0602004v5; Jacobson–Mattingly gr-qc/0402005v2; Yagi et al. 1311.7144v4.
- **Abstract or citation level only:** Zosso 2610.06528; Taherasghari & Will 2506.03843.
- **New this round** [derived]: the gauge-invariant tidal response with a general conserved source; the exact
  spherical-source exterior field; the κ ∝ 1/c_S² trade-off; the link budget.
- **Reproduced:** speeds, G_N, the GR quadrupole field, and Yagi's breathing I_kk term.
- **Unresolved:** a factor ½ on the O(c14²) n̂n̂Q̈ term relative to Yagi Eq. 106.
