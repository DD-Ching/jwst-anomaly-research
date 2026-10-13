# X2 S4: a matter-consistent axial gravitational equation in place of the test-field proxy

Owner request 2026-10-13. Question: can gravitational perturbation equations with an explicit supporting matter
and its perturbations replace X2's test-field proxy (premise iii / assumption P1) and narrow one model
uncertainty? Budget: stage 1 ≤ 2 h, stage 2 ≤ 3 h, from X2's remaining 10.5 h. S3's failed spacing gate
(4L + 4.1 %) stays on record. No rotating mergers, scans or multi-event work.

## Stage 1: literature check (done before any stage-2 computation)
Neither paper was assessed before (no hit in DECISIONS, SOURCES, CHANGELOG or the X2 / S0 files).

**Yang & Lan 2026, arXiv:2603.20594v1** ("Branch-dependent ringdown in black-bounce spacetimes").
- **Scope.** Axial (odd-parity) gravitational perturbations only. Static, non-rotating SV, ℓ = 2. Free
  evolution of a Gaussian packet on a fixed background, with QNMs from Prony fits. Not a merger. Polar
  perturbations are left to future work.
- **Two matter interpretations of the same metric:**
  - **(AF) anisotropic fluid.** δT_(1)(3) = δT_(2)(3) = 0 because the fluid's u and x vectors are kept
    unperturbed. That is a closure by assumption (no action; axial fluid degrees of freedom frozen). The result is
    one channel, V_AF = f μ²/Σ² + Σ (1/Σ)_{,r*r*}, with μ² = (ℓ−1)(ℓ+2).
  - **(NED) nonlinear electrodynamics + scalar.** Action in §2.2. The scalar decouples from the axial sector by
    parity. The axial sector is a closed two-channel system: H₂ (gravitational) coupled to H₁ (EM), with
    potentials V₁₁, V₁₂ = V₂₁/4π and V₂₂ = V_AF (Eq. 56). These depend only on G⁰₀ − G²₂, so the electric and
    magnetic sources give the same matrix.
- **Wormhole-branch results (a > 2M).** Published QNMs for a = 2.13–10 M. In 2.01–2.12 M the fits are
  "unidentified" or jump between modes. There are no values for 2 < a < 2.01, which is where X2's
  GW150914 cells lie (a − 2M ≈ 2e-4 to 7e-4).
- **Echo statements** (NED echoes suppressed; distinct delays) are qualitative only. Nothing is computed.

**Bronnikov & Walia 2022, PRD 105, 044039 (arXiv:2112.13198).**
- Any static, spherically symmetric SV metric is sourced by a **phantom** scalar (φ′² = a²/(x² + a²), so ε = −1
  is required) plus a **magnetic** NED field. For Schwarzschild-SV, L(F) = (12/5) M a²/(x² + a²)^{5/2}.
- An electric version needs fractional powers of F < 0 and is ill-defined. The magnetic model is the
  well-defined field theory, and by the duality above it carries the same axial equations.
- Stability "cannot be formulated until we know their precise dynamic context". It may differ between sources.

**Checks we ran** (`scratchpad` sympy / numpy; reproduced in `scripts/x2_axial_ned.py` tests):
- Background: an independent Einstein-tensor computation for ds² = f dt² − dr²/f − Σ² dΩ² reproduces the
  paper's G⁰₀ and G²₂ exactly.
- **G⁰₀ − G²₂ = 3Ma²/Σ⁵ > 0 everywhere.** Hence L_F > 0, the √ terms are real, and the EM channel is not a ghost.
  Our L_F ∝ 1/Σ agrees with Bronnikov–Walia's magnetic L_F = (3Ma²/2q²)/Σ up to a constant normalisation (F
  convention).
- Closed forms: V_AF (Eq. 55), V₁₁ and V₁₂ (Eq. 56) equal their definitions (V_AF symbolically; V₁₁ and V₁₂ to
  ≤ 2e-16 at 8 test points). V₁₂ = √(3/2π) μ a √M (Σ − 2M)/Σ^{9/2} is independent of the charge.
- **Schwarzschild limit:** V_AF(a → 0) = f[(μ² + 2)/r² − 6M/r³], the Regge–Wheeler potential. V₁₂ → 0
  (decoupling). V₁₁ does not tend to Maxwell's f ℓ(ℓ+1)/r², because the L_F-gradient term is scale-free in a. That
  is expected for an NED channel, not an error.
- **Dimensions:** every V scales as M⁻², including V₁₂ ∝ M^{1/2} a · M / M^{9/2}.
- **Throat:** all potentials are smooth, even functions of r at r = 0, so the throat is regular.
- **Boundaries:** all V → 0 as r → ±∞, and V₁₂ falls as |r|^{-3.5}. Both ends are radiative, matching the
  paper's outgoing / "ingoing" (outgoing into the other universe) conditions.
- **Barrier shapes (ℓ = 2):**
  - V_AF peaks at R = 3.196M with V = 0.1439, about 5 % below the X2 "RW-form" proxy (0.1513).
  - V_AF has a double barrier for 2 < a < 3, the same range as the test fields.
  - The larger local eigenvalue of the NED matrix peaks at R ≈ 2.77M.

**Doubts found (recorded, not resolved by the paper):**
- **D1.** "V₁₂ ≠ V₂₁ makes the system non-Hermitian" (§5.4). Because V₁₂V₂₁ > 0, replacing H₁ by √(4π)·H₁
  makes the potential real and symmetric (off-diagonal √(4π)·V₁₂). The asymmetry is a normalisation artifact.
  Openness comes only from the boundaries. The equations themselves are unaffected. We evolve the symmetric form;
  its conserved flux fixes the relative channel normalisation.
- **D2.** The eikonal formula r_c = √(9M² − a²) is stated "for a ≤ 2M"; it holds for a < 3M. Cosmetic.
- **D3.** The transition-region table (2.01–2.12 M) jumps between modes. For example the NED GR-led real part
  goes 0.336 → 0.223 between a = 2.14 and 2.15. Only the smooth sequences from about 2.15 M on are usable as
  reproduction targets.
- **D4.** The initial data are under-specified (centre, width, observer position; time derivatives
  "self-consistent"). A reproduction can match frequencies, not amplitudes.
- **D5.** The AF closure is an assumption. The NED magnetic model is the only one with an action and a closed,
  dynamical matter perturbation in this sector.

**Verdict: what it supplies.**
- It replaces P1 **for the axial sector only**, with a gravitational master equation derived from the field
  equations with explicit, perturbed matter.
  - In the NED magnetic + phantom-scalar model, gravitational waves couple to the EM field. Energy can leave the
    GW channel, which no single-channel proxy can do.
  - The AF equation gives the pure-gravity barrier V_AF.
- This is more physical information than the test-field proxies. Stage 2 proceeds, with the NED magnetic model.

**What is still missing:**
- **The polar sector.** Isospectrality is broken, and a merger excites both parities. The phantom scalar couples
  to the polar sector, where instabilities of phantom-supported throats are usually found; nothing here tests it.
- **The source excitation η.** A Gaussian or ringdown packet amplitude is not η.
- **Spin.**
- **Stability beyond linear axial modes.**
- **Any formation or merger context.**

## Stage 2 design (frozen before computing)
**One model, one parameter set.** NED magnetic + phantom scalar, axial ℓ = 2, symmetrised potential matrix
[[V₁₁, √(4π)V₁₂], [√(4π)V₁₂, V₂₂]]. Tool: `scripts/x2_axial_ned.py`, a two-channel time-domain leapfrog that
reuses `x2_scatter.Grid` for r(r*). A coupling-off run (V₁₂ = 0, which is AF / pure EM) is a control, not a
second model.

**Reproduction (published result):**
- Target: Yang & Lan supplementary Table "WH_GR_1", a = 2.20M, NED (GR-led) ω = 0.246140 − 0.000111i.
- Setup: Gaussian in H₁ (σ = 2, centred 10 M outside the right peak), H₂ = 0, zero initial time derivative.
  Observer at r* = +50. Modes are fitted by matrix pencil on H₂, t ∈ [350, 6000].
- Pass:
  - an extracted mode with |ΔRe ω|/Re ω ≤ 0.5 % and 0.67 ≤ |Im ω|/0.000111 ≤ 1.5;
  - dx = 0.1 and 0.05 agree to 0.1 % in Re ω and 10 % in Im ω.
- Coupling-off control: the AF mode at a = 2.20 should appear near the published 0.308238 − 0.000266i (same
  tolerances). The control is informative only.
- Fail → the tool is not used for the comparison; record why and stop stage 2.

**Comparison question (frozen):**
- Setup: at X2's applicable GW150914 cell a − 2M = 3.58e-4, an inward ringdown-shaped packet in the **GW
  channel** (same w₀ = 0.5307, τ = 12.2 M as S3) starts 25 M inside the right GW peak. The observer is at
  r* = 1500 (GW channel, and the EM channel for leakage).
- Measured in the time domain, so the full transmission and reflection phases are included and nothing assumes 4L:
  - **Q1 timing:** the first-echo arrival relative to the direct signal, and the pulse-1 → pulse-2 spacing
    (energy centroids). Each is compared with that model's own 4L, and with the same quantities for the X2 proxies
    s = 1 and "rw" computed by the same code.
  - **Q2 strength:** pulse-1 and pulse-2 h_rss in the GW channel relative to the packet, and the fraction of
    outgoing energy in the EM channel. Then η_max,NED = 0.21 / √(Σ_k p_k²) over pulses inside the PMW
    [0.215, 1.185] s, with times converted using M_det = 68.8 M☉. This uses the ringdown-packet weight, which S3
    showed agrees with the frequency-domain Lorentzian weight to 3–5 % for s = 1.
  - **Q3 decay:** the ratio p₂/p₁.
- **Numerics:** dx = 0.2 and 0.1 agree to 1 % in each p_k and 2 M in each arrival time.
- **Reading:**
  - If η_max,NED lies inside the proxy bracket 0.70–1.24, the proxies bracketed the matter-consistent axial
    model, and P1 for the axial sector narrows to that single value.
  - If it lies outside, the proxy bracket was inadequate for this sector.
  - Either way the result remains conditional on P2–P4, η and the missing polar sector.

## Stage 2 results (derived; `scripts/x2_axial_ned.py`, `tests/test_x2_axial_ned.py`, `results/x2/s4_axial_ned.json`)
**Reproduction: passes.** Yang & Lan, a = 2.20M, ℓ = 2.
- NED GR-led mode: 0.246142 − 0.000116i against the published 0.246140 − 0.000111i (Re 0.001 %; Im ratio 1.04).
- Converged: dx = 0.1 and 0.05 agree to 3e-6 in Re ω and 4e-5 in Im ω.
- Coupling-off (AF) control: 0.308249 − 0.000263i against 0.308238 − 0.000266i.
- The matrix pencil also finds faster-decaying modes with larger amplitude in [350, 6000] M, for example
  0.3117 − 0.0034i. The published mode is the long-lived one, as the paper's late-time Prony fit implies.
- A first run used a 60 M peak-search grid. At a − 2M = 3.6e-4 the peak lies about 330 M from the throat, so
  the echo windows were wrong. This was a code bug, fixed (guard and test) before the comparison below. The
  reproduction was unaffected: the peak at a = 2.20 is 13.7 M out.

**Comparison at GW150914's S2 cell** (a − 2M = 3.575e-4; dx = 0.1, converged against 0.2 within the frozen
1 % / 2 M criteria):

| Model | Own 4L (M) | First echo − free flight (M) | Spacing (M) | Spacing vs 4L | p₁ | p₂ | p₂/p₁ | GW-echo energy in EM channel | η_max (PMW) |
|---|---|---|---|---|---|---|---|---|---|
| NED magnetic + phantom (axial GW) | 1330.8 | 1312.8 | 1373.5 | +3.2 % | 0.122 | 0.063 | 0.52 | 67 % / 70 % (pulse 1 / 2) | **1.52** |
| s = 1 proxy (EM) | 1327.9 | 1305.9 | 1375.9 | +3.6 % | 0.239 | 0.116 | 0.48 | — | 0.79 |
| "rw" proxy | 1332.0 | 1308.9 | 1373.4 | +3.1 % | 0.168 | 0.082 | 0.49 | — | 1.13 |

- **Q1, timing.** It is set by the geometry, not by the matter model.
  - The three models tested (not AF; one a) agree on the spacing to 2.5 M (0.2 %) and on the first echo to 7 M.
  - Arrival is measured against free flight from the packet's start x₀ = x_peak − 25 M. That is a deviation
    from the frozen "relative to the direct signal" (post hoc label): an inward packet has no direct signal.
    - The first echo's geometric extra path is 2(x₀ + x_peak) = 4L − 50 M.
    - The measured 4L − 18 to 23 M therefore means the barriers add **+27 to +32 M** of delay, not an early
      arrival.
    - The spacing is 4L + 41 to 48 M.
  - In seconds: pulses at 0.445 s and 0.910 s are inside the PMW, and pulse 3 at 1.376 s is outside. The S2
    applicability of this cell stands.
  - This is the full-phase measurement the owner asked for. 4L understates the spacing by 3–4 %, consistent with
    S3's +4.1 % at a different a. S3's failed 1 % spacing gate remains a failure: the cause is the 4L reference,
    and it is quantified here.
- **Q2, strength.**
  - In the NED model about two-thirds of the outgoing echo energy leaves in the EM channel. The GW pulses are
    0.51× / 0.55× the s = 1 proxy and 0.73× / 0.78× the "rw" proxy (pulse 1 / pulse 2).
  - η_max,NED = 1.52 lies **outside the proxy bracket 0.70–1.24**.
  - By the frozen reading, the test-field bracket was inadequate for the axial sector: single-channel proxies
    cannot lose GW energy to the matter field. With the matter-consistent axial equations, GW150914 gives no
    constraint at η ~ 1 for this cell.
  - The EM fraction moves by 1.6 % (pulse 1) and 4 % (pulse 2) between dx = 0.2 and 0.1. It was not in the
    frozen convergence list.
- **Q3, decay.** p₂/p₁ is 0.52 (NED) against 0.48–0.49 (proxies). The per-round-trip loss is similar; the
  difference is mainly in the first transmission.
- **Weight check.** The same code gives s = 1: 0.79 and "rw": 1.13, against S2's frequency-domain 0.80 and 1.24.
  Ringdown-packet weighting agrees with the Lorentzian weighting to 1.4 % (s = 1) and 9 % ("rw").

**What this narrows, and what it does not.**
- It **narrows** P1 for the axial sector only, in one matter model and at one a.
  - The question was which barrier GWs obey. In the axial sector it is answered with a closed matter
    perturbation, and the answer falls outside the earlier bracket because of GW→EM conversion.
- It does **not** address:
  - the polar sector (non-isospectral; the phantom scalar enters there);
  - η (the packet amplitude is not η);
  - spin, stability beyond linear axial modes, or the electric/fluid alternatives (AF is a different
    closure; not compared);
  - other a values or other events.
- Not a statement about SV existence.

**Next executable step (if the owner resumes X2).**
- Same tool, AF (V₂₂ alone) at the same cell (≤ 0.3 h). This separates barrier-shape from EM-leakage effects.
- The polar-sector equations for the NED + phantom source: a literature check first. None is in our record, and
  Yang & Lan list polar modes as future work.
- Parameters for a rerun: a − 2M = 3.575e-4 (Δt = 0.45 s at M_det = 68.8 M☉); packet w₀ = 0.5307, τ = 12.2 M,
  started 25 M inside the GW peak; observer at r* = 1500; dx 0.2 / 0.1.

## S4b: physical strain and energy flux of the S4 readouts (owner request 2026-10-13; written before computing)
**Derivation (Yang & Lan definitions, checked against Nagar & Rezzolla gr-qc/0502064v5 §6 in the Schwarzschild limit):**
1. **H₂ is the Moncrief / Regge–Wheeler function, not the strain.**
   - Yang & Lan Eq. 34: H₂ = sin³θ Σ f (q_{1,θ} − q_{2,r}) / C^{-3/2}_{ℓ+2}.
   - In Regge–Wheeler gauge (q₂ = 0, so g_{rφ} = h₁ sinθ ∂_θY), q_{1,θ} = (h₁/Σ²) ∂_θ(∂_θY/sinθ).
   - sin³θ ∂_θ(∂_θY/sinθ) = c_ℓ C^{-3/2}_{ℓ+2} (checked: c = 8, 40, 120 for ℓ = 2, 3, 4).
   - Hence H₂ = c_ℓ f h₁/Σ. As a → 0 this is Nagar–Rezzolla's Q^(o) = (f/r) h₁ (their Eq. "mfop").
2. **Strain needs one time integral.**
   - Nagar–Rezzolla: h₊,× = −Ψ^(o)/r × angular factor (the CPM function), and ∂ₜΨ^(o) = −Q^(o) + (odd source)
     (Eq. "Q_vs_Psi").
   - In NED the corresponding equation is Yang & Lan's G13-SM. Using H₁'s definition and
     dC^{-3/2}_{ℓ+2}/dθ = 3 sinθ C^{-1/2}_{ℓ+1}, it becomes
     ∂ₜ(Σ⁴ Q₀₁ sin³θ) = −3 sinθ C^{-1/2}_{ℓ+1} [Σ H₂ + (κ² q √L_F / μ) H₁].
   - With L_F = 3Ma²/(8π q² Σ), the EM term relative to the GW term is 2.45 Σ^{-3/2} in the symmetrised H₁ (a ≈ 2M,
     μ = 2). That is 4e-5 at the observer (r* = 1500) and → 0 at infinity.
   - So far from the source, h ∝ (1/r) ∫ H₂ dt, i.e. h̃ = H̃₂ / (−iω), as in Schwarzschild. Near the throat
     (Σ ≈ 3M) the EM term is about 0.4 of the GW term, so H₂ is not the strain there.
3. **Energy flux.**
   - GW: Nagar–Rezzolla Eq. "power": dE/dt ∝ |Ψ̇^(o)|² = |Q^(o)|², i.e. ∝ H₂², not (∂ₜH₂)².
   - EM: the physical (Killing) energy and the symmetrised form's flux Σᵢ|∂ₜψᵢ|² are both conserved by the same
     scattering matrix. At fixed ω they can both be conserved while the channels mix only if their channel weights
     are equal.
   - Hence the physical flux is ∝ H₂² + H₁,sym² with the same constant. The symmetrised conserved quantity is the
     physical energy up to a factor ω² per frequency, i.e. its time-domain flux is ∫ψ̇² rather than ∫ψ².

**Consequences for the S4 readouts:**
- `hrss_ratio` = √∫H₂²dt / √∫p²du. Propagation is linear and time-invariant, and integration commutes with it. So the
  strain response to a source whose strain is p equals the H₂ response to an H₂ input p. The ratio is therefore the
  strain h_rss ratio for a ringdown-shaped strain source, as P2 intends.
  - No integration constant enters: the output is zero before arrival.
  - At low frequency the source strain p carries a nonzero ∫p du, but h_rss = √∫p² is finite. The physical H₂ input
    ṗ has zero mean (no memory).
  - **Check (written before running):** a smooth packet g and its derivative ġ, same cell, NED, dx = 0.1. The time
    integral of the ġ output (from t = 0, no detrending) must match the g output to ≤ 1 % in each pulse h_rss and
    ≤ 2 % L2 over the pulse windows. The integrated waveform's end value must be ≤ 1 % of its peak.
- `em_energy_fraction` used ∫(∂ₜψ)² dt. The physical flux uses ∫ψ² dt. It is recomputed from saved waveforms,
  and two quantities are reported separately:
  - (a) the receiver-side fraction: EM energy over GW + EM energy in each pulse window at the observer;
  - (b) the budget relative to the initial packet energy E₀ = ∫p² du, at r* = +1500 and −1500 up to t_end
    (5 round trips). The residual is energy still trapped or not yet arrived.
- **Diagnostic (≤ 0.3 h, not a self-consistent model):** the same cell, initial data and readout with the coupling
  removed and V₂₂ unchanged. This separates barrier shape from channel coupling.
- **Reruns:** only this cell, dx = 0.1 (convergence at this cell was shown in S4). The verified modes and solver are
  unchanged.

### S4b results (derived; `scripts/x2_s4b_units.py`, `results/x2/s4b_units.json`; waveforms kept in the scratchpad only)
**Correction to the pre-written item "consequences"** (post hoc; made after the budgets were seen; it keeps S4's
67–70 % share).
- The physical flux is ∝ ψ² **of the physical H₂ and H₁**. The question is which waveform is physical.
- P2 takes the source *strain* to be the ringdown packet p. The physical H₂ input is then ∓ṗ:
  - with u = x₀ − x − t, d/du = −d/dt;
  - by linearity, the physical outputs are ∓ d/dt of the run outputs.
- **Consequences:**
  - The strain is the run output itself, so `hrss_ratio` stands.
  - The physical energy flux is ∝ (d/dt run output)², so S4's ∫ψ̇² readout was the right one for this source.
  - The ∫ψ² readout instead belongs to a source whose *H₂* is ringdown-shaped. That source's strain ∫p dt
    carries a DC offset (memory). It is reported for comparison only.
- **The choice rests on P2 alone.** The budgets cannot decide between the readings: the linear evolution conserves
  both measures, since ∫ψ dt is itself a solution.
  - Budget closures, for the record: ψ̇² 96.4 % (NED) and 98.9 % (coupling off); ψ² 81.9 % and 84.7 %; smooth
    packet g 96.5 % / 93.7 %.
  - The ψ² shortfall reflects the free-wave E₀ = ∫p² du. It omits the ½∫Vψ² term (the packet starts 25 M inside the
    peak) and p's DC and low-frequency content, which leaves the cavity late.
  - An earlier version of this record said the budgets decided between the readings; that was withdrawn after
    review.

**Strain check (h = ∫H₂ dt, from t = 0, no detrending):**
- Pulse h_rss agree to 0.5 % and 0.7 %, and the end value is 1e-10 of the peak. Both pass as written.
- The L2 criterion **fails as written** (1.99): the integrated waveform is the negative of the reference, from the
  d/du = −d/dt sign above.
- With that sign (post hoc) and a trapezoid integral, the L2 difference is 0.75 % and 1.08 %, which passes. A
  rectangle rule gave 1.3 % / 1.4 % (half-step bias), replaced after review.
- **What this shows:** the solver is linear and time-invariant, the integral has no constant, and there is no
  low-frequency drift. The *physical* relation h ∝ ∫H₂ dt is not tested by it; that rests on the derivation above.

**Numbers** (one cell, a − 2M = 3.575e-4, dx = 0.1, t_end = 8820 M):

| Quantity | S4 | S4b (strain-shaped source, used) | H₂-shaped source (comparison) |
|---|---|---|---|
| p₁, p₂ (strain h_rss ratio) | 0.122, 0.063 | 0.122, 0.063 (unchanged) | — |
| η_max (PMW, pulses 1–2) | 1.52 | **1.52 (unchanged)** | — |
| EM share of energy in the receiver-side pulse windows (pulse 1 / 2) | "67 % / 70 %" | **67 % / 70 %** (justification corrected) | 53 % / 56 % |
| EM share of the **initial** packet energy, both sides, by t_end | not computed | **30.8 %** (our side 2.5 %, far side 28.3 %) | 22.0 % |
| GW energy out: our side / far side | not computed | 1.2 % / 64.4 % | 2.4 % / 57.5 % |
| Residual (trapped or not yet arrived) | — | 3.6 % | 18.1 % |

- The receiver-side fraction (67–70 %) is a ratio within the echo windows at our observer. The initial-energy
  fraction (31 %) is a budget over both universes. They answer different questions.
- Most of the inward packet (w₀ = 0.53, above the V₂₂ barrier top at ω = 0.38) passes through the far barrier.

**Diagnostic: coupling removed, V₂₂ unchanged** (same cell, initial data and readout; not a self-consistent model):
- p₁ = 0.152, p₂ = 0.075, η_max = 1.24. The budget closes to 98.9 % (strain reading).
- D's own arrivals are 1314.7 / 2692.9 M (spacing 1378.2 M), against NED's 1312.8 / 2686.3 M. Both pulses are in
  the PMW. An earlier draft reused NED's timing for D; η_max is unchanged.
- Time-domain η_max at this cell, same method throughout:

  | Potential | η_max |
  |---|---|
  | s = 1 proxy | 0.79 |
  | "rw" proxy | 1.13 |
  | V₂₂ alone | 1.24 |
  | NED (coupled) | 1.52 |

- **Barrier shape** (V₂₂ against the "rw" proxy) raises η_max by about 10 %. **Channel coupling** raises it by a
  further 23 %.
- The earlier frequency-domain proxy bracket was 0.70–1.24. Its upper end ("rw", Lorentzian weight) happens to
  equal the V₂₂-alone value here. The coupled NED value stays outside it.

**Remaining caveats:**
- The source is defined in free flight; near the throat (Σ ≈ 3M) H₂ is not the strain (the EM term is about 0.4).
  This is the inherited P2/P4 approximation.
- The equal-weight argument for the physical EM flux assumes the canonical (Killing) energy of the axial NED sector
  is conserved and positive. That holds because L_F > 0 and the scalar decouples.

## S4c: polar / radial stability literature (≤ 1 h; IDs and titles checked on arxiv.org)
- **No paper derives polar (ℓ ≥ 2) or radial (ℓ = 0) perturbation equations for SV**, or for any black bounce, with
  the phantom + magnetic NED source.
  - Black-bounce perturbation work found is axial or test-field: 2603.20594, and 2608.20081 (Alencar et al., de
    Sitter black bounces, axial).
- **Closest radial equation:** Bronnikov et al. 2026, arXiv:2607.10718. It treats NED in scalar-tensor gravity
  with a general L(ψ, F), magnetic charge and a possibly phantom scalar, and gives ℓ = 0 potentials.
  - Its gauge δβ = 0 "causes certain problems when applied to configurations with extremum values of r (such as
    wormhole throats) … irrelevant in the present study" (checked in the source).
  - So it does not cover SV with a > 2M as written.
- **Nearby evidence points to instability outside the axial sector** (suggestive, not shown for this source):
  - phantom-scalar wormholes are unstable at ℓ = 0 (Bronnikov, Konoplya & Zhidenko 1205.2224; Gonzalez, Guzmán &
    Sarbach 0806.0608; Cremona et al. 1911.13103);
  - charged Ellis–Bronnikov wormholes are unstable, radially (2510.11406), in the axial sector (2507.02258) and in
    the polar sector at ℓ = 2 (2607.18399).
- **Executable gap.** The smallest well-posed next step is the radial ℓ = 0 equation for SV with phantom +
  magnetic NED. Take 2607.10718's general system in its GR limit, re-derived in a throat-compatible gauge, as
  1205.2224 and 1911.13103 do for phantom-only wormholes. Then solve for unstable modes on one a.
  - It answers a concrete unknown: is the background used in S4 linearly unstable at ℓ = 0? If it is, echo
    bounds on it are moot.
  - The polar ℓ = 2 system (GW + EM + scalar coupled) is a larger derivation with no published template for
    black bounces.
