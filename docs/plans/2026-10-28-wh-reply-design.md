# From pre-loaded two-way propagation to a reply chosen after receipt: model and feasibility (2026-10-28)

Owner-requested bounded design round (≤ 60 min). **No evolution, scan or new computation was run.** This note
does not change any conclusion of #592 (`docs/results/wh_twoway.md`); it asks whether that model can be
extended, with a physically consistent minimal change, so that B chooses reply 0 or 1 only **after** receiving the
outbound pulse and A can tell which. It is a model-internal causal-reply question: no receiver noise or coding
model (so no bit rate), and the model joins two separate universes (so no faster-than-light statement).

## 1. What is reused and what was checked

- **From #592** (unchanged): Kain arXiv:2602.18231 double-null equations and scheme; the outbound χ₁ pulse
  (E_out = 10⁻⁴ x₀) and stations at r = 2 x₀; stored facts used below: χ₁ crosses B's sphere outward with Kodama
  energy 8.27 × 10⁻⁵ at τ_B = 4.468 (10 % at 4.280, 90 % by 4.711; `scan_n800.json`, E_ret = 10⁻⁸); over B's whole
  track the Misner–Sharp mass at B ranges from 0.196 to 0.2500022 (static value 0.25), the low end coming from the
  later collapse, so its drift *before* the outbound pulse arrives is not stored; certified arrivals for returns
  crossing B by τ_B ≈ 6.54 with E_ret ≤ 0.1 (§7 of the record).
- **In #592 nothing at B is physical**: B is a massless observer; the return pulse is pre-loaded initial data.
- **Published methods checked this round** (two targeted searches; not a general survey):
  - Thin-shell junctions: Israel (timelike shells) and Barrabès–Israel (null shells), the latter used for null-dust
    shells (e.g. arXiv:1703.07414); a Hamiltonian thin null-dust shell (gr-qc/9708012); a double-null scalar-collapse
    code with a matter shell (arXiv:1508.05578). None evolves a *timelike emitting* shell inside a double-null
    scalar code.
  - Self-consistent matter–scalar energy exchange with monopole scalar radiation: Novak, "Spherical neutron star
    collapse in tensor-scalar theory of gravity", gr-qc/9707041 (PRD 57, 4789, 1998): coupling α(φ) times the trace
    of the matter stress tensor, perfect fluid, Cauchy code (spectral scalar + shock-capturing hydrodynamics).
    Supported, but it needs a fluid solver that the #592 code does not have.
  - Not found in this round's search: any wormhole evolution with a triggered emitter or receiver.

## 2. The minimal problem

**Device D_B** (new): a smooth spherical matter layer around B's station, areal radius r ∈ [2 − w/2, 2 + w/2] on
the − side, width w = 0.1 x₀ (≥ 40 cells at n = 400), rest-energy density ρ(τ, r) > 0, held **rigidly** at constant
r (its worldlines are the static, Kodama-direction ones), with no radial stress and whatever hoop stress p_t the
radial force balance requires. Stress tensor T_D^{ab} = ρ u^a u^b + p_t q^{ab} (q = projector on the spheres). In
the double-null variables of #592: T_D,uu = ρ u_u², T_D,vv = ρ u_v², T_D,uv = ρ e^{2σ}/4, T_D,θθ = p_t r², so it enters
Kain's r_uv and σ_uv equations and both constraints.

1. **B's local measurable for "outbound received"**: the Misner–Sharp mass on the device's inner face,
   m_in(τ) = (r/2)(1 + 4 e^{−2σ} r_u r_v), obtainable locally from the expansions of inward and outward radial light
   rays at that sphere. The outbound pulse lowers it by ≈ 8.3 × 10⁻⁵ as it crosses outward. **Trigger** (in model):
   the first τ at which m_in(τ) − m_in(τ_ref) < −4.1 × 10⁻⁵ (half the expected drop; τ_ref = the face's first stored
   time). Whether the ghost-field drift of m_in stays well below this threshold until the outbound pulse arrives
   is not stored in #592; it is pass item 1 below.
2. **How B chooses after receipt**: an **external local control** b ∈ {+1, −1}, applied at B and entering the
   equations only after the trigger. The two runs of a pair have identical initial data (wormhole, outbound pulse,
   device mass) and are identical before τ_trig; they differ only in b. The model does not produce the choice; no
   claim is made that one closed-system initial state yields two outcomes.
3. **Energy source and exchange**: the device's own rest energy. After τ_trig + Δ (fixed latency Δ), the device
   sources the return field: □χ₂ = −Q with Q = b · Q₀ · S(r) · W(τ − τ_trig − Δ) (S the device's radial profile, W a
   fixed unit-width window). To keep the total stress tensor conserved, the device obeys
   ∇_a T_D^{ab} = +Q ∇^b χ₂, i.e.
   - energy: dρ/dτ + ρ θ_u = −Q dχ₂/dτ (θ_u the expansion of the device worldlines): the device pays for the work
     done on χ₂;
   - radial momentum: p_t = r [ρ (n·a) − Q n·∇χ₂] / (2 n·∇r): the recoil and the field force go into the hoop
     stress that holds the device at fixed r.
   In χ₂'s equation in #592's variables this is one extra term, χ₂,uv += (e^{2σ}/4) Q.
4. **A's observable, fixed in advance**: O_A = ∫ χ₂(A, τ) dτ over τ_A ∈ [τ_trig + 3, 12] on A's worldline
   (r = 2, + side); decision b̂ = sign(O_A). A is a test observer (no detector backreaction, no noise).

**Why the sign is the encoding.** The equations are invariant under (χ₂, Q) → (−χ₂, −Q): T_χ₂ and the device's
energy and momentum equations are even in this pair. So both replies have **the same** geometry, device energy loss
and received energy; only the sign of χ₂ differs. A cannot tell the bit from gravity or energy, only from the field.

## 3. Conservation, constraints, back-reaction, and what is idealised

- **Energy–momentum**: ∇_a(T_ghost + T_χ₁ + T_χ₂ + T_D)^{ab} = 0 exactly by construction (the device's equations are
  defined as the balance of χ₂'s source term).
- **Constraint propagation**: with the total stress tensor conserved, the Bianchi identity propagates Kain's null
  constraints; the double-null scheme still evolves only r_uv, σ_uv and the fields, so the null-constraint residual
  and the Misner–Sharp bookkeeping (now including T_D) remain the checks.
- **Device back-reaction**: the device's mass, hoop stress, energy loss and recoil all enter the Einstein equations;
  nothing is a test object except A's detector.
- **Idealisations and their validity conditions** (checked a posteriori, each a stop condition if violated):
  1. Rigid support (r = const, no radial stress, hoop stress as needed): admissible only while ρ > 0 and the
     dominant energy condition |p_t| ≤ ρ holds throughout the device; rigidity is non-causal across the layer, so
     w = 0.1 x₀ must stay ≪ the pulse duration (1 x₀/c).
  2. Phenomenological source Q (not derived from a material Lagrangian or an equation of state): its energy is
     balanced exactly, but which real material could do this is not addressed.
  3. External control b: outside the model, applied at B after τ_trig only.
  4. Spherical device (a whole sphere), s-waves only, as in #592.
  5. A static device must sit on the unstable wormhole in equilibrium: new static initial data (Ellis + ghost +
     layer with its hoop stress) are needed; the interior of a static spherical layer is unchanged up to a lapse
     constant, so the throat region is the same Ellis solution before the outbound pulse arrives.

## 4. Verdict: **B** — only a clearly labelled, idealised approximation is possible

- **Not A.** No published method plugs into the #592 code: the supported field–matter coupling with self-consistent
  monopole scalar radiation (Novak 1998) needs a self-gravitating fluid and a shock-capturing hydro solver, and the
  thin-shell junction literature has no evolved timelike emitter inside a double-null scalar code.
- **B.** The rigid-layer device of §2–3 can be added to the #592 code as a finite extension: one layer stress
  tensor in the r_uv and σ_uv updates and constraints, one source term in χ₂'s equation, the device energy equation
  along static worldlines, an algebraic hoop stress, a trigger, and static initial data with the layer. Exact
  energy–momentum balance and constraint propagation are kept; the idealisations in §3 are stated and testable.
  Estimated cost: one implementation-and-validation round (≈ 3–4 h), then the pair (minutes).
- **Not C for this question**, but C for a *physically complete* device: a material receiver/emitter with an
  equation of state, a Lagrangian coupling and an internal switch (e.g. a Novak-type fluid shell, or a boson shell
  around the throat coupled to χ₂) needs a new matter model and a new solver, and is not worth building now.

## 5. The single minimal paired test (not run)

- **Pair**: b = +1 and b = −1, same initial data; n = 400 and 800 for each (four runs; the 800 runs are the result,
  400 measures convergence). Device rest energy E_D = 0.1 x₀ (≥ 5 × the energy it must radiate); Q₀ set before the
  runs so that the device radiates ≈ 10⁻² x₀ in total (about half of it inward); Δ set so that the inward
  emission's 50 % point at B falls at τ_B ≈ 5.5, inside #592's certified window.
- **Pass** (all must hold, at n = 800, with the n = 400 → 800 change as the error):
  1. the trigger fires during the outbound pulse's 10–90 % crossing of the inner face, at the same τ in both runs,
     and the two runs agree to round-off before τ_trig;
  2. the device's energy loss equals the χ₂ energy that leaves it plus the χ₂ energy still near it, to ≤ 1 %, and
     |p_t| ≤ ρ, ρ > 0 everywhere in the device until the emission ends;
  3. the null-constraint residual in A's exact causal past (r > 0.5, not the over-inclusive rectangle of #592)
     converges at about second order;
  4. O_A(+1) = −O_A(−1) to ≤ 10⁻¹⁰ relative (the symmetry), |O_A| ≥ 10 × its resolution change, and sign(O_A) = b
     in both runs; at least 90 % of the χ₂ energy A receives by τ_A = 12 has arrived by τ_A = 11.
- **Stop** (end the line, report the reason, do not switch method): no static device-on-wormhole initial data with
  ρ > 0 and |p_t| ≤ ρ; static preservation with the device not second order; the trigger fires before the outbound
  pulse (on drift); constraint residuals near the layer do not converge; or any pass item fails at n = 800.
- **What it would add over #592**: a reply whose start is set by B's own local detection of the outbound pulse, whose
  energy is debited from a device that stays in the geometry (with its recoil), and whose content is chosen at B
  after receipt, recovered by A from a pre-fixed observable. **What it cannot add**: the timing window and arrival
  are already fixed by #592 (the inward emission is a return of ≈ 0.5 × 10⁻² x₀ inside the certified window), and the
  identification is guaranteed by the exact sign symmetry once anything arrives. The expected outcome is therefore
  largely set in advance; only the device-model items (static admissibility, energy-condition margin, trigger
  robustness) are genuinely open.

## 6. Plain answers

1. **Beyond #592**: whether the return can be a *reply*: started by B's own measurement, paid for by B's own
   energy, and chosen after receipt, with A able to read the choice.
2. **Still put in by hand**: the ghost field (as before); the device's rigidity and its coupling Q (no real material
   is specified); the choice b (an external control); A's noiseless test detector.
3. **What the minimal test could support**: that inside this idealised model a causally triggered, energy-paid,
   after-receipt choice reaches A and is identified — a consistency demonstration of the modelling, not new physics
   about wormhole communication, and nothing about rates or shortcuts.
4. **Worth a computation round?** Not now. Its likely result is predictable from #592 plus the sign symmetry, and it
   leaves every limit that matters (assumed ghost matter, two universes, no physical device, no noise model) in
   place. It becomes worth one round only if a physically specified device model (an equation of state and a
   Lagrangian coupling) is wanted, or if the owner wants the idealised demonstration on record.
