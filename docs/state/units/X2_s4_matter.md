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
- **D1.** "V₁₂ ≠ V₂₁ makes the system non-Hermitian" (§5.4). Because V₁₂V₂₁ > 0, rescaling H₁ → H₁/√(4π)
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
  - All three models agree on the spacing to 2.5 M (0.2 %) and on the first echo to 7 M.
  - In all of them the first echo arrives about 4L − 18 to 23 M after free flight, and the spacing is
    4L + 41 to 48 M.
  - In seconds: pulses at 0.445 s and 0.910 s are inside the PMW, and pulse 3 at 1.376 s is outside. The S2
    applicability of this cell stands.
  - This is the full-phase measurement the owner asked for. 4L understates the spacing by 3–4 %, consistent with
    S3's +4.1 % at a different a. S3's failed 1 % spacing gate remains a failure: the cause is the 4L reference,
    and it is quantified here.
- **Q2, strength.**
  - In the NED model about two-thirds of the outgoing echo energy leaves in the EM channel. The GW pulses are
    0.51× the s = 1 proxy and 0.73× the "rw" proxy.
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
