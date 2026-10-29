# Is the FGM wormhole's opening window compatible with a finite quantum wave packet? (2026-10-29)

Owner-requested bounded analytic round (≤ 90 min). Only analytic derivation, formula substitution and small
integrals; no evolution solver and no scan. This is not D-112's literature screen and does not touch #592. One
construction only: Fu, Grado-White & Marolf, arXiv:1908.03273v2 (FGM; CQG 36, 245018), read from its arXiv TeX
source. Labels: [lit] FGM's text or equations; [der] our arithmetic from them; [assumption]; [inference].

## 1. Has it been answered? Not found in this round's search

- Repository: D-112 (`docs/plans/2026-10-29-wh-support-screen.md`, on #592's branch) quotes FGM's fragility and
  footnote 6 but computes nothing for a wave packet; goal round 5 records only FGM's transit d + logs.
- FGM themselves treat signals as null geodesics ("test signals", §1) and leave the signal's back-reaction undone
  (fn 6).
- All 43 INSPIRE citers of 1908.03273 (titles; abstracts of the six closest: 2606.12528, 2010.15023, 2304.07356,
  2012.07821, 2103.01920, 1908.03998) contain no wave-packet traversal or self-back-reaction calculation for FGM.
  This describes the search, not the literature as a whole.

## 2. Frozen plan (written and committed before any computation)

**Observers and normalisation (all from FGM §3.1).** Near each mouth the throat is Reissner–Nordström (RN) with
r₊, r₋, κ₊. Killing time t of that RN patch (normalised at large r in the patch; for d ≫ r₊ this is the asymptotic
time up to O(r₊/d)); tortoise r* with r* = 0 at r = 0 (FGM (3.2)); u = t − r*, v = t + r*; Kruskal
U = ∓κ₊⁻¹e^{∓κ₊u}, V = ±κ₊⁻¹e^{±κ₊v} (FGM (3.3)), the normalisation in which FGM's (2.20) holds with
x*₀ = d/2 + O(ln d/r₊) (FGM (2.21)). Energies are Killing energies (conjugate to t), i.e. measured by static
observers at large r of the patch; frequencies ω are with respect to t. All redshifts are taken from these
definitions; no other normalisation is used.

**Direction and window.** A forward signal enters mouth A as an ingoing ray of advanced time v and moves along
increasing U near V = 0. It reaches the other mouth B iff V(v) < |ΔV| (FGM (3.8)–(3.9)), i.e. v < v_max(θ) =
κ₊⁻¹ ln(κ₊|ΔV(θ)|), with ΔV = C_RN · 8πG · (H(θ)/r₊²) · ∫T_UU dU, C_RN from FGM (3.9), H from (3.14) and ∫T_UU dU
from (2.20). The factor 1/r₊² comes from the string's transverse delta function; without it the argument of the
logarithm in FGM's (3.17) has dimension length² [inference, checked by dimensional analysis]. The emerging label is
v_out = −κ₊⁻¹ ln(κ₊(|ΔV| − V)). What opens the window is the counter-moving component T_UU along V = 0 (FGM (3.8),
fn 6).

**Parameter point (illustrative, fixed).** r₊ = 1 km; r₋ = 0.5 r₊ (κ₊r₊ = 0.25); d = 10⁴ r₊; N = 10 strings (c = 20);
Gμ = 10⁻⁷ (so μr₊² ≫ 1 and NGμ = 10⁻⁶ ≪ 1). Evaluated at θ = π (antipode of the compact string, smallest H, most
conservative) and reported also at θ = π/2.

**Packets (one quantum of a massless field, ingoing at mouth A at angle θ).**
- P1: Gaussian spectrum, ω₀r₊ = 500, σ_ω r₊ = 100, truncated to ω > 0.
- P2: compactly supported smooth spectrum on ωr₊ ∈ [300, 700]: |f(ω)|² ∝ exp(−1/(1 − x²)), x = (ωr₊ − 500)/200.
- Arrival density in advanced time at mouth A: P(v) = |F(v)|², F(v) = (2π)^{−1/2} ∫ f(ω) e^{−iωv} dω (narrow-band
  approximation to the Klein–Gordon flux) [assumption]. Mean energy ⟨E⟩ = ħ∫ω|f|² dω.
- Placement: centre v_c = v_max(θ) − 100 r₊.

**Requirements (each must hold).**
- R1 window: at least 0.9 of P(v) lies in (v_early, v_max), v_early = −t_merge with t_merge = (d³/r₊)^{1/2}
  [inference from FGM §1, "the time scale for the black holes to merge is ∼ d^{3/2}"].
- R2 geometric optics: at least 0.99 of the spectral weight at ωr₊ ≥ 5 [assumption: well above the low-ℓ barrier,
  so FGM's null-geodesic treatment applies].
- R3 first-order regime: G⟨E⟩/r₊ ≤ 10⁻³ [FGM's linear back-reaction].
- R4 co-moving self stress: in 1+1 a right-moving signal makes only right-moving stress and does not hinder itself
  [lit: fn 6]; no condition.
- R5 counter-moving self stress (4D coupling, no equation in FGM): FGM's covering-space argument (fn 6: the self-delay
  "is independent of when the signal is sent") implies ΔV_self = αV [inference]; α unknown, not set to zero. Report
  α_crit = e^{κ₊(v_max − v_c)} − 1 and the shape factor A = α_crit r₊/(G⟨E⟩) that would be needed to close the window.
- R6 positive energy into mouth B (fn 6: "prevent the signal from sending perturbations into the other (right)
  mouth at early times"): infall of Killing energy δE into B at retarded label u_B adds δE e^{−κ₊u_B} to ∫T_UU dU,
  so it closes the window if δE > δE_crit(u_B) = (cκ₊/16) e^{κ₊(u_B − x*₀)} [der from FGM (2.20), (3.3)]. The packet's
  sub-barrier part (ωr₊ < 1) is taken as fully reflected at A (worst case), reaching B with fraction f_geo = (r₊/d)²,
  arriving at u_B = v_c + d. Its energy must be below δE_crit. The above-barrier reflection is not in FGM: report
  the bound R(ω) must meet, with the standard estimate R ∼ e^{−πω/κ₊} [assumption, external] for comparison only.
- R7 receiver isolation (fn 7): a passive receiver at ≥ 10 d from B, so detection happens after ≥ 10 d [design
  choice]; FGM's criterion is qualitative.
- R8 FGM validity: d ≫ r₊, μr₊² ≫ ln n, Gμ ≪ 1, NGμ ≪ 1, κ₊d ≫ 1 (non-extremal, perturbative) [lit: §2, §4].

**Decision rule.** **A** if at least one packet passes R1–R3, R6 (sub-barrier part) and R8, and every unknown (R5's
shape factor, R6's above-barrier reflection, the O(ln d/r₊) term in x*₀) would have to lie more than ten orders of
magnitude beyond its natural or standard value to change the outcome. **B** if every packet fails an evaluable
requirement. **C** otherwise. Two-way consequences (a reply; a crossing back-signal) are computed separately with
FGM's same linear equations and do not enter the verdict.

## 3. Results [der unless marked] (`scripts/fgm_window.py` → `results/fgm_window/window.json`)

**3.1 Window.** At θ = π: v_max = −5695 r₊ (−19.0 ms for r₊ = 1 km); minimal transit 11,396 r₊ = 38.0 ms against an
exterior light time of 33.4 ms (t_min/d = 1.14). At θ = π/2 the edge moves by 1.3 r₊. Our ΔV reproduces FGM's
(3.17) exactly once H is read as H/r₊² (unit test). The window is open toward the past: every ray entering A before
v_max traverses, and it emerges at B no earlier than −v_max = +5695 r₊.

**3.2 Packets.** Both are one quantum with ⟨E⟩ = 500 ħc/r₊ = 9.9 × 10⁻⁸ eV (about 24 MHz); arrival widths
σ_v = 0.005 r₊ (P1) and 0.0089 r₊ (P2). R1: more than 0.99999998 of each lies in the window (centre 100 r₊ before
the edge). The packet energy is the mean of its spectrum; ΔE·Δt is not used as a bound.

**3.3 The signal's own back-reaction.**
- Co-moving stress: no condition in 1+1 [lit: fn 6].
- Counter-moving stress through 4D couplings: FGM give no coefficient; α is kept unknown, not set to zero. With
  ΔV_self = αV the window becomes V < |ΔV|/(1 + α). At the frozen margin closure needs α > 7.2 × 10¹⁰. Since α is
  first order in G⟨E⟩/r₊ = 1.3 × 10⁻⁷³, its dimensionless shape factor would have to exceed 5.5 × 10⁸³.
- The closing condition for arbitrary positive energy is not applied to the signal's own co-moving stress. It is
  applied only to energy the packet sends into mouth B (reflection, leakage), which is what fn 6 names.

**3.4 Energy into mouth B.** At the reflection arrival of the frozen packet (u_B − x*₀ = −795 r₊) the closure
threshold is δE_crit = 1.5 × 10⁻⁸⁷ ħc/r₊ = 2.9 × 10⁻⁹⁷ eV; one quantum at the signal's own frequency exceeds it by
3.4 × 10⁸⁹.
- P1's sub-barrier tail (weight 1.5 × 10⁻⁸, energy 7.7 × 10⁻⁹ ħc/r₊) reaching B with f_geo = 10⁻⁸ exceeds the
  threshold by 5.3 × 10⁷⁰: **P1 fails R6.**
- P2 has no sub-barrier weight. Its above-barrier reflection must stay below 2.9 × 10⁻⁸²; the standard estimate
  R ∼ e^{−πω/κ₊} gives ln R = −3770 at the band edge [assumption]. **P2 passes R6 as frozen** (centre only).

**3.5 Post hoc (outside the frozen rule).**
- *P2's past tail.* A compact spectrum cannot give a time profile that falls exponentially (Paley–Wiener). For P2 the
  saddle-point asymptote is |F|² ∼ e^{−√(2WΔt)} (W = 200/r₊), while δE_crit falls as e^{−κ₊Δt}. With the standard
  reflection estimate, the reflected tail reaches B above threshold for arrivals more than 2.8 × 10⁴ r₊ earlier,
  inside the background's lifetime of 10⁶ r₊. Any nonzero above-barrier reflection gives such a crossover, so **P2
  fails once its past tail is counted.** An earlier draft used a decay rate √2 too fast (crossover 3.6 × 10⁴ r₊);
  the corrected value is stored.
- *A narrower Gaussian* (ω₀/σ_ω = 20, σ_v = 0.02 r₊, not frozen) puts 1.9 × 10⁻⁹⁷ ħc/r₊ into B, 1.3 × 10⁻¹⁰ of the
  threshold, and its Gaussian time tails fall faster than the threshold. It passes every evaluable condition here,
  by ten orders, which is comparable to the factors not counted (geometric fraction, angular kernel, the
  O(ln d/r₊) offsets).

## 4. Condition table

| # | Condition | Source | P1 (Gaussian, ω₀/σ_ω = 5) | P2 (compact spectrum) | Inputs not in FGM |
|---|---|---|---|---|---|
| R1 | ≥ 0.9 of the arrival probability in (−t_merge, v_max) | FGM (3.8)–(3.9), (2.20), (3.13)–(3.14), §1 | pass | pass | edge known to O(ln d/r₊)/κ₊; margin 100 r₊ |
| R2 | ≥ 0.99 of the spectrum at ωr₊ ≥ 5 | assumption (geometric optics) | pass | pass | — |
| R3 | G⟨E⟩/r₊ ≤ 10⁻³ | FGM first order | pass (1.3 × 10⁻⁷³) | pass | — |
| R4 | co-moving self stress | fn 6 | no condition | no condition | 1+1 statement |
| R5 | counter-moving self stress, ΔV_self = αV | fn 6 + inference | closes only if shape factor > 5.5 × 10⁸³ | same | coefficient unknown, not set to zero |
| R6 | energy into B below δE_crit, at the centre | fn 6, (2.20), (3.3) | **fail** (× 5.3 × 10⁷⁰) | pass; needs R < 2.9 × 10⁻⁸² | reflection, f_geo |
| R6′ | same, over the past tail (post hoc) | as R6 | pass (Gaussian tails) | **fail** (crossover 2.8 × 10⁴ r₊ earlier) | any nonzero reflection; lifetime 10⁶ r₊ |
| R7 | receiver isolation | fn 7 | by design | by design | FGM criterion qualitative |
| R8 | FGM validity | §2, §4 | pass | pass | d/r₊ = 10⁴, μr₊² = 3.8 × 10⁶⁸, Gμ = 10⁻⁷, NGμ = 10⁻⁶, κ₊d = 2500 |
| — | transmitter leakage, ambient radiation into B | fn 6 | not evaluated | not evaluated | must stay below 2.9 × 10⁻⁹⁷ eV at that time |

## 5. Verdict

- **Frozen rule: A through P2** (R1–R3, R6 at the centre and R8 pass; the R5 and R6 unknowns lie far beyond ten
  orders of magnitude).
- **Post hoc, A is withdrawn.** The frozen R6 checked only the packet's centre. Over P2's past tail, R6 fails for any
  nonzero above-barrier reflection within the background's lifetime. With the tail counted, both frozen packets fail
  an evaluable requirement at this point: **B for the frozen packets** (P1 through its sub-barrier spectral tail,
  P2 through its past tail).
- **Not settled:** whether some other packet passes. A narrow-band Gaussian (ω₀/σ_ω ≥ 20) passes every evaluable
  condition, post hoc and by a margin comparable with the uncounted factors. It is a candidate for a pre-registered
  check, not a result.

## 6. Two-way consequences (same first-order equations; not part of the verdict)

- **No reply after receipt.** The earliest receipt at B is +5695 r₊ (+19 ms), while a reply must enter B before
  −5695 r₊ (−19 ms; the two mouths are mirror images [inference]), a shortfall of 11,390 r₊ ≈ 38 ms [der].
- **No crossing exchange.** A traversing quantum's own co-moving stress crosses the other horizon blueshifted. For a
  reverse ray it is an external positive perturbation, about 10¹¹⁷⁵ times the strings' negative energy [der; the
  reverse negative energy is taken equal to the forward one by time reversal, inference]. So the reverse direction
  closes, and a pre-arranged two-way exchange through the same window is excluded at first order.

## 7. What this round excludes and leaves

- **Excluded (this point, first order):** broadband Gaussian packets; compact-spectrum packets with long past tails;
  replies after receipt; crossing two-way exchange.
- **Left open:** narrow-band packets (post-hoc candidate); the 4D self-coupling coefficient α; the actual RN
  reflection and two-centre geometric fraction; apparatus and environment isolation (≈ 3 × 10⁻⁹⁷ eV into the
  receiving mouth at the critical time, smaller still earlier).
- **Not claimed:** a channel; a bit (one quantum without a receiver and noise model is not a reliable bit); a rate;
  anything faster than light (transit 1.14 d here, longer than the exterior light time).
- **Next dynamic calculation:** only one is worth naming. For the one-way question, a pre-registered test-field
  calculation on a fixed RN background (time-domain linear wave equation for a narrow-band packet, with the energy
  reflected toward the second mouth compared with δE_crit(u_B)) would settle the candidate cheaply. For the
  communication goal it is not worth running: replies are excluded by timing and crossing signals by blocking,
  whatever the packet.

**Deviations (post hoc, labelled):** the extension of R6 to the past tail (R6′); the narrow-band Gaussian; the
corrected tail decay rate.
