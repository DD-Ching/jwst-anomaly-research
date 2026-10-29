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
