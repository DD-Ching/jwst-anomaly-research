# Goal round 7: rescoped Krasnikov evidence, a decided part of G2, and the G1 method

Owner brief (2026-10-11). Continue from #158:
1. Separate "a theorem of the right type exists" from "this model satisfies all its premises".
2. Separate negative-energy states on a fixed background from fields that self-consistently source the tube.
3. Assess G2 first, on a fully specified geometry, without claiming achronality from a mid-tube profile or from
   the absence of a counter-example.
4. Assess G1 for method feasibility only.
5. Scope every stop condition.

No new candidates, no heavy computation, all PRs unmerged.

Code `scripts/goal_krasnikov_g2.py`; results `results/goal_round7/krasnikov_g2.json`; tests
`tests/test_goal_krasnikov_g2.py`. Round-6 text is corrected in place (`docs/goal_round6/README.md`, marked
"round 7"). Labels: [derived] = computed or proved here; [lit] = stated in the cited source;
[ASSUMPTION] = a choice made here.

## 1. Quantum limits: a theorem type exists, but its premises are only estimated for this model

**What exists** [lit]:
- Fewster–Smith (gr-qc/0702056) prove an absolute QEI for the free, minimally coupled Klein–Gordon field along any
  smooth timelike curve in a globally hyperbolic spacetime.
- Kontou–Olum (1410.0665, Eq. 130) make it explicit along a geodesic, to first order in the Ricci tensor and its
  derivatives. They need bounds R_max, R″_max and R‴_max on the Fermi-frame components over the region involved.

**What round 6 actually computed** is a **conditional estimate**, not a numerical upper bound:
- ρ-derivatives of the static-frame Ricci components stand in for Fermi-frame derivatives.
- The tails of the Gaussian sampling were ignored.
- The counterterm coefficients a, b were left as an unknown factor.
- The O(R²) terms beyond Kontou–Olum's first order were not bounded.

**What a numerical upper bound would need:**

| Item | Status | What closes it |
|---|---|---|
| Frame | Fixed where ∂t is Killing [derived] | That is the eternal tube, or after the front has passed. ∂t is then Killing and geodesic, so ∇_a ξ_b annihilates u: the static frame differs from the Fermi frame by a pure spatial *rotation* (no boost), at rate ~k′. Rotation-invariant norms are frame-independent. |
| Derivatives | Open | Kontou–Olum's R″ and R‴ are Fermi-coordinate partial derivatives, which equal ∇ⁿRic only to first order. Use sup-norms of ∇ⁿRic (n = 0, 2, 3) plus that correction, over the region above. |
| Region | Open | Kontou–Olum need the bounds at every point in both the causal future and the causal past of points of the sampled path (1410.0665, conclusion). A Gaussian has unbounded support, so this is the causal hull of the whole worldline. In practice that is the whole tube, including caps, wall–cap junctions and, for the built tube, the front. |
| a, b | Open (physics input) | Their values come from the gravitational action and are not fixed by the QEI. Without them the bound holds only for 24\|a\| + 11\|b\| ≲ 10⁴ (where that term stays ≤ 10 %). |
| Second order in curvature | **Not available** | Kontou–Olum give no remainder bound beyond first order. At R t0² ≈ 10⁻⁴ the truncation is plausible, but no published estimate controls it. |

**Conclusion.** The round-6 latency-row numbers (ε_max ≈ 1.0–1.45 × 10⁴ l_P) are a conditional estimate. They hold
if:
- the first-order Kontou–Olum expansion is accurate;
- 24|a| + 11|b| ≲ 10⁴;
- the derivative stand-in is within O(1).

They are **not** a verified upper bound. Upgrading them needs a second-order remainder bound that the literature
read here does not supply.

## 2. Negative-energy states versus a self-consistent source

- **[lit] Fewster–Osterbrink (0708.2450 §3).** For ξ > 0 (Minkowski space), "Given any bounded subset O of
  Minkowski space and an arbitrary constant ρ0 > 0 we will construct a Hadamard state Ψ in which the expected
  energy density ⟨ρ⟩_Ψ is less than −ρ0 throughout O". Their §6 adds that such states "require more energy to
  achieve than positive energy densities of the same magnitude".
- So on a **fixed background**, macroscopic regions of arbitrarily negative energy density do exist for
  non-minimally coupled fields. Round 6's "no known example ... sustained macroscopic negative density" was too
  broad and is corrected.
- **What has not been found** is a field and a state shown to **self-consistently** source a Krasnikov tube. That
  needs all of the following:
  - the full ⟨T_ab⟩ equal to G_ab/8π, not only the energy density;
  - for the duration of use;
  - including back-reaction.
- Fewster–Osterbrink's states are bounded in space and time and carry compensating positive energy. Whether such
  a configuration can match ER's wall profile is not known.

## 3. G2: ANEC on the tube's achronal and shortcut rays

**Geometry specified.** The complete ER tube (Eqs. 12–13): the front θ(t − x − ρ) laid at light speed, both end
caps, and walls of width ε.

k = 1 − η θ(ρ_max − ρ) θ(t − x − ρ)[θ(x) − θ(x + ε − D)]

The bracket lies in [0, 1] for a monotone step, and η < 1, so k ≥ 1 − η > 0. Two idealisations are covered:
- the **built** channel (with the front);
- the **eternal** tube (no front factor).

Two step functions are used:
- ER Eq. 7: compact, exactly 0 or 1 outside [0, ε]. The code uses a C∞ smoothstep.
- ER's own Eq. 35: a tanh step whose tails reach everywhere.

The metric factorises: ds² = −(dt − dx)(dt + k dx) + dρ² + ρ²dφ² [lit: ER Eq. 4]. With u = t − x it becomes
ds² = −du² − (1 + k) du dx + dρ² + ρ²dφ².

**Result A: outgoing rays.** Any k(u, x, ρ), any step, built or eternal [derived; symbolic].
- **Achronal.** For a future causal vector, (dt − dx)(dt + k dx) ≥ dρ² + ρ²dφ² ≥ 0. With k > 0 both factors are
  ≥ 0, so du ≥ 0, strictly for timelike vectors. u is a null coordinate (g^uu = 0), not a time function, but it
  never decreases along causal curves. So each line L = {u, ρ, φ fixed} is **achronal**.
  - Monte Carlo check: dx/dt stays below 1 for random causal vectors.
- **Complete null geodesics.** g_xx = 0, the only nonzero Γ^a_xx is Γ^x_xx = ∂_x k/(1 + k), and the affine
  parameter is dλ ∝ (1 + k)dx with 1 + k ∈ [2 − η, 2].
- **R_xx ≡ 0, and the screen tidal tensor R_AxBx (A, B ∈ ρ, φ) vanishes.** The rays are shear-free, carry exactly
  zero ANEC, and fail the generic condition. This is consistent with Hawking–Ellis: a complete achronal null
  geodesic cannot satisfy both NEC and the generic condition.

**Result B: core shortcut rays with the compact step** [derived].
- Wherever k depends on x alone, the spacetime is flat (symbolic). With the compact step this is the core
  {ρ < ρ_max − ε, t − x − ρ > ε}, both caps included, plus the exterior.
- **They are the fastest −x paths.** For a future causal vector with dx < 0, dt − dx > 0 forces
  dt ≥ k(t, x, ρ)|dx|. For the ER form, k(t, x, ρ) ≥ k_ax(x) := 1 − η[θ(x) − θ(x + ε − D)], with equality in the
  core.
- **They stay in the core.** An incoming ray at fixed (ρ, φ) has d(t − x)/dt = 1 + 1/k > 0. If it enters at x = D
  after completion (t − D − ρ > ε) with ρ < ρ_max − ε, it never meets the front or the walls.
- **Numerical check:** 64 rays (η = 0.999, D = 400ε, ρ_max = 50ε) find ∂_ρ k = ∂_u k = 0 along every ray, with
  transit time 11.4ε against 410ε for light outside.
  - Positive control: rays launched before completion register ∂k = 2.0/ε, with transit 214–238ε.
- **So their ANEC is exactly 0, but this is a marginal case.** It needs a step that is exactly 1 inside the core.
  Oblique rays can reach the wall and are not covered.

**Result C: the axial shortcut ray, for any step monotone in ρ** [derived; symbolic, with numerical checks].
Premises:
- (P1) η < 1.
- (P2) k = k(x, ρ) along the ray. This holds in the eternal tube; in the built tube it needs a compact front and
  entry after completion.
- (P3) The axis minimises k(x, ·) for every x. This is automatic for ER Eq. 12 with any step that is monotone in ρ.

Then:
- **Achronal.** w := t + ∫k(x, 0)dx never decreases along causal curves. For dx ≤ 0: dt + k(x, ρ)dx ≥ 0 and
  k(x, 0) ≤ k(x, ρ). For dx > 0 it holds trivially. w increases strictly along timelike curves.
  - Monte Carlo check over the whole tube: min dw/|v| = +0.002, for both steps.
  - The axial incoming ray keeps w constant. It is a null geodesic: the radial acceleration ∝ ∂_ρ k vanishes on
    the axis, and the ρ = 0 plane is fixed by the rotations. It is complete, so it is a **complete achronal null
    geodesic**.
- **Its null curvature:** R_ab l^a l^b = −(1 + k0) k2 / k0² (symbolic), with k0 = k(x, 0) and k2 = ∂²_ρ k(x, 0).
  - P3 gives k2 ≥ 0, so **R_ll ≤ 0 along the whole ray**. ANEC is < 0 unless k2 = 0 almost everywhere on the axis.
  - Compact step (ρ_max ≥ ε): k2 ≡ 0, ANEC = 0, marginal.
  - ER's tanh step (Eq. 35), mid-tube, in units of 1/ε²:

    | η | ρ_max = 3ε | 10ε | 50ε |
    |---|---|---|---|
    | 0.95 | −5.3 × 10⁻⁵ | −2.5 × 10⁻²⁹ | −2.7 × 10⁻¹⁶⁸ |
    | 0.999 | −0.13 | −6.3 × 10⁻²⁶ | −6.7 × 10⁻¹⁶⁵ |

    These values are exponentially small but strictly negative.

**What this decides** (conditional on the self-consistent AANEC conjecture of Graham–Olum, which is not a theorem):
- An eternal tube, or a built one with a compact front, in class P1–P3, can be self-consistent only if
  **T_ab l^a l^b = 0 along the entire axis** (k2 ≡ 0).
  - ER's own tanh profile fails this.
  - The compact-step profile meets it only marginally.
  - Any semiclassical source whose null energy on the axis is negative anywhere, while the axis stays the fastest
    line, would conflict with the conjecture.
- Round 6's inference ("negative ANEC on an achronal shortcut ray via Olum") is thus **right in generic cases,
  wrong only in the marginal flat-core case**. Round 6 is annotated accordingly.

**What it does not decide:**
- Profiles where the axis is not the fastest line for some x (k2 < 0 there, i.e. positive null energy on the axis).
  The fastest rays are then elsewhere.
- The built tube with a tailed front (P2 fails).
- Oblique wall- and front-crossing geodesics, whose achronality is not determined.

## 4. G1: is the source computation feasible, and what would it decide?

**Question.** Can a free quantum field's renormalised ⟨T_ab⟩ supply ER's wall?

**Scaling statement** [derived, dimensional analysis]. Take a massless free field on the thin-wall profile, whose
only length is ε (with ρ_max/ε → ∞). Then

⟨T_ab⟩ = (ħ/ε⁴)[A_ab(η) + B_ab(η) ln(με)].

Here μ is the renormalisation scale; changing μ shifts the counterterm coefficients a, b. The required source is
T_ab = C_ab(η)/(8πG ε²). Self-consistency therefore gives

**ε_sc = l_P √(8π N F(η))**, with F a pure number fixed by the field, η, the profile shape and a, b.

What this means:
- F contains ln(με). Unless μ is tied to ε, the equation for ε_sc is implicit.
- A wall much thicker than l_P needs N F ≫ 1.
- Round 6's one-loop estimate corresponds to F ~ 10⁻⁴.
- Large a, b raise F. That is a change to the gravitational action, not a matter source, and it is not pursued
  (no new candidates).
- Massive fields:
  - with mε ≫ 1 (where DeWitt–Schwinger is valid), ⟨T⟩ ~ N/(m²ε⁶). Self-consistency then gives ε⁴ ~ N l_P²/m²,
    which together with mε ≫ 1 needs m ≫ m_P/√N. So for fewer than ~(m_P/m)² fields there is no semiclassical
    self-consistent wall in that regime;
  - with mε ≪ 1 they reduce to the massless case plus O(m²ε²).
- **So the DeWitt–Schwinger expansion cannot rule out light or massless fields; they need F(η) itself.**

**Method** (feasible offline, medium cost; not started):

| Element | Choice | Why |
|---|---|---|
| Geometry | Planar thin-wall limit: ds² = −dt² + (1 − k(y))dt dx + k(y)dx² + dy² + dz² | Killing symmetries in t, x and z make the wave equation separate into one radial ODE per mode (ω, p, q). |
| Global state | The ground state for ∂t. ∂t is timelike everywhere (g_tt = −1), so there is no ergoregion and positive Killing frequency is well defined. | For free fields on a stationary background, the long-time average along stationary curves is lowest in the ground state [lit: Kontou–Sanders §4.3]. So the ground state gives a state-independent necessary condition. |
| Built tube | The physical state is the Minkowski vacuum before construction, evolved forward. | The ground-state bound is proved for stationary spacetimes, so it covers the eternal idealisation only. For the built tube it is an assumption that the late-time stationary region behaves the same way. |
| Boundary conditions | None material. Mode functions bounded (delta-normalisable) as \|y\| → ∞. | The thin-wall planar limit drops the axis, so corrections are O(ε/ρ_max). The axis value needed by §3 needs the cylindrical problem. |
| Renormalisation | Point-splitting along a Killing direction: the pragmatic mode-sum method, t-splitting (Levi–Ori, arXiv:1503.02810 [lit, abstract only]). Hadamard subtraction, with a, b and μ kept explicit. | Per the abstract it requires only stationarity, though their demonstration is Schwarzschild. Adapting it to continuous planar transverse momenta is not covered. |
| Controls | Flat (η = 0) gives 0. The small-η limit must agree with flat-space linear response (to be sourced). Changing μ must reproduce the a, b shift. | Validates the code before any η ≈ 1 number. |
| Output | F(η) for the minimally coupled scalar, all components of ⟨T_ab⟩ against C_ab, and the sign of ⟨T_ll⟩ on the axis (cylindrical version; §3). | The energy density alone is necessary, not sufficient. The conformal scalar (ξ = 1/6) is excluded: for ξ > 0 the ground state is not established as the long-time minimum (§2). |

**What it cannot decide.**
- **Finite use versus infinite averages.** The ground-state bound concerns infinite-time averages. For a channel
  used for a finite time T_use, a difference QEI allows a deficit below the ground state that shrinks as T_use
  grows. Its constant on this background is not known, so finite-use statements are conditional.
- **Back-reaction and stability.** Even if F were large, matching all of C_ab, back-reaction and stability would
  still be untested.
- **Other fields.** Interacting and non-minimally coupled fields are outside its reach (§2).

## 5. Scoped stop conditions

- **G1 failing** (N F ≲ 1 for the computed fields) would rule out only those free fields, under the stated state,
  renormalisation and planar-wall assumptions.
- **An AANEC conflict** (§3: the axial ray of any P1–P3 tube with k2 ≠ 0) limits the tube only if the
  self-consistent AANEC conjecture holds.
- **Pausing for low research value** is a resource decision. It is not a statement that the tube is physically
  impossible.

## 6. Answers

**New, testable physical judgement** [derived; symbolic plus numerical checks]:
- In the Everett–Roman tube with η < 1 (eternal, or built with a compact front), take any step that is monotone in
  ρ. Then the axial shortcut ray is a **complete achronal null geodesic**.
- Its null curvature is R_ll = −(1 + k0)k2/k0² ≤ 0, so its ANEC integral is negative unless ∂²_ρ k = 0 along the
  whole axis.
- The outgoing rays are achronal with ANEC exactly 0.
- **Consequence, conditional on the self-consistent AANEC conjecture:** a self-consistent source must make the
  null energy T_ll vanish along the entire axis, or make the axis not the fastest line. ER's own tanh profile
  violates this; a compact flat-core profile meets it only marginally.
- **How to test it:**
  - recompute R_ll on the axis for another profile;
  - find a timelike curve with Δw < 0 under P1–P3;
  - or show the conjecture fails.

**Does it change the research decision?** Yes, in two ways.
- G2 is no longer a cheap way to *exclude* the tube, but it is not idle either. It turns into a sharp necessary
  condition: on-axis T_ll = 0. That condition joins G1: any G1 computation must also report the sign of ⟨T_ll⟩ on
  the axis, which needs the cylindrical problem rather than the planar limit.
- A source passing G1 (large enough N F) would still face this axis condition, which looks like a fine-tuning
  requirement along 20–100 ly. That lowers the expected value of G1 further.

**Proposal for the owner.**
- Pause the Krasnikov line for research value. Its numbers (≥ 10³⁸ M☉, Planck-thin walls on the one-loop scale)
  and the on-axis AANEC condition make any engineering use remote.
- Keep G1, with the axis-sign output, as a well-posed unit.
- The pause is not a claim of physical impossibility. The AANEC condition is conditional on a conjecture.

## 7. Bookkeeping

- **Read this round:** ER §§2–3 (geometry; previously read), Fewster–Osterbrink §§3 and 6 (re-checked), and the
  Levi–Ori abstract (arXiv:1503.02810; method only, not relied on for any result).
- **Independent review** (21 findings in round 6; 12 in this round) applied. Its main finding: "exactly zero
  ANEC" held only for compact steps. That led to Result C.
- **Not covered:**
  - profiles where the axis is not fastest, tailed fronts, and oblique wall- and front-crossing geodesics (§3);
  - the second-order QEI remainder (§1);
  - the G1 mode sum itself (§4).
