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
| Frame | Fixed [derived] | ∂t is Killing and geodesic, so ∇_a ξ_b annihilates u: the static frame differs from the Fermi frame by a pure spatial *rotation* (no boost), at rate ~k′. Rotation-invariant component norms are frame-independent. |
| Derivatives | Open | Replace ∂_ρ of components by sup-norms of the covariant derivatives ∇ⁿRic (n = 0, 2, 3) on the wall, computed symbolically. |
| Gaussian tails | Partly fixed | In the eternal idealisation the static observer stays in the wall forever, so sup-norms over the wall bound all times. For the built tube, the sup must also cover the front (t − x − ρ ∈ [0, ε]) in the observer's past. |
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

**Geometry specified.** The complete ER tube (Eqs. 12–13), i.e. the built channel:
- the front θ(t − x − ρ), laid at light speed;
- both end caps;
- walls of width ε;
- k = 1 − η θ(ρ_max − ρ) θ(t − x − ρ)[θ(x) − θ(x + ε − D)], with 0 ≤ the bracket ≤ 1 and η < 1 (so k ≥ 1 − η > 0).

The eternal idealisation is this metric without the front factor. Both are covered below.

The metric factorises: ds² = −(dt − dx)(dt + k dx) + dρ² + ρ²dφ². With u = t − x this becomes
ds² = −du² − (1 + k) du dx + dρ² + ρ²dφ² [lit: ER Eq. 4; the rest is derived].

**Result A: outgoing rays are achronal and carry zero null energy** [derived; symbolic, any k(u, x, ρ)].
- **u is a time function.** For a future causal vector, (dt − dx)(dt + k dx) ≥ dρ² + ρ²dφ² ≥ 0. With k > 0 both
  factors are non-negative, so du ≥ 0, strictly for timelike vectors. Hence each line L = {u, ρ, φ fixed} is
  **achronal**: a timelike curve between two points of L would need Δu > 0.
  - Monte Carlo check: the maximum dx/dt over random causal vectors stays below 1.
- **They are complete null geodesics.** g_xx = 0, and the only nonzero Γ^a_xx is Γ^x_xx = ∂_x k/(1 + k), so L is
  a null pregeodesic. Its affine parameter is dλ ∝ (1 + k)dx, which covers (−∞, ∞).
- **R_xx ≡ 0 for every k(u, x, ρ).** So the ANEC integral along every outgoing ray is exactly 0. This holds through
  the walls, the end caps and the front, including rays that ride inside the front from x = 0 to D.

**Result B: the shortcut rays sent after completion lie entirely in flat spacetime** [derived].
- **Flatness.** Wherever k depends on x alone, the Riemann tensor vanishes (symbolic). In the ER tube this is the
  core {ρ < ρ_max − ε, t − x − ρ > ε}, end caps included, plus everything outside 0 < x < D.
- **The rays stay in the core.** An incoming ray has d(t − x)/dt = 1 + 1/k > 0. If it enters at x = D after the
  front has passed (t − D − ρ > ε) at ρ < ρ_max − ε, it never meets the front or the walls.
  - Numerical check: 64 rays through the ER profile (η = 0.999, D = 400ε, ρ_max = 50ε) find ∂_ρ k = ∂_u k = 0 along
    every ray.
  - Their transit time is 11.4ε against 410ε for exterior light, so they are the shortcut.
- **Consequences:**
  - The complete geodesic has Riemann = 0 everywhere, so its ANEC integral is exactly 0.
  - The generic condition fails.
  - Olum's Condition 1 fails: parallel neighbours in the flat core arrive together.
  - Their achronality is therefore irrelevant to AANEC.

**What this decides.**
- Two families are now settled, and neither violates ANEC, so neither conflicts with Graham–Olum's self-consistent
  AANEC conjecture:
  - the achronal outgoing rays;
  - the fastest shortcut rays used after completion.
- Round 6's inference (an achronal shortcut geodesic must carry negative ANEC by Olum's defocusing) was wrong for
  these rays and is corrected there.

**What it does not decide.** Some null geodesics cross the walls, the wall–cap junction or the front obliquely.
- Their R_ll takes both signs (round 6 §1).
- Their achronality is **not determined**. Neither a mid-tube profile nor finding no counter-example settles it.
  It needs the global causal structure of the tube.
- Whether any of them is complete, achronal and ANEC-negative is **undecided**.

An AANEC conflict, if one were found, would also be a limit conditional on the conjecture. It would not be a
theorem.

## 4. G1: is the source computation feasible, and what would it decide?

**Question.** Can a free quantum field's renormalised ⟨T_ab⟩ supply ER's wall?

**Scaling statement** [derived, dimensional analysis]. Take a massless free field on the thin-wall profile, whose
only length is ε (with ρ_max/ε → ∞). Then

⟨T_ab⟩ = (ħ/ε⁴)[A_ab(η) + B_ab(η) ln(με)].

Here μ is the renormalisation scale; changing μ shifts the counterterm coefficients a, b. The required source is
T_ab = C_ab(η)/(8πG ε²). Self-consistency therefore gives

**ε_sc = l_P √(8π N F(η))**, with F a pure number fixed by the field, η, the profile shape and a, b.

What this means:
- A wall much thicker than l_P needs N F ≫ 1.
- Round 6's one-loop estimate corresponds to F ~ 10⁻⁴.
- Large a, b raise F. That is a change to the gravitational action, not a matter source, and it is not pursued
  (no new candidates).
- Massive fields:
  - with mε ≫ 1 (where DeWitt–Schwinger is valid), ⟨T⟩ ~ 1/(m²ε⁶), so ε_sc is smaller still;
  - with mε ≪ 1 they reduce to the massless case plus O(m²ε²).
- **So the DeWitt–Schwinger expansion cannot rule out light or massless fields; they need F(η) itself.**

**Method** (feasible offline, medium cost; not started):

| Element | Choice | Why |
|---|---|---|
| Geometry | Planar thin-wall limit: ds² = −dt² + (1 − k(y))dt dx + k(y)dx² + dy² + dz² | Killing symmetries in t, x and z make the wave equation separate into one radial ODE per mode (ω, p, q). |
| Global state | The ground state for ∂t. ∂t is timelike everywhere (g_tt = −1), so there is no ergoregion and positive Killing frequency is well defined. | For free fields on a stationary background, the long-time average along stationary curves is lowest in the ground state [lit: Kontou–Sanders §4.3]. So the ground state gives a state-independent necessary condition. |
| Built tube | The physical state is the Minkowski vacuum before construction, evolved forward. | The ground-state bound is proved for stationary spacetimes, so it covers the eternal idealisation only. For the built tube it is an assumption that the late-time stationary region behaves the same way. |
| Boundary conditions | None material. Regular at the wall's two flat sides, outgoing/regular at \|y\| → ∞. | The thin-wall planar limit drops the axis, so corrections are O(ε/ρ_max). |
| Renormalisation | Point-splitting along a Killing direction: the pragmatic mode-sum method, t-splitting for stationary backgrounds (Levi–Ori, arXiv:1503.02810 [lit, abstract only]). Hadamard subtraction, with a, b and μ kept explicit. | It is designed for stationary, non-static backgrounds. |
| Controls | Flat (η = 0) gives 0. The small-η limit must agree with flat-space linear response (to be sourced). Changing μ must reproduce the a, b shift. | Validates the code before any η ≈ 1 number. |
| Output | F(η) for the minimal and conformal scalar, and all components of ⟨T_ab⟩ against C_ab. | The energy density alone is necessary, not sufficient. |

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
- **An AANEC conflict** (none found; §3) would limit the tube only if the self-consistent AANEC conjecture holds.
- **Pausing for low research value** is a resource decision. It is not a statement that the tube is physically
  impossible.

## 6. Answers

**New, testable physical judgement** [derived; checked symbolically and numerically]:
- In the complete Everett–Roman tube with η < 1, every outgoing null ray is a complete achronal null geodesic with
  R_ab l^a l^b ≡ 0.
- Every shortcut ray used after completion runs through exactly flat spacetime.
- So the ANEC integral is exactly zero on both families. Neither can conflict with the self-consistent AANEC
  conjecture, and Olum's superluminal-travel theorem does not apply to the core shortcut rays (Condition 1 and the
  generic condition fail).
- **How to test it:**
  - recompute R_xx for any other k(u, x, ρ);
  - find a timelike curve with Δu = 0;
  - find a core ray that meets ∂_ρ k ≠ 0 or ∂_u k ≠ 0.

**Does it change the research decision?** Yes.
- Round 6 ranked G2 as the cheapest possible field-independent limit. Its two natural targets are now settled with
  zero ANEC.
- What remains of G2 is a global-causal-structure problem for wall- and front-crossing geodesics. It is undecided
  and not cheap, so it is demoted.
- The only open discriminator is G1. It reduces to one pure number F(η) per field (§4), computable with a planar
  mode sum. It is feasible but medium cost, and a negative result would rule out only the fields computed.

**Proposal for the owner.**
- Pause the Krasnikov line for research value. Its numbers (≥ 10³⁸ M☉, Planck-thin walls on the one-loop scale)
  make any engineering use remote.
- Keep G1 as a well-posed, ready-to-run unit if the owner wants a decisive free-field answer.
- This pause is not a claim of physical impossibility.

## 7. Bookkeeping

- **Read this round:** ER §§2–3 (geometry; previously read), Fewster–Osterbrink §§3 and 6 (re-checked), and the
  Levi–Ori abstract (arXiv:1503.02810; method only, not relied on for any result).
- **Not covered:**
  - wall- and front-crossing geodesics (§3);
  - the second-order QEI remainder (§1);
  - the G1 mode sum itself (§4).
