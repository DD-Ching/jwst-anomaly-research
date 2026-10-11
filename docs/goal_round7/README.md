# Goal round 7: rescoped Krasnikov evidence, the axial-ray ANEC result for smooth profiles, and the G1 method

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

**Geometry specified.** The ER tube (Eqs. 12–13): the front θ(t − x − ρ) laid at light speed, both end caps, and
walls of width ε.

k = 1 − η θ(ρ_max − ρ) θ(t − x − ρ)[θ(x) − θ(x + ε − D)]

The bracket lies in [0, 1] for a monotone step, and η < 1, so k ≥ 1 − η > 0. Two idealisations:
- the **built** channel (with the front);
- the **eternal** tube (no front factor).

Steps:
- ER Eq. 7: compact, exactly 0 or 1 outside [0, ε]. The code uses a C∞ smoothstep.
- ER's own Eq. 35: a tanh step with tails.

The metric factorises: ds² = −(dt − dx)(dt + k dx) + dρ² + ρ²dφ² [lit: ER Eq. 4]. With u = t − x it becomes
ds² = −du² − (1 + k) du dx + dρ² + ρ²dφ².

**Axis regularity: the ER metric as written is not smooth on the axis** [derived; owner's check, revised in
round 7]. A function of ρ alone is differentiable at the axis, in Cartesian (y, z), only if its ρ-derivative
vanishes there. It is C^n only if it is a C^n function of ρ² = y² + z².
- **ER's tanh radial factor.** k_ρ(0) = η θ′(R) ≠ 0. For η = 0.999 and R = 10ε, k_ρ(0) = 7.87544 × 10⁻³³/ε, so
  ∂_y k jumps by 2k_ρ(0) across the axis.
  - The metric is only Lipschitz on the axis, and the general null curvature keeps a term
    −(1 + k)k_ρ/(2ρk²), which diverges like 1/ρ there.
  - The finite on-axis values reported before for the tanh profile (e.g. −6.30665 × 10⁻²⁶/ε²) are **withdrawn**:
    the axis expansion k = k0 + k2ρ²/2 + … does not apply to that profile.
- **The front factor θ(t − x − ρ).** Even with the compact step, it has ∂_ρ = −θ′(t − x) ≠ 0 on the axis while the
  front crosses it (0 < t − x < ε). The code finds −2.0/ε at t − x = ε/2 and 0 once the front has passed.
  - So the **built** ER geometry is not C¹ on the axis segment {ρ = 0, 0 < t − x < ε, 0 < x < D}.
  - Until that set is regularised, which would be a **model modification** and is not done here, claims about a
    "complete smooth geometry" for the built tube are withdrawn. Geodesic statements below are limited to
    geodesics that avoid that set.
- **What is smooth:**
  - the compact radial step with R > ε, which is identically 1 near the axis;
  - radial factors that are smooth functions of ρ², such as exp[−(ρ/R)²] and exp[−(ρ/R)⁴] (§3, Result C checks).

**Result A: outgoing rays** [derived; symbolic, any k(u, x, ρ)].
- **Achronal.** For a future causal vector, (dt − dx)(dt + k dx) ≥ dρ² + ρ²dφ² ≥ 0. With k > 0 both factors are
  ≥ 0, so du ≥ 0, strictly for timelike vectors. u is a null coordinate (g^uu = 0), not a time function.
  - The argument uses only the pointwise cone, so it holds for piecewise-C¹ timelike curves even where the metric
    is merely continuous.
  - Each line L = {u, ρ, φ fixed} is therefore **achronal** as a point set.
- **Complete null geodesics.** Where the metric is smooth along L, g_xx = 0, the only nonzero Γ^a_xx is
  Γ^x_xx = ∂_x k/(1 + k), and dλ ∝ (1 + k)dx with 1 + k ∈ [2 − η, 2].
  - Scope: lines with ρ > 0, and axis lines where k is smooth at the axis.
  - Excluded: built-tube axis lines with 0 < u < ε, which lie in the front–axis set.
- **R_xx ≡ 0, and the screen tidal tensor R_AxBx (A, B ∈ ρ, φ) vanishes** for any k(u, x, ρ). The rays are
  shear-free, carry ANEC = 0, and fail the generic condition. This is consistent with Hawking–Ellis.

**Result B: core shortcut rays with the compact step** [derived].
- Wherever k depends on x alone, the spacetime is flat. With the compact step and R > ε this is the core
  {ρ < R − ε, t − x − ρ > ε}, caps included, plus the exterior. It is smooth, including at the axis.
- **Fastest −x paths.** For a future causal vector with dx < 0, dt ≥ k(t, x, ρ)|dx| ≥ k_ax(x)|dx|, where
  k_ax(x) = 1 − η[θ(x) − θ(x + ε − D)]. Equality holds in the core.
- **They stay in the core.** An incoming ray at fixed (ρ, φ) has d(t − x)/dt = 1 + 1/k > 0. If it enters at x = D
  after completion (t − D − ρ > ε) with ρ < R − ε, it never meets the front or the walls.
- **Numerical check:** 64 rays (η = 0.999, D = 400ε, R = 50ε) find ∂_ρ k = ∂_u k = 0 along every ray, with
  transit 11.4ε against 410ε for light outside.
  - Positive control: rays launched before completion register ∂k = 2.0/ε, with transit 214–238ε.
- **Their ANEC is exactly 0, a marginal case** that depends on a step that is exactly 1 in the core.

**Result C: the axial shortcut ray, eternal tube, radial factor smooth at the axis** [derived; symbolic, with
numerical checks].

Premises (required regularity stated):
- (P1) η < 1.
- (P2) Eternal tube: k = k(x, ρ), so ∂t is Killing.
  - The built tube would also need its front regularised near the axis, which would be a model modification
    (not done).
- (P3) k = K(x, ρ²) with K ∈ C⁴ in both arguments: the metric is C⁴ and the curvature C², smooth at the axis.
  The exterior is asymptotically flat (k → 1 as |x| → ∞, with ∫|∂²_ρk(x, 0)|dx < ∞).
- (P4) The axis minimises k(x, ·) for every x. For k = 1 − η g(ρ) B(x) with B ≥ 0, this means g is non-increasing.

Then:
- **Achronal.** w := t + ∫k(x, 0)dx never decreases along causal curves (cone inequality plus P4). It increases
  strictly along timelike curves.
  - Monte Carlo check of the cone inequality: min dw/|v| = +0.002.
  - The axial incoming ray keeps w constant. It is a null geodesic: the radial acceleration ∝ ∂_ρ k = 0 on the
    axis by P3, and the ρ = 0 plane is the fixed set of the rotations. With P3 it is complete, so it is a
    **complete achronal null geodesic**.
- **Local curvature, not yet ANEC.** The general incoming R_ll, normalised with l^t = 1, is
  −[(1 + k)k_ρρ − k_ρ²]/(2k²) − (1 + k)k_ρ/(2ρk²). Under P3, k_ρ/ρ → k2 := ∂²_ρk(x, 0), so on the axis
  R_ll → −(1 + k0)k2/k0², with k0 = k(x, 0). P4 gives k2 ≥ 0.
- **Affine normalisation and ANEC.** Normalise to unit exterior Killing energy, E = −g(∂t, l) = 1:
  - l = (2k0/(1 + k0), −2/(1 + k0), 0, 0); checked symbolically: g(l, l) = 0 and E = 1.
  - This is the affine tangent, because ∂t is Killing (P2). Then R_ab l^a l^b = −4k2/(1 + k0) and
    dλ = −(1 + k0)dx/2.
  - So **ANEC = ∫R_ab l^a l^b dλ = −2∫k2(x)dx**. It converges under P3's integrability condition and is ≤ 0, with
    equality iff k2 = 0 almost everywhere on the axis.
- **Checks on two smooth radial factors** (mathematical examples, not source models; η = 0.999, R = 10ε,
  D = 400ε, compact box B(x)):

  | Radial factor g(ρ) | k2 (mid-tube) | −2∫k2 dx | Direct ∫R_ll dλ at ρ = 10⁻⁵R | R_ll at ρ = R/2 (l^t = 1) |
  |---|---|---|---|---|
  | exp[−(ρ/R)²] | 0.020/ε² > 0 | −15.944/ε | −15.944/ε (agrees to 2 × 10⁻⁷) | −0.23/ε² |
  | exp[−(ρ/R)⁴] | 0 | 0 | −6.4 × 10⁻⁹/ε (→ 0 as ρ²) | −4.6/ε² |

  The second example shows that k2 = 0 needs neither a compact step nor a flat core: the curvature off the axis
  is nonzero.

**What this decides** (three kinds of statement, kept separate):
- **Geometry, recomputable without any conjecture:**
  - Under P1–P4 the axial shortcut ray is complete and achronal, with ANEC = −2∫k2 dx ≤ 0.
  - Outgoing rays (Result A) have ANEC = 0.
  - The ER tanh profile and the built front are not smooth on the axis, so these statements do not apply to them
    as written.
- **Conditional on Graham–Olum's self-consistent AANEC conjecture** (not a theorem): a self-consistent eternal tube
  satisfying P1–P4 needs k2 = 0 almost everywhere on the axis, i.e. null energy T_ll = 0 on the axis. Equivalently,
  if a self-consistent source keeps the axis as the fastest line, it cannot have negative null energy anywhere on
  it.
- **Experimentally testable predictions:** none. Nothing here predicts a laboratory or astronomical observable.

**What it does not decide:**
- Profiles where the axis is not the fastest line.
- Any regularised built tube (a model modification).
- Profiles that are not smooth at the axis.
- Oblique wall- and front-crossing geodesics, whose achronality is not determined.

Round 6's inference ("negative ANEC on an achronal shortcut ray via Olum") holds under P1–P4 whenever k2 ≠ 0 on a
set of positive measure. It fails in the k2 ≡ 0 cases, which include, but are not limited to, flat cores.

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
- **An AANEC conflict** (§3: the axial ray of an eternal P1–P4 tube with k2 ≠ 0) limits the tube only if the
  self-consistent AANEC conjecture holds.
- **Pausing for low research value** is a resource decision. It is not a statement that the tube is physically
  impossible.

## 6. Answers

**What this round adds that can be checked** [derived; symbolic plus numerical checks]:
- **Geometry, valid without any conjecture:**
  - In the eternal Everett–Roman-type tube with η < 1, take a radial factor that is smooth at the axis (P3) and
    makes the axis the fastest line (P4). Then the axial shortcut ray is a complete achronal null geodesic.
  - With unit exterior Killing energy, its ANEC integral is exactly −2∫∂²_ρk(x, 0)dx ≤ 0.
  - Outgoing rays are achronal and shear-free, with ANEC = 0.
  - ER's own tanh profile, and its built front, are not smooth on the axis (k_ρ(0) = 7.87544 × 10⁻³³/ε for
    η = 0.999, R = 10ε). The axis statements do not apply to them until they are regularised, which would be a
    model modification.
- **Conditional on the self-consistent AANEC conjecture:** a self-consistent source of such a tube needs zero
  null energy along the axis, or an axis that is not the fastest line.
- **Experimentally testable predictions:** none.
- **How to check it:**
  - recompute ∂²_ρk on the axis for another smooth profile;
  - find a timelike curve with Δw < 0 under P1–P4;
  - or settle the conjecture.

**Does it change the research decision?** Partly.
- G2 becomes a necessary condition, conditional on AANEC: on-axis T_ll = 0. It is not an exclusion.
- Any future G1 computation should also report the sign of ⟨T_ll⟩ on the axis, which needs the cylindrical
  problem.
- **Retracted:** the earlier remark that this condition is a "fine-tuning requirement along 20–100 ly". No
  perturbation analysis supports it. The exp[−(ρ/R)⁴] example satisfies k2 = 0 with a curved core, and how a
  self-consistent source would respond near the axis has not been computed.

**Status (owner decision 2026-10-11).**
- This specific Krasnikov model (ER metric, η < 1) is **paused at "no self-consistent physical source found"**.
- The G1 method draft (§4) is kept. No large computation is started.
- The pause is not a claim of physical impossibility, and the AANEC condition is conditional on a conjecture.

## 7. Bookkeeping

- **Read this round:** ER §§2–3 (geometry; previously read), Fewster–Osterbrink §§3 and 6 (re-checked), and the
  Levi–Ori abstract (arXiv:1503.02810; method only, not relied on for any result).
- **Independent review** (21 findings in round 6; 12 in this round) applied. Its main finding: "exactly zero
  ANEC" held only for compact steps. That led to Result C.
- **Owner review (2026-10-11):** the axis regularity of the tanh profile and of the front. Result C is now
  restricted to smooth radial factors in the eternal tube, the finite tanh axis values are withdrawn, and the
  normalised ANEC cross-check and the two smooth examples are added.
- **Not covered:**
  - profiles where the axis is not fastest, tailed fronts, and oblique wall- and front-crossing geodesics (§3);
  - the second-order QEI remainder (§1);
  - the G1 mode sum itself (§4).
