# Goal round 6: can a concrete field and state source the Krasnikov tube of round 5?

Owner brief (2026-10-11). Continue from round 5 (#157). Do not add candidates. Check the physical source and the
quantum limits of the Krasnikov tube that was actually computed:
- compare Everett–Roman (gr-qc/9702049, "ER"), Krasnikov 2003 (gr-qc/0207057) and the follow-up literature;
- decide whether each counter-example and each limit applies to this geometry;
- check the sampling-time and curvature approximations, and the model, size and observer behind the negative-energy
  estimate;
- restate "no semiclassical wall exists" and "no matter model exists" as what the evidence supports;
- answer whether a concrete field and state can self-consistently support the channel, and if not, what the
  smallest computable gap is;
- distinguish "one second early" from "only one second of waiting".

Code `scripts/goal_krasnikov_source.py`; results `results/goal_round6/krasnikov_source.json`; tests
`tests/test_goal_krasnikov_source.py`. Offline; no data. Labels: [derived] = computed here; [lit] = stated in the
cited paper; [ASSUMPTION] = a choice made here.

## 0. "One second early" versus "one second of waiting"

One tube of opening η along a distance D. A message sent from the far end at t_s is read at t_s + D(1 − η); exterior
light arrives at t_s + D (round 5 §3.1). Both times are on the common clock of the two (static, mutually at rest)
ends.

| Goal | Quantity | η | Status |
|---|---|---|---|
| "1 s early" | advance ηD = 1 s | 1 s / D ≈ 1.6 × 10⁻⁹ (20 ly) | Round 5's advance rows. Nearly useless for communication: the latency is still ≈ D. |
| "wait only 1 s" | latency D(1 − η) = 1 s | 1 − 1 s / D ≈ 1 − 1.6 × 10⁻⁹ | The owner's target. This round's latency rows. |

- Two-way, low latency means two tubes, one per direction, each with η just below 1. Closed timelike curves need η > 1
  (round 5 §4), so the target sits just below that threshold, not beyond it.
- The tubes are laid at ≤ c. A channel over D is usable only after the build time, ≥ D/c, has passed: ≥ 20 yr or
  ≥ 100 yr [lit: ER §II; Krasnikov 2003 Remark]. The owner's goal concerns the channel after it is built, so this
  is a fixed upfront cost, not a defect.
- Nothing in the latency rows changes between 1 s, 1 min, 1 day and 1 yr of latency: η ≥ 0.95 throughout, and the
  wall numbers below differ by < 5 %. The physics is that of an η ≈ 1 tube.

## 1. Facts about the computed geometry [derived]

ER metric at mid-tube: ds² = −dt² + (1 − k)dx dt + k dx² + dρ² + ρ²dφ², with k = 1 − η θ_ε(ρ_max − ρ).
1. **Stationary, not static.** ∂t is a Killing vector, but the (t, x) cross term does not vanish in the wall. The
   literature subagent checked that no timelike Killing field is hypersurface-orthogonal there. Theorems stated for
   static spacetimes (Pfenning–Ford 1998; Fewster–Roman 2005 §8) therefore do not cover the tube directly.
2. **Static observers are geodesic.** g_tt = −1, so the worldlines along ∂t have zero acceleration. The QI
   observer of ER (and of this round) is a geodesic, as curved-spacetime QEIs require.
3. **Global hyperbolicity for η < 1.** g^tt = −4k/(1 + k)² < 0 wherever k > 0, and k ≥ 1 − η > 0 for η < 1. So t is a
   time function. The coordinate light speeds are +1 and −1/k, bounded by 1/(1 − η). So every inextendible causal
   curve crosses each t = const slice: the slices are Cauchy surfaces (mid-tube metric). For η ≈ 1 − 10⁻⁹ this
   still holds, with incoming light speeds up to ~10⁹ (in units of c). The full tube with its front and end caps was
   not checked.
4. **Where the null energy condition fails.** R_ab l^a l^b = 0 identically along the outgoing null direction
   l = ∂t + ∂x. Curves with that direction at fixed ρ are null geodesics (zero radial acceleration). Along the
   incoming (shortcut) direction l = ∂t − k⁻¹∂x, R_ab l^a l^b takes both signs in the wall. Its minimum, with
   l^t = 1, is −0.062/ε² at η = 0.01 and −145/ε² at η = 0.95. These values depend on the normalisation of l; only
   the sign is invariant. So NEC is violated only for light travelling in the shortcut direction.

## 2. Sampling time and curvature validity

ER bound the most negative static density by the flat-space QI, |T| ≤ 3/(32π²τ0⁴). The sampling time is
τ0 = σ r_c with σ = 0.01 [ASSUMPTION, ER's and Ford–Roman's choice].

**What the literature supports** [lit]:
- Ford–Roman 1996 §2 and Roman 2004: the flat QI should hold for τ0 small compared with the smallest local curvature
  radius and the distance to any boundary. This is a heuristic based on near-plane-wave modes.
- Fewster–Roman 2005 §2.2 call it "strictly speaking ... an assumption".
- Krasnikov 2003 and 2004 (gr-qc/0409007) call the curved-spacetime version "a conjecture". That was accurate then.
- Rigorous since then:
  - Fewster–Smith 2008 (gr-qc/0702056) prove an absolute QEI for the minimally coupled Klein–Gordon field on any
    4D globally hyperbolic spacetime. It needs no symmetry, so it covers stationary, non-static metrics. The
    curvature dependence is implicit.
  - Kontou–Olum (arXiv:1410.0665) make it explicit to first order in the Ricci tensor along a geodesic (their
    Eq. 130). The corrections are small when τ0 is below every curvature radius and below the scale on which the
    curvature varies.

**Applied to our wall** [derived]. The relative size of the Kontou–Olum correction terms, at t0 = σ r_c, against the
flat term 3.76:

| Row | η | t0 / ε | Ricci term | R″ term | R″ term per unit (24\|a\| + 11\|b\|) | R‴ term |
|---|---|---|---|---|---|---|
| wait ≤ 1 yr (20, 100 ly) | 0.95 – 1 | 0.005 | 7 × 10⁻⁵ | 2 × 10⁻⁷ | 1 × 10⁻⁵ | 3 × 10⁻⁸ |
| 1 yr early, 100 ly | 0.01 | 0.06 | 7 × 10⁻⁵ | 3 × 10⁻⁵ | 1.4 × 10⁻³ | 5 × 10⁻⁵ |
| 1 day early, 20 ly | 1.4 × 10⁻⁴ (row: 10⁻⁴) | 0.57 | 7 × 10⁻⁵ | 3 × 10⁻³ | 0.15 | 0.05 |
| 1 s early, 20 ly | 1.6 × 10⁻⁹ | 140 | 7 × 10⁻⁵ | 160 | 9 × 10³ | 7 × 10⁵ |

The derivatives are taken along ρ of the static-frame Ricci components: an order-of-magnitude stand-in for the
Fermi-frame derivatives. a and b are the undetermined coefficients of the curvature-squared counterterms.

**What it settles.**
- **Latency rows (the owner's target).** τ0 = σ r_c ≈ 0.005ε is far inside the curvature radius (r_c ≈ 0.48ε) and
  the variation scale (≈ ε). The first-order corrections are ≤ 10⁻⁴, unless 24|a| + 11|b| ≳ 10⁵. So, for a free,
  minimally coupled scalar, the flat-space bound used for these rows is backed by a theorem (Fewster–Smith, made
  explicit by Kontou–Olum), not by an extrapolation. The remaining assumptions are:
  - the first-order truncation, justified here because R t0² ≈ 10⁻⁴;
  - the free, minimally coupled field;
  - a Hadamard state;
  - Kontou–Olum's Gaussian-type sampling, whose flat coefficient 3.76/(16π³) ≈ 0.0076 is close to ER's Lorentzian
    3/(32π²) ≈ 0.0095.
- **Advance rows with small η (round 5).** τ0 = σ r_c grows as η^−1/2 and exceeds ε below η ≈ 10⁻⁴. There the
  derivative terms are O(1) or larger, so the first-order theorem does not apply. Round 5 already labelled those rows
  an extrapolation; this round confirms it quantitatively. Their wall bounds stay unproven, in either direction.
- **σ dependence.** The allowed thickness scales as ε_max ∝ σ⁻². σ = 0.01 is the conventional, generous choice.
  The table shows corrections of ~10⁻⁴ at σ = 0.01, so a larger σ (e.g. 0.1, corrections ~10⁻²) would also be
  defensible. That would tighten the bound to ε_max ≈ 10² l_P. The latency rows therefore read: **ε_max ≈ 10⁴ l_P
  at σ = 0.01, and smaller for any larger defensible σ.**

## 3. Krasnikov 2003 and later: does each loophole apply to this geometry?

| Loophole or limit | Source | Applies to our tube? | Why [derived unless marked] |
|---|---|---|---|
| The curved-space QI is only a conjecture | K2003; K2004 | **No longer, for latency rows** | §2: Fewster–Smith and Kontou–Olum cover a free, minimally coupled scalar in a globally hyperbolic spacetime at small t0; our latency rows meet the conditions. Advance rows: still open. |
| Box-sampled QI fails (2D de Sitter counter-example) | K2004 | No | Fewster (gr-qc/0409043) [lit]: box sampling is outside QEIs, and the smooth absolute QEI holds in that example. ER use smooth Lorentzian sampling. |
| Weyl-dominated regions decouple the density from the curvature | K2003 §III A | No | max\|C\| / max\|R_ab\| = 0.33 (η → 0), 0.30 (η = 1), 0.22 (η = 1.9); never ≫ 1. |
| The chain ρ ~ ℓ⁻⁴ ~ ρ² forces Planck density | K2003 | Yes, as a statement about this profile | Our ε_max gives \|T\| ≈ 1.7 × 10⁻⁹ Planck density (η ≈ 1, ε = 10⁴ l_P): sub-Planckian, but fixed by the QI, not by the chain. |
| Casimir / periodic identification | K2003; Ford–Roman 1996 | No | The tube is simply connected and has no boundary; ER's wall is a smooth field configuration. A Casimir-plate wall would be a different model. |
| "Total negative energy is meaningless" for sub-Planck layers | K2003 | Advance rows only | Latency rows have ε ≈ 10⁴ l_P and r_c ≈ 5 × 10³ l_P: formally semiclassical. |
| Portal / Van Den Broeck reduction of E_tot | K2003, for wormholes and bubbles | No | The tube needs its wall along the whole length D. Only ρ_max can shrink, and a message needs ρ_max ≳ the carrier wavelength. |
| Many fields: the bound becomes N·C/τ0⁴ | Roman 2004 | Yes, as a cost | ε_max grows only as N^1/2 at fixed η: \|T\| ∝ 1/ε² against N/τ0⁴ ∝ N/ε⁴. |
| Non-minimally coupled field: no state-independent QEI | Fewster–Osterbrink 2008 [lit] | Possibly | The flat C/τ0⁴ bound does not hold for ξ > 0. Only a state-dependent bound does, which makes large negative densities "require more energy" [lit]. No absolute QEI exists for ξ ≠ 0 [lit: Kontou–Sanders 2020 §3.2]. The tube with a ξ ≠ 0 source is neither excluded nor shown possible by QEIs. |
| Interacting fields evade QIs | Olum–Graham 2003 [lit] | Not shown to help | Their 2+1D static negative density is small, set by the Lagrangian, confined near a wall with positive energy, and the system obeys ANEC, which "is sufficient to rule out superluminal travel" [lit]. No 3+1D interacting state that sources a tube profile is known. Interacting QEIs remain an open problem [lit: Kontou–Sanders §6]. |
| Static ground state is the lowest constant density | Fewster–Roman 2005 §8 (static) [lit]; stationary curves: Kontou–Sanders §4.1, citing Ford–Roman 1995 and Fewster 2017 [lit] | Yes (free fields with a ground state; ∂t timelike since g_tt = −1) | See §5: this turns the question into one about the ground state of this background. |
| AANEC forbids compactly generated Cauchy horizons; with Olum 1998, warp drives | Graham–Olum 2007 via Kontou–Sanders §5.3 [lit] | Time-machine part: no (η < 1, no CTC). Shortcut part: not checked | §1 item 4: along the outgoing achronal null geodesics R_ll = 0, so they violate nothing. A violating achronal geodesic would have to cross the wall in the shortcut direction. Olum 1998 and Graham–Olum 2007 were not read in this round. |

## 4. The negative-energy estimate: model, size and observer

What the numbers assume [ASSUMPTION unless marked]:
- **Model:** ER's mid-tube profile (smooth tanh wall, thin-wall limit ρ_max = 10⁶ε), in classical GR.
- **Observer:** static geodesic observers along ∂t. The density is ER's T_t̂t̂ in their frame.
- **Size:** length D (the exterior distance), radius ρ_max from 1 m (human scale) through 1 mm (mm-wave) to 1 µm
  (an optical carrier, the smallest radius for which a message is plausible).
- **Integration:** the coordinate volume 2πρ_max dρ dx across the wall. ER's band estimate uses \|T_min\| · α ε with
  α = 0.01. This round integrates the negative part instead [derived], giving α_eff = 0.33.
  - Round 5's energy column used ER's α = 0.01, so it understates the negative energy at its stated width by ×33.
  - Weighting by the proper volume of the t-slices (√k) lowers the integral ×3 at η ≈ 1 (k is small where T is most
    negative). These slices are not orthogonal to the static observers, so neither choice is a conserved energy. The
    numbers below are order-of-magnitude.

| D | Latency | η | ε_max (σ = 0.01) | r_c | E_neg, ρ = 1 m | 1 mm | 1 µm |
|---|---|---|---|---|---|---|---|
| 20 ly | 1 s – 1 day | 1 − (2 × 10⁻⁹ – 1.4 × 10⁻⁴) | 1.02 × 10⁴ l_P | 4.9 × 10³ l_P | −2.8 × 10⁴⁴ M☉ | −2.8 × 10⁴¹ | −2.8 × 10³⁸ |
| 20 ly | 1 yr | 0.95 | 9.9 × 10³ l_P | 4.9 × 10³ l_P | −2.7 × 10⁴⁴ | −2.7 × 10⁴¹ | −2.7 × 10³⁸ |
| 100 ly | 1 s – 1 yr | ≥ 0.99 | 1.02 × 10⁴ l_P | 4.9 × 10³ l_P | −1.4 × 10⁴⁵ | −1.4 × 10⁴² | −1.4 × 10³⁹ |

- The wall also carries positive energy of about the same size (I_pos ≈ 0.89 \|I_neg\|; it exceeds \|I_neg\| with the
  √k weighting).
- E_neg ∝ ρ_max D / ε. A larger defensible σ, or N-field content (ε_max ∝ N^1/2, so E_neg ∝ N^−1/2),
  changes it by a few orders of magnitude for any plausible N. No choice in this table brings it within 20 orders
  of magnitude of the mass of the Galaxy (~10¹² M☉).

## 5. Restated scopes

- Round 5 said "**no semiclassical wall exists**" for the 1 s *advance* rows.
  - What the evidence supports: under ER's assumptions (free, minimally coupled massless scalar; flat QI at
    τ0 = σ r_c with σ = 0.01; static observer), the largest wall allowed for a 1 s advance is sub-Planckian (0.34 l_P
    at 20 ly). So within those assumptions no wall exists for which the semiclassical treatment that produced the
    bound is valid.
  - This round adds: for those rows τ0 ≫ ε, so even the bound itself is unproven (§2). The statement is therefore
    "**the semiclassical QI analysis gives no consistent wall for a 1 s advance; it neither proves nor excludes a
    wall outside its assumptions**" (non-minimal coupling, interacting fields, quantum gravity).
  - For the owner's latency target this phrase does not apply. A wall of ~10⁴ l_P is consistent with the QI and
    formally semiclassical.
- Round 5 said "**no self-consistent matter model exists**".
  - What the evidence supports: **none is known in the literature checked, and none was found here.**
    - ER give no source; they bound one.
    - Krasnikov 2003 shows only that specific QI-based prohibitions do not follow. He constructs no source.
    - Olum–Graham obey ANEC.
  - §6 adds a necessary condition, derived for free fields with a ground state, that the one-loop scale fails by a
    large factor. That is an estimate, not a theorem. The claim "no matter model can exist" is **not** established.

## 6. Is there a concrete field and state that self-consistently supports the channel?

**Answer: none is known.** For free fields the question reduces to one computable condition. It fails on our
order-of-magnitude estimate, but it has not been computed.

1. **Reduction** [lit + derived]. A channel in use is stationary: the required T_t̂t̂(ρ) is constant along each
   static geodesic for years. For a free field with a ground state on this stationary background, the long-time
   average of the energy density along a stationary curve is lowest in the ground state [lit: Kontou–Sanders §4.1;
   Fewster–Roman 2005 §8 for static]. Hence a necessary condition is
   **⟨T_t̂t̂⟩_ground(ρ) ≤ T_required(ρ) wherever T_required < 0**.
   - Excited states (coherent, squeezed, thermal) cannot do better on long-time average.
   - A pulsed source can do better only if it drops below the ground state for a while and rises above it at other
     times. The wall, though, has to be negative at all times.
   - So the source must be the background's own vacuum polarisation in the ground state, possibly plus classical
     positive matter elsewhere.
2. **Estimate** [derived; ASSUMPTION about the coefficient]. The vacuum polarisation of N free fields scales as
   N c R² with c ~ 1/(2880π²), the one-loop scale. Its sign and its O(1) to O(10) coefficient on this background are
   UNKNOWN.
   - Matching it to T_required = −0.175/ε² at η ≈ 1 gives ε ≈ 0.06 √N l_P.
   - At the QI-allowed ε ≈ 10⁴ l_P this needs N ≈ 2.7 × 10¹⁰ fields. With a coefficient 10× larger, it needs ~3 × 10⁹.
   - A single field supports only sub-Planck walls.
3. **What could change the answer.**
   - Non-minimal coupling: no state-independent QEI holds [lit: Fewster–Osterbrink]. Whether the stationary
     ground-state argument of step 1 covers ξ ≠ 0 was not checked.
   - Interacting fields: open.
   - Large curvature-squared counterterm coefficients a, b: they enter ⟨T⟩_ground directly, as ∝ ∇∇R terms. With
     ∇∇R ~ R/(ε/4)², they reach the required size if a, b ~ 10⁻³–10⁻² (ε/l_P)², i.e. ~10⁵–10⁶ at ε = 10⁴ l_P
     (compare Ford–Roman 1996's (l/l_P)² for wormholes [lit]). That is a term in the
     gravitational action (higher-derivative gravity), not a matter source, and its own consistency (ghosts) is
     separate. Its experimental bounds were not looked up in this round.

## 7. The smallest computable gap

**G1 (decisive for free fields).** Compute the renormalised ground-state ⟨T_t̂t̂⟩(ρ) of a free, minimally coupled
massless (or light) scalar on the stationary ER wall background at η ≈ 1, including its sign. The background is
globally hyperbolic, with ∂t timelike everywhere (§1). Compare it with T_required(ρ) = ER Eq. (14).
- **Inputs:** the metric; mode sums in the (t, x, φ) Killing directions plus a radial ODE across a wall of width ε;
  point-splitting or Hadamard subtraction; the a, b ambiguity kept explicit.
- **Cheaper first pass:** the local DeWitt–Schwinger (large-mass) or Anderson–Hiscock–Samuel-type approximation. It
  is local in curvature and gives the sign and coefficient in hours of offline compute. Its validity for massless
  fields across a wall of width ε must be checked against a mode sum.
- **Outcomes:**
  - The ground state is positive, or negative but ≥ N × T_required for every plausible N: no free-field source.
    The line stops for free fields.
  - It is negative and within the reach of a known field content: write down the self-consistent semiclassical
    problem.

**G2 (shortcut-specific energy condition).** Find the fastest null geodesics between points outside the tube that
use the shortcut. Check whether they are achronal, and integrate R_ab l^a l^b along them. A negative integral on an
achronal geodesic conflicts with AANEC. Kontou–Sanders §4.3 report the evidence for AANEC [lit]:
- proofs for the massless scalar in 2D globally hyperbolic spacetimes (Wald–Yurtsever);
- proofs for spacetimes with small curvature (Kontou–Olum);
- Flanagan–Wald: macroscopic violations, i.e. above the Planck scale transverse to the geodesic, are excluded at
  second order around Minkowski space.

The wall's transverse extent is macroscopic (2πρ_max), so a negative result would put the tube in direct conflict
with that evidence. Olum 1998 and Graham–Olum 2007 need to be read first.

Neither gap is a feasibility test. Passing G1 would give a candidate source, not a channel: G2, the tube's front,
end caps and build-out, and stability would still be untested.

## 8. Answers

- **Do Krasnikov 2003's counter-arguments apply to our geometry?** Mostly not.
  - The Weyl, periodic and portal loopholes do not apply (§3).
  - The "conjecture" objection has been met for free, minimally coupled fields at the owner's latency rows (§2).
  - The non-minimal and interacting loopholes remain formally open, but no known example uses them to source
    macroscopic, sustained negative density. The one explicit interacting example obeys ANEC.
- **Are the sampling-time and curvature approximations valid?** For "wait ≤ 1 yr" (η ≈ 1): yes, with corrections
  ≤ 10⁻⁴. For "1 s early" (η ~ 10⁻⁹): no. τ0 spans ~140 wall widths, so the round-5 advance rows are not established.
- **What do the negative-energy numbers assume?** ER's profile, static geodesic observers, coordinate volume, the
  exterior distance D, and radius 1 µm – 1 m. They give 10³⁸–10⁴⁵ M☉ for 20–100 ly. They are order-of-magnitude
  (×3 from the slicing; ×33 from ER's α guess, now replaced).
- **Is there a concrete field and state?** None known. For free fields the necessary condition is G1. Our estimate
  says it fails by ~10¹⁰ in field number at the QI-allowed thickness.
- **Inapplicability is not feasibility.** Showing that a prohibition does not apply (e.g. Weyl dominance) leaves the
  tube unexcluded. It does not make the tube possible.

**Decision (proposal for the owner).** Keep the Krasnikov tube as the reference trade-off, with the "wait ≤ 1 yr"
rows (§4) replacing the round-5 advance rows as the goal-relevant numbers. Before any further work on this line, run
G1 (offline, one script). Pause the line if G1 gives a non-negative or insufficient ground-state density.

## 9. Scope

- **Read in full this round:**
  - Krasnikov 2003;
  - Roman 2004 (gr-qc/0409090);
  - Krasnikov 2004 (gr-qc/0409007v3) and Fewster's reply (gr-qc/0409043);
  - Fewster–Roman 2005 (gr-qc/0507013, gr-qc/0510079);
  - Fewster–Smith (gr-qc/0702056);
  - Fewster–Osterbrink (0708.2450v2);
  - Olum–Graham (gr-qc/0205134);
  - Kontou–Sanders (2003.01815v2);
  - Pfenning–Ford (gr-qc/9710055);
  - Ford–Roman 1996 (gr-qc/9510071);
  - Kontou–Olum (1410.0665v2).

  A literature subagent extracted them, and the key quotations were checked against the PDFs. The search was not
  exhaustive. Citation lists of Krasnikov 2003 and ER were scanned on INSPIRE.
- **Not covered:**
  - full-tube geometry (front, end caps);
  - dynamical stability;
  - back-reaction beyond the necessary condition;
  - quantum-gravity effects;
  - non-scalar fields (spin ½ and 1 polarisation coefficients differ; same scaling).
- No new candidates were added, as the brief required.
