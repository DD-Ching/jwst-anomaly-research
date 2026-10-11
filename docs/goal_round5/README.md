# Goal round 5: which concrete candidates let a receiver read a freely chosen message before exterior light arrives?

Owner brief 2026-10-11. At most three not-yet-handled candidates that explicitly claim this effect are screened on
four points:
1. their action or self-consistent field equations, their global geometry, and how the ends connect to the users'
   space;
2. the time from the free choice to identification at the receiver, with the channel and the exterior path compared
   on **the same clocks**;
3. hidden reliance on ordinary communication, pre-set messages, non-local external couplings, or shortening only the
   traveller's own time;
4. the physical source, the range of validity, and what is known and unknown about stability.

One small calculation is then made for the most decisive unknown.

Labels: **[paper]** (full texts read 2026-10-11), **[derived]**, **[estimate]**, **[ASSUMPTION]**, **UNKNOWN**.
A candidate that fails here is not a statement that every mechanism is impossible.

Code: `scripts/goal_ftl_candidates.py`; tests: `tests/test_goal_ftl_candidates.py`; output:
`results/goal_round5/ftl_candidates.json` (about 5 s).

## 0. Reused conclusions: these are not recomputed

| Configuration | Conclusion already in the repository or in the benchmark | Where |
|---|---|---|
| Carried short wormhole | Bookkeeping under an assumed throat; one link has no loop; the premise is unsolved | round 1 (#152) §2 |
| Bancal and LC4/LC5 v-causal models | The theorem alone guarantees no low-latency channel; only the lab design LC4 v1 is kept; the line is paused | rounds 2–3 (#154, #155) |
| Beyond-Horndeski wormhole of 2212.05969 | The solution establishes no shortcut within one exterior; stability is not established | round 4 (#156) |
| Maldacena–Milekhin–Popov 4D wormholes (1807.04726) | "Long": going through takes longer than the exterior path | round 1 §1.4 [paper] |
| Gao–Jafferis–Wall (1608.05687) | Paths through it are chronal, so no faster travel over long distances; it also needs a boundary coupling between the two sides (a non-local external coupling under criterion 3) | round 1 §1.4 [paper] |
| **Benchmark: Fu, Grado-White & Marolf, CQG 36, 245018 (2019), 1908.03273v2** | Both mouths in the *same* asymptotic region, with distance d between them. A well-timed signal needs t_transit = d + logs, so t/d → 1 as d grows; MMP gives t/d > 2. The generalized second law forbids wormholes from being the fastest causal curves. In D ≥ 5, t ≥ d. **"Short transit" therefore means approaching the exterior light time, not beating it** (their §1, §4, footnote 8). | [paper] |
| Drummond–Hathrell photons in curved spacetime | The wavefront (high-frequency) velocity is exactly c; the low-frequency superluminal phase velocity does not carry a front (Hollowood & Shore 2008, 0707.2303) | [abstract] |
| Alcubierre bubble | The bubble's front is causally disconnected from its occupants, so it can neither be created nor steered from inside (Everett & Roman 1997 §1); NEC is violated | [paper] |

The benchmark fixes how claims are read. A short transit time, or a time compared against the coordinate length,
counts here only if it beats the **exterior light time between the same two events, on the same clocks**.

## 1. The three candidates

Each makes an explicit claim of the effect:
1. **Krasnikov tube.** Krasnikov, PRD 57, 4760 (1998), gr-qc/9511068; extended to 4D by Everett & Roman ("ER"),
   PRD 56, 2100 (1997), gr-qc/9702049. Claim: once a tube has been laid, a return trip takes arbitrarily little Earth
   time.
2. **Scharnhorst effect.** QED vacuum between perfect plates: Scharnhorst 1990, analysed by Liberati, Sonego & Visser
   ("LSV"), Ann. Phys. 298, 167 (2002), gr-qc/0107091. Claim: photons crossing the gap travel faster than c; LSV argue
   "it seems to be the actual signal speed that is modified".
3. **Lentz positive-energy hyper-fast soliton.** Lentz 2021, 2006.07125. Claim: superluminal solitons sourced by
   positive energy, from a plasma and electromagnetic fields.

## 2. Screening

| | Krasnikov tube | Scharnhorst effect | Lentz soliton |
|---|---|---|---|
| **1. Equations and geometry** | Geometry is specified first (ER Eq. 13); the stress tensor follows from Einstein's equations (ER Eq. 14, reproduced here). **No matter Lagrangian produces it.** Inside the tube space is flat with opened light cones (Fig. 2). Both ends are in the users' own Minkowski space. | QED via the Euler–Heisenberg action; flat space with plates; same space. The effective metric (LSV Eq. 3.4) applies to propagation **normal to the plates** only. | Einstein–Maxwell–plasma is claimed, but "it is now a matter of finding the right configuration" (§4 end). **No self-consistent source is given.** The soliton is eternal and moves uniformly in one space. |
| **2. Timing, same clocks** [derived] | Earth/star rest frame, distance D. The tube is laid by an outbound trip at ≤ c, so it exists only for t > x: **nothing is early before t = D**. Afterwards a message chosen at the star at t_s reaches Earth at t_s + D(1 − η), against t_s + D for light: **an advance of ηD (≤ 2D)**. Earth → star gets **no advance**. Two tubes for two-way use allow closed timelike curves (ER §4). | Across one gap a, the advance is a(1 − (1 + ξ)^−1/2)/c: **5 × 10⁻⁴⁷ s at a = 1 µm**, and only across the gap. Plates are opaque, so gaps cannot be chained. | A payload riding a superluminal soliton would beat light. But Lentz notes that superluminal speed forms horizons (§4), and no creation or loading is given. By Krasnikov's Prop. 1 (globally hyperbolic spacetimes), a traveller cannot hasten arrival without structure laid earlier. So the first free choice cannot beat light. |
| **3. Hidden reliance** | The infrastructure is pre-laid at ≤ c. The message itself is chosen freely afterwards, so this passes. Clocks are external, not the traveller's. | None. | Relies on a pre-existing soliton, so the deployment question is unanswered. |
| **4. Source, range of validity, stability** | NEC/WEC violated: the static energy density is negative on the inner wall (ER §5). Quantum-inequality bound: §3 below. Stability UNKNOWN. | Standard QED in the soft-photon regime ω ≪ m_e c²/ħ. The front velocity (ω → ∞) is not derived. LSV's Kramers–Kronig argument *suggests* it is enhanced, but Hollowood–Shore found Kramers–Kronig violated for the curved-space analogue. UNKNOWN. | Santiago, Schuster & Visser (2105.03079) [paper]: Lentz checked only Eulerian observers, which does not establish the WEC. "All physically reasonable warp drives will violate the null energy condition", and the null convergence condition even in modified gravity. |
| **Functional screen** | **Passes conditionally** (geometry level): a real advance for later messages, one way. | **Passes only formally.** The prediction holds inside a sub-µm gap, at a magnitude §3 shows to be unresolvable in its own regime. | **Fails**: deployment is not possible from free choice, and the positive-energy claim is not established. Not computed further. |

## 3. Calculations [derived]

### 3.1 Krasnikov tube: the advance against the wall thickness allowed by the quantum inequality

Only the Krasnikov tube scales with distance, so its most decisive unknown is the one checked: can its negative
energy be supplied for the owner's distances (20 and 100 ly) and some useful advance?

**Positive controls.**
- An independent sympy Einstein tensor of ER Eq. (13), with k = k(ρ), reproduces **ER Eq. (14) to 1.6 × 10⁻¹¹**.
- The end-cap value of ER Eq. (38) (δ = 10⁻³, ρ_max = 100ε) is reproduced to 0.1 %: −0.04203 against −0.04199.

**Small-η scaling.** η is the opening of the light cone; the advance is ηD.
- The most negative static energy density is T_min = −0.1228 η/ε².
- The largest orthonormal curvature component is R̂_max = 3.079 η/ε², so the curvature radius is
  r_c = ε/√(3.08 η) ≫ ε.

**Bound.** The Ford–Roman flat-space quantum inequality (ER Eq. 43; massless scalar, Lorentzian sampling) gives
|T| ≤ 3/(32π²τ0⁴).
- ER require the sampling time to be shorter than the curvature radius: τ0 = σ r_c with σ = 0.01 [ASSUMPTION].
  Applied consistently, this gives **ε ≤ √(0.73 η)/σ² Planck lengths**.
- So a *smaller* advance needs *thinner* walls. ER's remark that a "weak" tube (η ≈ 10⁻⁶⁶) could have centimetre
  walls used τ0 = ε, not σ r_c. With r_c ≫ ε that sampling time is far below the limit their own criterion allows,
  so it gives a much weaker bound than that criterion does.
- The looser prescription τ0 = σε is also listed below.

| D | Wanted advance | η | Thickest wall allowed (τ0 = σ r_c) | Looser (τ0 = σε) | Negative energy (1 m radius, at the thickest wall) |
|---|---|---|---|---|---|
| 20 ly | 1 s | 1.6 × 10⁻⁹ | **0.34 l_P** (5.5 × 10⁻³⁶ m) | 1.1 × 10⁻²⁷ m | −2.8 × 10³⁸ M☉ |
| 20 ly | 1 day | 1.4 × 10⁻⁴ | 100 l_P | 3.8 × 10⁻³⁰ m | −8.4 × 10⁴⁰ M☉ |
| 20 ly | 1 yr | 0.05 | 1.9 × 10³ l_P | 2.0 × 10⁻³¹ m | −1.6 × 10⁴² M☉ |
| 20 ly | 20 yr (= D/c) | 1 | 1.0 × 10⁴ l_P | 3.8 × 10⁻³² m | −8.5 × 10⁴² M☉ |
| 100 ly | 1 yr | 0.01 | 858 l_P | 4.5 × 10⁻³¹ m | −3.6 × 10⁴² M☉ |
| 100 ly | 100 yr | 1 | 1.0 × 10⁴ l_P | 3.8 × 10⁻³² m | −4.3 × 10⁴³ M☉ |

**What it settles** (under the flat-space QI, a free massless scalar and σ = 0.01; these are the ER assumptions):
- No useful advance over 20–100 ly comes with a wall thicker than about 10⁴ Planck lengths, or 10⁻²⁷ m with the
  looser sampling.
- The negative energy needed is of order 10³⁸–10⁴³ solar masses.
- Thinning the effect (small η) does not help; it makes the walls thinner still.

**What it leaves open.**
- Interacting or non-minimally coupled fields.
- Curved-space or non-flat versions of the quantum inequality.
- Sources that are not quantum fields in a Hadamard state.

### 3.2 Scharnhorst effect: is the advance resolvable inside its own regime?

The LSV Eq. (3.3) formula gives ξ = 11π²α²/(4050 (a m_e)⁴).
- That is 3.17 × 10⁻³² at a = 1 µm. LSV quote 4.36 × 10⁻³² beside it: a factor of 1.37 that we could not reconcile.
  The conclusion does not depend on it.
- The advance across one gap is Δt ≈ aξ/(2c).
- Detecting an advance Δt needs frequencies ω ≳ 1/Δt. ξ is derived only for ω ≪ m_e c²/ħ = 7.8 × 10²⁰ s⁻¹.
- The product Δt · m_e c²/ħ = (ξ/2)(a/ƛ_C) falls as a⁻³. Even at the very edge of the regime, a = ƛ_C, it is only
  **7 × 10⁻⁷**. At 1 µm it is 4 × 10⁻²⁶.
- So within the regime where the effect is derived, the advance is **never resolvable**. Whether a message arrives
  early depends on the front velocity, which no calculation found here supplies.

## 4. Answers

**Which candidates truly predict early reading of a freely chosen message?**
- **Krasnikov tube:** yes at the level of geometry, with conditions.
  - Only for messages sent *after* a tube has been laid at ≤ c.
  - Only in one direction, with an advance up to ηD.
  - Its source, under the stated QI assumptions, needs Planck-thin negative-energy walls and 10³⁸–10⁴³ M☉ of
    negative energy for the 20–100 ly task.
  - No self-consistent matter model exists.
- **Scharnhorst effect:** only formally. The predicted advance is confined to one sub-µm gap, and inside the regime
  where it is derived it is too small, by ≥ 10⁶, to be resolved in principle. The front velocity is UNKNOWN.
- **Lentz soliton:** no. It cannot be deployed from a free choice (horizons; Krasnikov Prop. 1), and the
  positive-energy claim is not established (Santiago–Schuster–Visser).

**Observable prediction.** None was found that current experiments could resolve. At µm gaps the Scharnhorst
advance is 10⁻⁴⁷ s.

**Next minimal verification step.** For the Krasnikov line: check whether any QI for interacting or non-minimally
coupled fields permits sustained negative energy density over macroscopic thickness. A literature check would
decide it; Kontou–Olum-type results are a starting point. If none is found, pause the line.

**Decision (proposal for the owner).**
- No candidate in this screen merits an experimental or engineering design.
- Keep the Krasnikov result as the reference trade-off.
- Do not extend the search unless a new candidate states (a) a self-consistent source and (b) a same-clock advance
  over the exterior light time.

## 5. Scope and hypotheses not covered

- **Searched:** the explicit claims named above and their direct critiques (arXiv, 2026-10-11). This was not an
  exhaustive search.
- **Not covered,** each a separate explicit-claim class:
  - Lorentz-violating gravity with superluminal modes, e.g. Einstein–aether or Hořava spin-0 and spin-1 modes. GW170817
    constrains only spin-2.
  - Rotating or Gödel-type global geometries and Gott-type cosmic-string pairs, which need infinite or cosmological
    structures.
  - Non-local quantum-gravity proposals.
  - Tachyonic fields.

  None is claimed impossible here.

## 6. Bookkeeping

- **Time.** 2026-10-11 06:50 → about 07:20 UTC, including the round-4 closing revision.
- **Compute.** About 5 s.
- **Data.** Six arXiv PDFs (1908.03273, gr-qc/9511068, gr-qc/9702049, gr-qc/0107091, 2006.07125, 2105.03079), kept in
  the scratchpad with sha256 recorded in SOURCES.md and not committed.
