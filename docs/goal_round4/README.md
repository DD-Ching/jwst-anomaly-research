# Goal round 4: does the beyond-Horndeski wormhole give a stable, usable channel?

Owner brief 2026-10-11. The question: can a concrete modified-gravity wormhole model provide a stable channel
that carries messages and gives both ends a time advantage? (Closing revision: see §4 and §8.)

Labels used below:
- **[paper]**: arXiv:2212.05969v1, Mironov, Rubakov & Volkova, PRD 107, 104061 (2023), and arXiv:2404.06297v2,
  Mironov & Volkova (v2 of 2025-01-17, IJMPA ICNFP2022 issue). Both were read in full; the formulas were read from
  page images.
- **[derived]**: this round.
- **[abstract only]**, **[ASSUMPTION]**, **UNKNOWN**.

Nothing here is a new-physics result. Reproducing the papers is not a discovery, and a failure of this one model
says nothing about other mechanisms.

Code: `scripts/goal_bh_wormhole.py`. Tests: `tests/test_goal_bh_wormhole.py`. Output:
`results/goal_round4/bh_wormhole.json` (about 37 s).

## 0. What the project already had, and what this round adds

| Earlier work | What it covered | Reused here |
|---|---|---|
| Goal round 1 (#152) §1.4, §2 | Energy-condition theorems; carried-wormhole bookkeeping under an *assumed* short rigid throat; "candidate C: a modified-gravity escape from Gao–Wald" left open | The bookkeeping and the theorem table. The open premise is the one tested here. |
| D-047, W2/W3 (Ellis lenses), docs/exotic_limits.md | Observational searches for Ellis-type (1/b² deflection) and negative-mass lenses | Only to check whether they apply to this model (§5). |
| Hypothesis rounds A/B (D-069, D-075, D-076), D1 (D-073), E1 (D-074) | New-physics hypotheses, distance self-consistency, event networks | None treats a modified-gravity wormhole, so nothing is reused. |

**The unsolved assumption this round tests.** Round 1 assumed a short traversable throat can be held open without
NEC-violating matter. The question here is whether one specific modified-gravity solution does that *stably*.

## 1. The fixed model [paper]

**Theory.** Quadratic beyond Horndeski, action (1) of 2212.05969, with the functions F(π, X), G4(π, X) and
F4(π, X).

**Solution** (Sec. 4, Eqs. 70, 72, 73, 80, 81 and 83):
- metric: A = 1, B = 1 + sech(r/τ), J = τ ln(1 + 2 cosh(r/τ)), for r ∈ (−∞, ∞);
- scalar: π = tanh(r/τ) − 1;
- parameters: τ = 100, w = 2, u = 0.

**Lagrangian.** It is reconstructed from the solution ("reverse engineering"):
- H = G = 1, Ξ = 0 and Γ = 0;
- H − 2F4B²π′⁴ = 1 − 2 sech(r/τ);
- c²_s2 = cosh/(1 + cosh).

The functions are given only as plots and asymptotics (Eqs. 84–86).

**Units, scaling and transmission are three separate questions.**
1. **Unit convention [paper + ASSUMPTION].**
   - The caption of 2212.05969 Fig. 1 says the choice u = 0, w = 2, τ = 100 "guarantees that the size of the
     wormhole throat safely exceeds the Planck length". That reads τ in Planck units.
   - The action fixes only 8πG = c = 1, so taking ħ = 1 as well is our reading of that caption. On it, the throat
     radius J(0) = τ ln 3 = 110 l_Pl ≈ 9 × 10⁻³³ m.
   - Every number below in metres, seconds or GeV rests on this convention.
2. **Scaling [derived].** The classical stability conditions are homogeneous in τ (§2, `scale_check`). A different τ
   is a different reconstructed Lagrangian with the same stability verdict.
3. **Wave transmission [not computed].** How much of a signal of a given frequency passes the throat was not
   calculated. Only an order-of-magnitude remark is given in §4.

**Consequence [derived].** Because Ξ = Γ = 0 and H = G = 1, every quantity in the 2404.06297 stability set is
fixed by the metric. That includes F, which follows from their Eq. (71): F = J²[B′J′/(2J) + B(J′²/J² + J″/J)].
So the complete set can be checked without the unpublished Lagrangian functions.

**Literature status** (INSPIRE citations and an arXiv search, 2026-10-11):
- **No erratum or second version of 2212.05969 exists.**
- **2404.06297 v2, by the same authors, states that this solution "does not satisfy the complete set of stability
  conditions and, hence, cannot claim to be fully stable"** (Conclusion). It does not say which condition fails.
- Of the five papers citing 2212.05969, none constructs a fully stable version. These include no-go and
  speed-of-gravity papers by the same group: 2411.05416, 2408.01480 [abstract only].
- The 2609.00664 obstruction concerns black holes with static scalar hair, not wormholes [abstract only].

## 2. Step by step: what holds, what fails, what is unknown

| Step | Status | Evidence |
|---|---|---|
| **Solution exists** | **Holds, as published** (not re-derived here) | Reconstruction in 2212.05969 §4.1–4.2. The background equations were not re-checked: the Lagrangian functions are not tabulated. |
| No ghosts, no radial or angular gradient instabilities (high ℓ) | **Holds** [paper, reproduced] | F_min = 1.0000001 > 1. K11, det K, G11, det G, M11^(ℓ²), det M^(ℓ²) > 0 on the grid, as in their Fig. 2. K11 ∝ Θ² and M11^(ℓ²) ∝ h² have isolated zeros, at Θ = 0 and at cosh(r/τ) = 2, which the grid misses. 2P1 − F > 0 and N > 0. |
| Odd-parity tachyons, condition (6d) | **Holds** [derived] | V + 6H/J² > 0 at every grid point for ℓ = 2–50. |
| Even-parity tachyons and low-ℓ angular modes, (39b) and (39c) | **Sufficient condition fails** [derived, new]. Whether a real tachyon exists is **UNKNOWN**. | (39b) M22 > 0 fails for ℓ = 2 at \|r\| < 0.455τ. (39c) fails near the throat or around the Θ = 0 crossing for every ℓ computed. It also fails in a far zone for every ℓ tested: beyond 6.3τ (ℓ = 2), 9.6τ (ℓ = 20), 12τ (ℓ = 50) and 16τ (ℓ = 200; `far_zone_probe`). Both kinds of failure may be artefacts (see §3). The set is sufficient, not necessary (2404.06297 §1, §4). |
| (Sub)luminality: a propagation property, not a stability criterion | **Not satisfied** [derived, quantitative; agrees with 2212.05969 §4.2]. Superluminal propagation is not an instability. It bears on causal structure and on whether a standard UV completion exists (EFT consistency, next row). | One even-parity mode has angular speed² up to **1.83** at \|r\| ≈ 1.4τ. It stays above 1 at every grid point and tends to 1 far away. The other mode lies in 0.16–0.82. Radial speeds² are 0.68–1 (tensor) and 0.50–1 (scalar). |
| Nonlinear stability | **UNKNOWN** | Not addressed by any source found. |
| Validity of the effective theory: consistency, separate from stability | **UNKNOWN**, with three warning signs: superluminal modes (in tension with positivity arguments for a standard UV completion; not evaluated here), and the two below | If the units are Planck units, the throat is 110 Planck lengths. The reconstructed functions grow as (−π)^(−7/2) and (−π)^(−11/2) as π → 0 at large r (Eq. 86), so the couplings of perturbations there are not shown to be weak. |
| Formation or control | **UNKNOWN** | No mechanism is given. The solution is static and eternal. |

**Scale.** The conditions are homogeneous in τ. A full rerun at τ = 10⁴ on the same grid (`scale_check`) gives the same angular speeds to 10⁻¹⁵ and identical failure segments for every ℓ.
So making the throat macroscopic does not cure the failures. That would in any case be a different reconstructed
Lagrangian.

## 3. The new calculation [derived; positive controls first]

**Unsolved question chosen.** Does the published solution satisfy the *complete* stability set of 2404.06297 that
followed it, and if not, which part fails? This decides whether the model can even be called linearly stable,
which is a prerequisite for any channel.

**Method.**
- The model's choices reduce every coefficient to the metric (§1).
- sympy builds Eqs. (5), (22), (24), (27), (33)–(37), (41)–(43) and (129)–(133) symbolically.
- numpy evaluates them on 4,000 points in −10 ≤ r/τ ≤ 10, avoiding the point r = 0, where J′ = 0.
- A 40-digit mpmath evaluation at four points agrees with float64 to ≤ 2 × 10⁻¹² relative (`mpmath_crosscheck`). Complex eigenvalues and non-finite points are counted; both are 0.

**Positive controls.**
- F > 1, matching their statement below Eq. 72.
- All six Fig. 2 quantities are positive (with the isolated zeros noted in §2).
- The ℓ → ∞ limit of the full-ℓ Eq. (33) matches Appendix B Eqs. (132) and (133) to 7 × 10⁻⁷ and 2 × 10⁻⁵
  (relative) at ℓ = 10⁴. M11 matches Eq. (41) except near its zeros.
- The tensor-like angular eigenvalue tends to 1 in the flat far zone, which confirms the normalisation of the speeds.

**Inconsistencies found in 2404.06297.**
- det K from Appendix B (129)–(130) is exactly 4 × Eq. (22) (`detK_appendixB_over_eq22_range`).
- det M^(ℓ²) from (41), (132) and (133) is 4F(r) × Eq. (42): the ratio runs from 4.0 far away to 5.9 at the
  throat (`detM_appendixB_over_eq42_range`). The independent review showed it analytically: the prefactor 1/(4J²F²)
  of Eq. (42) should read 1/(J²F). The sign, and so condition (44), is unaffected.
- Eqs. (132) and (133) print H − 2F4B²π′², where every other formula has π′⁴. This is probably a typo; the script
  uses π′⁴.

The speeds use the Appendix B components, because those agree with Eq. (33) and give a luminal tensor mode far
away. This choice changes no conclusion: Eqs. (22) and (42) also give a positive det K and det M^(ℓ²).

**Results.**
1. **Superluminality is quantified, consistent with 2212.05969 §4.2.** The maximum angular speed² is 1.83, at |r| ≈ 1.4τ.
2. **The low-ℓ sufficient conditions fail.** This is new: no earlier source evaluated them on this solution.
   - (39b) fails for ℓ = 2 around the throat.
   - Inner zone (|r| < 3τ): (39c) fails on 0.38τ < |r| < 1.79τ for ℓ = 2 (47 % of the zone) and on |r| < 1.6τ for
     ℓ = 3–4. For larger ℓ the failing band shrinks to 1.28τ < |r| < 1.35τ at ℓ = 50, around the Θ = 0 crossing
     (2 %).
   - (39c) also fails in a far zone: from 6.3τ outward for ℓ = 2, from 8.2τ for ℓ = 10 and from 9.6τ for ℓ = 20, all
     out to the grid edge. Beyond the grid, an mpmath probe finds it from 12τ for ℓ = 50 and 16τ for ℓ = 200, so every
     ℓ tested fails far enough out (`far_zone_probe`; independent review).
   - That zone is nearly flat (B − 1 ≲ 4 × 10⁻³). What drives the failure is 2P1 − F → 0 exponentially (det G ∝ N
     → 0) while Q12 stays finite, so the scalar's kinetic and gradient weights vanish asymptotically but the cross term
     does not. This is plausibly an artefact of that degenerate asymptotic scalar combined with the criterion treating
     v and v′ as independent. It is not proven: **UNKNOWN.**
   - The inner failures near Θ = 0 may also be a basis artefact, since K22, M22 and Q12 all carry 1/Θ. 2212.05969
     argues that the Θ = 0 crossing is not pathological.
3. **What it settles.** Three separate findings:
   - **Stability:** *not established*. The only published complete sufficient criterion fails at low ℓ. **It does
     not show that a tachyon exists**; only the low-ℓ eigenvalue problem can decide that (§6).
   - **Propagation:** one even-parity sector is superluminal, with c² up to 1.83. That is a causal-structure
     property, not an instability.
   - **Effective-theory consistency:** open (§2).

## 4. Clocks, signal time and the channel [derived]

**The two ends.** The ends are the two asymptotically flat regions r → +∞ and r → −∞. The solution connects two
separate asymptotic regions, so **it does not establish a shortcut comparison within one exterior space**: there is
no exterior path between the ends to compare against. This is a limit of what the solution provides, not a
demonstration that a shortcut is impossible in this theory.
- A shortcut inside one universe would need the two far regions glued into one exterior. No such solution of this
  theory was found in the literature searched; it is UNKNOWN.

**Clock comparison.**
- A = 1 everywhere, so static clocks at both ends tick at the Killing-time rate. Their ratio is exactly 1.
- They can be synchronised by Einstein exchange through the throat.
- There is no relative time shift between the ends. **That is not a failure.** Equal clock rates and no shift are
  what a usable two-way channel wants: neither end ages faster. Round 1's carried-mouth shift needs one mouth to
  move, which is a non-static configuration outside this solution.

**Signal time.**
- Radial light from r = −R to r = +R takes Δt = 2R − τI, with I = ∫(1 − 1/√(1 + sech x))dx = 1.1066, in the
  limit R ≫ τ.
- **The 1.1066τ is a difference against the coordinate length 2R. It is not an advantage over any exterior path,**
  because no exterior path exists in this solution. The light time equals the proper radial distance, since A = 1.
- The one-way time is about (proper distance)/c, and a round trip twice that.

**Deployment, continued use and round trips.**
- Deployment, i.e. creating a throat or moving its ends into place, is not described: UNKNOWN.
- Continued use needs the throat to survive finite-amplitude signals. That requires nonlinear stability: UNKNOWN.

**Message transmission [not computed; order of magnitude only].** Waves much longer than the throat are expected
to be mostly reflected, but the transmission coefficient was not calculated. On the Planck-unit reading (item 1 of
§1), a carrier that fits the throat would have a quantum energy of order 10¹⁷ GeV. For a larger-τ member of the
family the scale objection disappears classically, while the stability verdict is unchanged (item 2).

**Answer to the brief's channel question.** This solution does not establish a shortcut channel in one exterior
space:
- it links two separate asymptotic regions, with no exterior path to compare against;
- equal clock rates and no time shift are neutral, not a failure;
- its linear stability is not established (a separate question).

Superluminal propagation and EFT consistency are recorded separately (§2).

## 5. Observational limits: what applies to this model, background and parameters

- **Lensing.** The ADM mass is zero and B − 1 = sech(r/τ) decays exponentially. Deflection is exponentially small for
  impact parameters ≫ τ.
  - The project's Ellis-lens limits (W2/W3, D-047) assume a 1/b² deflection law. They **cannot be applied directly**:
    the model's deflection law and abundance would first have to be mapped onto those searches.
  - At τ = 100 on the Planck-unit reading, no astronomical observation is known to apply directly.
- **GW170817 (c_T = c).** It cannot be applied directly. That limit constrains tensor speeds on the *cosmological*
  scalar background of a theory, so applying it needs a cosmological completion of this Lagrangian and the map
  between the two backgrounds.
  - Here the Lagrangian is reconstructed only along this solution (π ∈ (−2, 0)), so its cosmological background is
    undefined.
  - The radial tensor speed² on *this* background, G/F = 1/F, is below 1 near the throat (down to 0.68). That is
    not excluded by GW170817, which says nothing about a local static background.
  - Whether a cosmological completion exists, and lies in the c_T = c class, is **UNKNOWN**. 2408.01480 [abstract
    only] relates the two for time-dependent backgrounds.
- **Lab and solar-system tests.** The modification is confined to r ≲ few τ around an object of size τ (Planck-size if the units are Planck units). None
  is known to apply directly.

Every item above is "cannot be applied directly until the model is mapped onto the background the limit refers
to". None says the theory is unconstrained by observation.

**No prediction found** that an observation or experiment could distinguish at the published parameters.

## 6. Answers

**Where it holds.** The solution exists, as published. It is free of ghosts and of high-ℓ gradient
instabilities, and of odd-parity tachyons.

**Where it fails or is unknown.**
- **Stability:** it fails the complete *sufficient* set at low ℓ (39b, 39c). Real tachyons and nonlinear stability
  are UNKNOWN.
- **Propagation:** superluminal in one sector (c_a² up to 1.83). This is not an instability.
- **EFT, formation, control:** UNKNOWN.
- **Channel:** no shortcut comparison is established within one exterior space (separate asymptotic regions). Equal
  clock rates and no time shift are neutral.

**What the new calculation settled.**
- The authors' statement that the solution fails their complete set is consistent with what we find. They do not say
  which condition fails, so we cannot know it is the same one.
- The failing conditions are located (low-ℓ even parity near the throat, superluminal angular modes) and quantified.
- Two internal inconsistencies of 2404.06297 are flagged.

**Distinguishable prediction.** None at the published parameters.

**Next minimal verification step.** Solve the ℓ = 2 even-parity eigenvalue problem, ω² from K⁻¹(G, Q, M), in a
formulation that is regular at the Θ = 0 crossing. That decides whether a real tachyon exists:
- a negative ω² would make the model linearly unstable (a FAIL);
- none would leave superluminality and the missing nonlinear and EFT analysis as the remaining obstacles.

**Decision (owner, 2026-10-11).** The communication use of *this particular solution* is paused. The results are kept.
This is not a conclusion about wormholes or modified gravity in general. An intra-universe solution of this or another
theory would be a new candidate.

## 7. Not done

- The low-ℓ eigenmode calculation (§6).
- Re-deriving the background equations, which needs the reconstructed functions.
- The 2018 predecessor (1812.07022). 2404.06297 also says it fails the complete set.
- Wave transmission through the throat (only an estimate is given above).

## 8. Bookkeeping

- **Closing revision (owner review 2026-10-11).**
  - "No time advantage" is narrowed to "no shortcut comparison is established within one exterior space".
  - Stability, superluminality and EFT consistency are now listed separately.
  - The Fig. 1 Planck-scale remark is cited, and units, scaling and transmission are separated.
  - Observational limits now read "cannot be applied directly".
  - The 1.1066τ difference is described as a coordinate difference, not an advantage.
- **Independent review** (separate agent, its own sympy/mpmath code built from the page images).
  - Confirmed: the reduction (Γ1 = Γ = 0 follows from A = 1), every transcribed equation, the positive controls, the
    speeds (max 1.832 at |r| = 1.40τ) and the low-ℓ results.
  - Fixed from it: the Eq. (42) factor (4F), the far-zone failure for every ℓ, the isolated zeros in Fig. 2, the
    unit caveat, and "consistent with" in place of "confirms".
  - Also: 2212.05969 Eq. (87)'s c²_a1·c²_a2 > 1 holds only near |r| ≈ 2τ. One eigenvalue exceeds 1 throughout.
- **Time.** 2026-10-11 06:15 → about 06:50 UTC.
- **Compute.** About 37 s for the script, including the τ = 10⁴ rerun and the mpmath cross-check.
- **Data.** Only the two papers (arXiv PDFs, sha256 `2a8ea802…` and `22eb223e…`, kept in the scratchpad and not
  committed), the INSPIRE citation lists and the arXiv search results.
