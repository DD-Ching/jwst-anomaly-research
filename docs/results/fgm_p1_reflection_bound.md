# Can P1's own reflection into mouth B be bounded from below? (2026-10-31)

Owner-requested follow-up to `fgm_wavepacket_window.md` (#619, D-113), at most 90 min. Only derivation, stored
numbers and small integrals. One question: with #619's FGM parameter point and packet P1 unchanged, is there a
physically grounded lower bound on the energy that P1 sends into mouth B by reflection at mouth A, within the
relevant time, that decides whether P1 closes the window? Numbers come from `results/fgm_p1_reflection/bound.json`
(`scripts/fgm_p1_reflection.py`) and `results/fgm_window/window.json`. Labels as in #619: [lit], [der],
[assumption], [inference], [estimate]. **Decision rule (the owner's, given before any computation):**
- a lower bound above the closure threshold supports excluding P1 by this mechanism;
- an upper bound below it means P1 passes this necessary condition;
- bounds that straddle the threshold, or a two-mouth transfer that lacks a needed input, give C, with the gap named.

## 1. Reproduction

From `window.json`: the threshold at the reflection arrival of the packet centre (u_B − x*₀ = −795.1 r₊) is
δE_crit = 1.46 × 10⁻⁸⁷ ħc/r₊. P1's energy below ωr₊ = 1 is 7.69 × 10⁻⁹ ħc/r₊. **Critical fraction: 1.90 × 10⁻⁷⁹**
of the low-frequency energy must reach B to close the window. #619 assumed a fraction (r₊/d)² = 10⁻⁸, giving its
"fail by 5.3 × 10⁷⁰".

## 2. What P1 fixes, and what each stage needs

P1 fixes the spectrum |f(ω)|² (Gaussian, ω₀r₊ = 500, σr₊ = 100, cut at ω = 0), one quantum, entry along the
horizon generator θ = π of A (the antipode of the compact string; also π/2), and a 1+1 radial arrival density. It
says "a massless field" without a spin, and gives no angular-momentum content (ℓ, m), transverse profile,
polarisation, or time and place of preparation.

#619's f_geo = (r₊/d)² equals isotropic (ℓ = 0) re-emission times a capture cross-section equal to B's horizon area
4πr₊². That area is the low-frequency absorption cross-section of a minimally coupled massless scalar (Das, Gibbons &
Mathur 1997). So the frozen R6 used a 4D completion that P1 does not state, together with full reflection, which is
an upper-side input.

Geometry: θ = π at A faces away from B, so along P1's approach B lies on the A–B axis behind A. The background (a
static axisymmetric two-centre solution with the strings on the axis) is symmetric about that axis, so a test
field's azimuthal number m about the axis is conserved [inference]. For m ≠ 0 every such field has a node on the
axis, where B sits.

| Stage | Input needed | Fixed by P1 | Missing | Can it push the bound below the threshold? |
|---|---|---|---|---|
| Reflection at A | barrier reflection R_{s,ℓ}(ω) near RN; P1's amplitude per (ℓ, m) at ωr₊ < 1 | spectrum only | ℓ, m, spin | Not alone: at ωr₊ ≪ 1 every partial wave is almost fully reflected |
| A → B transfer | far-field amplitude of the reflected wave toward B; two-centre geometry | d, θ | angular pattern (from ℓ, m) | Yes: m ≠ 0 about the axis puts a node where B sits (§5) |
| Capture by B | σ_B(ω; spin, ℓ_B) for the field arriving at B; deposit pattern on B's horizon | nothing | spin, field structure at B | Spin alone, no (factors ≈ 10⁻² and 10⁻⁴ [estimate]); with m, yes |
| Time support | arrival-time distribution at B within the lifetime | the 1+1 density, including power-law tails from the cut at ω = 0 | time and place of preparation; ω-dependence of (ℓ, m) | Yes for the early tail (§6) |
| FGM weights | time weight e^{−κ₊(u−x*₀)}, transverse kernel H | yes (FGM) | — | No: at θ = π any deposit weighs ≥ 1 relative to the string (§3) |

None of the missing inputs makes the bound exactly zero: the transfer stays positive for any finite |m|. They can make
it as small as wanted.

## 3. Closure functional with FGM's time and angle weights [der]

Positive Killing energy falling into B with flux F_B(u, Ω′) closes the window on the generator θ_s if

  J = ∫ du e^{−κ₊(u − x*₀)} ∫ dΩ′ F_B(u, Ω′) H(γ)/H(θ_s) > cκ₊/16,

where γ is the separation between θ_s and Ω′, and the strings sit at θ′ = 0. The time weight is from FGM (2.20) and
(3.3); the transverse kernel is from (3.14)–(3.15). #619's R6 is the special case of a point deposit at the string's
pole, arriving at the centre time. The integration runs over the background's lifetime (from −t_merge = −10⁶ r₊).
It stops at the signal's earliest emergence label at B (1,490 r₊ after the centre arrival u_c). Later arrivals meet
the signal far from B's horizon and do not close it.

- **Angle.** H decreases with separation (checked on a 400-point grid), so at θ_s = π the string sits at the
  kernel's minimum. Any deposit weighs at least 1.0 relative to the string (0.73 at θ_s = π/2). A uniform deposit
  (ℓ = 0 absorption at B) weighs ⟨H⟩/H(θ_s) = 1/(4πκ₊H(θ_s)/r₊²) = 1.60 at π and 1.16 at π/2. An equatorial
  deposit weighs 1.38 at π.
- **Time.** e^{−κ₊u} is convex, so Jensen's inequality gives J ≥ E_B · a · e^{−κ₊(⟨u⟩ − x*₀)}. Here ⟨u⟩ is the
  energy-weighted mean arrival label at B, ⟨u⟩ = u_c + ⟨τ⟩. The mean delay ⟨τ⟩ collects the reflection group delay,
  the inter-mouth Shapiro delay and FGM's O(ln d/r₊) offsets.
- **Tail after emergence.** The tail arriving after emergence holds at most 5.8 × 10⁻⁶ of the band energy, so
  dropping it changes nothing.

## 4. Bound in an added setting S (conclusion limited to S)

**S** [assumption]: a minimally coupled massless scalar whose sub-barrier part is in the ℓ = 0 wave about A,
launched within r < d. This is the completion behind #619's f_geo.

| Input | Value | Label |
|---|---|---|
| reflection at A | \|R₀\|² ≥ 1 − 8(ωr₊)², i.e. ≥ 0.92 at the band top | assumption: leading low-frequency \|T₀\|² = ω²σ₀/π → 4(ωr₊)², with a factor-2 allowance |
| capture by B | ≥ (r₊/d)²/2 = 5 × 10⁻⁹ of the energy A re-emits | assumption: isotropic far-zone wave, DGM area, factor-2 allowance |
| band used | 0.02 ≤ ωr₊ ≤ 0.1 (7.2 × 10⁻¹¹ ħc/r₊) | choice |
| mean delay | ⟨τ⟩ ≤ 200 r₊ | assumption: the delay terms are tens of r₊ |

**Result:**
- E_B ≥ 3.4 × 10⁻¹⁹ ħc/r₊. The transfer fraction of the sub-barrier energy is ≥ 4.5 × 10⁻¹¹, against a critical
  1.9 × 10⁻⁷⁹.
- J / threshold ≥ 3.7 × 10⁶⁸ at the centre and ≥ 7.2 × 10⁴⁶ with ⟨τ⟩ = 200 r₊. The bound stays above the threshold
  for any mean delay below 632 r₊.
- The early tail (§6) only adds to J.

**In S, the mean-field lower bound exceeds the threshold: this supports excluding P1 in S.**

## 5. Completions consistent with P1 that fall below the threshold [estimate]

Take a vortex about the A–B axis: the sub-barrier part carries azimuthal number m, with a sectoral pattern Y_{|m|m}.
B, centred on the axis, can then absorb only ℓ_B ≥ |m| [der]. Assume Γ_ℓ ≤ (ωr₊)^{2ℓ+2} for ωr₊ ≤ 1 [assumption].
Then the transfer is T_m ≤ [(2|m|+1)!!]² (r₊/d)^{2|m|+2}/4. With the equatorial deposit weight 1.38, at the centre
time:

| \|m\| | 1 | 8 | 12 | 13 | 20 | 50 |
|---|---|---|---|---|---|---|
| log₁₀ (J / threshold), upper | 63.2 | 21.3 | 0.05 | −5.1 | −39.5 | −168.9 |

P1's main part still enters A for |m| < ω₀b_c ≈ 1574 (b_c = 3.149 r₊ for r₋ = 0.5 r₊). So the lower side over
completions consistent with P1 lies below the threshold for the prompt reflection. It falls further as |m| grows.

## 6. Early tail and the finite lifetime

**Tail shape.** Because the spectrum is cut at ω = 0 with p(0) = |f(0)|² = 1.49 × 10⁻⁸ r₊, P1's frozen arrival
density has power-law tails. The energy flux far from the centre is p(0)/(8|Δ|³); a direct integral in the test
checks this. With FGM's weight, reflected tail arriving at B a lead Δ before u_c contributes about
T p(0) e^{κ₊Δ}/(8κ₊Δ³) against δE_crit. This holds for Δ < d, in the far zone.

**Where it crosses the threshold:**
- **S:** above at any lead, since the prompt part already is. At Δ = d the ratio is 10¹¹⁴⁴. Extrapolated to the
  background's start (Δ = 1.0 × 10⁶ r₊) it is 10¹⁰⁹⁰⁸³. This is order of magnitude only: below ω = 1/d the mouths
  are in each other's near zone, which is not modelled.
- **Vortex:** the crossing lead is 104 r₊ for |m| = 13, 1.6 × 10³ r₊ for |m| = 50 and 6.2 × 10³ r₊ for |m| = 200.
  From |m| = 367 the crossing lead exceeds d, where the estimate does not apply.

**What follows.** For every fixed-|m| completion up to |m| = 200, the packet's own early tail closes the window in
the mean-field sense within the background's lifetime. Avoiding this needs one of two things P1 does not specify:
- the tail's low-frequency content carries a larger |m| than the prompt part (an ω-dependent completion);
- the state is prepared later than the crossing lead.

Any one-quantum packet has non-exponential tails (Paley–Wiener), while the threshold falls exponentially into the
past. The same mechanism therefore applies to every packet in the mean-field sense [inference], as for #619's P2.

**Erratum to #619 §3.5.** The narrow-band Gaussian's tails are not Gaussian: the cut at ω = 0 leaves
p(0) = 2.2 × 10⁻⁸⁹ r₊. In the S completion its early tail crosses the threshold at a lead of 154 r₊. #619's "passes
by ten orders" holds at the packet centre only. Not re-evaluated further here.

## 7. Expectation flux versus a single detection

- **Mean field (what §§4–6 compute).** J uses the expectation value of the stress tensor in the one-quantum state.
  The reflected branch enters with its probability as a weight, as in FGM's linear semiclassical treatment.
- **Single run.** The quantum either enters A or it does not. Energy that P1 itself puts into B, by reflection or
  through its early tail, exists only in runs in which the quantum did not enter A. In every run in which P1
  traverses, its own contribution to B is zero, so this mechanism cannot close the window on it [inference: the
  linear shift acts branch by branch]. The mechanism lowers P1's success probability by at most the probability of
  not entering A: ≤ 1.5 × 10⁻⁸ for the sub-barrier part in S. In the closing runs (probability ≤ 3.0 × 10⁻¹⁶ in S),
  one quantum exceeds the threshold by ≥ 10⁸⁵, but that quantum is not the traversing signal.
- **Validity of the mean-field number.** Semiclassical gravity is a good approximation only when stress-tensor
  fluctuations are small (Kuo & Ford 1993). The relative fluctuation of the energy P1 puts into B is about
  p^{−1/2} ≈ 5.7 × 10⁷ in S. FGM's negative energy averages c = 20 string fields; the single quantum's branches are
  not averaged.
- **Energy present in every run** (apparatus, environment) is a different case. In a single run any one quantum
  entering B before the critical time exceeds the threshold, so the single-run condition is that no quantum enters
  B. That is stricter than the mean-energy condition. It is not evaluated here, as in #619.

## 8. Verdict under the owner's rule

- **Mean field, P1 as frozen: C.**
  - The lower side over consistent completions is below the threshold: the vortex with |m| ≥ 13 for the prompt part,
    and ω-dependent completions or later preparation for the tail.
  - The upper side is far above: in S, ≥ 7.2 × 10⁴⁶ (prompt) and 10¹¹⁴⁴ (tail at lead d). For fixed-|m|
    completions up to |m| = 200, the early tail alone exceeds the threshold within the lifetime.
  - **Gap:** the two-mouth transfer lacks three inputs: (i) P1's angular content about the A–B axis (m and ℓ as
    functions of ω); (ii) the field's spin; (iii) when and where the state is prepared, which sets the support of
    its early tail.
- **Mean field, added setting S:** the lower bound exceeds the threshold by ≥ 7.2 × 10⁴⁶. This supports excluding P1
  in S.
- **Single run:** P1 passes this necessary condition for every completion. Its own reflection cannot block the run
  in which it traverses (§7, inference).
- **P2 (labels unchanged, not re-evaluated):**
  - frozen rule: A (centre only);
  - post hoc: B (past tail; R ∼ e^{−πω/κ₊}; f_geo; mean field). The two caveats above (completion; mean field versus
    single run) bear on this post-hoc B.
- **#619 as a whole:** the post-hoc "B for the frozen packets" needed P1's failure. With P1 at C (mean field) and
  passing in single runs, the pair's post-hoc reading is C. The frozen-rule label (A through P2) is unchanged.

**Deviations and notes:**
- The setting S, the vortex family and the ω-band are this round's choices (labelled above).
- The early-tail analysis implements the owner's lifetime integral.
- The §3.5 erratum is a correction found in passing.
- Nothing here is a channel, a bit, a rate or a faster-than-light claim.

## 9. Plain summary

P1's verdict changed. #619 reported a failure by 5.3 × 10⁷⁰. That number assumed that P1's low-frequency part is
re-radiated evenly in all directions and absorbed by B over its whole horizon area. P1 never specified that, and
the number compares averages over the quantum state. The deciding quantity is the share of P1's low-frequency
energy that reaches B, compared with 1.9 × 10⁻⁷⁹, inside the inequality J ≥ E_B · a · e^{−κ₊(⟨u⟩−x*₀)} > cκ₊/16.
The share is at least 4.5 × 10⁻¹¹ for the spherical completion, but it can drop below 1.9 × 10⁻⁷⁹ if the field
circulates around the line joining the mouths (|m| ≥ 13). B sits on that line, where such a field vanishes. The
packet's power-law early tail again favours closure for simple completions. So the mean-field answer depends on
inputs P1 does not give (C). For one quantum, being reflected into B and traversing are exclusive outcomes. P1's
own reflection therefore cannot block the run in which P1 crosses.
