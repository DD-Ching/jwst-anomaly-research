# Can A's choice reach B through the FGM window? One-way bit, necessary conditions (2026-11-02)

Owner-requested analytic round, at most 90 min; derivation and small integrals only. This is a new model. It
does not reuse #619's packets P1 and P2, and #619's single-quantum branch argument is not an input. The FGM
parameter point is #619's (`results/fgm_window/window.json`). The closure functional is the one with FGM's time
and angle weights (`fgm_p1_reflection_bound.md` §3). Numbers come from `results/fgm_oneway_bit/bound.json`
(`scripts/fgm_oneway_bit.py`). Labels: [lit], [der], [assumption], [estimate]; "upper bound" and "lower bound"
are marked where a coefficient is bounded.

**Question:** at a pre-agreed time A chooses 0 or 1. Can B identify the choice inside a set read window with error
≤ 10 %, while the energy limit that keeps the window open is respected?

## 1. Model

**Background.** #619's point:
- r₊ = 1 km and r₋ = 0.5 r₊, so κ₊r₊ = 0.25;
- d = 10⁴ r₊; N = 10 strings (c = 20); Gμ = 10⁻⁷;
- FGM (3.3) normalisation; the signal travels along the generator θ = π; the window edge is v_max = −5695.1 r₊.

**Source S★ (fixed).**

| Item | Choice | Label |
|---|---|---|
| Field | minimally coupled massless real scalar test field, linear coupling to a classical current J; nothing else couples to it | assumption |
| Position | thin disc aperture, radius a = 5 r₊, on the A–B axis, L = 10 r₊ from A's centre on the side away from B; static, held by a support that never crosses either horizon | design |
| Angular distribution | amplitude bump(ρ/a), bump(x) = e^{−1/(1−x²)}, with a phase that focuses on A's centre (paraxial aperture a/L = 0.5); the beam enters A around θ = π | design |
| Time envelope | J = s · bump(ρ/a) · bump((t − t_c)/T) · cos ω₀(t − t_c), ω₀r₊ = 500, T = 20 r₊, compact support | design |
| The choice | s = ±1 is set at t_c − T; before that the two choices are identical (J = 0) | design |
| Timing | envelope centre enters A at v_c = v_max − 100 r₊ (#619's placement); the whole envelope enters 80–120 r₊ before the edge | design |
| Energy source | a local battery drives the aperture and supplies the radiated energy N_A ħω₀ (128 ħc/r₊ at N_A = 0.2554) | design |

**States and arrival** [der]:
- A classical current acting linearly displaces the background state. The two choices give coherent states |±α⟩,
  with α linear in J.
- The s-dependent field 2φ_ret[J] vanishes outside the causal future of the source (retarded response).
- **External path:** the earliest s-dependent field at B is the leak around A. From the envelope centre it arrives at
  u − x*₀ = −795.1 r₊ (light time d from A's plane; Shapiro terms of O(ln d/r₊) not counted).
- **Wormhole path:** the field emerges from B at labels ≥ −v_max = +5695.1 r₊, 1,490 r₊ after the external leak.
  The wormhole path is the slower one here, as in #619.

**Receiver and read window** [assumption]:
- The receiver sits on the axis beyond B, at ≥ 10 d (FGM fn 7 isolation; #619 R7).
- The read window is the emergence interval of the envelope, 8.2 × 10⁻⁹ r₊ long at B.
- The receiver mode is the image of the source mode under the emergence map: an ideal, mode-matched receiver.

**Phase reference:**
- **External path:** the static Killing time and the pre-agreed schedule give it.
- **Wormhole path:** at a fixed emergence label the received phase depends on the shift. The sensitivity is
  ∂φ/∂ln|ΔV| = ω₀κ₊⁻¹e^{κ₊·100 r₊} = 1.4 × 10¹⁴ rad [der]. A 0.1 rad reference therefore needs
  δln|ΔV| ≤ 6.9 × 10⁻¹⁶.
  - FGM fix |ΔV| only up to the O(ln d/r₊) term in x*₀. That term moves ln|ΔV| by O(κ₊ ln(d/r₊)) ≈ 2.3 [estimate].
    FGM give no fluctuations.
  - The reference must therefore be measured through the same window (a pilot), with the shift unchanged between
    pilot and signal at that precision [assumption, unverified].
  - Without a reference, the two signs average to the same state and the error is 0.5 at any strength [der].

## 2. Ideal benchmark and how FGM's background differs

**Benchmark** [der]:
- Fixed background, no noise, equal priors, Helstrom measurement: P_e = ½(1 − √(1 − e^{−4N})) ≤ 0.1 requires
  **N ≥ 0.2554** photons in the receiver mode (reproduced).
- With a pilot: a reference pulse plus a signed pulse of unknown common phase, interfered and counted. Then
  P_e = ½e^{−2n}, so n ≥ 0.805 per pulse, 1.61 in total.

| Item | FGM background | Effect here |
|---|---|---|
| Thermal occupation (Hartle–Hawking, T_H = κ₊/2π) | n̄ = 1/(e^{2πω/κ₊} − 1), log₁₀ n̄ = −5458 at ω₀ | negligible; it would matter only near ωr₊ ≈ 0.04 (n̄ = 0.58 raises the requirement to ≥ 0.55) [der] |
| Throat map | a translation in Kruskal V, exponential in Killing time | positive/negative-frequency mixing e^{−2πω₀/κ₊} (log₁₀ −5458): coherent states stay coherent, and B sees the Hartle–Hawking state again [der in 1+1; inference in 4D] |
| Blueshift | emerging carrier 2.4 × 10¹¹ – 5.3 × 10¹⁵ /r₊ (chirped), 8.2 × 10⁻⁹ r₊ long | the receiver must match it [assumption] |
| Shift fluctuations | not given by FGM | the phase reference (§1) |
| Closure | mean field: ⟨T⟩ of a coherent state is the classical stress of the displacement plus the background value | used as is [der] |

**Two paths, kept in separate accounts:**
- **External path:** the leak around A keeps the source carrier (500/r₊). It reaches B 1,490 r₊ before the
  wormhole emergence, arriving from A's side.
- **Wormhole path:** emerges radially from B along the θ = π generator, at ≥ 2.4 × 10¹¹/r₊.
- **Separation:** the read mode differs from the external field in frequency, time and direction, so the read counts
  only the wormhole path [estimate].

## 3. Two bounds on the same source

**Read-out (minimum strength):** N_A η_WH ≥ N_th. The coefficients of η_WH:

| Factor | Value | Label |
|---|---|---|
| A captures the beam (impact parameter b < b_c = 3.149 r₊) | 1 − h, h = 1.2 × 10⁻²⁶ | computed (paraxial focal field) |
| Envelope enters before the edge | 1 (compact envelope, 80–120 r₊ before v_max; curved-space tails of the retarded field neglected) | der |
| Transmission through the throat | ≈ 1 (geometric optics; mixing e^{−2πω₀/κ₊}) | assumption |
| Exit through B's barrier; mode match | ≈ 1 | assumption |
| η_WH | ≤ 1 | upper bound |

→ **N_A ≥ 0.2554** (1.61 with a pilot). This is a lower bound on the required strength, and it assumes the phase
reference.

**Closure (maximum strength):** J = N_A ħω₀ · h · τ · a · W ≤ cκ₊/16 = 0.3125 (in ħc/r₊² units, FGM (2.20)).

| Factor | Value | Label |
|---|---|---|
| h: beam energy share passing outside b_c | 1.2 × 10⁻²⁶ (Q_c = ω₀ a b_c/L = 787.25; at Q_c = 200 the tail formula is within 8 % of direct quadrature, ratio 0.93) | computed in the paraxial model |
| τ: share of that halo crossing B's horizon | ≤ 1 | upper bound |
| τ, isotropic redistribution by A's strong lensing | (b_c/2d)² = 2.5 × 10⁻⁸; B sits on the symmetry axis, a caustic, which would raise it | estimate; inference |
| Spectral leak below ωr₊ = 10 (unfocusable) | 10⁻⁸⁶·⁰ of the energy, 10⁻⁶⁰ of the halo term | der (saddle point) |
| A's above-barrier reflection | ∼ e^{−πω₀/κ₊} (10⁻²⁷²⁹); it heads back toward the source | assumption |
| a: FGM angular weight of the deposit | ≥ 1 for any deposit at θ_s = π; 1.60 for a uniform deposit | lower bound; computed |
| W: FGM time weight over the leak's arrival profile | 7.5 × 10⁸⁶; 10⁸⁴·⁷ with a 20 r₊ allowance for rays winding around A | computed; assumption |
| Signal's own 4D self-coupling (#619 R5) | matters only if its shape factor exceeds 2.2 × 10⁸⁴ | der; kept unknown, not zero |

→ N_A ≤ (cκ₊/16)/(ħω₀ h τ a W):
- with τ = τ_iso and the most favourable weights: **N_A ≤ 10⁻⁵⁴·⁴** [estimate];
- with τ = 1 and a uniform deposit (a = 1.60), the limit is 10⁻⁶⁴·⁴ [bound on τ, given h]. A deposit near B's θ = π
  pole would weigh more than 1.60 (H diverges logarithmically at zero separation), which only lowers the limit.

## 4. Do the two intervals overlap?

**No, for S★ (estimate).** B needs N_A ≥ 0.2554; the closure limit allows N_A ≤ 10⁻⁵⁴·⁴.

**Critical inequality**, in a form where numbers can be substituted:

  N_th · ħω₀ · h · τ · a · W ≤ (cκ₊/16) · η_WH.

For S★ it holds only if **τ ≤ 10⁻⁶¹·⁴**, even with the most favourable values (a = 1, the delay allowance,
η_WH = 1). Equivalently, the leak share must satisfy h·τ ≤ 10⁻⁸⁷·³.

**Missing coefficient: τ.** This is the two-mouth external transfer, from the part of the beam that misses A into
B's horizon. No existing calculation gives it; #619's f_geo is an isotropic s-wave estimate, not a bound for a beam.
- **Direction:** a lower bound τ ≥ 10⁻⁶¹·⁴ settles "no overlap". The isotropic estimate exceeds that by 10⁵³·⁸, so
  even a bound 53 orders below the estimate is enough. Only an upper bound below 10⁻⁶¹·⁴ would allow overlap.
- **Second missing input:** the in-situ phase reference (δln|ΔV| ≤ 6.9 × 10⁻¹⁶ between pilot and signal). If it
  fails, B cannot read the sign at any strength.

**Deciding quantity:** the energy share of the beam that misses A's capture radius and still reaches B's horizon,
h·τ. Against FGM's blueshift-weighted limit it must be ≤ 10⁻⁸⁷·³.

**What would change the verdict:** for this geometry and apodisation, the halo term reaches the critical level at
ω₀ ≈ 5.1 × 10³/r₊, 10.3 times the fixed carrier, with τ = τ_iso [estimate]. That would be a different source,
decided by terms at the critical leak share (≲ 10⁻⁸⁷) that nobody has computed.

**Not claimed:** a channel, a bit rate, or anything faster than light (the wormhole path is the slower one).

## 5. Plain answer

**Not through the wormhole, for this source.** To give B the 0.26 photons it needs to read the sign at 10 % error,
A must send at least that much. S★'s focused beam passes 1.2 × 10⁻²⁶ of its energy outside A's capture radius. If an isotropic share of that edge reaches B's horizon early, FGM's weighting makes it outweigh the strings'
negative energy by about 10⁵⁴, and the window closes.

**Which assumption it rests on:** the unverified transfer τ, the share of the beam's edge that ends up in B's
horizon. The answer would flip only if τ were below 4 × 10⁻⁶², which is 54 orders below the isotropic estimate. A
second, independent assumption is that B can get a phase reference through the window, which needs the shift stable
to 7 × 10⁻¹⁶.

**Via external space:** an ordinary signal would reach B sooner. S★ sends almost nothing that way, and the read
window excludes it.
