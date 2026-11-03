# Direct exclusion test for S★: a lower bound on its weighted leak into B under stated approximations (2026-11-02)

Owner-requested direct test, at most 90 min, limited to existing equations, derivations and small integrals.

**One question:** can an error-controlled lower bound K_low > 1.224 be established for K? K is the FGM-weighted
energy that the fixed source S★ sends into B's horizon, per unit mean emitted quantum.

**Setting:** units, carrier, aperture, pulse, spacetime and read window are those of `fgm_oneway_bit.md`.

**Sources and labels:**
- Numbers: `results/fgm_sstar_leak_bound/bound.json`, from `scripts/fgm_sstar_leak_bound.py`.
- The thresholds are read from `results/fgm_oneway_bit/bound.json` and `results/fgm_window/window.json`.
- Labels: [lit] FGM's text; [der] derivation; [computed] quadrature or root finding; [assumption]; [allowance] a
  margin taken in the direction that lowers the bound.

## 1. The two thresholds (reproduced)

- Read-out: N_A ≥ 0.2554, with an ideal phase reference and the receiving efficiency at its upper limit 1
  (Helstrom, Round 3).
- Mean-field closure: N_A K ≤ cκ₊/16 = 0.3125.
- The two hold together only if K ≤ 0.3125/0.2554 = **1.2235**. A provable K_low above that excludes S★.
- The planar sheet radiates half of its energy away from A, so η ≤ ½. That would halve the critical K. It is not
  used here; keeping 1.224 is the case most favourable to S★.

## 2. Method: why each step keeps a lower bound

1. **Receiving mode: B's s-wave channel.**
   - B's absorption is a sum of non-negative terms over its partial waves, so the s-wave share is a lower bound [der].
   - If this share is deposited uniformly over B's horizon, FGM's angular weight is ⟨H⟩/H(π) = 1.60 [assumption:
     uniform deposit; the kernel value is computed].
2. **Early region of B's horizon.** Only a window of labels u is kept, ending before any orbiting path arrives.
   The flux into the horizon is non-negative, and e^{−κ₊(u−x*₀)} is at least its value at the window end [der].
3. **All paths in one amplitude.** The static frequency-domain Green's function is symmetric (reciprocity; the radial
   Green's function with ingoing horizon conditions is symmetric).
   - The field at B's centre is therefore e^{iωd}/(4πd)·∫J Ψ, where Ψ is the field of a point source at B seen at the
     aperture, with B absent [der].
   - Every path enters this single overlap coherently.
4. **Propagation from A to B.**
   - Near A the throat is Reissner–Nordström with r₊ and r₋ [lit: FGM §3.1, (3.1)–(3.3)].
   - Ψ at the aperture is computed by geometric optics on RN null geodesics, by quadrature [computed], with
     kb = 2598 for family 2. The quadrature reproduces the second-order weak-field RN deflection, 4M/b + (15π/4)M²/b² −
     (3π/4)Q²/b², to < 3 × 10⁻⁴ at b = 50 r₊ (test).
   - Between the mouths propagation is flat. The two-centre corrections are O(r₊/d), and the logarithms go into
     the timing allowance [assumption].
5. **Paths from the aperture edge to B** [computed]:

| Path | Impact parameter b (r₊) | Arrival minus family 2 (r₊) |
|---|---|---|
| family 1 (same side) | 8.854 | −6.30 |
| family 2 (opposite side) | 5.196 | 0 |
| first orbit, same side | b_c + 0.0167 | +18.84 |
| first orbit, opposite side | b_c + 0.0080 | +21.77 |
| second orbit, same side | b_c + 0.0001 | +38.64 |

6. **Overlap.**
   - The overlap phase has no stationary point anywhere on the aperture. The family-2 mismatch is 0.43–0.62 ω,
     smallest at the edge [computed, six radii].
   - The overlap is therefore the edge saddle of the bump:
     |∫x·bump(x)e^{iqx}dx| ≈ e^{−1/4}√π 2^{−1/4} q^{−3/4} e^{−√q}, with q = 1079.5 [der].
   - Check on a simplified integral, the two-edge bump integral ∫_{−1}^{1} bump(x)e^{iqx}dx (not the full overlap
     with Ψ): its envelope peaks agree with the formula to ratios 1.011, 1.007 and 1.005 at q = 100, 300 and 600
     (double precision), and 1.004 (0.4 %) at q = 1079.5 (45-digit mpmath) [computed]. The full overlap was not
     integrated directly.
7. **Interference.**
   - Only families 1 and 2 are present in the window.
   - Family 1 is 7.7 × 10⁻⁶ of family 2 in amplitude (its edge parameter is q₁ = 2022). It can remove at most twice
     that share of the window energy [der].
   - The earliest orbit arrival over all aperture points is 18.84 r₊ after family 2, and the minimum sits at the edge
     [computed]. With the geometric-optics Green's function, each family contributes the source envelope delayed by
     its own path times, so the orbits are exactly zero in the window [der].
   - The spectral leak below ωr₊ = 10 (10⁻⁸⁶, Round 3) and A's above-barrier reflection are negligible.
   - The string's conical deficit on the axis (8πNGμ = 2.5 × 10⁻⁵ rad) shifts the convergence at B by about
     0.25 r₊ against a 9 r₊ central lobe [estimate, not stored].

## 3. The bound

| Item | Value | Label |
|---|---|---|
| s-wave energy into B per radiated energy | ≥ 1.7 × 10⁻⁴⁰ | computed, with allowances |
| window | 17.84 r₊ (the orbit gate minus 1 r₊); it holds 39 % of the envelope's energy | computed |
| family-2 delay relative to #619's reference u_c = v_c + d | 15.0 r₊ | computed |
| latest label in the window | u − x*₀ ≤ −682.3 r₊ (includes a 100 r₊ allowance) | allowance |
| **K_low** | **6.5 × 10³⁶** (log₁₀ 36.81), above 1.224 by 10³⁶·⁷ | computed under the stated approximations |

**Allowances** (each lowers the bound):
- B's s-wave transmission at ωr₊ = 500: Γ₀ ≥ 0.9.
- Curved-space source normalisation and coupling: ≥ 0.5 of the flat formula.
- Edge-saddle amplitude error: 30 %. The 0.4 % in §2 is the envelope-peak comparison on the simplified integral,
  not a measured error of the full overlap.
- Timing conventions: 100 r₊. This covers B's tortoise normalisation, the O(ln d/r₊) term of x*₀, the transmission
  delay into B, and the ≈ 1 r₊ forward smearing of the edge saddle.
- ⟨1/cos θ⟩ ≤ 1.118 for the sheet's radiation.

**Break-even.** The bound reaches 1.224 only in one of two cases:
- the timing offsets exceed the allowance by another 338 r₊;
- the amplitude factors are wrong by about 10¹⁸ (a factor 10³⁶·⁷ in energy).

**Isotropic comparison** (Round 3, kept as an estimate): K_iso ≈ 10⁵³·⁹. The bound is 17 orders lower because it
keeps only B's s-wave channel, one path and an early window.

## 4. Verdict

**Final conclusion:** Under the specified geometric-optics, aperture and mean-field approximations, the calculation supports
excluding the fixed source S★. A physical lower bound for the full two-mouth propagation remains to be verified.

Within these approximations (and #619's FGM point, an ideal phase reference, η ≤ 1) the calculation gives
K ≥ 6.5 × 10³⁶ against the critical 1.224.

**Error control:**
- Every factor is a lower bound or carries a stated allowance.
- The two asymptotic steps are approximations, not bounds: geometric optics at kb ≈ 2600 away from caustics, and the
  edge saddle. The edge-saddle formula was checked only on the simplified integral (envelope peaks, ≤ 1.1 %);
  neither step was checked against a full wave solution.
- This is not a theorem and not a verified physical lower bound. That needs the full two-mouth wave propagation.
  Within the approximations the margin is 10³⁶·⁷.

**Labels kept from earlier rounds:** P1 is C (mean field); P2 is frozen-rule A and post-hoc B. This round concerns
S★ only.

## 5. Calibration requirement and window fluctuations, kept separate

- **Calibration requirement** [computed, Round 3]: to hold the received phase to 0.1 rad, the emergence map must be
  known or stable to δln|ΔV| ≤ 6.9 × 10⁻¹⁶.
- **Model uncertainty:** FGM give x*₀ only up to an O(ln d/r₊) term, worth about 2.3 in ln|ΔV| [estimate]. This is
  a systematic of the model. Calibration through the window could remove it in principle.
- **Physical fluctuations of the window:** the quantum fluctuations of the shift (the strings' stress operator) are
  not computed by FGM, and their size is unknown. If they exceed 6.9 × 10⁻¹⁶ between calibration and signal, no
  phase reference exists.
- **This round** uses neither. It assumes an ideal reference, the case most favourable to S★.

## 6. Physical reading

S★ focuses its beam into A, but its aperture has a hard edge where the bump apodisation vanishes. That edge leaks a
tiny diffraction tail, e^{−√q} with q ≈ 1080. Part of the tail passes A on the far side at b = 5.2 r₊, is bent around A by its
gravity, and arrives at B on the axis together with the rest of the converging ring.

Even counting only B's s-wave channel, the energy reaching B is about 10⁻⁴⁰ of what S★ radiates. Its label lies at
least 682 r₊ before x*₀, even after the timing allowance. FGM's weighting multiplies positive energy arriving that
early by e^{κ₊·682} ≈ 10⁷⁴. In this calculation the leak therefore exceeds the strings' negative energy by about 10³⁶ at the smallest strength
B needs, and the window closes before B can read the sign.

Within these approximations, interference does not rescue it:
- The leak is fixed by S★'s own aperture and timing.
- The only paths that could cancel it arrive at least 17.8 r₊ later.

**Not claimed:** a channel, a bit rate, or anything faster than light. S★ now goes to human review.
