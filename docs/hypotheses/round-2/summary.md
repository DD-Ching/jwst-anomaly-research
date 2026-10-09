# Round 2 (2026-10-09): three System A worlds, three System B reviews

Everything here is a **hypothesis** or a reviewer's estimate ([B-est], [calc]). The System A texts were written
without literature or repository access; the System B reviews cite references their authors opened (lists at the end
of each B file).

## Brief for System A

Round 2 started from the three round-1 failure modes (D-069):
1. Lorentz invariance must be derived, with explicit bounds (|c_gw − c|/c < 10⁻¹⁵, laboratory 10⁻¹⁸), and the
   finite-valency obstruction must be avoided.
2. No rare long-range links (disordered locality).
3. No achromatic dimming (cosmic opacity).
4. Every amplitude must be derived from the axioms, not left free.

All three worlds satisfied rules 1–3 on paper. Each then put its whole risk into one fixed-amplitude prediction.

| Group | Seed | World | File | Review |
|---|---|---|---|---|
| R2-A | records + mutual consistency | "Co-Witness Consistency" (every two records share a future witness) | [R2A-consistency.md](R2A-consistency.md) | [B-on-R2A.md](B-on-R2A.md) |
| R2-B | a spectrum of rates, no positions | "Rate-Cone World" (rank-2 Jordan-algebra cone; position from stationary phase) | [R2B-spectrum.md](R2B-spectrum.md) | [B-on-R2B.md](B-on-R2B.md) |
| R2-C | minimal axiom count | "Two-Slot World" (N unlabelled C² spinors with the area form) | [R2C-minimal.md](R2C-minimal.md) | [B-on-R2C.md](B-on-R2C.md) |

## Verdicts

| World | Fixed-amplitude prediction | B verdict | Main reason |
|---|---|---|---|
| R2-A | Late GW energy equal to the GR horizon-absorbed energy (perfect reflector), echo delay ≈ 0.29 s for a GW150914-like remnant | **FAIL (excluded)** | Ergoregion instability at zero absorption (primary). Supporting: GW250114 reflectivity < 0.35 % (90 %; preprint arXiv:2610.12429, first return, coherent phase). GWTC-4.0 remnant tests are null in the predicted window (2603.19021). The idea is the perfectly reflecting exotic compact object (Mark et al. 1706.06155) |
| R2-B | Polarisation mirrored and circular handedness flipped in saddle (negative-parity) lensed images | **FAIL** for the linear mirror; conditional PASS (low yield) for the circular-polarisation flip and FRB echoes | Internally inconsistent: its own scalar-phase axiom already gives GR's answer, and the mirror is imposed. In B0218+357 the A/B polarisation-angle variations correlate positively at the lens delay (1802.10088), but the mirror requires anti-correlation |
| R2-C | Lightest neutrino massless ⇒ Σm_ν = 58.8 meV (NO) | conditional PASS, **not novel** (minimal seesaw) and not testable beyond collaboration analyses | 2.2–3.4σ tension with DESI DR2 + CMB [reviewer estimates: 2.2σ profile likelihood on 2503.14744, 3.4σ extrapolated with DES-Dovekie, 2605.21456; DESI quotes 3.0σ]. Not a formal kill under the world's rule (≥ 3σ with a ≤ 10 meV systematics budget; τ alone moves Σ by ~40 meV). Its element-number conservation law contradicts observed oscillations |

Failed in all three: "exact-zero" predictions (no dispersion, c_gw = c, η = 1, β = 0, w = −1). They coincide with
GR + ΛCDM and discriminate nothing.

## What this project can test

Nothing in round 2 earns search compute here:
- R2-A is excluded by published GW analyses.
- R2-B's only live parts need parity-resolved circular polarimetry of lensed radio quasars, which has not been
  published, or confirmed lensed FRBs, of which there are none.
- R2-C is tested by DESI, CMB and KATRIN analyses that this repository cannot improve on. It is tracked externally
  (DESI DR3, birefringence calibration, final KATRIN).

Optional, low priority: a sign test on public VLA monitoring of B0218+357 (programmes AH593, AB809). Correlate image B's
polarisation angle with ±image A's at the 11.3-day delay. That would formally close R2-B P1.

## Lessons for round 3

- Deriving Lorentz invariance from a symmetric primitive (a cone, an order, a spinor space) is easy and gives exact
  GR kinematics. The risk then moves entirely into one sector (horizons, polarisation transport, neutrino mass), and
  that sector already has precision data.
- A prediction with O(1) amplitude in a sector GR forbids is almost always already excluded. Round 3 should look for
  a derived small effect in a quantity measured with large N but never analysed for it, rather than an O(1) effect
  in a well-measured channel.
