# X2 S6 (screen): does rotation give a wormhole-stabilisation mechanism with explicit equations we can verify?

Owner request 2026-10-14. Literature and equation check only, ≤ 2 h (used about 0.8 h); no computation. The
model here is the phantom-scalar Ellis–Bronnikov (EB) family. It is **not** S4/S5's phantom + magnetic-NED SV
model, so nothing here transfers to it automatically. None of these papers had been assessed in the repository
before.

## What exists (IDs and titles checked on arxiv.org; sources read where marked)
**2301.05243** (PLB 848, 138349; v2 "results unchanged") and **2403.08387** (PRD 109, 124051; only v1).
- Source of 2403.08387 read, §III–VI.
- Method: radial-led modes (ℓ = 0 polar, coupled to ℓ = 1 axial) of 4D EB, to second order in slow rotation.
- **Main mode:** ω_I = ω_I⁽⁰⁾ + J²Δω⁽²⁾. For C = 0, r₀ = 1: −1.182 + 16.358 J².
  - Our S5 tool's static EB value is 1.18199. That independently matches their ω_I⁽⁰⁾.
  - The main mode reaches zero at J_c/A, which is ≲ 50 % of the extremal-Kerr bound J/A = 1/(8π).
- **Second mode:** it grows out of the static zero mode, with Δω⁽²⁾ = −28.185 for C = 0.
- The two branches **cross** near J/A ≈ 0.013, about one third of the bound.
  - The authors attribute the crossing to the quadratic truncation, which cannot produce a merger.
  - By analogy with 5D they *conjecture* the branches merge and vanish.
  - Their own conclusion: "a verification of this conjecture will need more sophisticated methods than
    perturbation theory".

**2510.11406** (PLB 875, 140335).
- Versions: v1 13 Oct 2025, titled "Charged wormholes can be long-lived"; v2 13 Mar 2026, retitled "Radial
  perturbations of charged wormholes" and the "long-lived" sentence removed. v1 and v2 abstracts compared; v2 source read.
- **Charged EB, closed-form background.** Two purely imaginary unstable modes merge at a critical charge. They
  then **continue as a complex pair**: equal imaginary parts, opposite real parts, still growing.
- The instability disappears only in the extremal Reissner–Nordström limit. Near it, the growth rate falls as
  (1 − M/r_T)^{3.1}.
- By analogy the authors now conjecture **the same for rotation**: the merged modes persist as an oscillating
  instability and fade only near extremality. This replaces the earlier "merge and disappear" picture.
- 2607.18399 (charged, static) adds a **new non-radial polar ℓ = 2 instability**.

**1309.2448** (5D, equal angular momenta, full rotation; source read, §V).
- This is the only non-perturbative rotating-stability calculation.
- Background and perturbations reduce to ODEs: three coupled G₀, G₁, G₂ with explicit potentials V_ij, after
  the gauge H₁ = H₃ − H₂ and H_φ = 0.
- Two unstable modes merge at (|J|/A)_cr = 0.10063 for symmetric wormholes and 0.1381 for |a₋∞| = 1.
- **Real Ω² only.** The authors solve with the auxiliary equation (Ω²)′ = 0 and real boundary conditions.
  "Beyond J_cr we do not find any further unstable modes **within the type of modes investigated**". So a
  complex pair after the merger, the 2510.11406 behaviour, was not searched for.

**4D full rotation.**
- Backgrounds exist: 1409.1503 and 1608.05253.
- Perturbations by spectral decomposition, 2401.02898 (PRD 109, 084013): sectors M_z = 2, 3 only, "no
  instabilities". The radial-led M_z = 0 sector has not been done.
- 2509.22118 (Universe 11, 325) derives the slow-rotation ℓ = 2, M_z = 2 equations, with no verdict yet.

**Other.**
- 2603.12748: 5D large-asymmetry backgrounds; no stability calculation.
- 2610.13206 (9 Oct 2026, not refereed): non-rotating 3+1 EB binaries; single-throat e-fold time within 1.9 %
  of the linear rate.
- No errata or comments found.

**Code.** No public rotating-wormhole background or perturbation code.
- FIDISOL/CADSOL: "free for scientific use" (1999 notice), no licence file, no confirmed download.
- Kadath: GPL, but no wormhole setups.
- A BSD-3 GRTeclyn fork (Zenodo 10.5281/zenodo.23257999): non-rotating binaries only.

## What the evidence supports (kept separate)

| Claim | Status |
|---|---|
| One growth mode gets slower with rotation | **Shown**, 4D, 2nd-order slow rotation (2403.08387); 5D full rotation (1309.2448) |
| All examined growth modes disappear | **Shown only in 5D for purely imaginary modes** (1309.2448). Not shown in 4D: the branches cross in the perturbative calculation |
| Complex (oscillating) growth modes after the merger | **Not searched** in either rotating case. **Found** in the charged analogue (2510.11406) |
| Full linear stability (all sectors) | **Not shown** anywhere. 4D M_z = 2, 3 is stable on full rotation (2401.02898); radial M_z = 0 on full rotation is open; charged analogues have a polar ℓ = 2 instability |
| Nonlinear stability, formation, phantom-field consistency | **Open**; no rotating nonlinear evolution exists |

**Validity of the approximations.** The second-order slow-rotation results are argued to hold up to about
50–60 % of the extremal bound, by analogy with Kerr. The two branches cross at about 33 %, inside that range,
but the merger itself needs higher orders.

## Proposed next unit (one): S6a, complex modes of 5D rotating Ellis wormholes after the merger
**Unknown.** In 1309.2448's 5D equal-J rotating EB wormholes, beyond (|J|/A)_cr = 0.10063 (symmetric case): do
unstable modes with complex ω (ω_R ≠ 0, ω_I > 0) persist, as in the charged analogue, or is the throat linearly
stable in this sector?
- Persisting modes would remove the only non-perturbative support for "rotation stabilises".
- Their absence would strengthen it, in 5D only.

**Why it is executable.** Cohomogeneity-1, so everything is ODEs:
- background equations in 1309.2448 §II–IV;
- perturbations G₀, G₁, G₂ with the V_ij given there (Eqs. eqv00–eqv22).

The S5 tooling already covers the parts needed: dense spectra including complex eigenvalues, constraint-respecting
time evolution, and an independent-equation residual.

**Steps and saveable stages:**
1. Re-derive the reduction for complex Ω. For H_φ the argument becomes Im(Ω²)∫|H_φ|² = 0, so H_φ = 0 also holds
   for complex unstable Ω. Check this, and check that the gauge H₁ = H₃ − H₂ stays consistent.
   - Saved: derivation and sympy check.
2. Background ODE solver (scipy solve_bvp).
   - Gate: reproduce published background quantities, and the static EB unstable mode, which S5's tool already
     has (1.18199 at b = 1 in 4D; the 5D value comes from the paper's Fig. 5 scaling).
   - Saved: background solutions on a J grid near and above J_cr.
3. Real-mode gate: the two purely imaginary branches must merge at (|J|/A)_cr = 0.10063 ± 1 %.
   - This is a reproduction, used only to validate the tool.
4. **Main result:** the full complex spectrum (dense generalised eigenproblem, plus time evolution) at several
   J/A between J_cr and the Myers–Perry bound, with residual checks.

**Budget and stops.** ≤ 6 h inside X2.
- Stop if step 2 cannot reproduce the background within 2 h.
- Stop if step 3 misses 0.10063 by more than 1 %; record why.
- Stop if step 1 shows the 3-ODE reduction fails for complex Ω; record the extra equations needed.
- No 4D PDE solver, no nonlinear evolution, no new echo search.

**What it would not answer.**
- 4D (the angular dependence does not factorise; this needs a 2D spectral PDE solver, which is not available).
- Our NED + SV model.
- Non-radial sectors.
- Phantom-field quantum consistency, formation, control or communication.

## Plain answer
- **Support:** in 4D slow rotation, one radial growth mode weakens, but a second appears and the two only cross.
  In 5D full rotation, the purely imaginary modes merge and vanish. Neither shows full stability.
- **Since then, the evidence leans against rotation as a stabiliser.** The charged analogue (2510.11406 v2)
  shows that merged modes can persist as an oscillating instability, and the same authors now conjecture this
  for rotation.
- **Most important unknown:** whether complex unstable modes exist after the merger. This is checkable only in
  5D with ODEs.
- **Worth a cross-session unit:** yes, one bounded unit (S6a above). It tests the single non-perturbative
  stabilisation claim directly. Either outcome narrows the question, though only in 5D.
- **Not a step toward observation:** there is no observational evidence for rotating phantom wormholes and no
  usable communication capability, and S6a would not change that.
