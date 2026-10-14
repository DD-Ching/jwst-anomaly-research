# X2 S6a: what happens to the two unstable modes of 5D equal-J rotating Ellis wormholes past their merger

Owner approval 2026-10-14.
- Budget: ≤ 6 h inside X2, including validation and review.
- Scope: only the 5D equal-angular-momenta rotating Ellis–Bronnikov wormholes of Dzhunushaliev et al.
  (arXiv:1309.2448), radial-type sector.
- Question: track the two published unstable modes across their merger at (|J|/A)_cr = 0.10063 (symmetric
  wormholes; scaled with the throat area A, not unscaled J, not a Kerr spin parameter).
- Out of scope: other 5D branches, nonlinear evolution, echo searches. No transfer to 4D or to NED + SV.

## Written before any computation
**Time convention.** 1309.2448 uses perturbations ∝ e^{iΩt}.
- A root Ω grows iff Im Ω < 0.
- The equations contain Ω only as Ω². Negative real Ω² is a purely imaginary unstable pair. Any non-real Ω²
  has one root with Im Ω < 0, so it grows and oscillates: growth rate |Im Ω|, frequency |Re Ω|.
- The charged-wormhole paper (2510.11406) uses the opposite sign convention; it is not mixed in here.

**H_φ = 0 for complex Ω.** Multiply the H_φ equation by the complex conjugate H̄_φ and integrate:
[ξ H̄_φ H_φ′] = ∫(ξ|H_φ′|² − Ω² w|H_φ|²) dl, with ξ = p h^{3/2} > 0 and w = p²h^{3/2}e^{−2a−q} > 0.
- For modes decaying at both ends the boundary term is zero.
- The imaginary part gives Im(Ω²)∫w|H_φ|² = 0, so H_φ = 0 for non-real Ω².
- The real part gives H_φ = 0 for real Ω² < 0.
- So the reduction to the gravitational sector holds for every unstable mode. Still to check: the gauge
  H₁ = H₃ − H₂ (step 1 below).

**Can complex Ω² be excluded a priori? No.**
- The reduced system is (ξG′)′ = ξVG − Ω² ξ w B G, with B = [[1,1,0],[0,1,0],[0,0,1]].
- An inner product making the problem self-adjoint would need a constant symmetric S with SB and SV symmetric.
  SB symmetric forces s₀₀ = s₀₂ = 0. Then SB has the principal minor [[0, s₀₁], [s₀₁, ·]], and s₀₁ ≠ 0 is
  required for S to be invertible.
- So SB is always indefinite. No positivity (energy) argument restricts Ω² to the real axis.
- This shows only that complex modes are *allowed*, not that they exist. That V₀₁ ≠ V₁₀ is not by itself an
  argument either.
- (Whether a constant symmetriser exists at all is checked in step 1 with the re-derived V.)

**Typo found in the printed V_ij.**
- In the static limit (c_ω = 0, a = q = 0) the paper's own zero mode (G₀ = −2G₂, G₁ = 3G₂, G₂ const)
  satisfies the G₀ and G₂ equations but **not** the G₁ equation with the printed V₁₀: the residual is −32/h.
- With 4e^{2a} − e^{2q} in V₁₀ it would vanish.
- The perturbation equations are therefore **re-derived from the metric**, not copied.

## Plan and gates (written before computing)
1. **Derivation** (sympy, 5D): background equations; linearised Einstein–phantom equations for the paper's
   perturbation ansatz.
   - Check that the gauge and H_φ reduction gives a closed G system.
   - Check against the paper's equations (record every difference) and against the static zero mode.
2. **Background** (≤ 2 h, including step 1):
   - solve the a, q ODEs for symmetric wormholes (r₀ = 0.5) at given c_ω;
   - gates: the constraint (R_ll) residual ≤ 1e-6; the static limit equals the closed form (a = q = 0); the
     Smarr relation (2/3)M = 2ω₀J holds to ≤ 1e-5 relative;
   - on failure: stop, save the reason; this is not a physics result.
3. **Known modes** (gate):
   - the static unstable mode;
   - two purely imaginary branches below J_cr, merging at (|J|/A)_cr = 0.10063 ± 1 %.
4. **Main computation**: track both branches across J_cr in steps of J/A, with complex Ω² allowed.
   - Mode tracking by continuity of eigenvector overlap, not by eigenvalue order.
   - Candidates are checked for grid and domain convergence, and for the residual of the original
     (un-reduced) linearised equations and constraints.
   - Cross-check at one representative point by an independent method (shooting with complex Ω²).
5. **Reading** (fixed now):
   - credible complex growth → "the instability continues, oscillating, in this 5D model and sector";
   - none found → state the searched region and sensitivity, no stability claim;
   - lost track numerically → "unresolved", not "vanished".
