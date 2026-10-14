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
[ξ H̄_φ H_φ′] = ∫(ξ|H_φ′|² − Ω² ξw|H_φ|²) dl, with ξ = p h^{3/2} > 0 and ξw = p²h^{3/2}e^{−2a−q} > 0
(w = p e^{−2a−q} as in the G system).
- For modes decaying at both ends the boundary term is zero.
- The imaginary part gives Im(Ω²)∫ξw|H_φ|² = 0, so H_φ = 0 for non-real Ω².
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

## Step 1: reduction and signs (derived; `scripts/x2_s6a_derive5d.py` → `results/x2/s6a_equations.json`)
**Setup.** The 5D Ricci tensor of Eq. (pertlineel) is computed directly. The angular part is rewritten
exactly in Hopf form (θ̃ = 2θ, χ = ψ − φ, η = ψ + φ), which removes most of the angle dependence.

**Background.**
- Every Ricci component except R_ll vanishes identically on Eqs. (eqa), (eqq), (omprime) and (solp).
- R_ll equals the paper's Eq. (Rll-expr) exactly.
- So the ansatz and the rewrite are verified against the paper.

**Gauge.** A pure-gauge vector field ξ = e^{iΩt}(T, Z, 0, 0, Y) keeps the ansatz (g_tl = g_lη = 0) when
Y′ = 2ωT′ and T′ = iΩ p e^{−2a−q} Z.
- Under it, H_φ shifts by φ′Z, and H₁ − H₃ + H₂ shifts by 2(Z′ − (ξ′/ξ)Z) − 2iΩT.
- Since φ′ = Q/ξ ≠ 0 everywhere, Z = −H_φ/φ′ reaches the gauge δΦ = 0 for any perturbation.
- In that gauge the linearised scalar equation reduces to (Q·(H₁ − H₃ + H₂))′ = 0. So H₁ − H₃ + H₂ is a
  constant, and a constant T removes it (Ω ≠ 0).
- Therefore H_φ = 0 and H₁ = H₃ − H₂ hold together for **any** Ω, complex included. This needs no positivity
  argument; the integral argument above agrees with it.
- For unstable modes (Ω² ∉ [0, ∞)) all fields decay exponentially at both ends. The gauge functions keep that
  decay, so every boundary term vanishes.
- The coordinate l runs over the whole line with h ≥ r₀² and p ≥ ½, so there is no inner boundary. The
  throat l = 0 is a regular point, where only parity enters.

**G system (derived).**
- The first-derivative structure is exactly (1/ξ)(ξG_i′)′.
- The mass matrix is the paper's B (G₀ equation: G₀ + G₁) with weight w = p e^{−2a−q}.
- The derived V_ij equal the printed ones except **V₁₀**: derived = printed − 16e^{2q−2a}/h, i.e.
  4(4e^{2a} − e^{2q}). This confirms the typo found a priori.
- Numerically, the re-derived equations reproduce the paper's published static eigenvalue (−21.80, below),
  while the printed V₁₀ gives −146.16 (and −5.62); see `printed_V10_static` in `results/x2/s6a_rot5d.json`. This supports reading the printed V₁₀ as a typesetting error. It does not
  establish what the authors' own code implemented.
- The printed Eq. (hop) for H_ω′ is twice the derived one.
  - The derived one comes from the linearised R_lη equation. It also agrees with perturbing the first integral
    of Eq. (omeq), written for the general F's as h^{5/2}p²e^{F₃−4F₀+(F₃−F₂−F₁)/2}F_ω′ = c_ω.
  - H_ω′ is eliminated before the G system, so this does not affect it.
- Constraints: the corotating-frame R_tl (∝ Ω) and R_ll are not used to build the G system. They are kept as
  independent residual checks on every mode.
- Running the script reproduces all of the above (about 2.5 min).

## Step 2: background (gates passed; about 0.3 h)
Symmetric solutions with r₀ = 1 are solved with scipy solve_bvp on l ∈ [0, 400], with asymptotic Robin
conditions (a ∝ l⁻², q − a ∝ l⁻⁴).
- The plan said r₀ = 0.5; all reported quantities are scale-free, so this changes nothing.

**Units.** 4πG = 1, as in the paper's figures. This is fixed by the extremal Myers–Perry dot in Fig. 1 (read
from the EPS vector file):
- J/A = 0.5 = 1/(8πG) and M/A^{2/3} = 3.218 = 3/(8Gπ^{1/3}) agree only for 4πG = 1.
- So (|J|/A)_cr = 0.10063 is 20 % of the extremal bound.
- It is a scaled ratio. It is not unscaled J and not a Kerr spin.

**Gates.**

| Gate | Result |
|---|---|
| Static limit | a = q = 0 to machine precision |
| Constraint: 8πG Q² = −h³p²R_ll constant in l | relative variation ≤ 3.4e-7 up to c_ω = 2 |
| Smarr (2/3)M = 2ω₀J | relative residual ≈ 2e-7 |
| Fig. 1 M/A^{2/3}(J/A), digitised (`fig1_comparison`), at J/A = 0.067 / 0.142 / 0.214 / 0.285 / 0.346 | 0.1431/0.1435, 0.5857/0.5862, 1.1772/1.1776, 1.7936/1.7935, 2.2878/2.2888 |

Digitised values come from the paper's gnuplot EPS paths, calibrated by the tick marks. They are data read
from a figure, with about 0.2 % of the axis range as uncertainty.

## Method revisions made after trial runs (post hoc; recorded before the verification run)
The plan above fixed the checks (grid, domain, residuals of the un-reduced equations, an independent method)
but not their numerical settings or thresholds. The settings below were changed **after** seeing trial
results. They were frozen in `CRITERIA` (`scripts/x2_s6a_rot5d.py`) before the verification run. They were
not set in advance.

| Item | Before (trial) | What the trial showed | After (frozen) |
|---|---|---|---|
| Spectrum filter | eigenvalue must recur in boxes S = 5 and 6 (l ≤ 74 / 202) at the same grid step | the box turns the continuum into many spurious complex pairs; same-step grid modes also recur; above J/A ≈ 0.36 the mode spreads beyond S = 5 and dropped out of the filter (a filter failure, not a disappearance) | different boxes **and** grid steps; tracking by eigenvector overlap; verification boxes S = 5 / 6 / 7.5 |
| Shooting outer radius | S = 6 (l = 202), secant stop at relative step 1e-9 | did not converge at c_ω = 0.5. The determinant is smooth there, but it loses digits to cancellation at large S | S = 3.5 and 4.0 compared. The secant stops at stagnation and fails if the *minimum* relative step in Ω² is > 1e-6 (the docs first said "last"; see "Shooting stop"). At S = 3.0 / 3.5 / 4.0 / 4.5 the roots were −2.31965, −2.318702, −2.318702, −2.318709 (real part, scaled), with minimum steps 6e-10, 6e-10, 4e-10, 1.5e-7 |
| Thresholds in `CRITERIA` | none numeric | trial at c_ω = 0.5: Richardson-extrapolated dense and shooting agreed to about 1e-7 | grid ≤ 1e-5, box ≤ 1e-6, shooting radius ≤ 1e-6, eigenpair ≤ 1e-8, constraints ≤ 1e-4 and falling ≥ 3× per halving, localisation ≤ 1e-6 |
| Representative points | — | chosen from the exploratory scan | c_ω = 0.4, 0.5, 1.0 (J/A ≈ 0.121, 0.142, 0.214) |
| Merger test | complex pair among the filtered eigenvalues | — | non-real eigenvalue with Re z < −1 in the S = 6 spectrum, bisection at two grid steps (n = 480, 960) |

The exploratory scan also followed the pair to J/A = 0.40. That run is **not** part of the verified result:
these points are candidates, and per the owner's instruction (2026-10-14) the search was not extended toward
extremal rotation.

## What each number measures
- **Eigenvalue difference between methods.** |z_dense − z_shoot| / |z_shoot|, where z = Ω²A^{2/3}.
  - z_dense is the Richardson extrapolation from grid steps ds and ds/2 in box S = 6.
  - z_shoot is the shooting root at S = 4.0.
  - This is the closest thing to an eigenvalue accuracy reported here.
- **Grid and outer-boundary uncertainty.** Each is the change in z when only one setting changes:
  - box S = 6 → 7.5 at a fixed step;
  - grid step ds → ds/2 (before extrapolation);
  - shooting radius sinh 3.5 → sinh 4.0.
- **Shooting stop.** Relative secant steps |Δ(Ω²)| / max(1, |Ω²|). They measure convergence of the
  iteration. They are **not** a determinant residual and not an eigenvalue accuracy. The determinant itself is
  not normalised and is not used as a residual.
  - As run on 2026-10-14, `Shooter.root` returned the iterate after the *smallest* step and that step's size. It
    accepted a root if that minimum was ≤ 1e-6. The docs said "last step", which did not match the code.
  - Audit (2026-10-14, writeup round): `root` now records every step. The shooting was re-run on the stored
    dense results (`python scripts/x2_s6a_rot5d.py --shooting-only`). The static gate, whose dense spectra take
    seconds, was recomputed whole; it came out identical.
  - At every representative point, both radii and the static mode, the iteration reached the 1e-9 tolerance:
    last step = minimum step, 2e-10 to 9e-10, after 10–27 iterations.
  - The roots are bit-identical to the earlier run. So the earlier "last_step" values were in fact last steps,
    and no status changes. The previous records are kept as `shooting_before_rerun`.
- **Eigenpair residual** (algebraic, of the discretised problem): ‖Kv − λv‖ / (‖Kv‖ + |λ|‖v‖), with the full
  complex λ and v. It shows that the linear algebra is solved, nothing about discretisation error.
- **Constraint residual** (physics, un-reduced equations). The corotating-frame R_tl/Ω and R_ll linearised
  equations are evaluated on the complex mode with finite-difference derivatives, for l ≤ sinh(S)/3.
  - Normalisation: max |Σ terms| / max Σ|terms|, i.e. relative to the size of the individual terms of each
    equation.
  - These two equations are not used to build the G system.
  - A genuine mode shows second-order decrease with the grid step. A box artefact does not: at c_ω = 0.5 a box
    eigenvalue near z = 5.1 + 0.6i gave 0.03 / 0.06, unchanged under halving.

## Step 3: known modes (gates passed; `scripts/x2_s6a_rot5d.py` → `results/x2/s6a_rot5d.json`)
**Static mode.** z = Ω²A^{2/3}:

| Method | z |
|---|---|
| Dense Richardson | −21.801758 |
| Shooting | −21.801758 (difference 8e-9 relative) |
| Paper Fig. 5 (digitised) | −21.80 |

**Branches below the merger.** Both branches are followed by eigenvector overlap (≥ 0.998 per step).
- At J/A = 0.0496 / 0.0745 / 0.0905, against the digitised Fig. 5:
  - static branch −18.552 / −14.483 / −10.732, against −18.563 / −14.487 / −10.727;
  - zero-mode branch −0.363 / −1.219 / −2.518, against −0.360 / −1.212 / −2.521.
- Differences are ≤ 0.011, within the figure's read-off precision.
- The digitised point at 0.0970 lies beyond the last tracked point (0.0960) and is not compared.
- Our values are linear interpolations between tracked points (`branch_vs_fig5` in the JSON).
- These branches are reproduction gates. They were not held to `CRITERIA`: the R_ll residual on the zero-mode
  branch reaches 5.6e-4 at J/A = 0.048, where that mode is still wide.

**Merger.** The first non-real pair appears at J/A ∈ [0.100622, 0.100633].
- At the lower end of the bracket (grid step ds) the spectrum has exactly two eigenvalues with Re z < −1,
  both real: −5.930 and −5.674. These are the two branches.
- At the upper end it has one conjugate pair, −5.8009 ± 0.0143i, at their midpoint. Nothing else has Re z < −1.
- The same holds at ds/2.
- The bracket is the same at grid steps ds and ds/2.
- Against the paper's 0.10063, the deviation is ≤ 1e-4 relative. **The ±1 % gate is passed.**
- An identical bracket does **not** mean the nearby eigenvalues have converged. At the bracket ends they differ
  between the two grids:

  | Grid | Lower end (two real) | Upper end (pair) | Estimated J/A at the merger |
  |---|---|---|---|
  | ds | −5.9298, −5.6737 | −5.80092 ± 0.01427i | 0.1006324 |
  | ds/2 | −5.9036, −5.6997 | −5.80081 ± 0.07885i | 0.1006288 |

  - Near an exceptional point the half-splitting goes as √(J − J_EP), so a small grid shift of J_EP changes the
    split a lot.
  - J_EP is estimated by fitting the squared half-splitting linearly in c_ω across the bracket. It moves by
    3.6e-6 in J/A (relative 4e-5) between the grids.
  - The midpoint −5.8008/−5.8009 agrees to 2e-5.
  - This does **not** affect the representative points. They lie 0.02–0.04 away in J/A, where eigenvalues change
    by only 2–5e-5 under grid halving and the extrapolated values agree with shooting to ≤ 1e-7.

## Step 4: past the merger (main result)
**Representative points.** All checks of `CRITERIA` were run on the mode tracked from the merger. Values are
scaled: z = Ω²A^{2/3}, growth rate |Im Ω|A^{1/3}, frequency |Re Ω|A^{1/3}.

| c_ω | J/A | z (shooting, S = 4) | growth | frequency | Status |
|---|---|---|---|---|---|
| 0.4 | 0.1207 | −4.104519 + 5.527961i | 2.344 | 1.179 | **verified mode** (all 6 checks) |
| 0.5 | 0.1420 | −2.318701 + 7.461983i | 2.251 | 1.658 | **verified mode** (all 6 checks) |
| 1.0 | 0.2144 | 2.454167 + 8.524229i | 1.791 | 2.380 | **candidate**: fails the shooting-radius check |

**Checks at each point.** Thresholds are in `CRITERIA`; definitions are under "What each number measures".

| Check | c_ω = 0.4 | c_ω = 0.5 | c_ω = 1.0 |
|---|---|---|---|
| Dense (Richardson) vs shooting, eigenvalue | 9.0e-8 | 3.4e-8 | 5.9e-8 |
| Box S = 6 → 7.5, eigenvalue | 8e-12 | 6e-12 | 3e-11 |
| Shooting radius sinh 3.5 → sinh 4.0, eigenvalue | 8.5e-8 | 3.5e-8 | **4.0e-6** (limit 1e-6) |
| Grid step halving, eigenvalue (before extrapolation) | 5e-5 | 2e-5 | 4e-5 |
| Complex eigenpair residual (algebraic) | 6e-11 | 5e-11 | 3e-11 |
| R_tl, R_ll residual at ds/2 | 2.6e-6, 2.8e-6 | 2.9e-6, 1.8e-6 | 4.0e-6, 2.4e-6 |
| Residual ratio ds → ds/2 | ≥ 3.97 (second order) | ≥ 3.97 | ≥ 4.00 |
| Norm fraction beyond sinh(S)/3 | 6e-15 | 6e-15 | 9e-15 |

**The c_ω = 1.0 point.**
- Its root moves by 4e-6 between the two shooting radii. That is more than the frozen limit, so the point stays a
  **candidate**.
- The cause is **unexplained**. Slower decay does not account for it: the spatial decay rate Re κ is 0.79 there
  against 1.16 at c_ω = 0.4, which predicts contamination of about 4e-12, far below the observed 4e-6.
  No further radius was run, by the scope decision.
- Every other check passes at this point.

**Tracking between the points.**
- The pair is tracked from the upper end of the merger bracket. No start value is chosen by hand.
- The first step away from the exceptional point cannot separate the two conjugate members: overlap 0.9963
  against 0.9962. That is expected, since a real-like vector overlaps both members equally.
- From J/A = 0.104 to 0.214 (11 steps between 12 points), the overlap is ≥ 0.992 per step, while the next-best
  overlap is ≤ 0.96.
- The R_tl and R_ll residuals stay ≤ 1.6e-5 at every step (grid step ds).
- At the exceptional point the two real modes become one complex-conjugate pair. Which real branch becomes
  which member is not defined there, because conjugate members have equal overlap with a real vector.

**Odd parity: screened only.**
- At the three points, no odd-parity eigenvalue with 0.3 < |z| < 60 recurs (within 2e-3) in two boxes that
  also have different grid steps: S = 6 with ds, and S = 7.5 with 2ds/3.
- A first version of the screen used the same grid step in both boxes. It let through a continuum-like grid
  sequence with Re z ≥ 2.7. Review caught this, and it was fixed before the final run.
- This is a screen, not a proof of absence. Odd modes with |z| < 0.3, or modes failing the box/grid recurrence,
  are not excluded, and no odd candidate was put through `CRITERIA`.

**Exploratory, not verified.** A scratch scan with the earlier filter followed the pair to J/A = 0.40: z = 1.99 + 1.09i
there, growth 0.37. Those points are candidates at most. The search was not extended toward the extremal
limit (J/A = 0.5).

## Reading (rules fixed before computing)
- **Credible complex growth found.** At J/A = 0.121 and 0.142 the merged modes continue as a complex-conjugate
  pair. Ω² is non-real, so one root has Im Ω < 0: an **oscillating growing mode**.
  - It passes the re-derived 5D equations, two un-reduced constraint equations, grid, box and outer-radius
    convergence, and an independent shooting method.
  - Growth |Im Ω| ≈ 2.3 A^{−1/3}, about half the static rate (4.67 A^{−1/3}), with frequency 1.2–1.7 A^{−1/3}.
- **Statement:** in this 5D equal-J model, in this perturbation sector (the paper's ansatz, even parity), the
  instability continues past the merger in oscillating form. The paper's "no unstable modes beyond J_cr" holds
  only for purely imaginary modes.
- **Not shown:**
  - stability or instability at J/A ≳ 0.21 (candidate at 0.214, exploratory to 0.40);
  - the near-extremal limit;
  - odd parity (screened only) and other sectors;
  - 4D, other rotating wormholes, NED + SV;
  - anything about communication.

## Writeup-round corrections (2026-10-14)
- **NaN.** The first point of each tracked branch has no previous vector, so "next-best overlap" does not apply.
  It was written as NaN, which strict JSON readers reject. It is now `null`, and the writer uses
  `allow_nan=False`. Both result files parse with a strict reader.
- **Shooting record:** see "Shooting stop" above.
- **Merger-grid differences:** see Step 3.
