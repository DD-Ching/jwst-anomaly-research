# X2 S5: radial (ℓ = 0) stability of the S4 background at a − 2M = 3.575e-4

Owner request 2026-10-13. Question: does the S4 background have a radial growing mode at the original
cell, and does it grow faster or slower than the echoes? The S4 background is SV (Simpson–Visser) supported by
a phantom scalar plus magnetic nonlinear electrodynamics (NED), i.e. Bronnikov & Walia's source.
- Budget: ≤ 6 h, counted inside X2's 16 h. Used 1.7 h; X2 total 9.0 / 16 h.
- Out of scope: the polar sector, spin, parameter scans and nonlinear evolution.

## Equations (derived; `scripts/x2_s5_radial_derive.py` → `results/x2/s5_radial_equations.json`)
**Framework.** Cremona, Pizzocchero & Sarbach (arXiv:1911.13103, CPS):
- metric −α²dt² + γ²dx² + r²dΩ², with β = 0;
- gauge δΦ = 0, available because Φ′ = √2 a/(x² + a²) ≠ 0 everywhere, including at the throat;
- **δr is not fixed**, so the throat may move. There is no reflecting wall and no fixed throat radius.
- The variables are CPS's gauge-invariant A = δγ/γ and C = δr/r.

**Background with the S4 action** (κ = 1). Solving the field equations for the SV metric gives
- Φ′² = 2a²/(x² + a²)² > 0: a phantom scalar is required;
- U(r) = 6Ma²/(5r⁵), the magnetic-NED energy density, ½ of Bronnikov–Walia's L = (12/5)Ma²/r⁵;
- V = 4Ma²/(5r⁵), Bronnikov–Walia's potential.

**Fixed versus perturbed.**
- Fixed theory constants: the functions V(Φ) and L(F) (with M and a entering them as coupling constants), and
  the magnetic charge q.
- Perturbed fields: α, γ, r and Φ.

**Role of NED in the radial sector.**
- F = q sinθ dθ∧dφ is fixed by the Bianchi identity, so NED carries no radial degree of freedom. F_tx stays
  zero at ℓ = 0, because the electric charge is zero and L_F ≠ 0.
- Its stress is a fixed function of the *perturbed* areal radius: T^t_t = T^x_x = −U(r), and
  T^θ_θ = −U − rU′/2.
- It therefore enters the perturbation equations through U′ and U″ evaluated on δr.

**Structure of the linearised equations.**
- The tx equation gives CPS's constraint σ = 0, C′ = (α′/α − r′/r)C + (r′/r)A, with no NED term.
- The Klein–Gordon equation gives D′ (D = δα/α), also with no NED term.
- The xx equation gives C̈ with no C″. The θθ equation gives Ä.
- The NED terms appear in xx, θθ and tt as C·rU′ and C·(r/2)(rU″ + 3U′).

**Checks of the derivation.**
- U = V = 0 on the Ellis–Bronnikov wormhole reproduces CPS's master equation χ̈ − χ″ − 3b²χ/r⁴ = 0 symbolically
  (residual 0).
- On SV every coefficient is finite at the throat. The A″ coefficient there is f² ≈ 3.2e-8; there is no 1/r′.
- The Hamiltonian constraint is implied identically by σ = 0: residual about 1e-15 in a scratch check, not
  reproduced by the committed scripts. It is therefore not an independent test.

## Numerics (derived; `scripts/x2_s5_radial.py` → `results/x2/s5_radial.json`)
**Formulation.** A evolves with its wave equation. C is rebuilt from σ = 0 at every stage, so the solution stays
on the constraint surface. The C (xx) equation is not used and serves as the independent check. The grid is
uniform in the tortoise coordinate, with Dirichlet ends far from the mode.

**Failed approaches** (recorded so they are not repeated):
1. *Evolving C with its own equation.* This admits a pointwise, constraint-violating mode that grows at
   √max(coefficient) ≈ 0.168 M⁻¹, with σ violated at O(1). This was a scratch run, not reproduced by the
   committed scripts.
2. *Single-sided shooting with the parity condition at the throat.* In the long, nearly potential-free cavity
   the inward-dominant branch hides the true condition. It found the right λ, but placed the mode at the throat
   with a kink, which is spurious.

**Reference case: Ellis–Bronnikov.**
- Time evolution 1.18199, dense spectrum √1.397110 = 1.18199, Schrödinger operator of CPS Eq. MasterEllis 1.18207.
- This matches CPS's "τ of the order of the throat radius".

**SV + NED at a − 2M = 3.575e-4 (M = 1):**

| Method | λ (M⁻¹) | Mode location (r*) | C-equation residual |
|---|---|---|---|
| Time domain, dy = 0.1, seeded at ±315 | 0.1011238 | ±315.0 | 3.4e-5 |
| Time domain, dy = 0.05, seeded at ±315 | **0.1011244** | ±315.0 | 8.5e-6 |
| Time domain, seeded at ±200 | 0.1011238 | moves to ±315 | 3.4e-5 |
| Time domain, noise everywhere | 0.1011238 | ±315 | 3.4e-5 |
| Dense spectrum, Y = 1000, dy = 0.5 | √0.0102237 = 0.10111 | ±315 (even/odd pair) | — |
| Dense spectrum, Y = 2000, dy = 1 | √0.0102165 = 0.10108 | ±315 | — |

- **Convergence.** The rate converges in grid step, domain size and initial data. The C-equation residual
  falls about 4× per halving of dy, so it is converging to zero.
- **Location.**
  - The mode peaks at r* ≈ ±315 M: radial coordinate x = 0.444 M, areal radius R = 2.049M (2.4 % above a),
    f = 0.024. It is close to the throat in proper terms. The tortoise coordinate pushes it out because of the
    deep redshift there (f ≈ 1.8e-4 at the throat).
  - The Ellis–Bronnikov mode likewise peaks at x = b.
  - The local proper-time rate λ/√f = 0.65 M⁻¹ is close to Ellis–Bronnikov's 1.18/b ≈ 0.59 for b ≈ 2.
  - So this reads as the **throat-scale phantom-scalar instability**, redshifted, not a new mechanism.
  - The two sides give an even/odd pair, split only by tunnelling.
- **Seeding.** The ±315 seeds were chosen after the location was known (post hoc). The noise run and the ±200
  run do not depend on that choice and give the same rate.
- **Spurious eigenvalues.** The scheme uses only the A equation plus σ, so its spectrum also contains spurious
  positive eigenvalues that fail the C (xx) equation. Every dense eigenvalue is now residual-checked (see
  `runs.*dense*.top[].c_equation_residual`).
  - The dominant pair passes: residual 8e-4 at dy = 0.5 and 3e-3 at dy = 1, falling with dy.
  - The small throat eigenvalue (λ² = 2.1e-5 at dy = 0.5, 1.5e-6 at dy = 1, not converged) **fails**, with
    residual 0.96 / 1.00. It is spurious, not a second mode.
  - Review also found λ² = 1.2e-4 at Y = 600, dy = 0.25, with residual 0.88: spurious.
  - Only eigenvalues that pass the check count as modes.

## Answer: growth time against the echo clock
- Clock: the coordinate time t of the static background, i.e. the asymptotic observer's time (α → 1). That is
  the clock of S4's echo times.
- e-folds in each interval are independent of where t = 0 is put.

| Interval (S4 cell, M_det = 68.8 M☉) | Length | e-folds of the radial mode |
|---|---|---|
| e-folding time | 9.89 M = 3.35 ms | 1 |
| first echo (after free flight) | 1312.8 M ≈ 0.445 s | 133 |
| echo spacing | 1373.5 M ≈ 0.465 s | 139 |
| PMW start / end | 0.2 s / 1.2 s | 60 / 358 |

**The radial mode grows about 140 times faster than the echo cycle.** Its e-folding time (3.4 ms) is shorter
than the ringdown damping time (about 4 ms).

**Growth time is not collapse time.**
- The time to depart from the background depends on the initial radial-perturbation amplitude ε, as about
  τ·ln(1/ε), and on nonlinear evolution. Neither is computed here.
- Linear extrapolation only: unless ε ≲ e^{-133} ≈ 10⁻⁵⁸, the linear perturbation would reach order unity before
  the first echo. A merger is not expected to produce such a quiet initial state, but that is not shown here.

**Decision.** This applies the owner's qualitative rule ("if the perturbation can break the fixed-background
approximation within the relevant window, stop expanding echo searches on this background"). No numerical
threshold was written in advance; the margin (133 e-folds before the first echo) is far from any plausible one.
On this background the perturbation grows far within the relevant window, so the
fixed-background echo calculation is not self-consistent there. **No further echo-parameter search on this
background.** The S2/S4 numbers remain as computed, under their stated assumptions.

## What this removed, and what is left
- **Removed:** whether the S4 background (phantom + magnetic NED SV at a − 2M = 3.575e-4) is linearly radially
  stable. It is not, with τ = 9.9 M.
- **Not shown:**
  - other values of a or other matter sources (the anisotropic fluid has no action and no radial dynamics of
    its own);
  - nonlinear fate and the time to collapse;
  - the polar ℓ ≥ 2 sector;
  - spin.
- **Separate and still open:** quantum consistency of phantom fields, formation of such a throat, any control
  of it, and any communication capability.
- **Worth more investment?** Not for this background. Any further SV-echo work would first need a matter model
  that is radially stable. None is known for SV in our record; the nearby phantom-supported wormholes in S4c are
  unstable too.

## Re-run
- `python scripts/x2_s5_radial_derive.py` (symbolic, about 15 min) writes the coefficients.
- `python scripts/x2_s5_radial.py` (about 15 min, 4 processes) writes `results/x2/s5_radial.json`.
- Parameters: a − 2M = 3.575e-4 (Δt = 0.45 s at M_det = 68.8 M☉); κ = 1; Y = 600, dy = 0.1 / 0.05, T = 300; dense
  spectra Y = 1000 / dy = 0.5 and Y = 2000 / dy = 1.
