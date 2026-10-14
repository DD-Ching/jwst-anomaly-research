# X2 S6a: reproduction package (equations, parameters, results, how to rerun)

Everything needed to recompute the S6a result without this repository's history. Code is BSD-3-Clause.
Status and limits: `X2_s6a_summary.md`.

## 1. Model and conventions
- Model: arXiv:1309.2448, Eq. (pertlineel), with a massless phantom scalar and R_μν = −8πG ∂_μΦ ∂_νΦ.
- Symmetric wormholes: a = q = 0 at both ends. Throat parameter r₀ = 1. Every reported quantity is scale-free.
- Time dependence e^{iΩt}. A mode is unstable iff Ω² ∉ [0, ∞), and grows with rate |Im Ω|.
- Units: 4πG = 1, as in the paper's figures (fixed by the extremal Myers–Perry point of their Fig. 1).
- Notation: h = l² + r₀², p = (l² + r₀²/2)/h, ξ = p h^{3/2}.

## 2. Background (paper Eqs. eqa, eqq, omprime, surface, Sm2)
- **Field equations.**
  - (ξa′)′/ξ = c_ω² e^{4a−q} / (2p³h⁴)
  - (ξq′)′/ξ = 4(e^{2(q−a)} − 1)/h
  - ω′ = c_ω e^{4a−q} / (h^{5/2}p²)
- **Boundary conditions.** a′(0) = q′(0) = 0. At l = L = 400: L a′ + 2a = 0 and L(q′ − a′) + 4(q − a) = 0 (from
  a ∝ l⁻², q − a ∝ l⁻⁴).
- **Global quantities.**
  - Throat area A = 2π² (h³p³e^{−(2a+q)})^{1/2} at l = 0.
  - J = −π c_ω/(16G).
  - M = 3π l³a′(∞)/(8G).
  - ω₀ = −c_ω ∫₀^∞ e^{4a−q}/(h^{5/2}p²) dl.
- **Checks.**
  - Smarr: (2/3)M = 2ω₀J.
  - Constraint: 8πG Q² = −h³p²R_ll must be constant in l.

## 3. Perturbations (re-derived; `results/x2/s6a_equations.json`)
**Ansatz.** The paper's perturbation ansatz in the gauge H₁ = H₃ − H₂ with H_φ = 0. These hold together for any
complex Ω (proof in `X2_s6a_5d.md`, Step 1). The variables are G defined by H₀ = G₂ + G₀, H₂ = −G₂,
H₃ = 2G₀ + G₁.

**System.**
- (1/ξ)(ξG_i′)′ = Σ_j V_ij G_j − Ω² w (B G)_i
- w = p e^{−2a−q}, B = [[1,1,0],[0,1,0],[0,0,1]]

**Potentials.** V_ij are as printed in 1309.2448, Eqs. (eqv00)–(eqv22), **except**
- V₁₀ = 2[4(4e^{2a} **−** e^{2q}) e^q h³p³ − c_ω² e^{6a}] / (h⁴p³) · e^{−2a−q} (printed: + e^{2q}).
- The symbolic expressions are in `G_system.V` of the JSON.

**Boundary conditions.**
- G → 0 at both ends.
- Even parity at the throat, G′(0) = 0: this is the paper's sector, normalised there by G₂(0) = 1.
- Odd parity, G(0) = 0, was screened only.

**Independent checks** (not used to build the system):
- the linearised R_tl in the corotating frame (∝ Ω);
- the linearised R_ll.
- Both are in `constraints` of the JSON.
- The printed H_ω′ (Eq. hop) is twice the derived one. H_ω is eliminated, so the G system does not depend on it.

## 4. Numerical methods and settings
**Dense solver.**
- Finite differences in s, with l = sinh s. Conservative 3-point stencil, a parity ghost point at s = 0,
  Dirichlet at s = S.
- Standard complex eigenproblem for Ω² (scipy.linalg.eig), with **full complex eigenpairs** kept.
- Boxes S = 5 / 6 / 7.5 (l ≤ 74 / 202 / 904). Grid n = 400 / 480 / 600 (ds = 0.0125) and n = 960 (ds/2).
- Richardson extrapolation from n = 480 and 960 at S = 6.

**Shooting (independent).**
- The (G, ξG′) system is integrated in s with DOP853, rtol 1e-10. It starts from l = sinh S with decaying data
  G = e_k, ξG′ = ξ(−κ − 3/(2l))e_k, where κ = √(−Ω²) and Re κ > 0. The potentials are taken from a cubic spline of
  the background.
- The determinant of ξG′ at the throat (3×3, even parity) is solved by complex secant iteration.
- S = 3.5 and 4.0 are compared. Iteration stops at a relative step < 1e-9, which was reached at every reported
  point.
- Acceptance: the returned root is the iterate that followed the smallest relative step, and it is accepted if
  that step is ≤ 1e-6. Every step is stored.

**Tracking.** Modes are followed in c_ω by the eigenvector overlap |v_k* · v_prev| / (‖v_k‖ ‖v_prev‖).

**Pass criteria** (`CRITERIA` in the script; frozen before the verification run, after trial runs):

| Check | Requirement |
|---|---|
| Dense (Richardson) vs shooting | ≤ 1e-5 relative |
| Box S = 6 vs 7.5 | ≤ 1e-6 |
| Shooting S = 3.5 vs 4.0 | ≤ 1e-6 |
| Complex eigenpair residual | ≤ 1e-8 |
| Both constraint residuals | ≤ 1e-4, falling ≥ 3× when the step halves |
| Eigenvector norm beyond l = sinh(S)/3 | ≤ 1e-6 |

## 5. Results (`results/x2/s6a_rot5d.json`)
Scaled eigenvalue z = Ω²A^{2/3}.

| c_ω (r₀ = 1) | \|J\|/A | A | z, dense (Richardson) | z, shooting S = 4 | Status |
|---|---|---|---|---|---|
| 0 | 0 | 6.978864 | −21.8017583 | −21.8017585 | gate (paper −21.80) |
| 0.4 | 0.120731 | 8.174903 | −4.1045193 + 5.5279618i | −4.1045190 + 5.5279613i | **verified** |
| 0.5 | 0.142040 | 8.685558 | −2.3187017 + 7.4619833i | −2.3187015 + 7.4619831i | **verified** |
| 1.0 | 0.214391 | 11.508893 | 2.4541671 + 8.5242290i | 2.4541674 + 8.5242295i | candidate (S = 3.5 vs 4.0: 4e-6) |

**Merger.** J/A ∈ [0.100622, 0.100633] at both grid steps. The paper gives 0.10063.
- Estimated exceptional point: 0.1006324 (ds) and 0.1006288 (ds/2) (`merger_ep_estimate`).
- Method: the squared half-splitting is fitted linearly in c_ω between the bracket ends, then mapped linearly
  to J/A.

**Below the merger.** The tracked branches match the digitised Fig. 5 to ≤ 0.011 in z (`branch_vs_fig5`).

## 6. How to rerun
Python 3.12.3, numpy 2.5.4, scipy 1.18.1, sympy 1.14.0. The project's `dev` extras pin the minimums.

```
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -e ".[dev]"
python scripts/x2_s6a_derive5d.py              # about 2.5 min: symbolic derivation -> s6a_equations.json
python scripts/x2_s6a_rot5d.py                 # about 25 min, 4 processes -> s6a_rot5d.json
python scripts/x2_s6a_rot5d.py --shooting-only # under 1 min: re-run the shooting on stored dense results
python -m pytest -q tests/test_x2_s6a_rot5d.py # about 10 s: tool and gate tests
```

The JSON files are strict JSON (no NaN; a check that does not apply is `null`). Digitised paper values are in
`PAPER_FIG1` and `PAPER_FIG5` of the script, read from the paper's gnuplot EPS sources (about 0.2 % of the axis
range).

## 7. Suitability: research article or Comment? (assessment, not a decision)
**A Comment on 1309.2448: poorly suited.**
- The paper limits its claim to "the type of modes investigated", and its numbers are reproduced. So there is no
  error in a published *result* to comment on.
- The V₁₀ misprint does not change their numbers. It would fit an erratum by the authors, or a footnote in a new
  paper.

**A short research article (brief report or letter): suitable in principle, after the steps below.**
- A related question is posed in print as a conjecture. 2510.11406 states it mainly for 4D rotation, so applying
  it to 5D is our reading. The answer here is direct, and the method and the checks are reproducible.
- As it stands, the result covers two points just past the merger. A referee would very likely ask for:
  - the full |J|/A range up to the extremal limit, including the candidate at 0.214 and what happens there;
  - the odd-parity sector;
  - an independent time-domain evolution;
  - an expert check of the derivation.
- The owner's scope decision did not extend the computation in this round.

**Before either: review by a human domain expert is needed.** Section 8 lists what is unchecked.

## 8. Not yet checked by a human professional
- **Derivation:** the Hopf-coordinate rewrite, the linearisation, the gauge and H_φ argument, and the V₁₀ typo.
- **Numerics:** the box-artefact explanation (Jordan-type mass matrix), the shooting asymptotics, and whether the
  thresholds suit the problem.
- **Data and literature:** the digitised reference values, and the reading of 1309.2448, 2510.11406 and the 4D
  papers.
- **Model limits:** how this sector relates to physically relevant perturbations (non-radial sectors, other
  matter models).
