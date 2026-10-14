# Instability of 5D equal-J rotating Ellis wormholes past the mode merger: research summary (draft)

Status: internal draft, 2026-10-14. **Not reviewed by a human expert. Not submitted. No author contacted.**
Details: `X2_s6a_5d.md` (method and checks), `X2_s6a_repro.md` (reproduction), `X2_s6a_prior_work.md`
(literature search).

## Context (published results and conjectures, not ours)
**Published results.**
- Dzhunushaliev et al. (arXiv:1309.2448) built 5D Ellis wormholes, supported by a phantom scalar, rotating with
  two equal angular momenta. They studied linear radial-type perturbations with the ansatz of their Sec. V.
  - Two unstable modes with purely imaginary frequency merge at (|J|/A)_cr = 0.10063 (symmetric wormholes).
  - They state: "Beyond J_cr we do not find any further unstable modes within the type of modes investigated."
  - Their method (an auxiliary equation (Ω²)′ = 0 with real BCs) targets real Ω². That this restricts the search
    to real Ω² is our reading of the method, not their statement.
- In 4D, second-order slow rotation (arXiv:2301.05243, 2403.08387) finds the analogous two branches only up to
  their crossing.
- For *charged* (static) wormholes, arXiv:2510.11406 finds that the two merging modes continue as a
  complex-conjugate pair that still grows.

**Published conjectures and readings.**
- 2301.05243 conjectures that the 4D branches "merge and disappear".
- 2510.11406 conjectures, by analogy with the charged case, that rotating wormholes behave like the charged
  ones. It says this "can only be confirmed by a non-perturbative study in the future".
- Several papers read the 5D result as stabilisation beyond J_cr.

## What this work adds (possible new direct verification)
Linear perturbations in the same 5D model and sector were solved without restricting Ω² to be real.

1. **Equations.** The perturbation equations were re-derived symbolically from the metric.
   - They agree with the printed ones except one potential, V₁₀: the sign of its e^{2q} term.
   - The re-derived equations reproduce the published static eigenvalue (Ω²A^{2/3} = −21.80) and the merger
     point ((|J|/A)_cr ≈ 0.10063). The printed V₁₀ does not reproduce them.
   - This supports reading the printed V₁₀ as a typesetting error. It says nothing about the authors' own code.
2. **Two verified points.** At |J|/A = 0.1207 and 0.1420, just past the merger, the two modes continue as a
   complex-conjugate pair (time dependence e^{iΩt}, growth for Im Ω < 0).

| \|J\|/A | Ω²A^{2/3} | Ω A^{1/3} (growing root) | growth / static growth |
|---|---|---|---|
| 0 (static) | −21.80176 | −4.66924i | 1 |
| 0.1207 | −4.104519 + 5.527961i | ±1.179119 − 2.344107i | 0.50 |
| 0.1420 | −2.318701 + 7.461983i | ±1.657594 − 2.250848i | 0.48 |

   These are oscillating growing modes. **Each of the two points passed all six checks fixed before the
   verification run:**
   - grid convergence;
   - box size;
   - shooting outer radius;
   - algebraic residual of the complex eigenpair;
   - residuals of two linearised field equations not used to build the system, falling at second order;
   - agreement of a finite-difference eigenvalue solver with an independent Evans-function shooting method
     (≤ 1e-7 relative in Ω²).
3. **The 5D reading.** Within the searched literature, the "stabilised beyond J_cr" reading of the 5D result does
   not hold just past J_cr. There the instability continues in oscillating form, as conjectured in 2510.11406.

## What is not claimed
- **Points not verified.** No statement for |J|/A > 0.142 as a verified result.
  - |J|/A = 0.214 is a candidate: its shooting-radius check missed the frozen limit, 4e-6 against 1e-6.
  - Exploratory tracking reached 0.40.
  - The near-extremal regime (|J|/A → 0.5) was not computed.
- **Other sectors.** Odd parity was only screened, and other perturbation sectors were not studied. There is no
  claim of a complete spectrum.
- **Beyond this model.** Nothing is transferred to 4D, to other rotating wormholes, or to other matter models.
  No astrophysical or observational implication is claimed.
- **No priority claim.** No "first" claim is made. The literature search covers INSPIRE- and arXiv-indexed work
  up to 2026-10-14 (`X2_s6a_prior_work.md`).

## Who did what
| Contributor | Actual work |
|---|---|
| Repository owner (DD-Ching) | Chose and scoped the question; set the gates, the reading rules and the corrections during the run (no fitting to results, keep full complex vectors, no extension toward extremality, wording). Has not reviewed the physics as a domain expert |
| Claude (Claude Code agents, this repository's sessions) | All derivation and numerical code, runs, documentation and state files. Review and literature search were also done by Claude subagents |
| GPT | No GPT output is recorded in this repository or in this work's sessions. If the owner used GPT, e.g. to draft or check instructions, that has to be added by the owner |
| Human domain experts | **None so far** |

**Not yet checked by a human professional:**
- the symbolic derivation (Hopf-coordinate rewrite, linearisation, gauge argument);
- the identification of the V₁₀ typo;
- the numerical methods and their thresholds;
- the reading of the published literature;
- whether the result is of interest to the field.
