---
title: "Unstable modes past the merger in 5D rotating Ellis wormholes (S6a)"
---

**Archival brief.** Version 2026-10-16, based on commit
[`a57b824`](https://github.com/DD-Ching/jwst-anomaly-research/commit/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c),
the last commit that touches the S6a files. Repository: DD-Ching/jwst-anomaly-research, licence BSD-3-Clause.
**No human domain expert has checked this work. It has not been submitted anywhere.**

> **What this is and is not.** This is a *linear-stability result in a theoretical model*. It is **not** an
> observation of a wormhole. It does **not** imply that wormholes exist, and it does **not** imply that they could
> be used to travel or to communicate. It is not a complete stability analysis of the model.

## Result

**Model and sector.** The 5D Ellis wormholes, supported by a phantom scalar and rotating with two equal angular
momenta, of Dzhunushaliev et al., arXiv:1309.2448 (PRD 88, 124028). Only: symmetric wormholes, that paper's
perturbation ansatz (its Sec. V), even parity.

**Published background.** In that paper two unstable modes with purely imaginary frequency merge at
$(\lvert J\rvert/A)_{\mathrm{cr}} = 0.10063$. Beyond it the authors "do not find any further unstable modes within the
type of modes investigated". Their method targets real $\Omega^2$; that this restricts the search is our reading.

**What we computed.** The same linear perturbation problem, solved without restricting $\Omega^2$ to be real.
At two verified points just past the merger, $\lvert J\rvert/A = 0.1207$ and $0.1420$, the two merged modes
continue as a **complex-conjugate pair that oscillates and still grows**. Time dependence is $e^{i\Omega t}$; a mode
grows when $\mathrm{Im}\,\Omega < 0$. The scaled eigenvalue is $z = \Omega^2 A^{2/3}$, with $A$ the throat area.

| $\lvert J\rvert/A$ | $z = \Omega^2 A^{2/3}$ (and its conjugate) | $\Omega A^{1/3}$, growing roots of the pair | growth / static growth |
|---|---|---|---|
| 0 (static) | $-21.80176$ | $-4.66924\,i$ | 1 |
| 0.1207 | $-4.104519 + 5.527961\,i$ | $\pm 1.179119 - 2.344107\,i$ | 0.50 |
| 0.1420 | $-2.318701 + 7.461983\,i$ | $\pm 1.657594 - 2.250848\,i$ | 0.48 |

The sign $\pm$ covers the pair: $-1.179 - 2.344\,i$ belongs to the listed $z$, $+1.179 - 2.344\,i$ to its
conjugate. The growth rate is about half the static one at both points.

**Relation to the literature (our reading).** Some papers read the 5D result as stabilisation beyond the merger.
In this sector, just past the merger, the instability continues in oscillating form. This is consistent with a
conjecture of arXiv:2510.11406, which is stated mainly for 4D rotation; applying it to 5D is our reading.
We found no prior computation of complex unstable modes past the merger for these 5D wormholes in the literature
indexed by INSPIRE and arXiv up to 14 Oct 2026 (search log in the prior-work file linked below).

## How it was verified

**Equations.** The perturbation equations were re-derived symbolically from the metric. They agree with the
printed ones except for one potential, $V_{10}$, where the sign of the $e^{2q}$ term differs. The re-derived
equations reproduce the published static eigenvalue ($z = -21.80$) and the merger point
($(\lvert J\rvert/A)_{\mathrm{cr}} \approx 0.10063$); the printed $V_{10}$ gives $z = -146.16$ for the static case.
Reading the printed sign as a typesetting error is **our own re-derivation**; no one else has confirmed it, and it
says nothing about the original authors' code.

**Two independent numerical methods.** A dense finite-difference eigenvalue solver (complex eigenpairs kept,
Richardson-extrapolated) and an independent Evans-function shooting method agree to $\le 10^{-7}$ relative in
$\Omega^2$.

**Six checks**, each passed at both points. The thresholds were fixed before the verification run, but after
trial runs:

1. eigenvector localised in the box ($\le 10^{-6}$ of the norm beyond $l = \sinh(S)/3$);
2. box size ($S = 6$ vs $7.5$, $\le 10^{-6}$);
3. shooting outer radius ($S = 3.5$ vs $4.0$, $\le 10^{-6}$);
4. algebraic residual of the complex eigenpair ($\le 10^{-8}$);
5. residuals of two linearised field equations not used to build the system ($\le 10^{-4}$, falling at second order);
6. dense solver vs shooting.

## Scope and checks not yet done

- **Larger rotation.** Nothing is claimed as verified for $\lvert J\rvert/A > 0.142$. The point
  $\lvert J\rvert/A = 0.214$ is only a candidate: its shooting-radius check gave $4\times10^{-6}$ against the
  frozen $10^{-6}$. Exploratory tracking reached 0.40.
- **Near-extremal regime** ($\lvert J\rvert/A \to 0.5$): not computed.
- **Other sectors.** Odd parity was only screened; other perturbation sectors were not studied. No complete
  spectrum is claimed.
- **No time-domain evolution** was run.
- **Beyond this model.** Nothing is transferred to 4D, to other rotating wormholes or to other matter models.
  No astrophysical or observational implication is claimed.
- **No human expert check** of the derivation (including the $V_{10}$ sign), the numerical methods and thresholds,
  the reading of the literature, or whether the result is of interest to the field.

## Who did what

The GPT row is the owner's cross-conversation record (supplied 2026-10-14). The repository holds no GPT output
of its own.

| Contributor | Actual work |
|---|---|
| Repository owner (DD-Ching) | Chose and scoped the question; set the gate rules (reproduce the known modes, merger within ±1 %, residuals of the original equations, an independent method, freeze rules before the run), the reading rules and the corrections during the run (no fitting to results, keep full complex vectors, no extension toward extremality, wording). Passed GPT's comments on to Claude. Has not reviewed the physics as a domain expert. |
| Claude (Claude Code agents, this repository's sessions) | Main derivation, implementation, computation and documentation, including all code, runs and state files. Chose the numerical settings and thresholds after trial runs. Review and literature search were also done by Claude subagents. |
| GPT (separate conversation; owner's record) | Read the original paper, the GitHub code and the results. Helped shape the research scope and the verification design. Pointed out three issues that the owner relayed to Claude (minimum vs last step in `Shooter.root`, the grid differences near the merger point, and the NaN formatting). Did **not** independently re-derive the equations or re-run the numerics. |
| Human domain experts | **None so far.** |

## Files (pinned to commit a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c)

- Research summary:
  [`docs/state/units/X2_s6a_summary.md`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/docs/state/units/X2_s6a_summary.md)
- Reproduction package (equations, settings, full results table, rerun commands):
  [`docs/state/units/X2_s6a_repro.md`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/docs/state/units/X2_s6a_repro.md)
- Method and checks in detail:
  [`docs/state/units/X2_s6a_5d.md`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/docs/state/units/X2_s6a_5d.md);
  prior-work search:
  [`docs/state/units/X2_s6a_prior_work.md`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/docs/state/units/X2_s6a_prior_work.md)
- Code:
  [`scripts/x2_s6a_derive5d.py`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/scripts/x2_s6a_derive5d.py),
  [`scripts/x2_s6a_rot5d.py`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/scripts/x2_s6a_rot5d.py),
  [`tests/test_x2_s6a_rot5d.py`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/tests/test_x2_s6a_rot5d.py)
- Results (JSON):
  [`results/x2/s6a_rot5d.json`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/results/x2/s6a_rot5d.json),
  [`results/x2/s6a_equations.json`](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/results/x2/s6a_equations.json)
- Licence: [BSD-3-Clause](https://github.com/DD-Ching/jwst-anomaly-research/blob/a57b824e85a1dda2633caa1e0b2f1f245d4f1c2c/LICENSE)

The S6a files are unchanged between `a57b824` and the main branch at the time of writing (`3a80e18`).
