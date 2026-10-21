# S3 hybrid copies in SN Refsdal light curves — result (2026-10-21)

Plan: `PLAN.md` (frozen before any fit). Numbers: `summary.json` (`scripts/s3_refsdal.py run`). Exploratory.

**Gates.** G1 passes: fitted delays relative to S1 are SX 385 d (F125W) / 370 d (F160W) against the
published 376.0 +5.6 −5.5 d (Kelly et al. 2023, arXiv:2305.06377 abstract), and S2–S4 −1 to 18 d (within ±30 d
of 0). G2 passes: r = 0.2 injected at δ = +150 d into SX F125W is recovered as 0.213 ± 0.057.

**Frozen decision.** 0 flagged lags (rule: r/σ ≥ 5 in both F125W and F160W at the same lag, with
|r₁₂₅ − r₁₆₀| < 2σ; `flagged()`, unit-tested). Per-target sensitivity r95 = max over δ of (r + 1.645σ), as
frozen: 0.21–1.69 over the eight target–filter scans (S2 0.30 / 0.30, S3 0.21 / 0.62, S4 1.69 / 0.57,
SX 1.20 / 1.05 for F125W / F160W). None reaches the plan's ≲ 0.05 branch; six of eight are ≥ 0.3. So the
Refsdal light curves do **not** give the r ~ 10⁻² reach that D-072's "Revisit if" anticipated. These are
formal sensitivities from one injection point, not limits.

**Single-band excesses (post hoc reading, not a frozen test).** F125W alone exceeds 5σ in S2 (9.3σ at
δ = +52 d, r 0.23), S3 (6.7σ at +64 d, r 0.17) and S4 (5.1σ at −336 d, r 0.71); F160W at the same lags gives
−0.5σ, 2.9σ and 1.2σ. The main-term fits in F125W have χ²/dof 3.85 (S2), 2.97 (S3), 2.11 (S4) against
1.40–1.47 in F160W; with σ rescaled by √(χ²/dof) the F125W peaks drop to 4.7σ, 3.9σ and 3.5σ. For S2 the
F125W r is positive at every δ > 0 and negative at every δ < 0 (`frac_lags_r_*` = 1.0 / 1.0), a broad
sign-flip pattern of a main-term shape mismatch (template time-stretch, differential microlensing or S1
template error), not a localised copy. No ordinary explanation is tested here; none of these is a candidate.

**Limitations.** Formal σ ignores S1 template noise and correlated residuals (χ²/dof up to 3.85); τ is fitted
on a 1-d grid without the copy term; the zero-flux-outside-span template is an ASSUMPTION; F105W (24 epochs)
not used; SN H0pe not attempted.

**Outcome.** S3 on SN Refsdal: no flag, sensitivity-limited at r ≳ 0.2–1.7 by the slow SN light curve and
main-term mismatch. S3 stays paused (D-072 addendum).
