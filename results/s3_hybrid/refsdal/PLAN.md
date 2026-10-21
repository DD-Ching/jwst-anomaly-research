# S3 hybrid copies in SN Refsdal light curves (D-072 "Revisit if": lensed transient) — frozen plan

Written 2026-10-21 15:25Z before any fit (cloud routine). Exploratory: the data are used once, no held-out set.

**Hypothesis (D-069 A1 P2, hypothesis only).** Image j of SN Refsdal carries a faint copy of the SN at a lag
δ relative to its own signal: F_j(t) = m F_i(t − τ_ij) + a F_i(t − τ_ij − δ), r = a / m.

**Ordinary explanations.** Microlensing of individual images (slow, achromatic to first order); photometric
covariance from the shared host-subtracted difference imaging; template error (S1 interpolated across gaps);
the SN's own slow (SN 1987A-like) light curve makes copies at |δ| ≲ the rise time collinear with the main term.

**Data (observed).** CDS J/ApJ/948/93 (Kelly et al. 2023, ApJ 948, 93): `table3` HST F105W / F125W / F160W
difference-image fluxes of S1–S4 and SX (805 rows), `images` positions and magnifications. Streamed via
TAPVizieR, not stored (manifest line only).

**Method (derived).** Template F_1 = S1 per filter, linearly interpolated, zero outside its sampled span
(ASSUMPTION: pre-explosion flux 0). For each target j ∈ {S2, S3, S4, SX} and filter F125W, F160W
(the two long-sampled bands): fit τ, m, c (constant) by least squares on a τ grid (step 1 d), then on a
δ grid (step 2 d, 30 ≤ |δ| ≤ 400 d) fit a with τ fixed; r(δ) and its formal σ from the weighted fit.

**Gates (frozen).**
- G1 (tool on a known case): the fitted SX − S1 delay lies within 376 ± 30 d (Kelly et al. 2023's
  headline delay is ≈ 376 d; value to be checked against the paper before any merge) and S2–S4 within ±30 d of 0.
- G2 (injection): an injected copy r = 0.2 at δ = +150 d into SX recovers r within 0.1.

**Decision rule (frozen, ASSUMPTION).** A lag is flagged if r(δ)/σ_r ≥ 5 in both F125W and F160W with
consistent r (difference < 2 σ); flags go to /vet-candidate, not interpretation. Otherwise report the
per-target 95 % sensitivity r95(δ) = max over δ of (r + 1.645 σ_r) as a sensitivity, not a limit (D-072
lesson: window maxima are not limits without injection calibration).

**What changes the plan.** G1 fail → tool not validated, stop and record. Flag → vet. Null with
r95 ≲ 0.05 → S3 gains a real constraint for a transient; r95 ≳ 0.3 → S3 stays sensitivity-limited, close.
