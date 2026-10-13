# X2 S2 design: one-event demonstration (GW150914), frozen before any comparison

Written 2026-10-13 before any echo amplitude was computed. Later changes are labelled post hoc.

## 1. The observational result used (Miani et al. 2023, arXiv:2302.12158v2)
- **Table II** (§IV.4), GW150914:
  - merger-ringdown h_rss h^MR = 16.0 × 10⁻²³ Hz^-1/2;
  - upper limit h^UL = 3.4 × 10⁻²³ Hz^-1/2;
  - ratio 0.21.
- h^UL is the **95 %-coverage upper limit on the h_rss injected inside the post-merger window (PMW)**, from a
  Feldman–Cousins-type confidence belt built with injections of their echo proxy.
- **The proxy:** two elliptically polarised sine-Gaussian pulses, f₀ = 140 Hz, τ = 20 ms, γ = 0.5, at
  t_echo = Δt_echo = 0.3 s, with the BBH's own inclination and sky position.
- **Table I's "h_rss at 50 % detection probability"** is a sensitivity figure, **not used** as a limit.
- **PMW for GW150914:** Table I predicts t_echo = 227 ms (> 200 ms), so §III.5 gives blind time 200 ms and a
  1 s window. The PMW is [0.2, 1.2] s after coalescence.
- **Validated range** (§IV.1):
  - the detection probability is unchanged when t_echo and Δt_echo vary, "as long as both pulses occur inside
    the analyzed time window";
  - Δt_PMW ∈ [0.2, 1.2] s and Δt_blind ∈ [0.05, 0.2] s;
  - morphology: sine-Gaussians at 80–400 Hz, a BBH-like single pulse and white-noise bursts give "almost
    identical" detection probability.
  - Not validated: trains of more than two pulses, and echoes much longer or broader than the proxies.
- **The h^UL/h^MR ratio** equals an upper limit on the relative amplitude A only "in the approximation that the
  merger-ringdown and each echo pulse share similar morphologies" (§IV.4).
- **Three quantities kept distinct:** h^MR is the integrated merger + ringdown h_rss; the ringdown-only h_rss
  is smaller; the peak strain h_max (their A = h₀/h_max) is a third quantity. Only h_rss ratios are compared
  here.

## 2. From the model to h₊, h× (what is known and what is assumed)
**Known (derived, S0/S1):**
- the geometry of the symmetric SV throat;
- the spacing Δt(a) between the two barriers;
- single-barrier |R(ω)|² and |T(ω)|² for **test fields**: scalar (s = 0) and electromagnetic (s = 1),
  ℓ = 2, from Eq. 6 of arXiv:1911.11823.

**Not known (no published result in our record):**
- the gravitational (spin-2) perturbation equations of SV, which depend on the matter that supports the
  throat (for example nonlinear electrodynamics plus a phantom scalar);
- that matter's own perturbations, and their coupling to h₊ and h×.

**Phenomenological assumptions P1–P4.** Every S2 number is a limit on this phenomenological model, **not** an
observational exclusion of the SV geometry.
- **P1.** ℓ = 2 gravitational waves obey a test-field barrier. Three proxies bracket it:
  - s = 1 (EM; the peak exactly at R = 3M);
  - s = 0 (scalar);
  - "RW-form", V = f[ℓ(ℓ+1)/R² − 6M/R³]: the Schwarzschild Regge–Wheeler form evaluated on the SV areal
    radius. This is ad hoc, with a barrier height near the true Schwarzschild spin-2 barrier.
- **P2. Source.** The wave sent inward from the light-ring region has the spectrum of the observed
  merger-ringdown times a **free amplitude η** (inward relative to outward). η is not fixed, so the result is a
  bound η ≤ η_max(a), not an exclusion of a.
- **P3. Spectral weight |h̃^MR(f)|².** Bracketed by:
  - (i) the ringdown Lorentzian of a remnant with M_det = 68.8 M☉ and χ = 0.69 (f_RD ≈ 251 Hz, τ ≈ 4 ms;
    GWTC-1, arXiv:1811.12907; source-frame M_f = 63.1 M☉, z = 0.09);
  - (ii) a flat band from 35 to 350 Hz.
  - Polarisation: the echoes inherit the ℓ = 2 ringdown polarisation (same inclination as the BBH), as in the
    proxy.
- **P4. Timing.**
  - Pulse k arrives at t_k ≈ k·Δt(a), counting from the light-ring excitation at about coalescence. The offset
    is ≲ 30 M ≈ 10 ms, which matters only within 10 ms of the PMW edges.
  - Each echo pulse takes one more round trip. Pulse k carries T·R^(2k−1) in amplitude: transmitted once, and
    reflected once at the far barrier plus twice per further round trip.
  - Pulses do not overlap (Δt ≥ 0.2 s ≫ pulse length), so their h_rss add in quadrature.

## 3. The comparison (computed after this file was committed)
F(a)² = Σ_{k: t_k ∈ PMW} ∫ W(f) |T|² |R|^(2(2k−1)) df / ∫ W(f) df, so that echo h_rss in the PMW =
η · F(a) · h^MR, and η_max(a) = 0.21 / F(a).

The comparison is applicable only where all four of these hold:
- (a) the first echo arrives at t₁ ≥ 0.2 s and the second inside the PMW. That means Δt ∈ [0.2, 0.6] s for two
  pulses, i.e. a − 2M within the range from `x2_echo_delay.a_for_delay`;
- (b) the echo-weighted mean frequency lies in 80–400 Hz;
- (c) at most about 2 pulses carry ≥ 90 % of the in-window h_rss (the validated case). More pulses are labelled
  "beyond the validated morphology";
- (d) the pulse duration is ≲ 0.1 s.

Outside these conditions no bound is quoted.

**Result classes:**
- η_max(a) < 1 for some a under all three proxies: the phenomenological SV echo model with an inward
  excitation comparable to the outward one is disfavoured there, conditional on P1–P4.
- η_max ≥ 1 everywhere: no constraint even at η ~ 1.
- Neither outcome is a statement about real remnants without P1–P4. Spin (χ = 0.69) is ignored (premise ii).
