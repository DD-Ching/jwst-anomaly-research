# Finite-energy two-way signals through an unstable wormhole (2026-10-28)

Owner-requested bounded research round (≤ 4 h wall clock, ≤ 4 cores, ≤ 1 GB download; about 50 MB of arXiv
sources were downloaded to the session scratchpad). Code: `scripts/wh_twoway.py`; tests: `tests/test_wh_twoway.py`;
every number below is in `results/wh_twoway/*.json` unless marked "not stored".

**Labels.** [lit] published and read in the paper body; [figure-read] read off a published figure; [sim] our
simulation (a model calculation, never an observation); [der] derived by arithmetic from stored numbers;
[assumption]; [hypothesis].

**Question.** In a published 4D wormhole model with explicit matter and dynamical equations: after an outbound
perturbation has passed through, can a *positive-energy return pulse that acts back on the geometry* still reach the
original sender? If so, over what verified range, and with what limits?

## 1. Has it been answered? No (checked this round)

Three different claims have to be kept apart:

| Claim | Status before this round |
|---|---|
| (a) A **test** light ray (no stress-energy) can return after an outbound pulse | Shown [lit]: CFK22 Fig. 3; Kain26 §3.4.2 (one ray, by eye); GGS09 Fig. 4 (a possible path with an unmodelled mirror) |
| (b) A **finite-energy, self-gravitating** return pulse reaches a fixed receiver | Not computed anywhere found (this round computes it) |
| (c) An **identifiable message** (receiver, noise, coding, bit rate) | Not computed anywhere found; not computed here |

- **CFK22** (Calhoun, Fay & Kain, arXiv:2210.04905v2, PRD 106, 104054) [lit]: return "light signals" are null
  geodesics $dr_\mathrm{null}/dt = \pm\alpha/\sqrt{A} - \beta^r$ (the last equation of their §II) integrated on the evolved 3+1 metric,
  launched from the peak of the matter pulse after it has crossed (their §V, Fig. 3). They carry no stress-energy,
  stop at the plot edge r = 5 r₀, and no receiver worldline or energy is computed.
- **Kain26** (Kain, "Dynamical wormholes", arXiv:2602.18231v1, CQG 43, 045014) [lit]: double-null review and
  evolutions of the Ellis–Bronnikov wormhole perturbed by **one** normal scalar pulse. Return is discussed only as
  a test null ray: a ray going up from (v,u) = (4,0) and turning right at (4,6) "will then cross back to the original
  side … never having crossed the event horizon at u ≈ 9.5" (§3.4.2). No second pulse.
- **Shinkai & Hayward** gr-qc/0205041 and the other papers in `docs/state/units/S0-2026-10-16.md` §1 were read on
  2026-10-16; none has a return pulse (SH02's maintenance pulse is a *ghost* pulse).
- **Citation check (this round)**: INSPIRE citers of 2210.04905, 2602.18231, 2503.07610 and 2507.00109 (14), citers
  of gr-qc/0205041 since 2020 (101, titles and abstracts) and full-text queries ("round trip", "two-way", "return
  signal/pulse") found no paper with a self-gravitating return pulse. Nearest: Hayward gr-qc/0110080 (2D dilaton
  gravity, one pulse); AdS bounds on information transfer (Maldacena et al. 1704.05333; Hirano et al. 1906.10715;
  Freivogel et al. 1907.13140), a different setting; Shirokov 2610.13206 and 2610.29941 (3+1 binary and rotating Ellis
  wormholes) without signals. INSPIRE full-text indexing may lag for October 2026.

## 2. The fixed problem (frozen before any return-pulse run, except the deviations listed at the end)

- **Model** [lit]: Kain26 §2–3 — 4D Einstein gravity, spherical symmetry, the symmetric massless Ellis–Bronnikov
  wormhole (throat radius x₀ = 1; G = c = 1; lengths, times and energies in units of x₀, x₀/c and x₀c⁴/G)
  supported by a massless ghost scalar. **[assumption] The ghost field is a premise of the model, not a known or
  manufacturable material.** The model joins **two separate asymptotically flat universes** (§6).
- **Signals**: two ordinary massless scalar fields, χ₁ (outbound) and χ₂ (return), minimally coupled; both enter
  $T_{\mu\nu}$, so the geometry includes the back-reaction of both legs (Kain26 eqs. 51, 53–54 with χ split in
  two; equation and figure numbers follow the arXiv v1 source).
- **Stations**: A (sender) on the + side and B (far side) on the − side, both on worldlines of constant areal radius
  **r = 2 x₀** (static observers before the perturbation; proper distance √3 x₀ from the throat).
- **Clocks**: each station's own proper time τ along its worldline, zero where the worldline meets the initial null
  surface. In the static era τ = t − √3 with t the Ellis Killing time, so the two clocks agree before the pulses.
- **Outbound** (fixed): χ₁ on u = 0 with $\chi_{1,v} = a_1\sin^2[\pi(v - v_1)]$ for v ∈ (v₁, v₁ + 1),
  v₁ = 2√3 + 0.5, a₁ set so that its initial-surface Misner–Sharp energy is **E_out = 1.0 × 10⁻⁴**.
- **Return** (fixed shape): χ₂ on v = 0 with $\chi_{2,u} = a_2\sin^2[\pi(u - u_1)]$ for u ∈ (u₁, u₁ + 1). Only
  its energy **E_ret ∈ {10⁻⁸ (test limit), 10⁻³, 10⁻², 10⁻¹}** and its timing (u₁) vary.
- **Energy**: the Kodama (Misner–Sharp) energy a field carries through a station's sphere, from
  $m_{,u} = -8\pi r^2 e^{-2\sigma} r_{,v} T_{uu}$ and $m_{,v} = -8\pi r^2 e^{-2\sigma} r_{,u} T_{vv}$, split per field
  and per direction; also the energy the station measures locally, $dE_\mathrm{loc} = dE_K/\sqrt{1 - 2m/R}$.
  In the static era $E_K = 0.866\,E_\mathrm{loc}$ at r = 2.
- **Success rule** (`passes`, unit-tested): the χ₂ Kodama energy leaving through A's sphere before A's proper time
  **τ_dead = 12 x₀/c** is at least **50 %** of the χ₂ Kodama energy that entered through B's sphere. A run is scored
  only if both stations stay untrapped with √(1 − 2m/R) > 0.5 up to the deadline, A's track reaches the deadline and
  B's track covers the whole return pulse (`run_valid`).
- **Timing labels**: B's emission time τ_B,emit and A's reception time τ_A,recv are the 50 % points of the
  cumulative flux. "Reply-compatible" means χ₂'s 10 % point at B comes after χ₁'s 90 % point at B. The return pulse
  is pre-loaded in the initial data, so every case is **two-way propagation**, not a reply (§5).

**Deviations from the frozen plan** (all labelled in the code):
1. **τ_dead 14 → 12 (post hoc).** Frozen at 14.0. Five test runs with return pulses (n = 200–800; u₁ = 9 all
   passed, u₁ = 10 gave 0.60 and 0.53 at n = 400 and 800 with τ_dead = 14) showed the Misner–Sharp bookkeeping on A's
   track failing beyond τ_A ≈ 13.4 at n = 200 (≈ 12.2 at n = 800 for a 10⁻⁶ tolerance): A's late track hugs the
   event horizon, where u-gradients grow exponentially and the uniform grid no longer resolves them. The deadline was
   moved to 12.0 before the scans, on that numerical check, but after those outcomes had been seen. Sensitivity to
   11.5 and 12.5 is reported; 14 is not reported because the solution there is not verified.
2. **E_ret = 0.3 and 0.5 and u₁ = 9.75, 10.25 added (post hoc)** after the first n = 800 scan showed that E_ret ≤ 0.1
   barely moves the boundary. Still only energy and timing vary.
3. **Validity rule added after the independent review** (finding 3): stations must stay untrapped up to the
   deadline and A must reach it. In the stored scans the station check covered the whole track; one run (n = 400,
   E_ret = 10⁻³, u₁ = 9.75) failed it only after the deadline (A's sphere becomes trapped at τ_A ≈ 17.2 in the
   under-resolved late track). Re-run with the deadline-limited check it is valid, √(1 − 2m/R) ≥ 0.866 up to
   τ_A = 12 (`validity_recheck.json`, which replaces its scan entry).

## 3. Method and checks

**Reuse.** CFK22, Kain26 and XCY25 do not state that their codes are public; Shinkai's `wh_dual` is public but has
no licence and was not downloaded (S0-2026-10-16 §6). So the round reimplements **only** Kain26's published
double-null equations (eqs. 3–4, 51, 53–54) and his second-order predictor–corrector cell update (§2.3), marched one
anti-diagonal (a spacelike slice u + v = const) at a time with numpy: `scripts/wh_twoway.py`, about 700 lines
including diagnostics. No general numerical-relativity framework was built. Radial null rays are the coordinate
lines u = const and v = const, so no geodesic integration is needed for the test-ray comparison.

**Checks passed before any return-pulse result** (`results/wh_twoway/checks.json`; tests in
`tests/test_wh_twoway.py`):

| Check | Result |
|---|---|
| Static wormhole kept (u, v ≤ 4) | max \|r − r_static\| = 7.39 × 10⁻⁴, 1.85 × 10⁻⁴, 4.62 × 10⁻⁵ at n = 40, 80, 160: second order |
| Linear instability vs GGS09 [lit] (e-folding 0.846 linear; 0.836 from their nonlinear code, as recorded in S0-2026-10-16 §6) | 0.842, 0.840, 0.840 at n = 100, 200, 400 (fit of 1 − r_throat over t ∈ [5, 6.5], Kain26 pulse) |
| Kain26 Fig. 9 case (A = 0.002, v ∈ (1, 3)) | proxy for the event horizon: the last u whose ray is beyond r = 3 at v = 20 (a lower bound on the horizon's u): 9.40, 9.20, 9.15 (n = 100, 200, 400). Kain's figure shows the horizon near u ≈ 9.2 [figure-read]; his text says "u ≈ 9.5". Contours, horizons and singularity agree with the figure by eye (`fig_benchmark.png`) |
| Null constraints, Misner–Sharp bookkeeping | converge at second order on test grids (unit tests); per-run values in the scan files |
| Resolution of the production runs | n = 400, 800 (all cases) and 1600 (boundary cases) grid points per x₀; §4 |

**Numerical hazard found (and why E_out = 10⁻⁴).** The static wormhole is unstable, so truncation error itself
seeds the unstable mode from t = 0, with the sign of *expansion*. Launched from r = 2 (it reaches the throat at
t ≈ 4.5 instead of 2), Kain's own pulse energy 5.4 × 10⁻⁵ gives throat radius 1.109, 1.012, 0.984 at t = 10 for
n = 100, 200, 400: the code's verdict flips from expansion to collapse with resolution [sim]. With E_out = 10⁻⁴
the collapse is kept at n = 200 and 400 (0.991, 0.962 at t = 10), so E_out = 10⁻⁴ was fixed before any
return-pulse run, and every decisive number below is checked at n = 1600. The same caveat may apply to published
runs at n ≤ 400 with weak pulses that reach the throat late (inferred from this code's truncation error only).

## 4. Results [sim]

All numbers are model outputs (`results/wh_twoway/`: `scan_n400.json`, `scan_n800.json`, `converge_n1600.json`,
`boundary.json`). Energies are in x₀c⁴/G, times in x₀/c. "K" = Kodama energy through the station sphere, "local" =
energy the station observer measures.

**4.1 Outbound leg alone (fixed).** χ₁ crosses A inward at τ_A = 1.010 carrying 9.46 × 10⁻⁵ K (1.09 × 10⁻⁴ local)
and reaches B at τ_B = 4.468 (50 %; 90 % by 4.711) with 8.27 × 10⁻⁵ K (n = 1600). It makes the wormhole collapse to
a black hole. The **no-return-load reference**: a test null ray leaving B inward reaches A by τ_A = 12 only if it
leaves before τ_B = 7.656, 7.556, 7.533 (n = 400, 800, 1600).

**4.2 Reproducible pass/fail set** (`boundary.json` → `cases`; reproduce one row with
`python scripts/wh_twoway.py scan --n 1600 --e-ret 0.5 --u1 9.5 10.0 --out x.json`):

| n | E_ret | u₁ | τ_B emit | E_emit K / local | E_recv K / local | reflected K | received fraction | τ_A recv | round trip τ_A | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| 1600 | 10⁻⁸ | 9.5 | 6.541 | 8.50e-09 / 9.81e-09 | 5.72e-09 / 6.42e-09 | 5.4e-10 | 0.673 | 10.187 | 9.177 | pass |
| 800 | 10⁻⁸ | 10.5 | 7.538 | 8.39e-09 / 9.66e-09 | 1.45e-09 / 1.59e-09 | 4.9e-10 | 0.173 | 11.611 | 10.601 | fail |
| 1600 | 0.1 | 9.5 | 6.538 | 0.0851 / 0.101 | 0.0548 / 0.0604 | 0.0039 | 0.643 | 10.182 | 9.172 | pass |
| 800 | 0.1 | 10.5 | 7.535 | 0.0839 / 0.0994 | 0.0141 / 0.0154 | 0.0039 | 0.168 | 11.605 | 10.596 | fail |
| 800 | 0.5 | 8.0 | 5.023 | 0.434 / 0.604 | 0.309 / 0.327 | 0.0053 | 0.712 | 8.523 | 7.514 | pass |
| 1600 | 0.5 | 9.5 | 6.522 | 0.427 / 0.591 | 0.230 / 0.244 | 0.0051 | 0.540 | 10.161 | 9.151 | pass |
| 1600 | 0.5 | 10.0 | 7.021 | 0.424 / 0.587 | 0.159 / 0.169 | 0.0050 | 0.375 | 10.807 | 9.797 | fail |

Every scanned case is reply-compatible (χ₂ reaches B after χ₁ has passed it). The round trip is measured on A's
clock from χ₁'s emission to χ₂'s reception. The part of the emitted energy that is neither received nor reflected
(0.026 K, 31 % of E_emit, in the pass case E_ret = 0.1, u₁ = 9.5 at n = 1600 [der]) is trapped, arrives after
τ_A = 12, or is exchanged with the time-dependent geometry; this run does not separate the three (not counted).

**4.3 Escape range against return energy.** The B emission time at which the received fraction crosses 50 %
(`boundary.json`; linear interpolation; a shape-preserving cubic differs by ≤ 0.007):

| E_ret | n = 400 | n = 800 | n = 1600 | h²-extrapolated | last change | shift vs. test limit, n = 800 / 1600 [der] |
|---|---|---|---|---|---|---|
| 10⁻⁸ (test) | 7.139 | 7.009 | 6.975 | 6.964 | 0.033 | 0 |
| 10⁻³ | 7.138 | 7.008 | — | — | — | −0.001 / — |
| 10⁻² | 7.134 | 7.003 | — | — | — | −0.006 / — |
| 0.1 | 7.087 | 6.948 | 6.918 | 6.907 | 0.031 | −0.061 / −0.058 |
| 0.3 (post hoc) | 6.962 | 6.830 | 6.803 | 6.794 | 0.028 | −0.179 / −0.173 |
| 0.5 (post hoc) | 6.822 | 6.688 | 6.657 | 6.647 | 0.031 | −0.321 / −0.319 |

- **The return pulse reaches A for every energy tried, if B sends it early enough.** With E_ret ≤ 0.5 x₀ (up to
  5 000 times the outbound energy), a pulse sent right after reception (τ_B ≈ 5.0) delivers 71–86 % of its Kodama
  energy (n = 800, u₁ = 8.0).
- **Its own energy shortens the escape range only modestly.** Measured from the outbound's 90 % point at B
  (τ_B = 4.711), the 50 % escape range is 2.26 x₀/c in the test limit and 1.95 x₀/c at E_ret = 0.5 (n = 1600)
  [der]. The shift grows roughly in proportion to E_ret: −0.006 at 10⁻², −0.06 at 0.1, −0.32 at 0.5.
- **Energy criterion vs. light ray.** Even a test pulse has an earlier limit (τ_B = 6.975 at n = 1600) than a test
  null ray (7.533): the pulse is spread over Δu = 1, and its later part is lost to reflection, trapping or late,
  redshifted arrival before the 50 % mark is reached.
- **Deadline sensitivity** (n = 1600): moving τ_dead to 11.5 or 12.5 moves every boundary by −0.05 / +0.05 (e.g. test
  limit 6.926 / 7.024; E_ret = 0.5: 6.602 / 6.703); the energy-induced shifts are unchanged within 0.01.

**4.4 An energetic return closes the door behind it.** The pulse's own energy reaches the throat after the pulse and
advances the event horizon. After a return sent at u₁ = 8.0 (n = 800), the last B time from which a later test ray
still reaches A by τ_A = 12 is 7.556 (E_ret = 10⁻⁸), 7.514 (10⁻³), 7.217 (10⁻²), 6.202 (0.1), 5.632 (0.3) and 5.395
(0.5), against the return's own emission at τ_B = 5.02–5.04. So a 0.5 x₀ return gets through, but anything B sends
more than ≈ 0.37 x₀/c after it does not (`scan_n800.json`, `test_ray`).

**4.5 Does numerical error change the verdict?**
- The boundary converges at about second order: the ratio of successive changes (400→800 over 800→1600) is
  3.9 (test), 4.5 (0.1), 4.8 (0.3) and 4.4 (0.5) against 4 expected. Remaining uncertainty ≈ 0.03, the last change
  (n = 1600 differs from the h² extrapolation by ≤ 0.011).
- Shifts caused by the return energy agree between n = 800 and 1600 to ≤ 0.006, so the energy effect at 0.1–0.5 is
  resolved. Below 10⁻² the shift (≤ 0.006) is within the uncertainty: no effect is claimed there.
- **Verdicts within ≈ 0.15 of the boundary depend on resolution at n = 400**: E_ret = 0.5, u₁ = 9.75 passes at
  n = 400 (0.516) and fails at n = 800 (0.476) and 1600 (0.466). Cases ≥ 0.4 x₀/c from the boundary (u₁ ≤ 9.5 pass,
  u₁ ≥ 10.5 fail) keep their verdict at every resolution.
- The dominant error is in the fixed outbound collapse, which the static wormhole's truncation error delays by
  seeding the expanding mode (§3): at n = 400 the boundary lies 0.175 x₀/c later than the extrapolation [der].

**4.6 Energy-accounting checks.**
- Misner–Sharp bookkeeping along A's track up to the deadline: |Σ pieces − Δm| ≤ 3.1 × 10⁻⁵ in all n = 1600 runs.
- Null-constraint residual in A's causal past up to τ_A = 11: ≤ 1 × 10⁻³ (n = 1600). In the last unit before the
  deadline A's track enters a strip within ≈ 0.1 in u of the event horizon where the uniform grid does not resolve
  the exponentially steepening u-gradients; there the residual reaches O(10⁵) for E_ret ≥ 0.3 (diagnosed on one
  n = 800 run; not stored). The energy received in that strip converges: the boundary change from τ_dead 11.5 to 12
  is 0.058 / 0.055 (E_ret = 0.5, n = 800 / 1600).
- Unperturbed wormhole (no outbound pulse, `static_return.json`): transmitted 0.877, reflected 0.136, re-entering at
  A 0.025 of the emitted energy (n = 800); the balance closes to 1.2 % (0.9 % at n = 400). The residual does not
  shrink with resolution and is not explained here; it bounds the accounting at about 1 % of E_emit.
- Pre-loading imprint (review item 2, `imprint.json`): the χ₂ pulse disturbs the ghost field on its way to B, which
  changes the ghost inflow through B during the χ₂ crossing by +0.9 % (E_ret = 0.1) and +0.5 % (0.5) of the χ₂
  energy (n = 400 and 800 agree). A transmitter at B would not cause this; the effect is small at these energies.

## 5. What this is and is not

- **Two-way propagation, not a reply.** The return pulse is pre-loaded on the initial null surface outside B and
  only *crosses* B's sphere after B has received the outbound pulse ("reply-compatible" timing in every scanned
  case). For everything inside B's sphere this is the same as B emitting that influx at those times (spherical
  symmetry and causality), but no trigger, detector, transmitter or energy reservoir at B is modelled. Calling it a
  "reply" would need a causally triggered source with its own energy budget.
- **No communication rate.** There is no receiver, noise or coding model, so no bit rate or "successful call" is
  claimed. The received energy only bounds what a receiver could have to work with.
- **One model, one configuration.** Symmetric massless Ellis–Bronnikov, stations at r = 2 x₀, E_out = 10⁻⁴ x₀, unit
  sin² pulses, s-waves only. Other station radii, outbound energies, pulse shapes or other wormhole models were not
  computed and are not implied. A failure here is not a statement about all wormholes.
- **The window is not the e-folding time.** The 0.846 x₀/c e-folding time is a growth rate; the computed window is
  the result of the full nonlinear evolution for this outbound energy and these station positions.

## 6. Distance from the final goal

- **Not the same universe.** This spherically symmetric wormhole joins two separate asymptotically flat regions.
  There is no external path between A and B to compare with, so **no faster-than-light or shortcut statement is
  possible** from this model. An intra-universe version needs a non-spherical global construction that no paper
  we found evolves dynamically.
- **The supporting matter is assumed.** The ghost scalar that holds the throat open is a premise. Nothing here
  says it exists or could be made.
- **The window is short and closes.** In this configuration the usable window at B is of order x₀/c and the
  wormhole then becomes a black hole; nothing here keeps it open.
