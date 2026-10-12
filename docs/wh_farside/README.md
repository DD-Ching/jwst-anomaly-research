# Far-side light through a wormhole throat: a period–astrometry test on Sgr A* flares

> **Revised 2026-10-12 after owner review.** Two claims are withdrawn: that the flares are excluded as far-side
> sources "for every allowed a", and that astrometric systematics cannot enlarge the loop. Gaps are now labelled
> as summary-level standardised differences. Geometry and numbers are unchanged; §4–§7 were rewritten.

Owner brief 2026-10-12: a bounded, high-risk exploration unit. Find one new prediction that existing data can test
and that separates a specific channel model from ordinary explanations.

**Files:**
- Code: `scripts/wh_farside.py`, which runs in under a second.
- Tests: `tests/test_wh_farside.py`.
- Output: `results/wh_farside/check.json`.

**Labels:**
- [lit]: a published number.
- [pred]: a model prediction of a conditional, unestablished geometry.
- [ours]: our calculation.

Nothing here is evidence that any wormhole exists.

## 1. What is new relative to the record (nothing below is recomputed)
| Earlier work | What it tested | Difference from this unit |
|---|---|---|
| W1–W6 (D-047, `docs/exotic_lensing.md`), W3 MOA-II | Ellis or negative-mass objects **as lenses** of sources on our side | Here the source is on the **other side** of the throat, and its light crosses the throat |
| S1 (D-071; #167, #169), S3 (D-072) | A second, deformed copy of a source in our universe | Here there is no copy. The one image comes from the other region |
| E1 (D-074), D1 (D-073) | Cross-messenger coincidences; distance self-consistency | These do not involve horizon-scale imaging |
| Hypothesis rounds A/B (D-069, D-075, D-076); goal rounds 1, 4 and 10 (#152, #156, #166) | Theory: energy-condition theorems; published traversable constructions give no shortcut (round 1, Gao–Wald, Maldacena et al.); beyond-Horndeski wormhole stability fails (round 4) | Those concern channel feasibility. This unit asks only an observational question: do we already see light from the other side? A "yes" would still say nothing about a shortcut |
| B-on-A2 (round 1), B-on-R2A | The EHT shadow size against horizonless objects | Here the EHT bound is applied to the Simpson–Visser critical curve, together with flare astrometry |

## 2. Models compared (two), and the one chosen
**Model 1: Simpson–Visser (SV) traversable wormhole** (arXiv:2408.07350, Chen et al. 2024).
- Metric: ds² = −f dt² + dr²/f + R² dΩ², with f = 1 − 2M/R and R = √(r² + a²). It is a wormhole for a > 2M.
- Geodesics are explicit. 2408.07350 ray-traces hot spots on the far side: the images lie inside the critical
  curve, the centroid motion is smaller, and the light curve can have extra peaks.
- No public code was found in the paper. We do not reproduce its figures.

**Model 2: Dai & Stojkovic flux leakage** (arXiv:1910.00429). A mass on the far side perturbs orbits on our side
(S2, binaries). Simonetti et al. (arXiv:2007.12184) already set limits this way. It concerns gravity, not light
from the other side, and an orbit residual is degenerate with ordinary extended mass. **Not chosen.**

**Premises that are assumed, not established:**
- (i) Sgr A* is a traversable SV throat. SV is an ad hoc metric: the matter that supports it, how it forms and
  whether it is stable are not established.
- (ii) An accretion flow on the far side produces NIR flares.
- (iii) The wormhole is reflection-symmetric (the same M on both sides).
- (iv) Far-side spots move on Keplerian orbits.
- (vi) The flare polarisation period equals the orbital period, also for light that crossed the throat.
  Polarisation transport through the throat was not computed.
- (v) The "other side" is another region, and nothing says it is a shortcut in our universe. Goal round 1 found no
  published construction that gives one. **This test is unrelated to usable communication.**

## 3. Joint prediction (shared parameters M, a and R0; parameters are not tuned per observable)
- **Inputs used to fix parameters:**
  - M = (4.297 ± 0.012) × 10⁶ M☉ from stellar orbits [lit];
  - R0 = 8.277 kpc [lit];
  - the flare **polarisation** period P = 60 ± 3 min [lit, GRAVITY 2023, arXiv:2307.11821].
- **Step 1, period → orbit:** Kepler's law holds exactly in areal radius for SV (Ω² = M/R³, [ours]).
  - If the polarisation period equals the orbital period (premise vi), P fixes the far-side orbit at
    R_e = 9.02 ± 0.30 M.
  - A far-side Kepler orbit then needs a ≤ R_e.
- **Step 2, prediction for data not used above:** the astrometric loop. Three quantities must be kept apart.
  - **(a) Geometric bound on the true far-side images** [pred, ours]: every photon that crosses the throat has
    b < b_c(a). That is √27 M for a ≤ 3M, and a/√(1 − 2M/a) for a > 3M; a numerical check is in the tests. All
    far-side image positions therefore lie inside b_c·θ_g, with θ_g = 5.12 µas: **26.6 µas for a ≤ 3M**.
  - **(b) Centroid of the far-side images alone:** an intensity-weighted mean of points inside that disk, so it
    is also inside b_c·θ_g.
  - **(c) The observed loop radius:** a pipeline estimate from the centroid of **all** light, after phase
    wrapping and fitting. It is not bounded by (a) unless the conditions in §5 hold.
- **Ordinary same-side control:** loop radius ≈ (R_e + 1)·θ_g = **51.3 ± 1.5 µas** [ours, GRAVITY's Appendix A
  approximation].

## 4. Result against published summary data [lit numbers; ours comparison]
**What the 61 ± 9 µas is** (GRAVITY 2023, Sect. 3 and Fig. 3):
- One "mean astrometric loop size" from four flares (2018–2022).
- The flares were time-shifted to the solution that best represents a sinusoid in both sky coordinates. All
  points were wrapped at a 60-min period and binned in 5-min bins; the error bars are the standard deviation of
  the binned data.
- The ± 9 is not defined further in the text.
- The phase alignment, the averaging, any correlation between flares and the hot-spot model are **not**
  propagated into our comparison. Per-flare astrometry is not tabulated publicly.

**Every gap below is a standardised difference (observed − value)/9 µas on this one summary number. None is a
full exclusion significance.**

**4a. Without the EHT calibration** (astrometry plus the period):

| a / M | 2.1–3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| far-side bound (µas) | 26.6 | 29.0 | 33.1 | 37.7 | 42.4 | 47.3 | 52.3 |
| summary gap, observed above the bound | 3.8 | 3.6 | 3.1 | 2.6 | 2.1 | 1.5 | 1.0 |

- **Same-side control:** gap 1.1, consistent.
- **a ≤ 7.1M:** gap ≥ 2.
- **a ≈ 7.1–9.0M:** the gap is < 2, so **this summary does not distinguish a far-side origin there** (1.5 at
  8M, 1.0 at 9M).
- **a > 9.0M:** a far-side Kepler orbit with a 60-min period cannot lie outside the throat. That limit depends on
  premise (vi) and on the 3-min period error. It is not unconditional.

**4b. If EHT's shadow calibration transfers to this case:**
- The EHT shadow-size bound (δ = −0.04 +0.09/−0.10 Keck, −0.08 ± 0.09 VLTI; arXiv:2311.09484) would limit the SV
  metric to a ≲ 4.1–4.4M at 2σ.
- In that range the summary gap is ≥ 3.4.
- **What EHT covered:**
  - its ring-to-shadow calibration library includes Kerr and non-Kerr GRMHD and semi-analytic images: JP and
    Kerr-Sen metrics, a boson star, a dilaton black hole, and one semi-analytic flow in a **spinning traversable
    (Teo) wormhole** whose throat equals the Kerr horizon radius;
  - its metric tests include a **non-spinning Morris–Thorne wormhole** with one specific redshift function
    (δ ≈ −0.48, ruled out).
- **What EHT did not cover:**
  - the Simpson–Visser metric;
  - emission that reaches us **through** the throat from the other side;
  - whether the ring still traces b_c in that configuration.
- So 4b is conditional on the calibration transferring, which is not verified.

**Withdrawn:** "the four flares are not far-side sources for any allowed a". Without the EHT transfer, a ≈ 7–9M
is not distinguished by this summary.

## 5. Conditions, mimics and systematics
**When can (c) be compared with (a)?** If all emitting light (flare and static) lies inside a disk of radius ρ
around the hole, the total-light centroid stays in that disk, and any fitted loop radius is ≤ ρ up to noise. This
fails in three ways:
- **Off-centre static emission or confusion.** If the static component's centroid sits a distance d from the
  hole, the total-light centroid moves inside the convex hull of the far-side disk and that point.
  - Flux changes then drive centroid motion along that line, with half-extent up to (ρ + d)/2.
  - Reaching 61 µas with ρ = 26.6 µas needs d ≳ 95 µas, about 19 R_g, with matching flux modulation. That motion
    is roughly linear rather than a loop, but after phase wrapping and fitting this was not checked.
  - Static flux that is **not** central is therefore a real loophole. It is excluded only if the static emission
    is centred within the disk, or its offset and flux are measured.
- **Astrometric errors.** Phase-referencing errors, aberrations, noise and the wrapping/fitting step can bias a
  fitted loop radius, including upward. A fitted radius of a noisy trajectory is biased high.
  - **Withdrawn:** "systematics cannot enlarge the loop".
  - The bound is unaffected only if the astrometric error model is known and propagated, which needs per-flare
    data.
- **Mimics of a small loop** (dilution by central static flux, sub-Keplerian motion, pattern motion): these
  shrink the observed loop. A small loop would therefore not identify a far-side source either.
- **Spectrum and redshift:** the local gravitational redshift and orbital speed are the same on both sides of a
  symmetric SV metric. The observed spectrum also depends on emission direction, Doppler factors and image
  weights along the different light paths, which were not computed. **The discriminating power of redshift is
  not established**, which is not the same as having none.
- **Magnification:** not computed.

## 6. Observability and budget
- **Data:**
  - published GRAVITY flare summaries (2018, arXiv:1810.12641; 2023);
  - EHT δ (arXiv:2311.09484);
  - four paper PDFs read, in scratch only.
  - No archive data.
- **Required:** per-flare loop radii to about ±10 µas, plus the static-flux fraction and offset. GRAVITY reaches
  this precision for bright flares, but the reduced per-flare data are not public.
- **Search accounting:** 2 models considered, 1 observable, 1 published summary. No search was run and no choice
  was tuned after seeing the result.

## 7. Decision
**Paused**, in a limited sense: **this observable is not yet sufficient to identify, or reliably exclude, a
far-side origin.**
- At summary level, the far-side origin is disfavoured for a ≲ 7M, or for all a if EHT's calibration transfers.
- These conclusions rest on one averaged number, a central-static-flux condition and the premise P_pol = P_orb.
- Not excluded is not support. The pause concerns this observable only, not the wormhole direction as a whole.

Inputs that would make the comparison decisive against a stated ordinary model (for example GRAVITY's same-side
orbiting hot spot plus a measured static component):
- (a) Per-flare astrometry with the static-flux fraction and offset, and the astrometric error model. This needs
  GRAVITY reduced data, or days of ESO raw-data reduction.
- (b) A ray-traced far-side prediction for the SV metric:
  - the image-weighted centroid;
  - the light curve (2408.07350's extra peaks);
  - the spectrum (Doppler weights);
  - polarisation transport through the throat, which would test premise (vi).
- (c) A check that EHT's ring-to-shadow calibration holds for through-throat emission.

Not done this round is not the same as physically impossible.
