# Far-side light through a wormhole throat: a period–astrometry test on Sgr A* flares

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
- (v) The "other side" is another region, and nothing says it is a shortcut in our universe. Goal round 1 found no
  published construction that gives one. **This test is unrelated to usable communication.**

## 3. Joint prediction (shared parameters M, a and R0; parameters are not tuned per observable)
- **Inputs used to fix parameters:**
  - M = (4.297 ± 0.012) × 10⁶ M☉ from stellar orbits [lit];
  - R0 = 8.277 kpc [lit];
  - the flare **polarisation** period P = 60 ± 3 min [lit, GRAVITY 2023, arXiv:2307.11821].
- **Step 1, period → orbit:** Kepler's law holds exactly in areal radius for SV (Ω² = M/R³, [ours]). P therefore
  fixes the far-side orbit at R_e = 9.02 ± 0.30 M. This also requires a ≤ R_e, since the orbit must lie outside
  the throat.
- **Step 2, prediction for data not used above:** the **astrometric** loop radius.
  - Every photon that crosses the throat has impact parameter b < b_c(a): √27 M for a ≤ 3M, and
    a/√(1 − 2M/a) for a > 3M. This is the potential-barrier bound [ours]; a numerical check is in the tests.
  - So the far-side image, and hence its centroid, lies inside b_c·θ_g, with θ_g = 5.12 µas. That is
    **≤ 26.6 µas for a ≤ 3M, whatever the orbit radius** [pred].
- **Ordinary same-side control:** loop radius ≈ (R_e + 1)·θ_g = **51.3 ± 1.5 µas** [ours, GRAVITY's Appendix A
  approximation].

## 4. Result against published data [lit numbers; ours comparison]
- **Observed:** mean loop radius of the four astrometric flares = **61 ± 9 µas** (GRAVITY 2023).
- **Ordinary control:** consistent, 1.1σ above the same-side prediction.
- **Far side** (σ = the published 9 µas):

| a / M | 2.1–3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| far-side bound (µas) | 26.6 | 29.0 | 33.1 | 37.7 | 42.4 | 47.3 | 52.3 |
| observed above the bound | 3.8σ | 3.6σ | 3.1σ | 2.6σ | 2.1σ | 1.5σ | 1.0σ |

- **The constraints on a:**
  - Astrometry excludes a far-side origin of these flares for **a ≤ 7.1M at 2σ** (a ≤ 5.2M at 3σ).
  - The period excludes **a > 9.0M**: no Keplerian 60-min orbit lies outside the throat.
  - The EHT shadow-size bound (δ = −0.04 +0.09/−0.10 Keck, −0.08 ± 0.09 VLTI; arXiv:2311.09484) independently
    limits the SV metric itself to **a ≲ 4.1–4.4M at 2σ**. That is our application: EHT's ring-to-shadow
    calibration was not made for this wormhole.
- **Combined:** in the symmetric SV model, **the four GRAVITY flares are not far-side sources for any allowed a**.
  For a ≤ 3M, where the SV shadow equals Schwarzschild's, the exclusion is at 3.8σ.
- **For asymmetric wormholes** (different mass on the far side): the period no longer fixes the orbit, but the
  image bound b < √27 M still holds whenever the near side has its photon sphere at 3M. The 3.8σ exclusion holds
  there too.
- **Caveat:**
  - The test uses one averaged loop radius from four flares. Per-flare astrometry is not tabulated publicly.
  - Our implementation of GRAVITY's Appendix A mass estimator returns 6.6 × 10⁶ M☉ for 61 µas, against the
    paper's quoted 5.1 ± 1.6. The bound test does not use this estimator.

## 5. Can ordinary physics mimic the far-side signature?
The far-side signature is "loop ≤ b_c·θ_g at a 60-min period". It **can be mimicked**:
- **Dilution** by the static quiescent flux shrinks the centroid motion. It would need about a 50 % static fraction.
- **Sub-Keplerian motion:** a spot at R ≈ 4M moving at about 0.3 of the Keplerian speed has P ≈ 60 min.
- **Non-orbital pattern motion.**
- **Astrometric systematics** (phase-referencing).

So a **small** loop would **not** identify a far-side source, and as a detection test this observable has **no
unique discriminating power**.

In the other direction it is robust. None of these mimics can make a far-side image **larger** than b_c·θ_g. A
large observed loop therefore excludes a far-side origin regardless of the mimics. That is the only conclusion
drawn here.

Redshift and magnification: a Keplerian far-side orbit in a symmetric SV metric has the same Doppler and
gravitational shifts as the same-side orbit, so the redshift adds no discrimination. Magnification was not
computed.

## 6. Observability and budget
- **Data:** published GRAVITY flare summaries (2018, arXiv:1810.12641; 2023) and EHT δ. Downloads: three paper
  PDFs, about 10 MB. No archive data.
- **Required precision:** a per-flare loop radius to about ±10 µas, which GRAVITY reaches for bright flares. The
  expected far-side signal is ≤ 27 µas against about 50 µas for the same side.
- **Search accounting:** 2 models considered, 1 observable and 1 published data summary. No search was run, so
  there was no unseen validation set to keep. No threshold or choice was tuned after seeing the result.

## 7. Decision
**Pause this line.** The existing data already exclude a far-side origin for the observed flares within this
model. A positive far-side identification would not be unique, because of dilution, sub-Keplerian motion and
astrometric systematics.

Minimum inputs needed to reopen it:
- (a) A second, model-specific observable that ordinary flows cannot produce, computed by ray tracing. Candidates
  are the extra light-curve peaks of 2408.07350 or a transmitted-light polarisation signature. It must be shown
  not to be mimicked before any data step.
- (b) Per-flare astrometry with a measured static-flux fraction: GRAVITY reduced data from the collaboration, or a
  reduction of ESO raw data, which takes days of work.

"Not done this round" is not the same as "physically impossible".
