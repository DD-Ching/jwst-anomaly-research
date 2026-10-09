# R3-A — Record-Correlation Ontology and the Frozen Pair-Correlation Floor

System A, round 3, group R3-A. **Everything below is a hypothesis.** No measurement or citation
is claimed. Every number marked "derived" depends only on the axioms plus two measured inputs:
the comoving particle-horizon radius and the scale where the matter field becomes nonlinear.

---

## 1. Axioms

**A1 (Records).** There is a finite set R of records. The only primitive is a real, symmetric,
positive-semidefinite matrix C on R with C_ii = 1. Equivalently, every record i is a unit vector
psi_i in a real inner-product space H, and C_ij = <psi_i, psi_j>. Nothing else exists: no events,
no manifold, no fields, no clock.

**A2 (Split).** H = H_S (+) H_N splits into two orthogonal parts. In H_S, C has *structured*
correlations: chains that compose. In H_N, C is *typical*. Write psi_i = cos(theta_i) s_i + sin(theta_i) n_i,
with s_i in H_S and n_i in H_N. The weight sin^2(theta_i) is the fraction of record i that its
own structured past cannot predict. Operationally, this is the variance fraction of its
**innovation** (see Section 4).

**A3 (Composition).** On H_S the correlation of two records is fixed by the best chain between
them: C^S_ik >= C^S_ij C^S_jk, with equality along a unique maximal chain. Define the
correlation length d_ij = -ln C^S_ij. Records with C^S_ij = 0 have no chain between them.

**A4 (Exchangeability).** The joint law of C is invariant under every permutation of R that
preserves C^S. No record is privileged, so there is no preferred origin, direction or rest frame.

**A5 (Typicality of novelty).** For every new record i, n_i is uniformly distributed (Haar measure)
on the unit sphere of H_N, independently of n_j for j != i. Correlations between the n_i come only
from the geometry of finite dimension.

**A6 (Dimension of novelty).** dim H_N = D, where D is the number of mutually orthogonal
*self-records* in the record's causal past. A self-record is a maximal block of H_S whose internal
correlation exceeds its correlation with its complement. In emergent language, a self-record is a
region whose own fluctuation is of order unity, which means the variance of the coarse-grained
density contrast is sigma^2(R_*) = 1.

These six statements are the whole theory. A1–A4 produce spacetime and ordinary physics. A5 and A6
produce the effect.

---

## 2. Lorentz invariance (derived)

1. **Causal order.** A3 makes "there is a maximal chain from i to k" a partial order: it is
   transitive and antisymmetric, because a chain cannot return to its start with C < 1. Records
   joined by no chain are incomparable, which is the emergent notion of spacelike separation.
2. **Universal cone.** A chain's correlation is multiplicative, so d is additive along chains.
   Coarse-grain R into cells of n >> 1 records each. A4 forces the largest d accumulated per cell
   step to be the same for every chain, whatever kind of record the chain carries. This number is
   the unique limiting speed c. Photons, gravitons and every other excitation are patterns that
   propagate along H_S chains, so they all share one cone. This gives c_gw = c_EM exactly, not merely
   to some precision, and so meets the 10^-15 bound from GW170817.
3. **Group.** A4 provides homogeneity, isotropy and the relativity principle (no preferred frame).
   Item 2 provides a finite invariant speed. The standard Ignatowsky argument then leaves only the
   Lorentz group as the symmetry of coarse-grained H_S.
4. **The noise sector cannot break it.** H_N carries an O(D) Haar measure (A5). That is an internal
   symmetry with no spacetime index, and H_N is orthogonal to H_S (A2). So the noise sector has
   no vector or tensor that could define a preferred frame. Every operator that violates Lorentz
   invariance needs a background object with a spacetime index, and A1–A6 contain none. Laboratory
   tests at 10^-18 and clock-comparison tests therefore predict exactly zero: there is no
   dispersion, no birefringence and no direction dependence.
5. **Locality.** The effect in Section 4 is a correlation between records formed **at the same
   detector event**: the same observer at the same proper-time interval. It is not a correlation
   between the distant sources. It adds no H_S chain, so it creates no shortcut and carries no
   signal. No experimenter can choose n_i (A5), so the no-signalling condition holds.

---

## 3. Recovery of standard physics (sketch)

- **Fields.** Coarse-grained H_S correlations between cells define two-point functions. A4 together
  with the Lorentz cone of Section 2 forces Lorentz-covariant, positive (reflection-positive)
  correlators. Reconstruction then gives a Hilbert space and local fields. Free QFT is the case
  where the correlators are Gaussian. Interactions are the connected higher correlations.
- **Born rule.** The probability of finding record i given record j is |C_ij|^2. This is the only
  non-contextual additive measure on a real or complex inner-product space of dimension 3 or more.
- **Gravity.** The emergent metric is the coarse-grained correlation length d. Requiring the cell
  count (an entropy) to be stationary under small deformations of chains gives Einstein's equations
  with Lambda as the integration constant. **This step is asserted, not derived here** (see
  weakness 3).
- **Cosmology.** Nothing in A1–A6 changes the mean of any H_S observable. Background expansion,
  distances, the CMB and growth are those of LCDM. The noise sector has **zero mean**, so it cannot
  dim or brighten anything. Cosmic opacity is exactly unchanged, which satisfies lesson 3.
- **Single-object statistics.** Each n_i is a unit vector, so the noise variance of any single
  record is exactly the standard one. Poisson shot noise, DRW quasar variance and pulsar
  white noise are all unchanged. Any per-object test therefore sees standard physics.

---

## 4. The derived small effect: the frozen pair-correlation floor

### 4.1 Formula

Take two distinct records i and j: two sources measured in the **same exposure** (the same detector
event). Let their innovations (the parts their own pasts cannot predict) be the noise components.
Then A5 implies:

    rho_ij = <n_i, n_j> sin(theta_i) sin(theta_j)

    E[<n_i,n_j>] = 0,   E[<n_i,n_j>^2] = 1/D   (exact for Haar vectors in R^D)

For a series that has been whitened and contains only innovations, sin(theta) is close to 1. The
whitening step is: fit each light curve's own DRW/ARMA model, subtract the one-step prediction and
normalise by the predicted error. This includes photon noise, which is pure innovation. So the
**frozen pair correlation** has

- mean 0 over pairs, so there is no common mode and no dimming;
- rms epsilon = D^(-1/2), the same for every pair;
- a **sign and value fixed per pair** that persist across epochs, seasons and instruments, because
  the coefficient belongs to the records, not to the detector;
- **no dependence** on angular separation, distance, redshift, colour, magnitude or object class.
  It depends only on the derived invariant D.

Standard physics predicts rho_ij = 0 exactly for the innovations of unrelated sources once
systematics are removed. A higher-order prediction follows for the eigenvalues of the
per-exposure covariance of N sources: the bulk should be the Marchenko–Pastur law convolved with a
semicircle of radius 2 sqrt(N/D). This is far too small to see in a single exposure.

### 4.2 Amplitude (derived)

By A6, D is the number of self-records, that is nonlinear cells, in the observer's causal past:

- The comoving particle horizon is R_H ≈ 14.2 Gpc, which gives V_H ≈ 1.2 × 10^13 Mpc^3.
- The nonlinear scale today is sigma(R_*) = 1 at R_* ≈ 8 h^-1 Mpc ≈ 11.8 Mpc. A top-hat cell then
  has V_* ≈ 6.9 × 10^3 Mpc^3.
- **D ≈ 1.7 × 10^9.** Different conventions for the cell change this: a Gaussian rather than
  top-hat window, counting along the light cone rather than on the t_0 slice, and the redshift
  evolution of R_*. Together they give the band **D ∈ [5 × 10^8, 5 × 10^9]**.

So:

    epsilon = D^(-1/2) ≈ 2.4 × 10^-5        (band 1.4–4.5 × 10^-5)
    S ≡ E[rho_ij^2] = 1/D ≈ 6 × 10^-10       (band 2 × 10^-10 – 2 × 10^-9)

This sits far below per-object precision. A light curve with T ≈ 10^3 epochs measures any r_ij
only to about T^(-1/2) ≈ 0.03, roughly 1,300 times larger than epsilon.

There is no free parameter. D is computed from two measured quantities (R_H and sigma_8 via R_*).
The scheme band of a factor of about 3 is the only latitude.

### 4.3 Estimator

The signal has random sign, so averaging r_ij gives nothing. Use split-sample persistence instead:

    Ŝ = (1/N_p) Σ_pairs r_ij^(A) r_ij^(B)

Here A and B are **disjoint, independent** subsets of the simultaneous epochs of the pair. They can
be two time halves, the ZTF primary versus secondary field grid, or two instruments. Under the null,
E[Ŝ] = 0. Under the hypothesis, E[Ŝ] = 1/D. With T_A ≈ T_B ≈ T effective independent epochs per
subset, σ(Ŝ) ≈ 1/(T sqrt(N_p)).

### 4.4 Datasets and N

Only pairs recorded in the same exposure count. That makes wide-field, high-cadence photometric
surveys the natural data:

- **ZTF public data releases.** There are of order 10^9 light curves, many with several hundred
  r-band epochs, and each exposure covers 47 deg². Pairs within one field are ~5 × 10^10
  (high latitude) to ~2 × 10^12 (plane), summed over ~1,700 primary and secondary fields. That gives
  **N_p ≈ 10^14 – 10^15**, with T ≈ 150–300 per half.
- **Gaia DR4 epoch photometry.** The field of view is ~0.5 deg² and there are ~35 transits per half,
  giving σ(Ŝ) ≈ 6 × 10^-9. This is not competitive.
- **Rubin LSST** (from about 2027 onward). Each visit covers 9.6 deg², with ~800 visits over ten
  years. This gives N_p ≈ 10^15 and T ≈ 400.

---

## 5. N for 5σ, comparison with data, kill criterion

The requirement is 5σ(Ŝ) ≤ 1/D, so

    N_p(5σ) = (5 D / T)^2

| D | T per half | N_p(5σ) |
|---|---|---|
| 1.7 × 10^9 (nominal) | 300 | 8 × 10^14 |
| 1.7 × 10^9 | 150 | 3 × 10^15 |
| 5 × 10^9 (pessimistic) | 300 | 7 × 10^15 |
| 5 × 10^8 (optimistic) | 300 | 7 × 10^13 |

**Now.** The ZTF public releases give an expected significance of roughly 2–5σ at nominal D. That
is enough to detect the optimistic edge and to test the nominal value, but not enough to rule out
the pessimistic edge. Adding ATLAS or Pan-STARRS cross-instrument pairs raises the effective N_p by a
factor of 2–3. **Rubin DR3–DR5 crosses 5σ across the whole band.** The data exist today, and as far
as this author knows nobody has computed this cross-split persistence statistic for unrelated
sources.

**Kill criteria** (any one of these suffices):

1. After the systematics regression in Section 6, ZTF + ATLAS gives Ŝ < 2 × 10^-10 at 95% CL. A
   null at the optimistic and nominal amplitudes kills the nominal hypothesis. Rubin must then reach
   Ŝ < 2 × 10^-10 at 5σ to kill the whole band.
2. A detection with Ŝ outside [2 × 10^-10, 2 × 10^-9]. The amplitude is derived, so the wrong
   amplitude means the axioms are wrong.
3. A detected Ŝ that depends on separation, colour, magnitude, chip pair, airmass or Galactic
   latitude. The hypothesis forbids any such dependence.
4. A cross-instrument covariance E[r^ZTF r^ATLAS] that is consistent with zero while the
   within-instrument Ŝ is not zero. This means the signal belongs to the instrument, not to the
   records.
5. **Laboratory test** (cheaper, and possibly already decisive). Two independent physical random
   sources, such as thermal or shot-noise sources, sampled by one acquisition system on a common
   clock. Raw samples must be used, before any randomness extractor. The hypothesis predicts a
   frozen |r| ≈ 2.4 × 10^-5 per device pair, which needs T ≳ 4 × 10^10 samples for 5σ per pair. If
   several device pairs give rms |r| < 5 × 10^-6, the hypothesis is dead.

---

## 6. Ordinary systematics that mimic it, and how to separate them

| Mimic | Why it persists across splits | Separation |
|---|---|---|
| Per-exposure zero point / transparency | Common positive mode | Subtract the per-exposure median innovation. This removes 1 of D dimensions, about 6 × 10^-10 of the signal. The signal has zero mean; this mimic does not. |
| Colour-dependent extinction and colour terms | The sign is set by the pair's colours, so it persists | Regress r_ij^(A) r_ij^(B) on a flexible function of (colour_i, colour_j, airmass) and keep the residual. Use colour-matched pairs. |
| Flat-field and PSF errors at fixed pixel positions | ZTF uses a fixed field grid, so a pair lands on the same pixels | Split by **primary versus secondary grid** or dither, so the pair lands on different pixels. Real records are unchanged; pixel systematics change. |
| Crosstalk and amplifier ghosts between CCDs | Fixed geometric coupling | Mask known crosstalk pairs. Test for dependence on (amplifier_i, amplifier_j). |
| Blending and shared PSF wings | Close pairs share pixels | Exclude separations < 30″. The hypothesis predicts no dependence on separation. |
| Difference-imaging reference reuse | A shared reference image correlates the residuals | Use different reference builds for A and B, or use forced photometry. |
| Scintillation / seeing | Correlated over the isoplanatic patch | This depends on separation (arcmin). Fit and remove it; it vanishes beyond ~1°. |
| Truly correlated astrophysics (lensed images, binaries, clustered variables) | A physical common cause | Remove known associations. They are too few to matter at N_p ~ 10^14. |
| Whitening model error | An imperfect DRW fit leaves autocorrelated residuals, which inflate Ŝ only through pair overlap in time structure | Use a null with **time-shifted** pairs: correlate i at epoch k with j at epoch k+Δ, where Δ is larger than any model memory. The hypothesis predicts 0 for these (the innovations are white), so this null calibrates any bias. |

The cleanest discriminant is **instrument transfer**. Systematics belong to a detector, and the
hypothesised coefficient belongs to the pair of records. A second strong check: the hypothesis
predicts the **same** Ŝ for star–star, star–quasar and quasar–quasar pairs, and for faint pairs
dominated by photon noise as well as bright pairs dominated by intrinsic variability. A systematic
would essentially never be blind to both class and flux regime.

---

## 7. Three weakest points

1. **The laboratory may already exclude it.** The effect is universal: every pair of
   simultaneously recorded innovations correlates at about 2.4 × 10^-5. Cross-correlation checks of
   independent physical noise sources and of raw QRNG streams at ≳10^11 samples would therefore
   already see it, if anyone has run them without extractors. I cannot check this without the
   literature. Restricting the effect to astronomy would require an ad hoc exemption. I have not
   added one, so this is the most likely way the world dies. The brief asks for an effect that
   appears only in large-N data; here that is true for astronomical objects but not for the lab.
2. **A6 is the load-bearing step, and it is the least derived.** Why should the dimension of
   novelty equal the number of nonlinear cells in the causal past, rather than the number of Planck
   cells (D ~ 10^180, invisible) or of CMB modes (D ~ 10^7, already huge)? I chose the scale where
   the structured correlation becomes block-diagonal. That choice is motivated, not forced. The
   factor-3 scheme band is honest, but the choice of scale behind it could be wrong by orders of
   magnitude.
3. **The recovery is partial.** Lorentz invariance and the shared cone follow cleanly. Einstein's
   equations, the interacting Standard Model, and why the coarse-grained H_S has 3+1 dimensions are
   asserted (Section 3), not derived. The "same detector event" coherence window is also
   underspecified: exposure time, readout or decoherence time. A longer window would add pairs from
   different exposures and change N_p and the systematics budget.

Smaller issues: the time-dependent D (dD/dt ~ 3H D) is negligible over a survey and therefore gives
no extra handle. The whitening step could leak intrinsic variability into "innovations".
Superdeterminism-style worries apply to the lab test (see below).

---

## Resemblances I noticed

- **Welch bound / Johnson–Lindenstrauss / random-matrix theory.** The 1/D mean square overlap of
  random unit vectors is standard high-dimensional geometry. The spectral corollary is a
  Marchenko–Pastur plus semicircle law.
- **Measurement-independence violation / superdeterminism.** Shared latent randomness between
  independent measurement outcomes is a weak, frozen, zero-mean version of it.
- **Quantum Darwinism, relational QM, consistent histories.** All take "records" or relations as
  primary.
- **Emergent spacetime from correlation or entanglement structure** (entanglement-builds-geometry
  arguments; Jacobson's thermodynamic derivation of Einstein's equations, which my gravity step
  borrows in spirit).
- **Causal sets.** The order from chains and the 1/sqrt(N) fluctuation logic (Sorkin's
  Lambda ~ N^-1/2) resemble my D^-1/2.
- **Ignatowsky / von Ignatowsky–Frank–Rothe.** My Lorentz step is that derivation applied to
  chains.
- **Osterwalder–Schrader / Wightman reconstruction** for recovering fields; **Gleason** for the
  Born rule.
- **Holographic-noise proposals.** These also predict a universal correlated noise, but theirs is
  geometric and tied to Planck-scale position noise. Mine has no spacetime index and is random-sign.
