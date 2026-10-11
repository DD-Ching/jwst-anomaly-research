# B on R3-A — Destruction of the "frozen pair-correlation floor"

Target: `R3A-correlation.md`. **[lit]** = verified source, **[calc]** = my own
arithmetic, **[judg]** = judgement. Accessed 2026-10-09.

**Bottom line.** As written (A5), the theory predicts E[Ŝ] = 0: the "frozen" sign needs an unstated
axiom that A1 forbids. The axioms give a metric, not a causal order, so the Lorentz step fails. D is
chosen: the author's own light-cone convention gives 5e8 to above 1e15. The Holometer cross-spectrum
already excludes the nominal amplitude with about 67 % probability. ZTF fits on a laptop but reaches
about 1.8σ at best, and every persistent systematic mimics the signal with the same sign.

---

## 1. Internal consistency

**1a. The persistence the prediction needs contradicts A5 (fatal).** A1 says that only records exist.
The detection of star a in exposure k and the detection of a in exposure k' are therefore different
records. A5 makes their n vectors independent and Haar-distributed. So ⟨n_{a,k}, n_{b,k}⟩ is drawn
again at every epoch, with mean zero, and the per-pair coefficient averaged over T epochs is zero-mean
noise of variance 1/(DT). Then E[r^A_ab r^B_ab] = 0, because the A and B epochs carry independent draws:
**E[Ŝ] = 0 under the theory as written** [calc]. Section 4.1's claim that the coefficient "belongs to the
records, not to the detector" and persists "across epochs, seasons and instruments" needs a persistent
loading vector per *object*. That is an unstated axiom A5′, and it reintroduces objects, which A1
excludes. Literal A5 leaves only a fourth-moment signature: excess kurtosis 2/D ≈ 1e-9 in the summed noise of
each exposure, which needs about D²/4 ≈ 7e17 exposures [calc]. Kill criterion 4 then fires by
construction.

**1b. Lorentz invariance does not follow (fatal for Section 2).**
- C is symmetric (A1), so "a chain from i to k" is the same thing as "a chain from k to i". The relation
  cannot be antisymmetric, so it is not a partial order. The argument that "a chain cannot return with
  C < 1" misreads A3, which allows C_ii = 1 ≥ C_ij² [calc].
- A3, C_ik ≥ C_ij C_jk, means that d_ik ≤ d_ij + d_jk. That is the *ordinary* triangle inequality, so d
  is a Riemannian (Euclidean-signature) metric. Timelike proper time obeys the *reverse* triangle
  inequality [calc].
- If C^S_ij > 0 and C^S_jk > 0, A3 forces C^S_ik > 0. "Connected" is therefore an equivalence relation.
  In Minkowski space, two spacelike events can share a common past; here they would have to be
  connected, which contradicts the definition "spacelike = no chain" [calc].
- The only rescue, Osterwalder–Schrader rotation of H_S as Euclidean correlators, needs E(4) invariance
  and analyticity, which A4 does not give. Discrete Lorentz invariance needs random sprinkling
  (Bombelli–Henson–Sorkin) [lit], and Ignatowsky needs continuous frame maps.
- The "same detector event" window is a detector-frame simultaneity choice, so the noise sector does
  carry a frame [judg].

**1c. D is chosen, not derived.**
- A6 equates "dimension of novelty" with a count of σ = 1 cells. The author's own point 7.2 lists
  alternatives spanning 1e7 to 1e180.
- The arithmetic, given the inputs, checks out: V_H = 1.19e13 Mpc³ (Planck18, astropy) [calc]. But σ8 ≈
  0.81 means σ = 1 is *not* at 8 h⁻¹ Mpc. Linear σ(R) = 1 at R_* ≈ 5.8 h⁻¹ Mpc = 8.6 Mpc
  (Eisenstein–Hu no-wiggle, σ8 = 0.811), which gives V_* ≈ 2.6e3 Mpc³ and **D ≈ 4.5e9** on the t0
  slice, at the edge of the stated band [calc].
- The t0 slice is spacelike to the observer, not the "causal past".
- Counting cells *on the past light cone* with R_*(z) from linear growth gives [calc]:

| z_max | 1 | 2 | 3 | 5 | 8 |
|---|---|---|---|---|---|
| D_lightcone | 4.8e8 | 1.3e10 | 1.6e11 | 8e12 | 1.2e15 |

  It diverges as R_*(z) shrinks, so "no free parameter" is false.

**1d. Born rule and no-signalling.**
- *Born.* Section 3 says that free QFT is recovered. For independent systems in a product state, QM
  predicts that joint outcome statistics factorize exactly, so Cov = 0. A universal ρ ≈ 2.4e-5 between
  independent shot-noise sources is a violation of quantum theory's tensor-product rule, not a corollary
  of it [judg].
- The proposed rule P(i|j) = |C_ij|² also fails to normalise. Every unrelated record contributes about
  sin²θ_i sin²θ_j / D, so Σ_i P(i|j) ≈ 1 + N_rec/D. With N_rec far above D = 1.7e9 (every photon
  detection is a record), this sum is far above 1 [calc]. Gleason's theorem gives a trace-class measure
  on a projector lattice, not |C_ij|² between non-orthogonal records.
- *No-signalling* holds if the latent variable is uncontrollable: A5′ is a classical common cause [calc].
- *Bell tests.* ρ ~ 1e-5 carries about 3e-10 nats of mutual information. Hall's model needs about 14 %
  relaxed measurement independence [lit], so cosmic Bell tests (Handsteiner 2017; Rauch 2018) [lit]
  neither see this effect nor are threatened by it.

## 2. Already excluded?

**Lab: the Fermilab Holometer, the closest existing test (strains the effect, does not kill it).**
- Setup [lit]: two co-located, independent 39 m interferometers whose output is dominated by
  photodiode shot noise, which the hypothesis treats as pure innovation. Each output is digitised at
  100 MS/s, downsampled to 50 MHz, and cross-correlated in real time. The 145 h of 2015 data average
  noise down by √(145 h × 700 kHz) ≈ 6e5. The weighted broadband cross-spectrum is "consistent at 1.1σ
  with zero broadband correlation". Correlated electronic pickup is below 1 % of the statistical
  sensitivity.
- Mapping to R3-A [calc]: a frozen white correlation r adds a flat real cross-spectrum r·√(S₁S₂). The
  normalised 1σ is therefore about 1/6e5 ≈ 1.7e-6, and |r_eff| ≲ 5e-6 at roughly 3σ. The cross-spectrum
  combines 2 × 2 photodiode pairs, so r_eff ≈ the mean of 4 independent frozen r's, with rms ≈ 1.2e-5
  at nominal D. P(|r_eff| < 5e-6) ≈ 0.33 at D = 1.7e9 and ≈ 0.52 at D = 5e9. **One instrument pair
  excludes the nominal amplitude with probability ≈ 67 %, not decisively.**
- Escape route: the two digitisers are "isolated and independently synchronized to GPS" [lit], so this is
  arguably not "one acquisition system". That escape only works because "detector event" is undefined
  (1b), and using it makes the lab sector unfalsifiable.
- The author's kill test 5 is cheap: two noise diodes, one dual-channel ADC, and 4e10 samples per pair,
  which is about 400 s at 100 MS/s; repeat for many pairs. It would settle the question; nobody needs
  ZTF [calc].
- Johnson-noise thermometry (k_B to 3.5e-6 with a switched cross-correlator [lit]) may constrain this
  too, but its ratio method partly cancels the effect. No bound is claimed.

**Astronomy: not excluded, because nobody measures the right statistic.**
- Kepler PDC/CBV-type cotrending [lit] and ZTF per-quadrant zeropoints [lit] remove *common* modes.
  These are mean-sensitive and blind to a zero-mean, random-sign ρ_ab.
- PTA Hellings–Downs (NANOGrav 15-yr: 67 pulsars, 2,211 pairs [lit; pair count calc]) has a per-pair
  correlation uncertainty that is orders of magnitude above 1e-5. The TOAs are also not same-event.
  **No constraint.**
- LIGO Schumann correlations (Thrane et al. 2013 [lit]) involve separate sites, so they fall outside the
  theory's "same event". **No constraint.**
- Gaia: the author already concedes it is uncompetitive.
- I found no published analysis that bounds per-pair *squared* persistent correlations of same-exposure
  residuals at 1e-5. This is "not found", not "does not exist".

## 3. Renaming check

- **The mechanism** (A5′): every same-event innovation is a projection of one shared latent Gaussian
  vector onto fixed Haar loadings. That is a classical latent-factor common-cause model, a hidden
  variable. In foundations language it is a weak, zero-mean relaxation of measurement and outcome
  independence (the superdeterminism family; Hossenfelder & Palmer [lit]) at about 1e-9 bits.
- **The 1/D overlap and the MP-plus-semicircle spectrum** are textbook random-vector and random-matrix
  results (Welch, JL), as the author says.
- **Holographic noise** (Hogan 2008 [lit]) is a parameter-free correlated noise between co-located
  instruments, tested by the same Holometer cross-spectrum. R3-A differs in sign structure and scale,
  not in experimental class.

## 4. Novelty verdict

The mechanism is not novel. It is a common-cause hidden variable plus standard high-dimensional geometry.
The only new element is the numerology D = N_cells(σ = 1), which is unmotivated (1c) and inconsistent
with A5 (1a). The *statistic* Ŝ is a legitimate and apparently unpublished survey null test. Its value is
as a systematics diagnostic, not as a test of this theory [judg].

## 5. Triage

| # | Prediction | Verdict | Why |
|---|---|---|---|
| P1 | ZTF split-sample Ŝ = 1/D ≈ 6e-10 | **Conditional PASS** (only for repaired A5′; FAIL for R3-A as written) | A5 ⇒ E[Ŝ] = 0; ≤ 1.8σ ideal at nominal D; positive-definite mimics |
| P2 | Ŝ independent of separation, colour, class, flux | FAIL | Needs ≥ 5σ per bin, so N_p × n_bins beyond ZTF; the only signature that separates signal from systematics is untestable |
| P3 | Cross-instrument ZTF × ATLAS persistence | FAIL | Different records ⇒ 0 under A5; kill criterion 4 fires by construction |
| P4 | Lab: frozen \|r\| ≈ 2.4e-5 per device pair | **PASS** (cheapest decisive test; outside this repo's scope) | Holometer already ~67 % excludes nominal; dual-ADC test needs about 400 s per pair |
| P5 | MP ⊞ semicircle eigenvalue bulk | FAIL | Width 2√(N/D) ≈ 0.05 at N = 1e6, below the finite-T MP width 2√(N/T) ≈ 115 |
| P6 | c_gw = c_EM, zero LV | FAIL (not discriminating) | GR predicts the same; the derivation is broken (1b) |
| P7 | Gaia DR4, Rubin | FAIL / defer | Gaia conceded; Rubin inherits P1's problems |

**P1 details (the only astronomy item worth a decision).**

- *Data* [lit]: ZTF DR light curves are public bulk Parquet. DR18 has 174,033 files (one per field,
  CCD-quadrant and filter) across 1,174 fields, about 6.0 TB in total, with about 2.3e9 r-band light
  curves. DR23 is on AWS (`s3://ipac-irsa-ztf/contributed/dr23/lc/hats`, no account needed); DR24 is
  the latest release.
- *Columns:* the files have per-epoch hmjd, mag, magerr, clrcoeff and catflags. **There is no exposure
  ID, zeropoint or airmass.**
- *Same-exposure matching* has to be rebuilt from hmjd, field and rcid. The heliocentric correction
  varies by up to about ±60 s across a 7° field [calc], which is comparable to the exposure spacing in
  high-cadence fields. So hmjd must be converted back to the observed time per object. Airmass and
  zeropoint come from IRSA image metadata, with one batched query per field.
- *Primary and secondary grid objects have different objectids*, so they must be cross-matched before the
  grid split.
- *Disk:* stream one field at a time (about 2–5 GB, files about 35 MB each, so within the 200 MB rule),
  then delete. Peak disk is under 10 GB. The total r-band download is about 3 TB, roughly 3 days at
  100 Mbit/s. **Bandwidth is the binding constraint.**
- *CPU* [calc]: do not loop over 1e14 pairs. Use Σ_{a≠b} r^A_ab r^B_ab = ‖X_Aᵀ X_B‖_F² / (T_A T_B) −
  Σ_a r^A_aa r^B_aa, with X the N × T standardised innovations and zero-fill for missing epochs. That is
  unbiased under the null. Cost is N·T_A·T_B ≈ 1e6 · 300² ≈ 1e11 flop per field, about 1e14 flop in
  total. That takes about an hour of laptop BLAS, plus per-object ARMA/DRW whitening (O(N·T)).
  **Feasible.**
- *Sensitivity* [calc]: σ(Ŝ) = 1/(T√N_p) = 3.3e-10 at N_p = 1e14 and T = 300. The signal-to-noise is
  1.8σ at nominal D, 0.6σ pessimistic and 6σ optimistic, before cuts (catflags: about 11 %
  of images bad [lit]) and the dof lost to whitening, which cost a further factor of about 1.5–2.
- *Mimics:* Ŝ estimates E[ρ²], which is positive for **every** persistent correlation of either sign.
  ZTF repeatability is 8–25 mmag at
  the bright end [lit], dominated by systematics, so the residual structure in flat-field, colour and PSF
  must be suppressed to an rms pair correlation below 1.4e-5. That means modelling the systematics to
  about 1e-5 of the innovation variance. No ZTF work I found demonstrates that.
- *Verdict:* run it only as a cheap systematics null. A positive Ŝ cannot be attributed to R3-A (P2 is
  untestable), and a null kills only the unwritten A5′.

---

## Verified references

1. A. S. Chou et al. (Holometer), "First measurements of high frequency cross-spectra from a pair of large Michelson interferometers", PRL 117, 111102 (2016). arXiv:1512.01216.
2. A. Chou et al. (Holometer), "Interferometric constraints on quantum geometrical shear noise correlations", CQG 34, 165005 (2017). arXiv:1703.08503 (abstract only checked).
3. C. J. Hogan, "Measurement of quantum fluctuations in geometry", PRD 77, 104031 (2008). arXiv:0712.3419.
4. M. J. W. Hall, "Local deterministic model of singlet state correlations based on relaxing measurement independence", PRL 105, 250404 (2010). arXiv:1007.5518.
5. J. Handsteiner et al., "Cosmic Bell test: measurement settings from Milky Way stars", PRL 118, 060401 (2017). arXiv:1611.06985.
6. D. Rauch et al., "Cosmic Bell test using random measurement settings from high-redshift quasars", PRL 121, 080403 (2018). arXiv:1808.05966.
7. S. Hossenfelder, T. N. Palmer, "Rethinking superdeterminism", Front. Phys. 8, 139 (2020). arXiv:1912.06462.
8. L. Bombelli, J. Henson, R. D. Sorkin, "Discreteness without symmetry breaking: a theorem", Mod. Phys. Lett. A 24, 2579 (2009). arXiv:gr-qc/0605006.
9. NANOGrav Collaboration (G. Agazie et al.), "The NANOGrav 15-year data set: evidence for a gravitational-wave background", ApJL (2023). arXiv:2306.16213.
10. E. Thrane, N. Christensen, R. Schofield, "Correlated magnetic noise in global networks of gravitational-wave interferometers", PRD 87, 123009 (2013). arXiv:1303.2613.
11. M. C. Stumpe et al., "Kepler Presearch Data Conditioning I", PASP (2012), DOI 10.1086/667698. arXiv:1203.1382.
12. F. J. Masci et al., "The Zwicky Transient Facility: data processing, products, and archive", PASP (2019). arXiv:1902.01872.
13. J. Qu et al., "Improved electronic measurement of the Boltzmann constant by Johnson noise thermometry", Metrologia 52, S242 (2015). arXiv:1501.00195.
14. D. J. Eisenstein, W. Hu, "Baryonic features in the matter transfer function", ApJ 496, 605 (1998). arXiv:astro-ph/9709112 (no-wiggle fit used in 1c).
15. ZTF DR18 release notes (IRSA), https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/dr18/ztf_release_notes_dr18.pdf (Table 1; §12c columns). No arXiv ID.
16. NERSC ZTF DR18 Parquet README, https://portal.nersc.gov/cfs/cosmo/data/ZTF/lc/lc_dr18/README.txt (file count, 6.0 TB). AWS Open Data ZTF entry, https://registry.opendata.aws/ztf/ (DR23 HATS bucket). No arXiv IDs.

Calculations: Planck18 distances via astropy (repo `.venv`). The σ(R) and light-cone cell counts used an
inline Eisenstein–Hu no-wiggle script (not saved), with linear growth from the Heath integral.
