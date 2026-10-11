# Goal round 10: pre-registration of the NIST 2015 no-signalling check

Written 2026-10-11 20:46 UTC. The run-level counts had been opened only for 19_45 grouping "5", to identify the file
format; the inference below has not been run yet.

**Question.** In the only public loophole-free Bell dataset with spacelike-separated setting choices (NIST 2015, Shalm
et al., PRL 115, 250402), does one party's free setting choice change the other party's local click statistics?
This is the elementary "send / receive" operation that any entanglement-based channel would need.

**Data.** NIST belltestdata, `Peter's Data Calculations/Diagnostic Statistics/*.xlsx`. The first row of each
pulse-grouping tab holds the 16 whole-run counts N(oA, oB | x, y). Order: ++, +0, 0+, 00 for setting pairs ab, ab', a'b,
a'b'. Six runs: 19_45, 02_54, 03_43, ClassicalXOR 1–3. Groupings: 5, 456, 34567, 2345678 (reported in the paper) and
123456789 (not reported; spacelike separation less certain, so secondary only).

**Hypotheses.**
- H0 (quantum mechanics, no-signalling): P_A(+ | x, y) does not depend on y, and P_B(+ | x, y) does not depend on x.
- H1 (signalling): a setting-dependent marginal shift δ.

**Ordinary explanations that could fake H1:**
- electrical or optical crosstalk between a setting's RNG / Pockels-cell drive and the detector at the same site;
- slow drifts (alignment, temperature, fridge warming) correlated with setting frequencies;
- unequal setting probabilities combined with drifts;
- the analysis bug NIST flags for "00" entries in other files. This is checked here by
  N(x,y) = N++ + N+0 + N0+ + N00 against the tab's trial totals.

**Method.**
- Independent of NIST's spreadsheet formulas: binomial likelihood-ratio (G) tests and an exact conditional
  (Fisher-type) test per comparison.
- Then the reported NIST two-proportion z for the same cells, as a reproduction check.
- **Primary family:** 6 runs × 4 comparisons (A given x = a and a'; B given y = b and b') at grouping 2345678, which
  is the largest grouping reported as spacelike. That is 24 tests. Bonferroni threshold α = 0.05/24.
- Secondary: the other groupings. These are nested and correlated, so they are descriptive only.

**Decision rules.**
- **Continue** (treat as an anomaly worth chasing): a family-wise significant shift (p < 0.05/24) that has the same
  sign in ≥ 2 independent runs and survives the setting-frequency check. Next step would be to rule out crosstalk
  with the raw time tags (hdf5) and the lightcone-shift control runs.
- **Stop** (null): no family-wise significant shift. Record the 95 % upper bound on |δ| and the implied channel
  capacity in bits per trial, as the quantitative statement of "no usable signal" for this setup.

**Expected result.** Null (H0), with |z| ≲ 2. NIST's own diagnostic tabs already show |z| < 1.5 for 19_45.
