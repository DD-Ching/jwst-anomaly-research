# Round 1 (2026-10-09): three System A worlds, three System B reviews

Everything here is a **hypothesis** or a reviewer's estimate. The System A texts were written without literature,
and the System B reviews were written with references their authors opened (lists at the end of each B file).
Numbers marked [calc], [est] or [B-est] in the B files are reviewer estimates, not measurements.

## Origin

The owner's ideas (2026-10-09, in Chinese), condensed:
1. Distance self-consistency: luminosity, angular-diameter, time-delay and dynamical distances to one region that no
   ordinary model can reconcile.
2. Missing light with a paired excess elsewhere.
3. Transient connectivity whose rate depends on local conditions, not on sky position.
4. A causal event network: dependence between events that no admissible path explains. This is the owner's bet.
5. A "spacetime response material".

Plus: split into System A / System B; at least 30 % of exploration on quantities current searches never measure; an
"Ontological Reset".

Overlap with earlier work. Idea 1 was never tested here (D-001–D-067 are all morphology, light-curve or count
screens). Idea 2's deficit half is W3/W5 (D-052, D-057–D-063, D-066/D-067); its paired-excess half is new here.

## Worlds

| Group | Seed | File | Review |
|---|---|---|---|
| A1 | only states and allowed transitions; distance = transition count; rare long-range "b" transitions | [A1-relational.md](A1-relational.md) | [B-on-A1.md](B-on-A1.md) |
| A2 | no fundamental identity or locality; one quantity X ("carry", mutual information per update) | [A2-no-identity.md](A2-no-identity.md) | [B-on-A2.md](B-on-A2.md) |
| A3 | a positive measure over configurations; sparse "bridges"; distance measures disagree in a fixed pattern | [A3-inconsistent-distances.md](A3-inconsistent-distances.md) | [B-on-A3.md](B-on-A3.md) |

**Common verdict on the foundations.** None of the three derives Lorentz invariance or gravity. All three long-range
link mechanisms are, in B's reading, Markopoulou–Smolin "disordered locality" (gr-qc/0702044; cosmology in
0903.5303) under new names. Each achromatic-dimming prediction is the existing cosmic-opacity / distance-duality
test, already limited to Δτ ≲ 0.01–0.05 (1004.2053).

## Survivor table

| # | Prediction (source) | Verdict | Why it survives | Data | Reach of existing data |
|---|---|---|---|---|---|
| S1 | **Burst twins**: pairs of GRBs or FRBs far apart on the sky, with matching light curves at any delay (A2 P1, made model-agnostic; s fixed at 1, free-s secondary) | conditional PASS | Lensed-burst searches cut on sky position first (Fermi GBM, 2006.07095) or used delays ≲ 100 ms (CHIME, 2204.06014); wide separations at arbitrary delay are untested. Also an instance of the owner's idea 4. | Fermi GBM burst catalogue + TTE/CTIME light curves; CHIME/FRB Catalog 2 | ~6 × 10⁻³ twins per burst from the 500 brightest multi-pulse GRBs [calc] |
| S2 | **Flat-kernel SN residuals**: SN Ia Hubble residuals regressed on a uniformly weighted foreground column versus the lensing-kernel column (A3 P2b) | conditional PASS | Lensing magnification is detected (DES-5YR, 2406.05047), but a flat-kernel term has not been fitted; a kernel-shape test, not an amplitude test | DES-SN5YR or Pantheon+ public tables + DES Y3 / Legacy Surveys galaxies; DESI DR2 BAO | σ_γ ≈ 0.007 mag, i.e. a new ~0.015 mag (2σ) limit; cannot reach A3's own kill threshold of 0.005 [calc] |
| S3 | **Hybrid images**: faint copies in lensed transients at model-fixed lags s_i − s_j, some before the leading image (A1 P2) | conditional PASS, limit only | Geometrically sound; precursor flux at these lags is unmeasured | COSMOGRAIL light curves; SN Refsdal / H0pe pre-explosion imaging | ~10⁻² stacked, against a predicted amplitude ≪ 10⁻¹⁰ under A1's own scaling [B-est] |

**Failed, with the main reason:**
- A1 P1 (consistency relation untestable at current EBL precision).
- A1 P3 (Bullet Cluster offset).
- A1 P4 (Lorentz-invariance limits put the broadening ≲ 0.2 ns).
- A1 P5 (internal error; laboratory bound 10⁻¹⁸).
- A2 P2 (internally inconsistent; no FRB below the Macquart floor).
- A2 P3 (cosmic opacity).
- A2 P4 (free function, unfalsifiable).
- A2 transfer covariance T (wrong sign under its own physics).
- A3 Δm = βΔκ (inconsistent with α ≳ 2; about 2,000 quads needed).
- A3 universal Q (undefined).
- A3 GW-vs-EM duality (sirens ~10 %).
- A3 aperture² mass deficit (gravitational slip / non-thermal pressure).
- A3 bridge dips (undefined duration, superluminal signalling).

## Correction to the causal-network idea

A1 claimed that blind timing searches between distinct extragalactic events cannot test light-cone violation even
in principle. B showed this is false. For a pair 1° apart at z = 0.5 the transverse separation (~108 Mly) exceeds
the radial error (~1–14 Mly), so a dependence within a 10-yr window is provably spacelike. The argument holds only
for pairs closer than about 0.1° at Gpc distances. The real obstacles are chance coincidences and common causes at
the observer (exposure, triggers, follow-up chains), which S1's scrambled-catalogue null must control.

## Next

S1 first: it is model-agnostic, it touches the owner's ideas 2 and 4, and existing catalogues can test it. Then S2,
then S3 as a limit. Ordinary twins (gravitational lensing, catalogue duplicates) are the positive control, and the
null is time-scrambled catalogues. Round 2 of System A should start from the failure modes above: no Lorentz
derivation, the disordered-locality renaming, and cosmic-opacity redundancy.
