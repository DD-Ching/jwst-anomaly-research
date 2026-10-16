# Plan: `arcs` on RX J2129.7+0005 (frozen before any `arcs` count; committed 2026-10-16 20:13Z, `bbe0306`; the run started 20:21Z)

**Why.** The Abell S1063 background anti excess (4 vs 0.52; P 0.031 clump-aware, #273, #279) was found post hoc in
one field. MACS0417 (#280) did not replicate it; MACS0647 (#281) gave no result; MACS1423 (#282) replicated at screen
level (P 0.016) but the clump-pair census explains it (P 0.259). RXJ2129 is the next independent field with
JWST-updated models (Dhanasingham et al. 2026, Zenodo 21460864; SOURCES.md). It has never been used to tune or choose
any `arcs` rule.

**Hypothesis tested (H0).** Background arcs in strong shear (g_min ≥ 0.2) align with the stretch the model predicts;
the number of `anti` background rows is consistent with the shape-noise expectation (σ_ε 0.26).

**Ordinary explanations for any `anti` row:** critical-curve flips (g_min ≥ 1), intrinsic shapes in weak shear,
deblending fragments and clump pairs, member/ICL photo-z blends, diffraction spikes, detector artefacts, lens-model
error.

**Inputs (sha256 in RESULT.md and `lensmodel.py`):**
- JWST program 2767 F200W pipeline catalogue `jw02767-o002_t001_nircam_clear-f200w_cat.ecsv` (2,061 s; MAST level 3,
  read from `s3://stpubdata/jwst/public/jw02767/L3/t/o002/`). The MAST query (322.4165, +0.0892, ±0.06°) found no
  other public F200W level-3 pointing.
- DJA grizli v7.0 `rxj2129-grizli-v7.0-fix.eazypy.zout.fits` (from the `.photoz.tar.gz`).
- Primary model `rxj2129-lenstool-maps`: the Lenstool maps of Zenodo 21460864 (`LENSTOOL/alphaX_DLSDS1`, `alphaY_…`,
  `kappa_…`; D_LS/D_S = 1 by name and README). Second model `rxj2129-glafic-maps`: GLAFIC `x/y-arcsec-deflect` and
  `kappa`, whose header says ZS = 3: used as published (normalisation at z_s 3, not D_LS/D_S = 1; this scales the
  deflection by one constant and leaves the predicted stretch axes of g < 1 rows unchanged only to first order, an
  ASSUMPTION stated here). z_lens 0.234 (CLASH, Postman et al. 2012), flat ΛCDM H0 70, Ωm 0.3. Each model's
  `validate` must reproduce its published κ from its deflections (median relative difference ≤ 10⁻³), else that
  model gets no `arcs` run.

**Rules (D-024 defaults, ASSUMPTIONs, unchanged):** e ≥ 0.5, a ≥ 2 px, S/N ≥ 10, r ≤ 50″, g_min ≥ 0.2, aligned ≤ 30°,
anti ≥ 60°, z-margin 0.1, σ_ε 0.26. No image list (map models): the control is the background alignment fraction.

**Outcome rule (frozen; `arcs_macs0647.outcome(primary, second)`, already unit-tested, reused unchanged):**
1. the primary decides if its control passes (≥ 70 % of strong-shear background rows aligned; n = 0 fails);
   otherwise the second decides if its control passes; otherwise **no result**;
2. the deciding run **replicates** if P(≥ n_anti) < 0.05 under the Poisson shape-noise expectation, else
   **does not replicate**.

**Pre-registered secondary checks (do not change the outcome):** background counts by g_min < 1 and ≥ 1; anti rows
of the deciding run get an F200W cutout looked at by eye; the other model's classes of the same rows.

**What changes the plan:** replication sends the anti rows to the S1063 vetting chain (ICL subtraction, measurement
noise, clump-pair census) before any reading; non-replication adds an independent non-replication; no result
leaves the tally unchanged.

**Data:** map members are read from the 4.30 GB Zenodo zip by HTTP range requests (no full download): three
Lenstool maps (100 MB each) and three GLAFIC maps (23 MB each), plus the 24 MB photo-z tarball and the 1 MB
catalogue, cached in the gitignored data root and deleted after the run. Nothing over 200 MB per file.
