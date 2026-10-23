# W2: one quasar seen twice, or two quasars? DESI DR1 spectra of rejected lensed-quasar pairs (bounded trial)

Owner brief 2026-10-23 (after #471):
- Budget: 4 h wall clock (start 23:47 UTC), ≤ 2 GB, ≤ 4 cores.
- One route, chosen within 45 min (§0).
- No expert contact, no account, no schedule.
- Frozen at commit, before any spectrum is retrieved or read. The only DESI information seen before freezing is the
  `desi_dr1.zpix` positions, redshifts and spectral types used for the availability count (§0).

## 0. Route choice (≤ 45 min, 23:47–00:15 UTC)
Three options were compared:
- **(a) W2 HST lens light with a joint lens-light + point-source model (D-060 "Revisit if").** Rejected: 5 of the 7
  targets are quads that STRIDES already modelled with their lens galaxies in HST (Schmidt et al. 2023, MNRAS 518,
  1260, arXiv:2206.04696; DESJ0405-3308 has a published deflector). Little new discriminating information, and the PSF
  problem risks tool development.
- **(b) A new JWST epoch for Sunrise n0153 (#61).** Rejected for this round:
  - the only newer NIRCam epoch of the field is VENUS 6882 (2025-07-08);
  - the question (AGN flare vs SN) is ordinary-transient vetting on an owner-review PR.
- **(c) This plan.** D-064 left 11 colour-matched rejected pairs "untestable": telling a binary from a lens needs
  spectra of both images.
  - DESI DR1 (public, 2025) has QSO spectra of both components for 20 rejected systems, 12 known lenses and 17
    Hennawi et al. 2006 binaries with θ ≤ 10″.
  - Availability count (`desi_dr1.zpix`, NOIRLab Data Lab TAP, 2026-10-23, positions and redshifts only) is stored as
    `availability.json`.
  - The known lenses and binaries let the spectral test be calibrated before any target is read.
- **Prior work** (does not settle this): the DESI DR1 quasar-pair catalogue (arXiv:2505.03103) classifies pairs by eye
  (spectral similarity, morphology).
  - J0728+2607: Lemon et al. 2023 (arXiv:2206.07714) report strongly different NOT spectra.
  - J0130+0725: arXiv:2607.00093 calls its spectra nearly identical with mild line-ratio differences and favours a
    dual quasar without excluding lensing.
  - Neither gives a control-calibrated test on the rejected list.

## 1. Question and scope
- **Question:** for each pair, are the two components' DESI spectra consistent with one quasar seen twice (lensing,
  which an ordinary *or* a dark deflector produces) or with two different quasars (binary or projected pair)?
- **Not tested:** whether a deflector is dark. "Consistent with one source" is necessary for a W2 dark lens, not
  sufficient. Ordinary lensing by a faint galaxy and binaries with near-identical spectra stay open.
- **ASSUMPTIONS:**
  - **Lensed images:** they have the same intrinsic spectrum. They may differ by a smooth factor (flux ratio,
    differential extinction, continuum microlensing) and by a changed broad-line equivalent width (microlensing or
    time-delay variability); both are absorbed by the free scale s below.
  - **Binaries:** independent quasars, differing in line profile and velocity.
- **Fibre crosstalk:** DESI fibres are 1.5″, so close pairs cross-contaminate, which makes binaries look alike. The
  separation of every system is reported.
- **Labels:**
  - spectra: observed;
  - statistics: derived;
  - catalogue classes: catalogued;
  - thresholds: ASSUMPTION.

## 2. Sample (frozen rules; `scripts/w2_desi_pairs.py`)
- **Systems:**
  - the D-064 sample from the pinned inputs (`scripts/w12_niq.py` `build_sample`; Lemon et al. 2023 table 1, SQLS
    DR3 / DR5 / DR7);
  - plus Hennawi et al. 2006 binaries (`J/AJ/131/1/binqso`, pinned) with θ ≤ 10″ and both catalogue redshifts with
    |Δv| < 2000 km/s, i.e. physical binaries; projected pairs are left out.
- **Roles:**
  - positive controls = D-064 "control" (catalogued lenses);
  - negative controls = the Hennawi binaries above;
  - targets = D-064 "rejected".
- **DESI spectra:**
  - `desi_dr1.zpix` rows with spectype QSO and zwarn = 0 within (separation + 1.5″) of the system position;
  - spectra from SPARCL, data release DESI-DR1 (anonymous public access; flux, ivar, wavelength of the brz coadd).
  - Per TARGETID one spectrum: the one with the most good pixels.
- **Pair rule:**
  - D-064 systems: two spectra whose separation equals the catalogued separation within 0.5″; if several, the closest
    match.
  - Hennawi: one spectrum within 1″ of each catalogued quasar position.
  - Systems without such a pair are "unknown (no DESI pair)".

## 3. Statistic (frozen; derived)
- **Redshift:** z = the mean of the two DESI redshifts.
- **Lines:** C IV 1549.06, C III] 1908.73, Mg II 2798.75, Hβ 4862.68 and Hα 6564.61 Å (rest). A line is used if its
  window λ0 (1 + z)(1 ± 15000 km/s / c) lies inside 3650–9800 Å, in both spectra.
- **Per line and spectrum:**
  - continuum C = an inverse-variance-weighted straight line through the pixels at |v| ∈ [9000, 15000] km/s;
  - profile p = F / C − 1 with σ_p = σ_F / C, on pixels with |v| ≤ 6000 km/s and ivar > 0 in both spectra (DESI
    coadds share one wavelength grid).
- **Per line, both spectra:**
  - fit p_A ≈ s · p_B with s ≥ 0 free;
  - χ²_ν = Σ (p_A − s p_B)² / (σ_A² + s² σ_B²) / (N − 1);
  - a line counts only with N ≥ 50 pixels and a line S/N of Σ p_B / √Σ σ_B² ≥ 10 in both spectra.
- **D:** the mean χ²_ν over the counted lines. No counted line → "unknown (no usable line)".
- **Δv:** the cross-correlation peak of p_A and p_B over ±3000 km/s, reported only.

## 4. Decision rule and validation (frozen)
- **Threshold T = 3.**
  - Noise alone gives χ²_ν ≈ 1. The factor 3 allows for underestimated coadd errors and the mild line-shape
    microlensing the free scale cannot absorb (ASSUMPTION).
  - D > T reads "two different spectra"; D ≤ T reads "consistent with one source at DESI precision".
- **Gate (controls only, computed and committed before any target spectrum is retrieved):**
  - G1: at least 6 usable positive and 6 usable negative controls, else stop (insufficient controls).
  - G2: at least 80 % of usable positive controls have D ≤ T. False "two sources" on real lenses ≤ 20 %.
  - G3: at least 70 % of usable negative controls have D > T. Power on binaries ≥ 70 %.
  - Wilson 95 % intervals are reported.
- **Gate fails:** the targets are not retrieved or read. The tool result is the deliverable. No threshold, line list
  or statistic change, and no second statistic.
- **Gate passes:** the rejected targets are retrieved and classified with the same code.
  - Pairs with D > T are recorded as "two different spectra (binary or projected)", which closes the dark-lens
    reading for them.
  - Pairs with D ≤ T are listed as "consistent with one source", with their separation, Δv, the D-064 colour match
    and any HST lens-light limit (D-064 addendum).
  - None is a candidate until ordinary explanations are vetted (`/vet-candidate`).
- **Fibre crosstalk check (reported, not gating):** controls and targets with separation ≤ 2″ are listed separately.

## 5. Budget and stop
- **Budget:**
  - downloads ≤ 100 MB (expected ~10 MB);
  - kept data are compact per-system line windows (`spectra/`, each < 1 MB) with SPARCL ids;
  - time cap 03:30 UTC for the analysis, then review and merge.
- **Stop:**
  - G1 or the gate fails;
  - SPARCL refuses or changes the data;
  - the work turns into tool development (e.g. needing a deblended or PSF-corrected spectrum).
- **Out of scope:** no expansion beyond this list if the gate passes, and no new search.

## Errata (added at review, 2026-10-23; the frozen text above is unchanged)
- **Header:** one SPARCL `find` + `retrieve` probe of control lens J0440-0905 (record format only; no flux
  examined) ran at about 00:07 UTC, before this plan was committed (RESULT, Deviation).
- **§0, Sunrise n0153:** "the only newer NIRCam epoch of the field is VENUS 6882" means no public NIRCam epoch of
  WHL0137-08 after VENUS 6882 (2025-07-08) turned up in a MAST query of 2026-10-23 (not stored). VENUS is the epoch
  that showed the brightening.
- **§0, J0130+0725:** the cited paper (Chen et al., arXiv:2607.00093) used the same DESI DR1 spectra, so it is not
  independent of this route.
- **§2:** two rejected pairs share their DESI TARGETIDs with two Hennawi controls; the overlap was not foreseen
  (RESULT §1).
- **§4:** the close-pair list was missing from the first RESULT and was added at review.
