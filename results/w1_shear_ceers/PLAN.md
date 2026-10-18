# W1 shear screen in CEERS (EGS, program 1345, no lens model): plan (frozen 2026-10-18, before any run)

**Question.** The sixteen-field W1 shear volume limit (`results/w1_shear_goodsn/RESULT.md`) is 1.82× weaker than
Takahashi & Asada (2013) at 2 × 10¹³ M☉ over 119.8 arcmin². CEERS has eleven F200W level-3 pointings in MAST
(2,835 s each) and a DJA v7.4 photo-z catalogue of the same area. They are the next blank-field W1 area (GOODS-N
handoff, TASKS "JWST focus" W1).

**Hypothesis being limited.** As `results/w1_volume/PLAN.md` (`hypothesis`: W1 lenses uniform in comoving volume).
A screen plus a limit; a peak that passes the E/B rule is a flag to vet, not a candidate.

**What would change the plan.** If the combined V_eff brings n₉₅ at 2 × 10¹³ M☉ below 10⁻⁴ h³ Mpc⁻³ (the Takahashi &
Asada value), W1 has matched the published bound and the next step is lower masses, not more area. If the CEERS
pointings add < 20 % to V_eff at every mass, blank-field additions stop (T4: cost per gain).

**Inputs (frozen).**
- Pointings (program 1345, F200W level-3 `_cat.ecsv`, MAST): footprints from the MAST `s_region` polygons
  (`footprints.json`, rasterised at 1.2″; queried 2026-10-18). o052_t022 is left out (it overlaps o002 by 6.2 and
  o064 by 2.7 arcmin²). The ten others are kept: o001_t021, o002_t022, o003_t023, o004_t024, o061_t043, o064_t044,
  o066_t055, o068_t048, o069_t046, o071_t047. Their pairwise overlaps are ≤ 0.58 arcmin² (the GOODS-N rule:
  ≤ 0.6) and sum to 3.1 arcmin², which the headline double-counts (ASSUMPTION; ≈ 3 % of the CEERS area).
- Photo-z: DJA `ceers-full-grizli-v7.4-fix.eazypy.zout.fits`, streamed out of the 397 MB
  `ceers-full-grizli-v7.4-fix.photoz.tar.gz`, which is never stored (stated reason: the only photo-z catalogue of
  these mosaics in the DJA format the loader reads; cloud disk rule; D-051 precedent). Tarball and member sha256 are
  recorded at the first fetch. Matched to the catalogue within 0.3″ (`attach_photoz`).
- Gaia spike veto (D-043) at each pointing's footprint centroid; streak masks from `w1_streak_mask.py` on each
  F200W `_i2d.fits` (frozen rule, `results/w1_streak_veto/PLAN.md`).

**Gate P (computable; before the screen).** `inject_shear.passes_gate_p`: a pointing is screened only if ≥ 50 % of
its catalogue rows that pass the screen's S/N and resolution cuts have a photo-z match. Fail → dropped and reported.

**Responsivity (ASSUMPTION).** R = 0.414 as in D-099 (Abell 370 F200W fit; program 1345 has no cluster pointing).
CEERS F200W is 2.8 ks per pointing, shallower than JADES 1181 (3.1–5.7 ks). Sensitivity (reported, not used in the
headline): the injections at Z = 0.3 repeated with R = 0.231.

**Screen (computable).** As `results/w1_shear_goodsn/PLAN.md`: `NullLensModel` at the pointing centroid; Schirmer
10″, 1″ grid, rotation null, E/B rule, streak veto, grid half-width 200″. **Screen rule:** unchanged E/B rule
(p_random < 0.05 and S > the field's max |S_×|). A flagged pointing gets a contact sheet and the frozen vetting rule
(`w1_shear_more` addendum A, V1 ∪ V2), and stays out of the headline either way (D-098).

**Injections (frozen).** `inject_shear.py --lens-z Z` on each pointing that passes gate P, Z ∈ {0.1, 0.3, 0.6, 1.0,
1.5}, masses 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉, 100 lenses per cell, seed 50. **Combined limit:** `w1_volume.py`
over the sixteen headline summaries plus the unflagged CEERS pointings, reported beside the sixteen-field values.

**ASSUMPTIONs.** Every threshold above; single-plane lenses as in `results/w1_volume`; R from a cluster fit in
another program applies to CEERS (same band and pipeline; shallower depth, see above).
