# W1 shear screen in GOODS-N (JADES 1181, no lens model): plan (frozen 2026-10-18, before any run)

**Question.** The fourteen-field W1 shear volume limit (`results/w1_shear_ncf/RESULT.md`) is 2.35× weaker than
Takahashi & Asada (2013) at 2 × 10¹³ M☉ over 97.5 arcmin². JADES GOODS-N has six F200W pipeline mosaics with small
mutual overlap (≈ 62 arcmin² of MAST footprint) and DJA v7.3 photo-z already pinned in the repository
(`orphan_pairs.DEEP_FIELDS["goodsn-dja"]`, D-051). They are the next blank-field W1 area (TASKS "JWST focus" W1).

**Hypothesis being limited.** As `results/w1_volume/PLAN.md` (`hypothesis`: W1 lenses uniform in comoving volume).
A screen plus a limit; a peak that passes the E/B rule is a flag to vet, not a candidate.

**What would change the plan.** If the twenty-field V_eff brings n₉₅ at 2 × 10¹³ M☉ below 10⁻⁴ h³ Mpc⁻³ (the
Takahashi & Asada value), W1 has matched the published bound and the next step is lower masses, not more area. If
the six pointings add < 20 % to V_eff at every mass, blank-field additions stop (T4: cost per gain).

**Inputs (frozen).**
- Pointings (program 1181, F200W level-3 `_cat.ecsv`, MAST, chosen from the s_region footprints on 2026-10-18 so
  that pairwise overlaps are ≤ 0.6 arcmin²): o001_t001, o003_t007, o004_t008, o005_t008, o009_t009, o010_t009.
  Left out because they overlap a chosen pointing by ≥ 1.1 arcmin²: o002, o006, o007, o011, o098. The one
  remaining overlap (o004/o005, ≈ 0.5 arcmin² of footprint) is double-counted area (ASSUMPTION; < 1 % of the total).
- Photo-z: DJA `gdn-grizli-v7.3-fix.eazypy.zout.fits` (member sha256 pinned), streamed out of the 371 MB tarball,
  which is never stored (D-051 stated reason; cloud disk rule). Matched to the catalogue within 0.3″
  (`attach_photoz`).
- Gaia spike veto (D-043) at each pointing's footprint centroid; streak masks from `w1_streak_mask.py` on each
  F200W `_i2d.fits` (frozen rule, `results/w1_streak_veto/PLAN.md`).

**Gate P (computable; before the screen).** A pointing is screened only if ≥ 50 % of its catalogue rows that pass
the screen's S/N and resolution cuts have a photo-z match (MACS0717 had 17 % and was nearly blind). Fail → that
pointing is dropped and the drop is reported.

**Responsivity (ASSUMPTION).** No shear to fit R in a blank field. R = 0.414: the smallest F200W R among the cluster
fits in `results/w1_volume/summaries` (Abell 370), leaving out MACS0647 (R 0.23, unexplained and nearly blind). A
smaller R paints a weaker W1 signal. Program 1181 is 3.1–5.7 ks per pointing in F200W, comparable to the CANUCS
pointings. Sensitivity (reported, not used in the headline): the injections at Z = 0.3 repeated with R = 0.231.

**Screen (computable).** `NullLensModel` at the pointing centroid; every other setting the D-053 default as in
`results/w1_shear_ncf` (Schirmer 10″, 1″ grid, rotation null, E/B rule, streak veto, grid half-width 200″).
**Screen rule:** unchanged E/B rule (p_random < 0.05 and S > the field's max |S_×|). A flagged pointing needs a
contact sheet and vetting and is kept out of the headline unless the frozen vetting rule (`w1_shear_more`
addendum A, V1 ∪ V2) explains it.

**Injections (frozen).** `inject_shear.py --lens-z Z` on each pointing, Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5}, masses
2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉, 100 lenses per cell, seed 50 (as `results/w1_shear_ncf`). **Combined limit:**
`w1_volume.py` over the fourteen stored headline summaries plus the unflagged GOODS-N pointings, reported beside the
fourteen-field values.

**ASSUMPTIONs.** Every threshold above; single-plane lenses as in `results/w1_volume`; R from cluster fits applies to
GOODS-N (same band, pipeline and comparable depth).
