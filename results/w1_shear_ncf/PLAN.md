# W1 shear screen in the CANUCS NIRCam flanking fields (no cluster model): plan (frozen 2026-10-18, before any run)

**Question.** The eleven-field W1 shear volume limit (`results/w1_volume/RESULT.md`) is 3.6× / 8.5× / 65× weaker
than Takahashi & Asada (2013) at 2 × 10¹³ / 8 × 10¹² / 2 × 10¹² M☉; at similar efficiency, reaching the bound at
2 × 10¹³ M☉ needs ~3.6× the area. The five CANUCS NIRCam flanking fields (NCF; program 1208, 3–6′ from the cluster
cores, catalogue μ 1.0–1.4; D-051) are the cheapest extra area with photo-z already pinned in the repository. They
have no strong-lensing model, so the screen runs without one: the measured ellipticities are used as they are.

**Hypothesis being limited.** As `results/w1_volume/PLAN.md` (`hypothesis`: W1 lenses uniform in comoving volume).
This is a screen plus a limit; a peak that passes the E/B rule is a flag to vet, not a candidate.

**What would change the plan.** If the sixteen-field V_eff brings n₉₅ at 2 × 10¹³ M☉ within 2× of 10⁻⁴ h³ Mpc⁻³,
blank deep fields (PRIMER, CEERS, COSMOS-Web) become the next W1 step. If the five NCF add < 20 % to V_eff at every
mass, blank-field W1 additions stop (T4: cost per gain).

**Inputs (frozen).**
- Fields and MAST level-3 catalogues (program 1208): macs0416-ncf `jw01208-o023_t002`, macs1149-ncf
  `jw01208-o029_t004`, abell370-ncf `jw01208-o020_t001`, macs0417-ncf `jw01208-o026_t003`, macs1423-ncf
  `jw01208-o032_t005` (the observations `orphan_pairs.DEEP_FIELDS` uses).
- **Band: F150W** (`*_nircam_clear-f150w_cat.ecsv`). **Deviation from the cluster runs (F200W):** the NCF
  observations have no F200W (S3 listing 2026-10-18: F090W, F115W, F150W, F182M, F210M and long-wavelength bands
  only). F150W is the deepest short-wavelength wide band there.
- Photo-z: the CANUCS DR1 NCF photometry catalogues (sha256 pinned in `orphan_pairs.DEEP_FIELDS`), converted with
  `inject_radial.canucs_zout`.
- Gaia spike veto (D-043) at the catalogue's median position; streak masks from `w1_streak_mask.py` on each F150W
  `_i2d.fits` (same frozen rule, `results/w1_streak_veto/PLAN.md`).

**Screen without a model (computable).** A `NullLensModel` (tangent plane at the catalogue's median RA/Dec,
κ = 0, g = 0 everywhere, `z_lens` 0) replaces the cluster model. Every other screen setting is the D-053 default as
in `results/w1_volume` (Schirmer 10″, 1″ grid, rotation null, E/B rule, streak veto). Grid half-width 200″ (covers
the NCF footprint; centres without enough sources drop out as in every field). The screen rows are all lensable
rows (z160 > 0.1 or no photo-z). The responsivity R is not fitted (no shear to fit against); it enters only the
injected images.

**Gate R (computable; before any NCF run).** R for F150W isophotal moments is fitted with the unchanged
`shear_responsivity` on the CANUCS cluster pointings of the same program in F150W: macs0417 `jw01208-o006_t003`
and macs1423 `jw01208-o010_t005`, with their `results/w1_volume` models and photo-z. Pass: both fits return R with
R > 3σ (the function's own rule) and |R_F150W − R_F200W| ≤ 0.15 in each field (R_F200W from
`results/w1_volume/summaries/z0.3_<field>.json`: 0.488, 0.604). Adopted R = the smaller of the two F150W values
(conservative: a smaller R paints a weaker W1 signal). Fail → stop; no NCF screen.

**Screen rule.** Unchanged E/B rule (p_random < 0.05 and S > the field's max |S_×|). A field passing it is flagged;
it then needs a contact sheet and vetting before any limit is quoted from it.

**Injections (frozen).** `inject_shear.py --lens-z Z` on each NCF, Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5}, masses
2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉, 100 lenses per cell, seed 50 — as `results/w1_volume` addendum A, so the cells
combine. **Combined limit:** `w1_volume.py` over the eleven stored cluster summaries plus the five NCF summaries
(sixteen fields), same integral, reported beside the eleven-field values.

**ASSUMPTIONs.** No cluster shear removal in the NCF (the residual cluster shear at 3–6′ is smooth across a 10″
aperture; the aperture mass of a constant shear is zero); R from the cluster pointings applies to the NCF (same
program, band and depth); single-plane as in `results/w1_volume`; every threshold above.
