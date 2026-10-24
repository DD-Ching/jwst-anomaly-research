# Abell S1063: WSLAP+ JWST 2025 extra-image test with the WLS source and the DJA v7 reference (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any DJA row, F200W pixel or counter-image prediction of Abell
S1063 was read for this test. Follows `../s1063_wslap/RESULT.md` (D-034 addendum al): gate K (mean source) fails as
written (rms 4.197″, entirely system 10); the WLS statistic declared before that run gives rms 0.581″, 100 % within
3″. Script `scripts/s1063_extra.py`; tests `tests/test_s1063_extra.py`. Route: El Gordo stage 2b
(`../elgordo_djaref/PLAN_2.md`, `PLAN_2b.md`: DJA v7 grizli total fluxes, pipeline-mosaic local noise, gate C2
shared) with the WLS source and frame check F of the Bullet Cluster (`../bullet_extra/PLAN.md`).

**Deviation (labelled, decided before any value of this stage is read).** Gating on the WLS source after the
mean-source gate failed is a post-failure deviation (the addendum k precedent; Bullet `PLAN_1b.md`). The WLS number
was declared in `../s1063_wslap/PLAN.md` rule 4 before that run, and this plan cites it unchanged. Gate K-WLS here is
a fit-reproduction check: WSLAP+ used BUFFALO identifications that include systems at z 0.73 and 1.261, so the 36
images are not held out.

**Disclosed before freezing (positions and file metadata only).** MAST: two NIRCam level-3 F200W products lie within
0.02° of the cluster centre, GO 1840 `o010_t003` (1.2 ks) and GO 3293 `o001_t001` (TARGPROP ACO-S-1063, 70.1 ks);
the GO 3293 SCI grid (10089 × 4775, `MJy/sr`) contains all 36 z_spec images (array bounds of the header WCS). DJA v7
`abells1063-grizli-v7.0-fix_phot.fits` exists (4,389,120 B, Last-Modified 2023-09-26); no header card or row read.
`../s1063_wslap/summary.json` and its residual tables were read (the gate-K run).

## Hypothesis, ordinary explanations, what would change the plan
- Question 0 (frame): do the Chan+2020 positions and the DJA v7 catalogue share a frame (frame check F)? The WSLAP+
  maps reproduce the Chan positions (gate K-WLS), so the model is taken in the Chan frame with no offset.
- Question 1 (method): does the DJA v7 reference with grizli total fluxes and local noise on the GO 3293 F200W
  mosaic pass gate C2 on Abell S1063?
- Question 2 (science, only if F, K-WLS and C2 pass): do the WSLAP+ maps (WLS source) predict a magnified, unmatched
  counter-image of a z_spec system, in the footprint, where F200W shows no source at the predicted flux (flag)?
- Ordinary explanations for a flag: model μ or position error (WSLAP+ is free-form; system 10 straddles a caustic
  under the mean source), low surface brightness or deblending, cluster-member / ICL light (the pipeline mosaic is
  not bCG-subtracted), a reference flux including unrelated light, the point-source `tot_corr`, an identification
  that differs from WSLAP+'s BUFFALO list. Each goes to `/vet-candidate` before any wording beyond "flag". Exotic
  readings are not considered at this stage.
- Expected: F passes or fails by ≲ 0.2″ (HFF-frame list vs a Gaia-aligned DJA mosaic; the S1063 CATS entry carries no
  frame offset > 0.1″, D-034); C2 passes or fails as on El Gordo (passed) / SMACS (failed); 0–3 flags.
- F fails → stop, nothing else read. K-WLS not reproduced → stop. Rule 4a / 4b fails → stop. C2 fails → the route
  is not validated on S1063; no flag is read.

## Inputs
- model_prediction: WSLAP+ JWST 2025 Case3 (Zenodo 10.5281/zenodo.17726321), the two α members only, re-fetched
  by range reads with `s1063_wslap.fetch_member` (CRC-checked), reading A (arcsec at D_ds / D_s = 1), z_lens 0.348,
  flat H0 70, Ωm 0.3 as in `../s1063_wslap/PLAN.md`. Transfer ≈ 262 MB (> 200 MB): the reason stated there holds
  unchanged (the deflection maps are the only released form of the model); session scratchpad only, deleted after
  the run.
- observed: the 36 Chan+2020 z_spec images (`s1063_wslap.spec_images`), z = z_spec. **ASSUMPTION:** z fixed.
- observed: DJA v7 `abells1063-grizli-v7.0-fix_phot.fits`, read into memory, sha256 recorded, not stored.
- observed: GO 3293 `jw03293-o001_t001_nircam_clear-f200w` `_cat.ecsv` (shapes; session cache, deleted) and
  `_i2d.fits` (stamps by S3 range reads; the SCI extension once by a range stream, ≈ 193 MB in memory).

## Frozen rules (computable; functions unit-tested before the run)
1. **Frame check F** (`bullet_extra.frame_check_f`, `min_n` = `s1063_extra.F_MIN_N`): nearest DJA row (finite
   `ra`, `dec`) within 0.5″ of each of the 36 images; passes iff ≥ 18 images are matched (ASSUMPTION: half of the
   list; the Bullet's 30 was for 135 images) **and** |median dx| ≤ 0.1″ **and** |median dy| ≤ 0.1″ (dx along
   RA cos δ). Fails → stop.
2. **Gate K-WLS repeated** (`s1063_wslap.k_block`, 0.25″ grid; `s1063_extra.k_repeat`): passes iff the repeat
   reproduces `../s1063_wslap/summary.json` `reported.K_wls_not_a_gate` (same `n_images`, rms within 1e-6″, same
   matched fraction) **and** passes (rms ≤ 1.0″ **and** ≥ 80 % within 3.0″). Else stop. The WLS image-plane table is
   written as `imageplane_abells1063_abells1063-wslap.ecsv` (its `z_used` feeds rule 3).
3. **Frame tables:** `m0417_extra.run_field` with `bt_hook = a2744_uncover.wls_backtrace`, GO 3293 `o001_t001`,
   rules otherwise unchanged; the model enters through `s1063_extra.install_model` (no frame offset). Its pipeline
   forced test and gate C are reported, not used.
4. **Reference catalogue** (El Gordo `PLAN_2b.md` rule 2 unchanged; `elgordo_djaref`):
   - 4a unit: `flux_unit_jy(with_tot_corr=True)`; 4b grid: `grid_rotation_deg` (|rotation| ≤ 1°, determinant < 0).
     Either fails → stop.
   - flux = `f200w_flux_aper_1 × flux_auto / flux_aper_1 × f200w_tot_corr` × unit (`to_cat_dja`), masks and shapes
     as there; ext = true for every row (ASSUMPTION). Depth by `m0417_canucsref.depth_mag`.
5. **Stage 3:** `a2744_uncoverref._reference_run` (= `m0417_arcref.run` with frame offset (0, 0), then
   `m1423_canucsref.apply_local` on the GO 3293 F200W SCI in nJy per pixel, constants unchanged: 0.3″ apertures,
   200 draws at 1.5–4″, 0.5″ exclusion around DJA positions, seed 20261022, 5 σ / 10 σ), containing-source rule on.
6. **Gate C2** (`m0417_arcref.gate_c2_shared`): ≥ 8 testable controls, ≥ 80 % recovered / confused overall **and**
   ≥ 80 % on the arc-referenced controls with ≥ 3 of them.
7. **Flag:** a predicted target with `class_local` = `absent` while F, K-WLS and C2 pass.
8. Reported, not rules: the same without the containing-source rule (`without_contain/`), plain `gate_c2`, the
   classes without the local-noise rule, and the pipeline gate C of rule 3.

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged. Any flag →
contact sheet first, then `/vet-candidate` in a later cycle; no interpretation here.
