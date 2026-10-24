# Bullet Cluster: CANUCS DR1 extra-image test with the WLS source (plan, frozen)

Frozen 2026-10-24, before any counter-image prediction, CANUCS flux, CANUCS pixel or JWST pixel of the Bullet Cluster
was read for this test (the downloads may start before the freeze; no value is read before it). Follows
`../bullet_gatek/RESULT.md`: gate K fails as written (`PLAN.md`), gate K-WLS passes as a post-failure deviation
(`PLAN_1b.md`, rms 0.391″, 135 / 135 within 3″; D-034 addendum y). Script: `scripts/bullet_extra.py`; tests
`tests/test_bullet_extra.py`. Route: `../abell370_buffalo_extra/PLAN.md` (CANUCS reference, local noise on the
CANUCS bCG-subtracted F200W mosaic) with the WLS source and the shared-pair C2 rule of
`../abell2744_uncoverref/PLAN_2b.md`.

**Disclosed before freezing.** `../bullet_gatek/summary_wls.json` lists 18 families with more predicted than
catalogued images (positions not read for this plan). The 135 gold images are the Lenstool fit's own constraints, so
gate K-WLS validated the maps on their own constraints only. Footprint (positions only): pipeline tile
GO 4598 `o001_t001` F200W (TARGPROP BulletCluster) contains all 135 gold positions in its pixel grid; `o002_t010`
(TARGPROP TOPLRD) contains none.

## Hypothesis, ordinary explanations, what would change the plan
- Question 0 (frame): do the gold list and the CANUCS DR1 catalogue share a frame (frame check F)?
- Question 1 (method): does the CANUCS DR1 reference route with the local-noise and containing-source rules pass
  gate C2 (shared-pair rule) on the Bullet Cluster?
- Question 2 (science, only if F and C2 pass): do the CANUCS best-fit maps (WLS source) predict a magnified,
  unmatched counter-image of a gold system, in the footprint, where CANUCS F200W shows no source at the predicted
  flux (flag)?
- Ordinary explanations for a flag: model μ or position error near critical curves (two gold families, K1b / K1e,
  straddle a caustic: `../bullet_gatek/RESULT.md`), low surface brightness or deblending, cluster-member / bCG
  residual light, a reference image whose catalogue flux includes unrelated light, a wrong system z. Each goes to
  `/vet-candidate` before any wording beyond "flag". Exotic readings are not considered at this stage.
- Expected: F passes (both lists come from the same team's JWST mosaics); C2 passes or fails on arc-referenced
  controls as on MACS1423 / Abell 370; 0–3 flags, each with an ordinary reading.
- F fails → no stage 2 / 3 (no reading). C2 fails → the route is not validated on the Bullet Cluster; no flag is read.

## Inputs
- model_prediction: `bullet-canucs-maps` (Lenstool best-fit 50 mas deflection maps, sha256 pinned in `lensmodel`).
  Download over 200 MB per file (219 / 205 MiB gzip): the reason stated in `../bullet_gatek/PLAN.md` holds unchanged
  (only these maps exist at 50 mas; gate K-WLS validated them); cloud session cache only, deleted after the run.
- observed: the 135 gold images (`../bullet_inputs/bullet_gold_images.ecsv`), z = their system z as in gate K
  (`z_used` of `../bullet_gatek/imageplane_bullet_bullet-canucs-maps.ecsv`). **ASSUMPTION:** z fixed.
- observed: CANUCS Bullet DR1 photometry catalogue `hlsp_canucs_jwst-hst_multi_bullet-clu_multi_v1_photometry-cat.fits`
  and F200W 40 mas bCG-subtracted mosaic `hlsp_canucs_jwst_nircam_bullet-clu-40mas_f200w_v1_bgsub-sci.fits.gz`
  (130,587,001 B), CANFAR vault `canucs/BULLET_DR1`. sha256 pinned in the script (`CAT_SHA`, `IMAGE_SHA`) from the
  file bytes only, before the run; the script refuses to run while they are unset. Not committed; deleted after the run.
- observed: JWST pipeline level-3 GO 4598 `jw04598-o001_t001_nircam_clear-f200w` `_cat.ecsv` (shapes) and `_i2d.fits`
  (S3 byte-range stamps only).

## Frozen rules (computable; functions unit-tested before the run)
1. **Frame check F** (`bullet_extra.frame_check_f`): for each gold image, the nearest CANUCS catalogue row (any row
   with finite RA, Dec) within 0.5″ (ASSUMPTION); offsets dx = (gold − CANUCS) RA·cos δ, dy = gold − CANUCS Dec, in
   arcsec. F passes iff ≥ 30 gold images are matched **and** |median dx| ≤ 0.1″ **and** |median dy| ≤ 0.1″
   (ASSUMPTIONs; 0.1″ is below the 0.18″ of MACS0417 gate F2). Fails → stop.
2. **Gate K-WLS repeated** (`a2744_uncover.imageplane` via `bullet_gatek.run_wls`-equivalent code, 0.25″ grid):
   rms ≤ 1.0″ **and** ≥ 80 % within 3.0″. A different outcome from `../bullet_gatek/summary_wls.json` (map bytes are
   sha-pinned, so none is expected) → stop.
3. **Frame tables:** `m0417_extra.run_field` with `bt_hook = a2744_uncover.wls_backtrace` (every traced row's source
   is its family's WLS source), pipeline tile `o001_t001`, rules otherwise unchanged. Its pipeline forced test and
   gate C are reported, not used.
4. **Reference catalogue:** `m0417_canucsref.to_cat` (frozen CANUCS DR1 mapping, `FLUX_KRON_TOTAL_F200W` in nJy);
   depth `m0417_canucsref.depth_mag`.
5. **Stage 3:** `a2744_uncoverref._reference_run` (= `m0417_arcref.run` with frame offset (0, 0), then
   `m1423_canucsref.apply_local` on the CANUCS F200W bCG-subtracted mosaic in nJy per pixel, pixel scale 0.04″
   checked against its WCS to 1e-4″, constants unchanged: 0.3″ apertures, 200 draws at 1.5–4″, 0.5″ exclusion
   around CANUCS positions, seed 20261022, 5 σ / 10 σ), with the containing-source rule **on**.
6. **Gate C2** (`m0417_arcref.gate_c2_shared`, D-034 addendum i): ≥ 8 testable controls, ≥ 80 % recovered /
   confused overall **and** ≥ 80 % on the arc-referenced controls with ≥ 3 of them.
7. **Flag:** a predicted target (not `observed` / skipped classes) with `class_local` = `absent` while F and C2 pass.
8. Reported, not rules: the same without the containing-source rule (`without_contain/`), plain `gate_c2`, and the
   classes without the local-noise rule.

## Budget and stop
One cycle (≤ 1 h). If the run does not finish, the plan stays frozen and the next cycle runs it unchanged. Any flag →
contact sheet first, then `/vet-candidate` in a later cycle; no interpretation here.
