# MACS1423 `arcs` anti rows, vetting step 2: clump-pair census of the 55 strong-shear background rows (plan, 2026-10-16)

Written and committed before the census runs. Parent: `RESULT.md` (step 1: 2197, 2322, 4234 stay anti after ICL
subtraction; 1818 becomes `mixed`). Model of the test: the S1063 census, `results/arcs/abells1063-cats/vetting/
PLAN_census.md` on #273, whose primary test gave P 0.031 there.

- **Code:** `scripts/arcs_macs1423_census.py`. It imports #273's `arcs_vet_census.py` (and through it
  `arcs_vet_deblend.py`, `arcs_vet_clumps.py`, `arcs_vet_noise.py`, `lens_consistency.py`) **unchanged at commit
  `43fbbbd`**, extracted to a directory given by `--s1063-scripts`, and reuses its frozen functions
  `deblend_row`, `row_p_anti`, `clopper_pearson_upper` and `outcome`. Only the image and the anti sets differ.
- **Image:** `s3://stpubdata/jwst/public/jw01208/L3/t/o010/jw01208-o010_t005_nircam_clear-f200w_i2d.fits` (the
  screen's F200W image), 8″ S3 byte-range cutouts.
- **Rows:** all 55 rows of `../macs1423-canucs-maps/bgsub_shapes.ecsv`; g_min from `arcs.ecsv` (the screen's
  photo-z interval, `z_basis = photo-z` for all 55).
- **Anti classes at each row's photo-z (the deviation in step 1 is not repeated):** primary k = 4, the screen anti
  rows (1818, 2197, 2322, 4234; `orientation_class`, photo-z basis). Secondary k = 3, the screen anti rows that
  stay anti after ICL subtraction (2197, 2322, 4234). 1758 (anti only at z_s = 2) is not counted.
- **Rules (frozen, #273's, unchanged; all ASSUMPTIONs):** deblend `two-knot` if ΔBIC ≥ 10, both knot S/N ≥ 3,
  separation ≥ 0.10″; own-flux gate with 20 + 20 injections, ≥ 80 % correct each; per-row P_anti: P_cp(g_min) for
  gate-passed `two-knot`, P_sn(g_min; σ_ε 0.26) for gate-passed others, the f̂ mixture for gate-failed rows;
  Poisson-binomial P(≥ k). Fit radius 0.4″, 1″ median filter, seed 20261016 + row index.
- **Outcome (frozen, `outcome()` of #273):** fewer than 12 gate-passed rows → `undetermined`; primary P ≥ 0.05 →
  `anti excess explained by clump orientation (exploratory)`; P < 0.05 → `excess remains under the clump-aware
  null`. Sensitivity (f̂ and its 95 % upper mixture) reported, not used.
- **Pixel scale check (before reading the outcome):** #273's gate and fit radius are in arcsec; the script records
  the cutout pixel scale. If the deblend's known case fails here (the gate pass rate below 50 % of rows), the
  census is reported as `undetermined` for this image, whatever `outcome()` says.
- **What changes the plan:** explained → the MACS1423 replication is closed as a clump-orientation effect (two
  independent fields then fail to replicate S1063 under the clump-aware null). Excess remains → the second band
  (F150W/F277W) and photo-z checks of 2197 and 2322 (z_phot 3.42 both, 6″ apart) run next.
- **Labels:** pixels observed; fits and classes derived; injections simulated; P values model_prediction. The
  sample test is exploratory (the census idea is post hoc on S1063); MACS1423 never tuned it.
