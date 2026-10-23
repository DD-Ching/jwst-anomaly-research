# RXJ2129 reference-catalogue extra-image test

## Stage 1 (2026-10-23): gate K passes; the DJA v7.0 reference carries every stage-2b input

Frozen plan: `PLAN_1.md` (committed and pushed before the run). Script `scripts/rxj2129_djaref.py`; numbers from
`summary.json`; tables `dhanasingham26_images.ecsv` (observed), `imageplane_rxj2129_rxj2129-lenstool-maps.ecsv`
(model_prediction), `readings_rxj2129.ecsv`.

- **Inputs:** the Dhanasingham+2026 JWST-frame image table (arXiv:2604.11882 source, Appendix B), 31 images in 10
  families, all at their family's spectroscopic z. The transcription matches the TeX table row by row (checked in
  the session against the arXiv source; not stored) and is unit-tested (`tests/test_rxj2129_djaref.py`). The
  paper's text says 11 systems; the table has 10 families (not resolved).
- **Rule 1, gate K** (`rxj2129-lenstool-maps`, mean back-traced source, 0.25″ grid): **passes**, rms 0.449″,
  31 / 31 within 3″ (thresholds ≤ 1.0″ and ≥ 80 %). Per-family rms 0.19–0.78″, largest for the S3-1 / S3-2 knots
  (0.78″, 0.66″) in `per_family_rms_arcsec`. Own-constraint validation only: the model was fitted to these
  positions (the paper quotes rms 0.32″ for its own estimator, which is not this one).
- **Rule 2, reference inputs:** DJA v7.0 `rxj2129-grizli-v7.0-fix_phot.fits` (9,152,640 bytes, 3,295 rows, 391
  columns; headers only, HTTP range) has every column and card SMACS rule 2 named; `APER_1` 12.5 px. **Passes.**
- **Stage-2 inputs, not read** (`extra_predicted`, catalogued / predicted images per family): 2 (3 / 5), 3 (4 / 7),
  4 (3 / 5), 5 (3 / 4), S3-1 (3 / 5), S3-2 (3 / 5). No family has fewer predicted than catalogued images.
- **Data:** the three pinned Lenstool members (≈ 300 MB in the cache, sha256 verified by `stage_members`) were
  range-read out of the 4.30 GB Zenodo zip in 250 s (not stored); kept in the session cache for stage 2 and deleted
  before the session ends (never committed). No archive tar was stored.

**Next (done, stage 2 below):** stage 2 under a frozen plan with El Gordo stage 2b's rules on this reference (`flux_auto × tot_corr`
unit check, grid check, frame tables with `m0417_extra.run_field` on an F200W tile containing all 31 images
(program 2767 `o002_t001` first; footprint check before any read), gate C2 with shared pairs, local noise). The S3-1 / S3-2 knots are separate families here, as in the paper's constraint list, but they lie
inside system 3's arcs: the stage-2 plan must state how their extra predictions and system 3's are counted (one
target per predicted region, or the shared-pair rule) before any read (/code-review finding).

## Stage 2 (2026-10-23): gate C2 fails; no flag is read

Frozen plan: `PLAN_2.md` (pushed in 63db430 before any read). Run `scripts/rxj2129_djaref.py --stage 2`; numbers
from `summary.json` `stage2`; tables `control_rxj2129.ecsv`, `predicted_rxj2129.ecsv` (rule 3), `without_contain/`
(rule 7), `frame/` (rule 1).

- **Rules 2a / 2b pass:** u = 23.9 (µJy, 3,255 rows); grid rotation −2.0 × 10⁻⁵°, determinant negative, pixel
  scale 0.040″. DJA table sha256 `0a391f6f…`, 9,152,640 bytes, 3,295 rows (in memory, not stored).
- **Rule 1, frame tables** (reported, not used): pipeline forced gate C passes, 10 / 11.
- **Gate C2 (rule 5) fails:** 4 / 7 (0.57) testable controls and 2 / 5 (0.40) arc-referenced, against 0.8 and the
  minimum counts 8 and 3. Control classes: 3 absent, 3 recovered, 1 confused, 8 undetectable_local,
  8 inconsistent_reference, 4 no_reference. Also fails without the local-noise rule (8 / 15, arc 3 / 8) and without
  the containing-source rule (3 / 7, arc 2 / 5).
- **No flag is read** (rule 6 needs C2). Target classes are reported only: 4 absent, 2 recovered, 1 confused,
  2 inconsistent_reference, 2 undetectable_local; D1 joined two knot-group sets (`knot_sets`), not read further.
- **Post hoc (not rules):** the 3 `absent` controls (systems 1, 2, 7) all take one reference row, DJA 1506
  (`ref_mag` 16.38, `ref_a3_arcsec` 8.7), so their predicted magnitudes are 16.7, 18.1 and 16.8
  (`control_rxj2129.ecsv` `m_pred`). Contact sheet `sheet_2.png` (`--stage sheet`, the absent / ambiguous /
  confused controls) looked at: a faint source sits in the circle of the system 1 and 2 controls; the system 7
  circle is beside a bright galaxy. Without row 1506, C2 would be 4 / 4 (arc 2 / 2; computed in the session, not
  stored), still below the minimum counts, so the route cannot validate itself here.
- **Data:** F200W SCI extension (≈ 82 MB) and the DJA table read into memory; pipeline `_cat.ecsv` and the Lenstool
  members (≈ 300 MB, re-staged by range reads) in the session cache, deleted after the run.

**Reading:** the reference route is closed for RXJ2129 (D-034 addendum u). Across four fields it validated on
Abell 2744 and El Gordo and failed on SMACS 0723 and RXJ2129.
