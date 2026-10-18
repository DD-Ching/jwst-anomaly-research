# Same-source pairing: feasibility test on Abell 2744 (plan, fixed before any score is computed)

Owner instruction 2026-10-18 (after #323). Written 2026-10-18 03:07Z, after loading the inputs (row counts only:
66,172 catalogue rows, 149 published images) and before any pair score, match or rank was computed.

## Question
Can a cheap, model-free step **propose** which catalogue sources are images of the same background source (task 1),
at cluster separations, with a false-match burden small enough to review? Testing a given pair against ordinary
lensing (task 2) is a separate step; the repository already does it with published models (`lensmodel.find_images`,
`lens_consistency.py images`), and the morphological test of task 2 (`parity_test.py`) stays paused (D-032).

## Inputs (reused, unchanged)
- Catalogue: DJA grizli v7.2-fix Abell 2744 photometry + eazy-py zout, via `orphan_pairs.load_catalogue("abell2744")`
  (sha256-pinned; 233.5 MB streamed to the session scratch and deleted after use).
- Identities for scoring only: the Bergamini+23 multiple-image list (`lens_consistency.image_list
  ("abell2744-bergamini23")` + `shift_images`), 149 images. The retrieval step never reads it.
- Eligibility of every catalogue source (identity-blind): `orphan_pairs.select_sources(cat, 0.308)` unchanged
  (use flag, not BCG, summed LW S/N ≥ 10, S/N ≥ 10 in ≥ 8 bands, z025 > 0.408).

## Sample
- A published image is matched to the nearest catalogue row within 0.3″; it is in the sample only if that row is
  eligible and no other published image matches the same row. Fraction surviving is reported (completeness of the
  eligibility step, not of the method).
- Positive pairs: in-sample images of the same system at 2–60″ (pairs < 2″ are excluded: same-galaxy knots dominate
  there; their number is reported).
- Split by system: holdout if `zlib.crc32(system_id) % 2 == 0`, else development. Any adjustment is made on the
  development half only; the holdout is evaluated once.

## Retrieval (task 1), identity-blind
- Query = each in-sample image's catalogue row; candidates = every other eligible row at 2–60″.
- Allowed transformation: achromatic magnification only (free flux scale) plus the 3 % error floor; no dust, no
  differential magnification, no colour shift, no time evolution, no redshift freedom beyond the photo-z intervals.
- Score: `orphan_pairs.match_table` (SED χ² with free scale; p_sed; ≥ 8 shared bands; photo-z 16–84 overlap).
  Rank = candidates ordered by `match` first, then p_sed descending.
- Metrics per query with ≥ 1 positive counterpart: rank of the best true counterpart; recall@1/5/10; number of
  non-counterpart candidates passing `match` (false-match burden; candidates that are published images of another
  system are counted separately, since they are real images but not this source).

## Negative control
- For each positive pair (A, B): a decoy C, an eligible row at 2–60″ from A, > 1″ from every published image, with
  |Δm| ≤ 0.5 mag (summed F277W+F356W+F444W flux), |Δ colour| ≤ 0.15 dex (log10 LW/SW flux, `orphan_pairs.source_colours`)
  and aperture radius within ±30 % of B; the decoy minimising (Δm/0.5)² + (Δc/0.15)² + (Δr/0.3)² is taken. No decoy →
  pair not used in the control.
- Report the `match` pass rate of decoys (A, C) vs positives (A, B).

## Morphology (optional add-on, cost-capped)
- `parity_test.compare` on the UNCOVER F200W mosaic (as in `results/parity/`), plain statistic, 1″ halo filter and
  segmentation masks: score = min over both handedness families of χ²/n_pix. Allowed: affine maps with singular values
  0.2–5, either handedness, shift ≤ 0.3″, free flux scale and constant.
- Run on development positives and decoys first. Only if the development AUC (positive lower than decoy) ≥ 0.75 is it
  run on the holdout; otherwise it is recorded as not discriminating and not used.
- Cost cap 30 min wall on 4 processes; stopped work is recorded as not run.

## Pass conditions (holdout; all fixed now)
- P1: recall@10 ≥ 0.5.
- P2: median false-match burden per query ≤ 10, and decoy `match` rate ≤ 0.2 while positive `match` rate ≥ 0.6.
- Capability "passes in this limited scope" only if P1 and P2 pass. A pass is a feasibility statement for bright
  (S/N-selected) images in one cluster; it is not a measured false-pair rate for any other field or the whole sky.
- Fail → branch B of the owner instruction (record why; return to the fixed transient tools).

## Cost cap
- Retrieval and control: 15 min CPU. Total for this test including morphology: 60 min wall.
