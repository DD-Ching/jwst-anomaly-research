# SPT-CL J0546-5345 LTM input check, unit check U and free-z gate K (2026-10-24)

Next field of the fitted-model route (RXCJ2211 handoff). Model: Allingham+2025 LTM (arXiv:2507.08949, ApJL 990, 1),
Zenodo 10.5281/zenodo.16534337 (CC-BY-4.0): α maps at z_s 3.5, a Lenstool-format image list, a spec-z catalogue.
Plan frozen in aa0a8ef, functions and tests in 4da5b19, frozen run in 69e014c. Script `scripts/spt0546_inputs.py`;
numbers from `summary.json` (frozen) and `posthoc.json` (post hoc); image list `spt0546_images.ecsv`; tests
`tests/test_spt0546_inputs.py`.

| Check (PLAN rule) | Found | Status |
|---|---|---|
| Image list | 59 rows: 35 main images in systems 1–10, 24 `c` candidate rows; no redshift column values | parsed |
| Map header | 2600 × 2600, 0.04″ / px, RA---TAN / DEC--TAN, no BUNIT, no source-z keyword | read |
| C: coverage | 35 / 35 main images inside the map | **passes** |
| Z: spec-z systems | 14 conf 3 / 4 entries behind the lens; none ≤ 0.5″ from a main image | 0 systems: no spec-z gate possible |
| U: unit, per-system z free | median of the 10 per-system minimum rms: A (arcsec) 27.72″, P (pixels × 0.04″) 1.38″; threshold ≤ 0.5″ | **fails** (no reading accepted) |
| K: free-z gate | not run (U failed) | — |
| R: DJA v7 field | 35 fields listed, none named for the cluster | **fails** |

**Plan count correction.** The plan's listing line says 45 main / 14 candidate rows; the file has 35 / 24 (a miscount
when the plan was written). The rule itself ("the non-`c` rows of systems 1–10") is unchanged and the script applies it.

**Reading.** As frozen, U fails, so gate K does not run and the fitted-model route stops on this field. Reading A is
off by more than an order of magnitude (median 27.7″, systems 3–9 at 18–85″), reading P is not: 3 of 10 systems have
minima ≤ 0.5″ (1, 2, 10: 0.13–0.50″), the other 7 lie at 1.12–2.13″. No lensed image has a spectroscopic redshift
(the paper's redshifts are photometric), and the JWST imaging (GO-5594, SLICE) has no DJA v7 catalogue, so neither a
spec-z gate nor gate C2 could run here even if U had passed.

**Post hoc (added after U failed; no verdict depends on it; `posthoc.json`).**
- Exact image-plane residuals with reading P, each system at its fitted z (the frozen K computation): rms 1.66″,
  91.4 % of the 35 images within 3″; K's rms condition (≤ 1.0″) would fail. Per system 0.13–2.66″.
- The authors report rms ≈ 1.2″ over 39 constraints (paper Sect. III, not a released file; stored as a constant in
  `posthoc.json`). Their z per system and constraint set differ from ours, so the two numbers are not directly
  comparable, but a model with a reported rms above 1″ cannot meet the 1.0″ gate used on the other fields.
- Under reading P, system 3 (the paper's anchor, z_s 3.5) has its fitted z at 3.50, the grid value nearest 3.5
  (derived), which agrees with the record's normalisation statement.
- `c` candidates against the image predicted for their system's mean source (reported only, no flag): c5.6, c6.21,
  c6.22 at 1.10–1.13″, c6.23 1.80″, c4.11 / c4.12 2.90–3.65″, c2.3 and c7.6–8 at 9.7–19.4″.

**Data:** the two α maps (27 MB each) were held in the session scratchpad and deleted after each run (md5 matched the
record, sha256 in `summary.json`); the paper's HTML was read for the rms only. Nothing else stored. No image (cutout or
mosaic) was read; nothing here is a candidate.
