# Orphan-pair (dark-deflector) screen on Abell 2744: result (2026-10-16)

- **Rules:** `PLAN.md`, committed before the run. Every number below is derived from observed catalogues; injected
  images are simulated.
- **Exploratory:** this is the first run of this statistic on this field. Abell 2744 had already been studied with
  other tests, so it is not a blind sample.

## Checks
1. **Alignment.** 163 of the 258 published images (CATS v4.1 and Bergamini+23) have a DJA source within 1″. Their
   median offset is (−0.014″, +0.020″) and the median separation 0.12″, inside the ≤ 0.3″ limit. **Passed.** This was
   computed in the session; the numbers are not stored in `summary.json`.
2. **Positive control.** None of the 174 SED-matched close pairs lies within 1″ of a published image (class
   `published` = 0). **This field has no positive control at 0.3–3″.**
3. **Injection recovery** (`inject_pairs.py --fields abell2744`). The CLI was extended by the field name
   (`CLUSTER`). The cluster lens's magnification of the injected images is not modelled. Cluster fields write only
   per-field efficiencies: the pooled `limits.json` (masses at `Z_LENS_REF`, pooled area) is not valid for this field.

   | type | θ_E 0.15″ | 0.3″ | 0.7″ | 1.5″ |
   |---|---|---|---|---|
   | point: per deflector | 0.00 % | 0.30 % | 0.95 % | 0 (pair separation > 3″) |
   | point: recovered / both selected | 1 / 50 | 17 / 45 | 3 / 27 | 0 / 45 |
   | W2: per deflector | 0.00 % | 0.70 % | 0.75 % | 0 (pair separation > 3″) |
   | W2: recovered / both selected | 0 / 48 | 12 / 43 | 3 / 38 | 0 / 37 |
   | W1: per deflector | 0.05 % | 0.10 % | 2.55 % | 15.2 % |
   | W1: recovered / both selected | 0 / 4 | 2 / 16 | 3 / 11 | 18 / 40 |

   - **These match the deep fields within noise** (D-051: point 0.54 % / 0.84 % at 0.3″ / 0.7″; W2 0.64 % / 0.83 %;
     W1 3.6 % / 18.6 % at 0.7″ / 1.5″).
   - **Completeness is low in every field, not specially here.** Only 1–12.5 % of lensed sources have both images
     selected. Most selected, matched pairs are then classed `same_galaxy` (θ_E 0.15–0.3″) or `visible_lens`
     (0.7″; point 17 of 27), as in the deep fields.
   - **A null in this one field therefore constrains only a small part of a dark-lens population.**
4. **Artifact burden in the top 15.**
   - #12 sits on a mosaic edge with a blank region.
   - #1 lies next to a bright star's spike or halo; one member is barely visible in F277W/F444W.
   - #8 is outside the UNCOVER cutout footprint, so it cannot be judged here.

## Count test (frozen rule: an excess needs P < 0.05 under both null (e) and null (f))
| quantity | Abell 2744 |
|---|---|
| catalogue rows / sources kept | 66,172 / 9,385 |
| pairs at 0.3–3″ / SED-matched | 9,328 / 174 |
| SED-matched expected, null (a) / null (c) | 94.3 / 181.3 |
| matched pairs: published / same_galaxy / visible_lens / orphan | 0 / 67 / 56 / 51 |
| orphans expected, null (e) / null (f) | 43.8 / 42.4 |
| P(≥ 51), null (e) / null (f) | 0.16 / 0.11 |

**Null under the frozen rule.** How the excess fares under each null:
- **Random-pair null (a):** the excess is large (P = 4.5e-6).
- **Redshift-clustering nulls (c) and (d):** it shrinks but persists, at P = 0.016 and 0.021.
- **Conditioned nulls (e) and (f), which the frozen rule uses:** it disappears (P 0.16 / 0.11), as D-055 found in
  the deep fields.

So the count is consistent with chance plus physical (same-redshift) companions, not with chance alone.

## Top 15 orphans (contact sheet `contact_sheet.jpg`; CATS v4.1 μ and parity)
None has a published spectroscopic redshift in the DJA table.

| # | DJA ids | sep ″ | z_a / z_b | CATS μ_a / μ_b | sheet reading | frozen escalation rule |
|---|---|---|---|---|---|---|
| 1 | 11987 / 12301 | 1.18 | 1.04 / 0.94 | 1.37 / 1.37 | faint pair beside a star spike; one member barely visible in LW | meets it (artefact-affected) |
| 2 | 10730 / 10800 | 0.57 | 0.72 / 0.84 | 1.43 / 1.43 | two faint blobs, possibly knots of one faint galaxy | meets it (ambiguous) |
| 3 | 19685 / 19722 | 0.63 | 5.79 / 5.77 | 1.52 / 1.52 | two compact clumps at z ≈ 5.8: one clumpy galaxy or a close merger | meets it (ambiguous) |
| 4 | 25350 / 25726 | 2.13 | 4.09 / 4.01 | 8.0 / 8.9 | faint pair in a high-μ region | excluded (|μ| ≥ 3) |
| 5 | 10801 / 10835 | 1.18 | 1.11 / 1.04 | 1.54 / 1.53 | part of a same-z group of ≥ 3 (shares 10835 with #6) | meets it |
| 6 | 10738 / 10835 | 0.77 | 1.04 / 1.04 | 1.52 / 1.52 | same group as #5 | meets it |
| 7 | 26610 / 27980 | 1.30 | 0.96 / 0.99 | 1.55 / 1.60 | bright compact source plus a faint companion (χ²_ν 1.66) | meets it |
| 8 | 36434 / 36589 | 1.81 | 0.57 / 0.54 | 1.25 / 1.24 | outside the cutout footprint | not judgeable here |
| 9 | 6559 / 6608 | 1.52 | 1.03 / 0.96 | 1.32 / 1.32 | faint compact pair, nothing between | meets it |
| 10 | 32288 / 32591 | 1.56 | 1.11 / 1.10 | 1.31 / 1.32 | faint pair on a bright galaxy's halo | meets it |
| 11 | 19177 / 19225 | 1.82 | 1.53 / 1.72 | 1.39 / 1.39 | faint compact pair | meets it |
| 12 | 39191 / 39387 | 2.01 | 1.98 / 1.86 | 1.48 / 1.49 | on a mosaic edge (blank region) | excluded (artifact) |
| 13 | 6200 / 6231 | 1.25 | 2.48 / 1.68 | 1.41 / 1.41 | faint pair | meets it |
| 14 | 22548 / 22662 | 1.02 | 1.63 / 1.73 | 1.33 / 1.33 | faint pair near a bright galaxy | meets it |
| 15 | 35387 / 35397 | 0.81 | 2.62 / 2.50 | −7.9 / −6.4 | faint pair in a high-|μ|, negative-parity region | excluded (|μ| ≥ 3) |

- **Escalation.** Eleven top orphans meet the frozen escalation rule (#1–3, #5–7, #9–11, #13, #14). **They are not
  candidates until vetted.** The field-level count test cannot classify individual pairs.
- **Deviation from the plan.** The per-pair `/vet-candidate` step of the frozen rule is **postponed, not done** in
  this unit.
- **First ordinary pass** on the sheet, with the D-051/D-055 deep-field vetting categories (applied after the run):

  | category | top orphans |
  |---|---|
  | group, satellite or companion | #5, #6 (same-z group sharing 10835); #7 (flux ratio 16, bright compact member); #10 and #14 (bright neighbour) |
  | possibly knots of one galaxy | #2 (0.57″), #3 (0.63″, z ≈ 5.8) |
  | artefact-affected | #1 (one member barely visible in F277W/F444W, next to a star spike) |
  | featureless faint pairs, nothing between | #9, #11, #13 |

  In the deep fields, featureless pairs of this kind were not sent to `/vet-candidate` (D-051).
- **Size of a deflector.** A deflector making the 11 separations would need σ ≈ 110–210 km/s (SIS, `model_prediction`).
  A normal galaxy of that mass would be easily visible, and none is seen. This defines the dark-lens hypothesis; it is
  not evidence for it.
- **SED match.** The escalated pairs' SED-match probabilities are p = 0.02–0.56.
- **Thresholds.** The escalation |μ| < 3 is this plan's own threshold. The `cluster_explains` column in
  `top_orphans.ecsv` uses the code's `MAX_ORDINARY_MU = 10`, so #4 shows False there.

## Limits
- One cluster, DJA v7.2-fix photometry (apertures not PSF-matched; D-051) and its eazy photo-z.
- The injections omit the cluster magnification.
- `summary.json` keeps the original absolute scratchpad path of `contact_sheet.png`; the committed copy is
  `contact_sheet.jpg`.
- There is no positive control at 0.3–3″.
- The sheet readings are visual judgements.
