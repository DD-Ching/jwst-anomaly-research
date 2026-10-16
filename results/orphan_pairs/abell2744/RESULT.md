# Orphan-pair (dark-deflector) screen on Abell 2744: result (2026-10-16)

- **Rules:** `PLAN.md`, committed before the run. Every number below is derived from observed catalogues; injected
  images are simulated.
- **Exploratory:** this is the first run of this statistic on this field. Abell 2744 had already been studied with
  other tests, so it is not a blind sample.

## Checks
1. **Alignment.** 163 of the 258 published images (CATS v4.1 and Bergamini+23) have a DJA source within 1″. Their
   median offset is (−0.014″, +0.020″) and the median separation 0.12″, inside the ≤ 0.3″ limit. **Passed.**
2. **Positive control.** None of the 174 SED-matched close pairs lies within 1″ of a published image (class
   `published` = 0). Published multiple images here are mostly > 3″ apart. **This field has no positive control at
   0.3–3″.**
3. **Injection recovery** (`inject_pairs.py`, default counts, 400 lensed sources per lens type and θ_E). The cluster
   lens's magnification of the injected images is not modelled.
   - Few lensed sources reach selection. Both images pass the S/N and background cuts for only 4–50 of 400 sources.
   - Of the pairs that are selected and SED-matched, few are classed `orphan`:

     | lens type | θ_E (″) | recovered / selected and matched |
     |---|---|---|
     | point | 0.3 | 17 / 45 |
     | point | 0.7 | 3 / 27 |
     | W2 | 0.3 | 12 / 43 |
     | W2 | 0.7 | 3 / 38 |
     | W1 | 1.5 | 18 / 38 |
     | any | 1.5 (point, W2) | no SED-matched pairs |

   - The rest are classed `visible_lens` or `same_galaxy`, because the cluster field is crowded: point θ_E 0.7 gives
     17 `visible_lens` and 7 `same_galaxy`.
   - **So this screen has low completeness in Abell 2744.** A null here constrains only a small part of the dark-lens
     population. The deep-field efficiencies (D-051, D-055) do not apply to this field.
4. **Artifact burden in the top 15.**
   - #12 sits on a mosaic edge with a blank region.
   - #1 lies next to a bright star's spike or halo, though both members look like real faint sources.
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

**Null.** The orphan count is consistent with chance plus redshift clustering. As in the other clusters, the excess
over the random-pair null (a) (P = 4.5e-6) disappears under the redshift-clustering nulls (c)–(f).

## Top 15 orphans (contact sheet `contact_sheet.jpg`; CATS v4.1 μ and parity)
None has a published spectroscopic redshift in the DJA table.

| # | DJA ids | sep ″ | z_a / z_b | CATS μ_a / μ_b | sheet reading | frozen escalation rule |
|---|---|---|---|---|---|---|
| 1 | 11987 / 12301 | 1.18 | 1.04 / 0.94 | 1.37 / 1.37 | two faint compact sources beside a star spike, nothing between | meets it |
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

- **Eleven top orphans meet the frozen escalation rule** (#1–3, #5–7, #9–11, #13, #14). None is vetted in this unit.
- **Why they are not counted as candidates.** The field's orphan count is what chance and clustering predict (51
  observed against about 43 expected). Their SED-match probabilities (p 0.02–0.86) and appearance on the sheet do not
  set them apart from chance pairs.
- **Size of a deflector.** A deflector making their separations would need σ ≈ 110–280 km/s (SIS, `model_prediction`).
  A normal galaxy of that mass would be easily visible, and none is seen. This is the dark-lens hypothesis, not
  evidence for it.
- **Still to do:** a per-pair ordinary check beyond the sheet (`/vet-candidate`: segmentation, morphology, spectra).

## Limits
- One cluster, DJA v7.2-fix photometry (apertures not PSF-matched; D-051) and its eazy photo-z.
- The injections omit the cluster magnification.
- There is no positive control at 0.3–3″.
- The sheet readings are visual judgements.
