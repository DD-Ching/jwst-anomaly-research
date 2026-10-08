# Orphan image pairs: a blind, lens-model-independent search (MACS0416, MACS1149, Abell 370)

Run of 2026-10-08, `scripts/orphan_pairs.py`. **Result: null.** The search finds no more SED-matched close pairs
than chance plus ordinary redshift clustering would give. The pairs it ranks highest all have ordinary readings, or
are faint pairs that the null fully accounts for. Nothing here is evidence for a dark deflector.

**Hypothesis tested.** An unseen compact deflector (a "dark lens") would split one background source into two
images 0.3–3″ apart. The images would have identical SEDs (lensing is achromatic) and no visible galaxy between
them. The search looks for such *orphan pairs*: pairs that no published multiple-image system explains. This is a
hypothesis. An SED match is not evidence of lensing (charter).

## Data (observed) and models (model_prediction)

- CANUCS DR1 PSF-matched photometry and EAzY photo-z catalogues (`FLUX_COLOR03_TOTAL_<band>`, 0.3″ colour
  apertures; HST + NIRCam, 20–22 bands). Readme:
  <https://archive.stsci.edu/hlsps/canucs/dr1/webpage/hlsp_canucs_jwst-hst_multi_v1_photometry-cat_readme.txt>.
- Published multiple images that the search excludes:
  - HFF CATS `arcs.txt` for each cluster, shifted to the JWST frame (D-034/D-038);
  - CANUCS `allmultim-cat` for MACS0416;
  - CANUCS Lenstool `multim` for Abell 370.
  - CANUCS DR1 has no MACS1149 image list; its readme says the list comes with v2.
- Lens-model check: the CATS v4/v4.1 deflection maps (`lc.load_model('<c>-cats')`). Outside those maps the check
  falls back to the CANUCS catalogue `MU`.
- Cutouts: MAST pipeline `_i2d` F150W/F277W/F444W of program 1208, read from S3 by byte range (no full download).
  Visits: macs0416 o004_t002, macs1149 o008_t004, abell370 o002_t001.

## Method (every threshold is an ASSUMPTION)

1. **Sources.** Each source must pass all of these:
   - `USE_PHOT_APER03` is set and `FLAG_BCG` is not;
   - summed F277W+F356W+F444W colour-aperture S/N ≥ 10;
   - S/N ≥ 10 in at least 8 bands;
   - behind the cluster: `Z025` > z_cluster + 0.1, or `Z_SPEC` > z_cluster + 0.1 where a spectroscopic redshift
     exists.

   The 8-band rule matters. Without it, 25 % of random far-apart faint pairs pass the SED test in MACS0416, so
   the test has no discriminating power.
2. **Pairs and the SED match.** Take all pairs 0.3–3.0″ apart. Fit SED A against a free multiple of SED B by
   iterated weighted least squares, using every band valid in both. Errors get a 3 % floor in quadrature. A pair
   matches when:
   - at least 8 bands are valid in both;
   - the χ² survival probability is ≥ 0.01, so 99 % of truly identical SEDs pass;
   - the 16–84 % photo-z intervals overlap.
3. **Classes.** Each pair gets the first class that applies:
   - **published**: a member lies within 1″ of a published image;
   - **same_galaxy**: either
     - the separation is less than 0.5 × the sum of the two Kron aperture radii (2.5 × `KRON_RADIUS` × `A`, or the
       isophotal radius where `KRON_RADIUS` = 0, which is about half the catalogue), or
     - both members are `FLAG_DEBLEND` and their segment boxes overlap;
   - **visible_lens**: another catalogued source of any redshift, BCGs included, lies at least 0.3″ from both
     members and also:
     - within 0.3″ of the midpoint, or
     - within 0.3″ of the joining segment, or
     - inside the circle that has the pair as its diameter.

     The last two conditions widen the brief's midpoint-only rule. A lens need not sit at the midpoint when the two
     fluxes differ, and the first contact sheet showed galaxies sitting between "orphans".
   - **orphan**: none of the above.
4. **Null.** Three estimates of the chance SED-match rate:
   - (a) the match fraction of cross pairs between the source list and copies shifted by 10–60″ (20 shifts);
   - (b) the match fraction of real pairs 10–30″ apart;
   - (c) the match fraction of 10–30″ pairs whose photo-z intervals already overlap. This null keeps redshift
     clustering: physical neighbours at one redshift share SEDs.

   Expected counts are each fraction times the number of close pairs. For a class, they are the fraction times
   the close pairs that the geometric class rules put in that class. Those rules ignore the SED, so they can be
   applied to every pair.
5. **Follow-up of the top ≤ 15 orphans per field.** Rank by the number of bands that agree within 2σ, then by
   reduced χ². For each pair:
   - make a contact sheet;
   - look up the CATS magnification (signed) at both members at the pair's mean photo-z. The pair counts as an
     ordinary cluster-scale image pair if |μ| > 10 or the parity flips between the members (a critical curve
     between them).
   - measure mirror symmetry: |Δφ_a + Δφ_b|, where Δφ is a member's major-axis angle relative to the joining line.
     A fold pair or a tangential pair gives a value near 0°.
   - compute the dark-deflector numbers (HYPOTHESIS): θ_E = sep/2, the mass inside θ_E, and the SIS velocity
     dispersion σ, all at z_l = z_cluster.

Reproduce with `python scripts/orphan_pairs.py --field {macs0416,macs1149,abell370} --cutouts`. Outputs go to
`outputs/orphan_pairs/<field>/`: `summary.json`, `matched_pairs.ecsv`, `top_orphans.ecsv` and
`contact_sheet.png`. Each field takes about 5 s without cutouts and about 30 s with them. The catalogues and image lists are
downloaded on first use and verified against the sha256 values in `FIELDS` (SOURCES.md); the CANUCS Abell 370 list
is shifted to the JWST frame by the `abell370-canucs` offset (D-044). `summary.json` writes NaN as `null`.

## Null comparison per field (derived)

| | MACS0416 | MACS1149 | Abell 370 |
|---|---|---|---|
| z_cluster | 0.396 | 0.543 | 0.375 |
| Catalogue rows / sources kept | 14149 / 1942 | 12851 / 1711 | 13567 / 1391 |
| Pairs at 0.3–3″ | 1664 | 1418 | 1028 |
| SED-matched (observed) | 62 | 98 | 35 |
| Expected, null (a) shift / (b) 10–30″ | 42.3 / 44.6 | 44.8 / 42.6 | 27.8 / 28.3 |
| Random-pair match rate (a) / (b) | 2.5 % / 2.7 % | 3.2 % / 3.0 % | 2.7 % / 2.7 % |
| Expected, null (c) z-clustered | 69.7 | 83.3 | 48.0 |
| Matched: published / same_galaxy / visible_lens / orphan | 12 / 11 / 28 / 11 | 8 / 19 / 53 / 18 | 0 / 11 / 15 / 9 |
| Orphans expected, (a) / (c) | 11.0 / 13.5 | 10.5 / 12.9 | 6.5 / 8.1 |
| P(≥ observed orphans), (a) / (c) | 0.54 / 0.79 | 0.021 / 0.10 | 0.20 / 0.43 |

Reading the table:
- **Null (a)/(b).** The SED-matched excess over (a)/(b) is significant in MACS1149 (98 against 45) and MACS0416
  (62 against 42). Null (c) removes it: the same match rate among photo-z-overlapping far pairs predicts 83 and 70.
  The excess is the ordinary clustering of galaxies at one redshift, such as the z ≈ 1 group in MACS1149 (pairs
  5113705/5113706/5113747), not lensing.
- **Orphans.** Orphan counts agree with the z-clustered null in every field. Summed over the three fields, 38 orphans
  are observed against 34.5 expected.
- **Match fraction against separation.** It is higher at 0.3–1.0″ (4–10 %) than at 1–3″ (2.4–7.7 %). That fits
  blending: the 0.3″ apertures of very close neighbours share light. It also fits clumps of one galaxy. Neither
  null models blending, so pairs under 1″ are least reliable.

## Top orphans and verdicts (contact sheets inspected)

These are visual verdicts on the contact sheets, plus the numbers in `top_orphans.ecsv`. "Chance" means a faint
pair with no feature beyond an SED match. The field-level null predicts this many such pairs. None of the pairs has
CATS |μ| > 10 or a parity flip: every pair lies where the models give |μ| ≈ 1.0–2.5, so a cluster-scale image pair
is not the explanation. Mirror angles are spread over 0–87°, as random orientation would give. The dark-lens column
lists θ_E, the mass inside θ_E and the SIS σ at z_l = z_cluster (HYPOTHESIS).

MACS0416 (contact sheet: `outputs/orphan_pairs/macs0416/contact_sheet.png`):

| # | Pair (CANUCS ids) | sep ″ | χ²_ν / bands | z_a / z_b | Verdict |
|---|---|---|---|---|---|
| 1 | 3116708 / 3116713 | 1.39 | 1.02 / 19 | 2.24 / 2.11 | chance or physical neighbours; one member barely visible, other elongated |
| 2 | 3112403 / 3112453 | 2.03 | 1.08 / 15 | 1.94 / 1.84 | faint compact pair, nothing between; consistent with chance |
| 3 | 3104853 / 3104963 | 2.58 | 1.14 / 13 | 1.41 / 1.11 | two resolved galaxies, different orientations; chance or neighbours |
| 4 | 3116012 / 3116013 | 0.63 | 1.01 / 10 | 2.80 / 2.80 | knots or interacting pair beside a larger elongated galaxy |
| 5 | 3100346 / 3106007 | 1.90 | 1.07 / 10 | 1.04 / 0.98 | big bright galaxy plus a compact companion (flux ratio 51); satellite |
| 6 | 3103892 / 3116048 | 2.12 | 1.36 / 10 | 4.07 / 4.19 | faint z≈4 pair, nothing between (a bad-pixel patch between them in F277W only); consistent with chance |
| 7 | 3105588 / 3105614 | 0.45 | 1.42 / 10 | 5.38 / 5.39 (z_spec 5.367) | clumps of one z≈5.4 galaxy or a close merger |
| 8 | 3100243 / 3119128 | 0.99 | 1.93 / 9 | 1.52 / 1.41 | knot or satellite of an extended galaxy with a tail |
| 9 | 3107134 / 3107170 | 1.28 | 2.60 / 8 | 1.08 / 1.04 | very faint; chance |
| 10 | 3102515 / 3102535 | 0.98 | 1.81 / 8 | 2.25 / 1.93 | faint; chance |
| 11 | 3100605 / 3107837 | 1.78 | 2.02 / 8 | 0.66 / 0.86 | bright compact galaxy plus a faint one (ratio 44); chance or satellite |

MACS1149 (contact sheet: `outputs/orphan_pairs/macs1149/contact_sheet.png`):

| # | Pair | sep ″ | χ²_ν / bands | z_a / z_b | Verdict |
|---|---|---|---|---|---|
| 1 | 5108605 / 5115231 | 0.84 | 0.70 / 20 | 1.31 / 1.55 | faint source plus an elongated one, morphologies differ; chance |
| 2 | 5106502 / 5114210 | 2.28 | 1.30 / 20 | 2.27 / 2.29 (z_spec 2.278) | compact source in a small group plus a faint isolated one (ratio 0.11); neighbours |
| 3 | 5115775 / 5127584 | 1.27 | 1.78 / 20 | 1.52 / 1.70 | knots along one diagonal streak; same galaxy, which the Kron rule missed |
| 4 | 5113706 / 5113747 | 1.80 | 1.66 / 20 | 0.86 / 1.03 | z≈1 group members |
| 5 | 5113705 / 5113747 | 1.75 | 1.90 / 20 | 0.96 / 1.03 | z≈1 group members (same group as #4) |
| 6 | 5103988 / 5104005 | 0.72 | 0.73 / 14 | 1.09 / 1.03 | two compact members of a compact group |
| 7 | 5103319 / 5103342 | 1.58 | 0.83 / 13 | 1.47 / 1.42 | extended galaxy plus a compact companion (ratio 8); satellite |
| 8 | 5101771 / 5101813 | 1.68 | 1.33 / 13 | 1.10 / 0.97 | two compact sources, nothing between; consistent with chance |
| 9 | 5104186 / 5104200 | 2.17 | 1.93 / 12 | 0.99 / 0.84 | faint pair next to a bright compact source; chance |
| 10 | 5112716 / 5112717 | 0.93 | 1.10 / 11 | 1.15 / 1.18 | two faint companions about 1″ beside a bright galaxy (not between them); satellites |
| 11 | 5101623 / 5101642 | 0.93 | 1.59 / 11 | 2.81 / 2.72 | compact plus a fainter source; merger or chance |
| 12 | 5101304 / 5101328 | 1.69 | 0.70 / 10 | 1.07 / 1.05 | faint; chance |
| 13 | 5112345 / 5112356 | 1.27 | 1.20 / 10 | 2.29 / 1.68 | one member on the edge of a bright irregular galaxy; knot or companion |
| 14 | 5101328 / 5101407 | 2.64 | 0.92 / 9 | 1.05 / 1.01 | faint compact plus an elongated galaxy; chance |
| 15 | 5101391 / 5112572 | 1.37 | 1.75 / 9 | 3.44 / 3.74 | one member is part of a chain galaxy; chance |

Abell 370 (contact sheet: `outputs/orphan_pairs/abell370/contact_sheet.png`):

| # | Pair | sep ″ | χ²_ν / bands | z_a / z_b | Verdict |
|---|---|---|---|---|---|
| 1 | 2111474 / 2111510 | 0.90 | 0.79 / 21 | 1.29 / 1.32 | faint plus compact (ratio 0.36); consistent with chance |
| 2 | 2100971 / 2100981 | 1.22 | 1.06 / 21 | 0.86 / 0.73 (z_spec 0.733) | compact pair in a z≈0.75 group |
| 3 | 2111456 / 2111616 | 2.94 | 1.27 / 21 | 1.01 / 1.14 | both faint; chance |
| 4 | 2101589 / 2124148 | 2.10 | 1.73 / 14 | 1.07 / 1.15 | bright knot of a large curved galaxy plus an isolated compact source (ratio 22) |
| 5 | 2104046 / 2104071 | 0.72 | 1.18 / 10 | 1.92 / 2.30 | two compact sources, nothing between; consistent with chance or a close pair |
| 6 | 2118938 / 2124544 | 1.15 | 2.32 / 9 | 1.51 / 1.29 | beside a bright elongated galaxy (about 0.7″ off the line); group or knots |
| 7 | 2100080 / 2102849 | 2.94 | 1.22 / 8 | 1.89 / 1.77 | faint; chance |
| 8 | 2124616 / 2124617 | 0.88 | 1.61 / 8 | 1.56 / 1.52 | both inside one large irregular galaxy; knots, which the Kron rule missed |
| 9 | 2102589 / 2102603 | 2.00 | 1.72 / 6 | 1.71 / 1.75 | faint; chance |

### What survives every ordinary test

Taken one at a time, no pair survives as an anomaly. Eleven pairs are faint pairs with no visible galaxy between them
and no other distinguishing feature:
- MACS0416 #2, #6, #9, #10;
- MACS1149 #8, #12;
- Abell 370 #1, #3, #5, #7, #9.

They are "unexplained" only in that the data cannot tell a chance SED match from a lensed pair. The field-level null
predicts this many (orphans: 11 against 13.5, 18 against 12.9 and 9 against 8.1 under null (c)), so they are what
chance gives. Their dark-deflector requirements are large:
- θ_E = 0.36–1.47″, so M(<θ_E) ≈ 2×10¹⁰–5×10¹¹ M☉ and σ_SIS ≈ 130–300 km/s at z_l = z_cluster (HYPOTHESIS).
- Under a Faber–Jackson anchor of σ* = 180 km/s at m*(F160W) = 19.0 for z ≈ 0.4 cluster members (ASSUMPTION, not
  fitted here), a luminous deflector of that σ would be m ≈ 16.8–20.4.
- That is about 9–12 mag brighter than a ~29 AB detection depth (ASSUMPTION: typical CANUCS depth, not measured here).

So a deflector of that mass would have to be essentially dark, with M/L orders of magnitude above any galaxy. Before
that becomes worth considering, the far likelier reading is a chance SED match. No pair is a candidate for
`/vet-candidate`.

## Limits

- **Detection power.** The SED test only works at high S/N. Requiring S/N ≥ 10 in 8 bands removes most faint
  sources (with the foreground and flag cuts, 86–90 % of each catalogue), so faint, highly magnified image pairs fall outside the search.
- **Photo-z.** A 16–84 % photo-z overlap is permissive for sources with broad posteriors.
- **Blending.** The 0.3″ colour apertures of pairs under about 0.6″ share light. This pushes their SEDs toward
  agreement, and the nulls do not model it.
- **same_galaxy rule.** The Kron/segment rule misses knots of large irregular galaxies (MACS1149 #3,
  Abell 370 #8). Visual inspection caught them. A segmentation-map adjacency test would be stricter. The CANUCS
  segmentation maps are not used here.
- **visible_lens rule.** Only catalogued sources count. A lens hidden in a BCG halo or the ICL would be missed. The
  wider "between" circle removes real neighbours along with possible lenses.
- **Coverage.** Every orphan lies at |μ| ≈ 1–2.5. Cluster cores, where `USE_PHOT_APER03` and the BCG cuts remove
  sources, are not searched.
- **Cutout depth.** The MAST `_i2d` cutouts come from single visits and are shallower than the CANUCS mosaics.
  Members at the S/N floor are barely visible in them.
- **Not tested.** A per-pair time-delay or flux-ratio test, spectroscopy, and pairs whose two members fall in
  different photo-z solutions (a catastrophic photo-z on one image).

## Deep fields: CANUCS flanking fields and GOODS-North (D-051)

Run of 2026-10-08, same script. **Result: null.** Six non-cluster fields (108 arcmin² searched) give 355 orphans
against 315 expected from chance SED matches at one redshift (null (e) below). No single pair survives the
ordinary explanations, and no pair goes to `/vet-candidate`. The injection-recovery that turns this into a limit on
dark deflectors is in [exotic_limits.md](exotic_limits.md), section "W2 / dark-deflector pairs".

### What changed in the script

- **One column layout.** `as_standard` turns each catalogue into the same columns (`STANDARD_COLUMNS`). CANUCS
  (`canucs_standard`) and DJA grizli (`dja_standard`, catalogue plus eazy-py zout) both run through
  `search` unchanged. The three D-048 cluster runs reproduce every count, every pair and every column of
  `matched_pairs.ecsv` / `top_orphans.ecsv`. The only change is the column rename `mu_canucs` → `mu_cat`.
- **Deep fields** (`DEEP_FIELDS`). They have no published images and no cluster model. Two rules replace the
  cluster ones (ASSUMPTIONs):
  - the redshift cut is z_low > `Z_LENS_REF` + 0.1 = 0.5, as behind a z = 0.4 cluster;
  - the ordinary-lensing test becomes "no lens model needed": |μ| is the catalogue `MU` (CANUCS model, 1.00–1.4
    in the flanking fields) or 1 (DJA).
- **S/N bands.** These are the bands of F277W/F356W/F444W that the catalogue has. The MACS0417 and MACS1423
  flanking fields have no F356W.
- **DJA specifics** (ASSUMPTIONs):
  - the 0.36″ apertures (`aper_0`, closest to CANUCS's 0.3″) are converted from µJy to nJy;
  - SEP aperture flags in `photometry.BAD_FLAGS` invalidate a band, and MIRI and `<band>u` duplicates are left out;
  - `flag & 1` (SEP `OBJ_MERGED`) is "deblended";
  - the same_galaxy radius is 3.3 × `flux_radius`. The DJA Kron apertures follow another convention:
    `kron_radius` lies between 2.4 and 3.8, which gives 8.7 × the half-light radius against CANUCS's 3.2–3.5.
    With them, the same_galaxy rule swallowed every injected pair under 1.5″. The first GOODS-N run (kept below
    as a failed approach) found 69 orphans this way, against 109 with the fix.
  - The DJA photometry is **not PSF-matched**: the SW bands are measured on 0.02″ images. Both members of a
    lensed pair of compact images share the aperture losses, so the free-scale χ² is unaffected. Resolved pairs of
    different sizes are not. The null (a) match rate is 1.5 %, against 2.3–3.7 % in CANUCS.
- **Two more nulls**, needed because the deep-field orphan counts exceed null (c) (P = 0.04–0.08 in three of the
  flanking fields):
  - (d) the z-overlap match rate of real pairs at 3–6″;
  - (e) null (c) conditioned on what makes two SEDs easy to match. The cells are the fainter member's summed
    S/N, the larger member's aperture radius and the LW/SW colour (`pair_cells`). Each close pair gets its cell's
    far-pair match rate.

  `zoverlap_match_fraction_by_sep` shows the z-overlap match rate is flat from 0.3″ to 10″ in every field.
- **Other changes.**
  - `searched_footprint` measures the searched area: 1″ grid points with a catalogue row within 4″ that has valid
    S/N bands and ≥ 8 valid bands.
  - Cutouts of deep fields use the MAST level-3 `_i2d` whose footprint (`s_region`) contains the pair: one
    query per band. They retry intermittent S3 errors.
  - `fetch_tar_member` streams the DJA photo-z tarball and keeps only the zout. Both checksums are verified, and
    the 371 MB archive never touches disk.

Reproduce:

```
python scripts/orphan_pairs.py --field {macs0416,macs1149,abell370,macs0417,macs1423}-ncf --cutouts
python scripts/orphan_pairs.py --field goodsn-dja --cutouts
```

Wall time per field (4-core cloud VM):
- flanking fields: 3–6 s without cutouts, 45–60 s with them;
- GOODS-N: 20 s without cutouts and about 3 min with them, plus a one-time 224 MB catalogue download and 371 MB
  tarball stream.

### Null comparison (derived)

| | M0416-NCF | M1149-NCF | A370-NCF | M0417-NCF | M1423-NCF | GOODS-N (DJA) |
|---|---|---|---|---|---|---|
| Catalogue rows / sources kept | 10804 / 2362 | 11657 / 2716 | 10267 / 2676 | 11027 / 2344 | 10412 / 2183 | 70421 / 13741 |
| Bands / searched area (arcmin²) | 29 / 9.92 | 27 / 9.87 | 29 / 9.90 | 16 / 9.86 | 18 / 9.82 | 21 / 58.72 |
| Pairs at 0.3–3″ / SED-matched | 3303 / 148 | 3920 / 143 | 3746 / 135 | 2985 / 152 | 2626 / 151 | 14789 / 366 |
| Random-pair match rate (a) / (b) | 2.7 / 2.7 % | 2.3 / 2.4 % | 2.3 / 2.3 % | 3.3 / 3.4 % | 3.7 / 3.7 % | 1.5 / 1.6 % |
| SED-matched expected (a) / (c) | 89.9 / 162.8 | 91.1 / 172.7 | 86.6 / 181.6 | 97.9 / 169.1 | 97.2 / 166.5 | 223.6 / 395.6 |
| Matched: same_galaxy / visible_lens / orphan | 28 / 69 / 51 | 27 / 66 / 50 | 22 / 70 / 43 | 32 / 71 / 49 | 18 / 80 / 53 | 94 / 163 / 109 |
| Orphans expected (c) / (d) / (e) | 39.2 / 38.2 / 41.4 | 40.1 / 40.2 / 40.6 | 43.0 / 38.9 / 43.1 | 37.3 / 37.5 / 38.0 | 43.3 / 45.9 / 47.9 | 99.4 / 99.3 / 104.4 |
| P(≥ observed orphans) (c) / (e) | 0.040 / 0.083 | 0.074 / 0.085 | 0.52 / 0.53 | 0.038 / 0.049 | 0.083 / 0.25 | 0.18 / 0.34 |

Reading the table:
- **SED-matched pairs.** As in the clusters, they exceed the random-pair nulls (a)/(b) and fall below the
  same-redshift null (c).
- **Orphans in GOODS-N** match every null.
- **Superseded by D-054** (section below): with a symmetric colour cell the excess is gone (246 / 229.5,
  P = 0.15). Original D-051 reading kept for the record:
- **Orphans in the five flanking fields** are 246 against 211 under null (e), P = 0.010. With GOODS-N, the total
  is 355 against 315, P = 0.015. This is a 10–15 % excess at about 2.3σ.
  - It is spread over 1–3″ and over fainter-member S/N 20–40. No single bin carries it (a one-off diagnostic
    split; per-field bins lie within about +2σ).
  - Pairs under 1″ show no excess. Lensing by the smallest deflectors would put pairs there, but the
    same_galaxy rule removes most of them anyway (see the injections).
- **Lensing is an implausible reading of the excess.** The per-deflector recovery efficiency measured by injection
  is about 0.5 % at θ_E = 0.3–0.7″. The ~35 excess pairs would then need ~7000 dark deflectors in 50 arcmin²:
  - that is 140 arcmin⁻², or 5 × 10⁵ deg⁻², each with M ≈ 2 × 10¹⁰–10¹¹ M☉ inside θ_E (model_prediction);
  - that is one invisible galaxy-mass deflector for every ~7 catalogued galaxies (about 1050 rows per arcmin² in
    the flanking fields).

  An ordinary explanation is far likelier.
- **The likelier reading** is that null (e) still under-models physical companions. Pairs at 1–3″ at one redshift
  (satellites, interacting pairs) share stellar populations more closely than z-overlapping pairs 10–30″ apart.
  The contact sheets show many such companions.

### Null (e) fixed and a companion-aware null (f) (D-054)

Re-run of 2026-10-08 (cloud). D-051's null (e) took the colour cell from member `i` of the pair only (pair order
is arbitrary) and put a NaN colour in the 0–0.3 bin. Now the cell holds both members' colour bins (unordered),
and NaN has its own bin. New null (f) is (e) with the 3–6″ z-overlapping pairs of null (d) as reference: physical
companions are common there and galaxy-scale lensing is not. Orphan counts are unchanged (bit for bit).

| Field | Orphans | (e) D-051 | (e) fixed, P | (f), P | (f) cells without reference |
|---|---|---|---|---|---|
| M0416-NCF | 51 | 41.4 | 44.4, 0.18 | 44.1, 0.17 | 14 |
| M1149-NCF | 50 | 40.6 | 45.4, 0.27 | 49.0, 0.46 | 9 |
| A370-NCF | 43 | 43.1 | 45.0, 0.64 | 43.9, 0.57 | 15 |
| M0417-NCF | 49 | 38.0 | 43.1, 0.20 | 45.7, 0.33 | 29 |
| M1423-NCF | 53 | 47.9 | 51.6, 0.44 | 54.5, 0.60 | 10 |
| GOODS-N | 109 | 104.4 | 103.9, 0.32 | 106.5, 0.42 | 4 |
| Five NCFs | 246 | 211 (P 0.010) | 229.5, 0.15 | 237.2, 0.29 | |
| Total | 355 | 315.4 (P 0.015) | 333.4, 0.12 | 343.7, 0.28 | |

- **Reading.** The flanking-field excess was a null-model artefact: no orphan excess remains under either null.
  The two changes (symmetry, NaN bin) were made together; their separate shares were not measured.
- "(f) cells without reference": close z-overlapping pairs whose cell has no 3–6″ reference pair; they take the
  global 3–6″ rate.
- **D-048 cluster fields**, re-run with the same code: orphans 11 / 18 / 9 (unchanged) against (e) 14.3 / 14.1 /
  11.1 (P = 0.84 / 0.18 / 0.77) and (f) 15.4 / 12.5 / 10.8 (P = 0.90 / 0.085 / 0.75) in MACS0416 / MACS1149 /
  Abell 370. Null.
- The unconditioned 3–6″ null (d) gives 300.0 (P = 0.001): conditioning on S/N, size and colour matters more than
  the reference annulus.

### Top orphans and verdicts (all six contact sheets inspected)

The top 15 orphans per field come from `outputs/orphan_pairs/<field>/contact_sheet.png` and `top_orphans.ecsv`.
They are vetted cheapest-first:
1. knots of one galaxy or blends;
2. group members, satellites and companions;
3. image artefacts;
4. chance SED matches.

Every pair has catalogue |μ| ≤ 1.3 and no cluster model to invoke. The cutouts are single-visit MAST `_i2d`
images, so the faintest members are barely visible.

| Field | Knot / same galaxy | Group, satellite or companion (flux ratio > 5 or a bright neighbour) | Artefact-affected | Faint pair, nothing between, no other feature ("chance") |
|---|---|---|---|---|
| M0416-NCF | #11, #12 (two knots in the disk of a z = 0.56 spiral) | #5, #6, #8 (ratio 11), #14, #15 (ratio 54) | #10 (detector stripe in F150W) | #1, #2, #3, #4, #7, #9, #13 |
| M1149-NCF | #11 (0.6″) | #2, #5, #6, #9, #13 (ratio 14), #14, #15 | — | #1, #3, #4, #7, #8, #10, #12 |
| A370-NCF | #1, #4 (≤ 0.6″ in one blob), #9/#10 (chain sharing 2222566) | #5, #6, #8, #12 (ratio 16), #14 (ratio 13) | — | #2 (z 4.74/4.75), #3/#7 (share an edge-on disk with different partners), #11, #13, #15 |
| M0417-NCF | #3, #10, #14 (both inside a face-on spiral) | #6, #8, #9, #11, #15 (ratio 13) | — | #1, #2, #4, #5, #7, #12, #13 |
| M1423-NCF | #4, #9/#10 (knots of one irregular) | #2, #6, #11, #12, #14 | #7 (on a scattered-light stripe; the CANUCS readme lists "dragon's breath" in this field), #8, #15 (beside a bright star) | #1, #3, #5, #13 |
| GOODS-N | — | #2, #15 (beside a spiral), #6 (ratio 28), #7, #9, #14 | — | #1, #3, #4, #5, #8, #10, #11, #12, #13 |

(GOODS-N ranks refer to the final run with the half-light radius rule.)

### What survives

**No pair survives as an anomaly.** Forty of the 90 inspected pairs are faint pairs with nothing visible between
them. The data cannot tell them apart from a lensed pair, and they are what the nulls predict. Their dark-deflector
requirements (`hypothesis` numbers in `top_orphans.ecsv`, z_l = 0.4) are:
- θ_E = 0.4–1.4″;
- M(<θ_E) ≈ 4 × 10¹⁰–7 × 10¹¹ M☉;
- σ_SIS ≈ 145–350 km/s.

As in the clusters, a luminous galaxy of that mass would be many magnitudes above the detection limit.

### Limits (deep fields)

- **Same rules as D-048.** The same S/N, Kron-rule and visible-lens limitations apply. In a deep field the
  visible-lens rule removes about 73 % of injected pairs at θ_E = 0.7″, because a catalogued source falls inside
  the pair's circle by chance at these densities.
- **The deep-field excess** (resolved by D-054: a null-model artefact). Null (e) does not model physical companions at 1–3″, so the deep-field orphan excess
  stays unexplained at the 2σ level. A companion-aware null would be one built from spectroscopic pairs or from
  pairs matched in redshift *and* environment.
- **DJA photometry** is not PSF-matched, and its same_galaxy radius is calibrated on CANUCS (3.3 × half-light
  radius; ASSUMPTION).
- **Failed approach (recorded so it is not repeated).** The DJA Kron apertures (`2.5 × kron_radius × a_image`)
  are ~3× the CANUCS ones. With them, the same_galaxy rule classed every injected GOODS-N pair under 1.5″ as one
  galaxy: 0 of 2000 point lenses at θ_E = 0.3″ were recovered.
