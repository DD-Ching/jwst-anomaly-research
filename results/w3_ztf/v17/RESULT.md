# W3 ZTF v17 — evidence check of the two v16 final survivors: both undetermined; 1624 has strong data-problem support, 1625 a likely catalogue blend

- **Post hoc:** this check was done after the v16 images had been seen.
  - `PLAN.md` was committed before any recomputation (baff0fa).
  - Two steps depart from it and are marked **deviation** below.
- **Labels:** IBE `_sciimg`, `_mskimg`, `_psfcat` and the DR24 light-curve API are `observed`; ratios and classes are
  `derived`.
- **Downloads (measured bytes, `check.json`):** 10.34 MB in total, made up of 40 cutouts (20 science, 20 mask) and
  20 psfcat files.
- **Scope:** no search, no threshold change, no new veto rule. No limit is derived, and nothing here bears on
  wormholes.

## 1. Evidence (`evidence.json`, copied from committed v16 files)

| | 1624: 708205400009336 | 1625: 709212400017577 |
|---|---|---|
| members | 5 objectids (`evidence.json`) | 5 objectids |
| position | 127.0853788, +39.3086507 | 128.6204723, +40.9935282 |
| run | r, HMJD 59815.509–59817.508, 9 NDs | g, HMJD 59896.382–59900.472, 3 NDs |
| final rule | frozen v13 W3, `--edge time`; all vetting tests true | same |
| achromatic test | no other-band epochs (`other_filter_epochs` 0): no information | same |
| v16 automated classes of the checked NDs | 4 × class 2 (poor image), 4 × class 5 (undetermined); 1 ND unchecked | 3 × class 2 (poor image) |
| v16 visual reading (separate) | star seen beside or partly on a black column in all 8 NDs | a source seen at the position in all 3 wide-seeing NDs |
| v16 2″ bracket detections | 2,446–2,560 µJy (4 values, median 2,496) | 131–156 µJy (4 values) |
| v16 2″ ND fluxes | −1,322 to 1,696 µJy | 101–147 µJy |
| official mask at the source (ND) | bit 2 LOW RESPONSIVITY in all 8 | bit 1 (detection) only |
| ND maglimit / seeing | 19.08–19.94 / 1.92–2.19″ | 19.06–19.92 / 4.05–4.64″ |

- The v16 "≈ 2,450 µJy" figure is the low end of the four 2″ bracket detections. Their median is 2,496 µJy.
- The DR24 reference magnitude is r 15.18. That is a PSF total, and it is not compared with 2″ aperture fluxes.

## 2. Planned checks (`check.json`, `check.png`)

1. **Consistency.**
   - All 20 science files are `_sciimg` (header IMGTYPE "object"), not difference images.
   - Flux is in µJy from each header MAGZP (25.92–26.07 for 1624).
   - The position lies inside every image.
   - The centroid within 3″ lies 0.24–0.81″ from the DR24 position (1624) in every exposure.
   - The DR24 records (LC API, by expid) agree with v16's ND / detection labels for every checked exposure. The 9th
     1624 ND (206350236) was found and checked.
2. **Aperture and mask** (radius = max(2″, seeing)).
   - **1624:** bad-bit pixels in the target aperture are 2–8 of 11–15 in all 9 NDs and 0 in all 4 detections.
   - **1624:** the two LOW RESPONSIVITY columns cover every row of the 178-row cutout in every exposure. Away from the
     star they lie 400–3,080 DN below their neighbours in the NDs and 434–1,111 DN below in the detections.
   - **1625:** no bad bits.
3. **Comparison stars (planned, image-based).** The 180″ cutouts held 0 (1624) and 1 (1625) stars in the planned
   brightness range, below the required 3. **The planned rule therefore gives "undetermined" for both** (every ND
   unknown).
4. **psfcat** (same images and photometry that feed DR24, so this is **not** an independent observation).
   - The target is **missing** (no source within 1.5″) in all 9 NDs of 1624 and all 3 NDs of 1625.
   - It is matched with flags 0 in every detection: PSF magnitude spread −0.014 to +0.005 mag for 1624 and
     −0.054 to +0.053 mag for 1625.
   - A missing row is not taken as zero flux.

## 3. Deviations (post hoc; `followup.json`, `blend_1625.json`; no new download)

- **D1 — comparison stars from the whole-quadrant psfcat**, replacing the failed image-based set.
  - Stars within 0.5 mag of the target's PSF magnitude, flags 0 in every bracket exposure, more than 15″ away:
    164 stars (1624) and 243 stars (1625).
  - **1624:** all 164 are present in every exposure, NDs included. Median offset is −0.003 to +0.007 mag, so zero
    points and transparency are consistent.
  - **1625:** median offset is −0.020 to +0.026 mag.
  - The scale enters only as a ratio; PSF magnitudes are never compared with aperture fluxes.
- **D2 — 1624 column-interpolated ratio.**
  - **Method** (the PLAN 2.2 rule with the D1 scale): bad pixels are replaced by the row's nearest good pixels. The
    same pixel offsets are interpolated in the 4 detections.
  - **R_i per ND** (`verdict.json` → `nd_R_i_by_expid`): 1.09, 1.14, 0.85, 0.90, 1.02, 0.93, **0.68**, 0.96 and 0.93.
    The same-offset bracket spreads are 4–31 %.
  - **Classes:** 8 NDs are "data problem supported" (R_i ≥ 0.75 and psfcat gap). 206350717 is **undetermined**:
    R_i 0.68, a deficit of 0.32 against a same-offset bracket spread of 0.215 (no σ stored).
- **D2b — position against the columns** (psfcat x of the columns, 2nd-order fit rms ≤ 0.07 px).
  - The target lies **0.00–1.70 px** from a LOW RESPONSIVITY column in all 9 NDs and **2.11–11.06 px** away in all 4
    detections. In 13 of 13 exposures, a position near the column goes with a catalogue miss.
  - **Not informative:** only one psfcat star (snr ≥ 20) falls on the column in each frame, the target itself.
- **D3 — 1625 neighbour test.**
  - A psfcat source of g 14.37–14.39 lies 6.55–6.71″ away in every exposure.
  - In the NDs (seeing 4.05–4.64″) the seeing-scaled aperture is contaminated: raw R 3.5–6.1.
  - **Wing estimate:** the same aperture placed at the target's position rotated 90°, 180° and 270° about the
    neighbour, without a PSF fit.
  - **Wing-subtracted ratios:**
    - 1.46 ± 0.69 (214237995), 2.01 ± 0.31 (214543142), 0.76 ± 0.41 (214647008);
    - in the detections 0.88–1.11.
    - All 3 NDs are **unknown (too shallow):** a 1-mag drop would be under 3σ.
  - **Same-exposure blend mechanism:** the one other faint star (g 18.55) with a 4.3-mag-brighter neighbour at 5–8″ is
    missing in both covered ND catalogues. Isolated similar stars are missing at 9/162 and 2/162. That is n = 1 pair:
    suggestive, not decisive.

## 4. Answers per source

### 1624 / 708205400009336: **undetermined**, with strong support for a data problem
- **Brightness drop confirmed in the images?** No. All 4 detections are normal (R 0.97–1.06, σ 0.04). No ND
  exposure is judgeable, because all 9 carry LOW RESPONSIVITY pixels in the aperture. The raw 2″ / seeing fluxes
  (R −0.52 to 0.60) are not used as a brightness, since they include the depressed column.
- **A specific data problem that explains the flag?** Yes, supported:
  - a detector column flagged LOW RESPONSIVITY by the official mask in every exposure;
  - the target on or within 1.7 px of it in exactly the 9 NDs, and ≥ 2.1 px away in the 4 detections;
  - the psfcat row missing in exactly those 9;
  - 164 comparison stars present in all of them;
  - column-interpolated R_i ≥ 0.75 in 8 of 9.
- **Which exposures are judgeable?**
  - judgeable: the 4 detections;
  - unknown: all 9 NDs (bad pixels);
  - data problem supported: 8 NDs (deviation D2);
  - undetermined: 206350717 (R_i 0.68).
- **What is missing to close it:** a flux measurement that does not use the bad columns, either a mask-aware PSF
  fit or forced photometry on difference images. Neither was built or requested here (owner rule; ZTF forced
  photometry needs an account).

### 1625 / 709212400017577: **undetermined**; likely a catalogue blend in wide seeing, not shown
- **Brightness drop confirmed?** No. None of the 3 NDs is judgeable. The wing-subtracted ratios (0.76–2.01) are not
  low, but their errors (0.31–0.69) cannot separate a 1-mag drop.
- **A specific data problem?** Plausible, not proven:
  - seeing 4.05–4.64″ in all 3 NDs against 1.86–2.82″ in the detections;
  - a psfcat star 3.9 mag (about 35×) brighter, 6.6″ away; psfcat keeps the neighbour and loses the target exactly in these
    3 exposures;
  - one comparable pair behaves the same way.
- **Which exposures are judgeable?** Judgeable: the 4 detections (normal). Unknown: the 3 NDs (too shallow after
  wing subtraction).
- **What is missing:** a good-seeing exposure inside HMJD 59896.38–59900.47, or a neighbour-subtracted PSF /
  difference-image measurement. The window was not searched further; the run has no other-band epochs.

- **Neither source** shows real-dimming evidence, so no variable-star, occultation or known-event check was opened.
  Nothing is escalated.

## 5. #385 batch, for completeness (`../v16/summary.json`)
- **Queued and started:** 8 areas, 500 calibration groups each (4,000).
- **Passed calibration:** 3 areas (1627, 1624, 1625).
- **Searched:** 5,995 groups, 11,919 series:
  - 1627: 2,064 groups, 9 runs, 0 survivors;
  - 1624: 1,887 groups, 16 runs, 1 survivor;
  - 1625: 2,044 groups, 5 runs, 1 survivor.
- **Failed areas:** 5 (187, 1626, 2274, 1358, 1357). Their 5,914 search groups were never read.
- **Cumulative searched with v13:** 17,773 groups. None gives an image-confirmed target event; that is not the same
  as "no star dimmed". The 30 flagged runs removed by vetting were not image-checked.

## Deviations
- **D1, D2 / D2b and D3 above:** added after the planned comparison set failed. Each is labelled where used.
- **`blend_script.py.txt`:** the pair test was first run with the target counted as one of the pairs (2 pairs). That
  was fixed before recording; the stored result has 1 pair.
- **The column test (D2b):** has no independent on-column stars, so it shows only the target's own pattern.
