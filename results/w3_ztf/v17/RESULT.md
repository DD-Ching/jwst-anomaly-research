# W3 ZTF v17 — evidence check of the two v16 final survivors: both undetermined, no real-dimming evidence

- **Post hoc:** done after the v16 images had been seen.
  - `PLAN.md` was committed before any recomputation (baff0fa).
  - Deviations D1–D3 are marked wherever they are used.
- **Labels:** IBE `_sciimg`, `_mskimg`, `_psfcat` and the DR24 light-curve API are `observed`; ratios and classes are
  `derived`.
- **Downloads** (measured bytes, `check.json`): 10.34 MB, made up of 40 cutouts (20 science, 20 mask) and 20 psfcat
  files.
- **Scope:** no search, no threshold change, no new veto rule. No limit is derived, and nothing here bears on
  wormholes.

## 1. Evidence (`evidence.json`, copied from committed v16 files; the 9th 1624 ND from `verdict.json`)

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
| official mask, 8 v16 NDs | bit 2 LOW RESPONSIVITY in the 2″ aperture in all 8; at the centre pixel in 5 | bit 1 (detection) only |
| ND maglimit / seeing | 19.08–19.82 for the 8 v16 NDs (19.94 with the 9th) / 1.92–2.19″ | 19.06–19.92 / 4.05–4.64″ |

- The v16 "≈ 2,450 µJy" figure is the low end of the four 2″ bracket detections. Their median is 2,496 µJy.
- The DR24 reference magnitude is r 15.18. That is a PSF total, and it is not compared with aperture fluxes.

## 2. Planned checks (`check.json`, `check.png`): **both sources undetermined, every exposure unknown**

1. **Consistency.**
   - All 20 science files are `_sciimg` (header IMGTYPE "object"), not difference images.
   - Flux is in µJy from each header MAGZP (25.92–26.07 for 1624).
   - The position lies inside every image.
   - The centroid within 3″ lies 0.24–0.81″ from the DR24 position (1624) in every exposure.
   - The DR24 records (LC API, by expid) agree with v16's ND / detection labels for every checked exposure. The 9th
     1624 ND (206350236) was found and checked.
2. **Aperture and mask** (radius = max(2″, seeing)).
   - **1624:** bad-bit pixels in the target aperture are 2–8 of 11–15 in all 9 NDs and 0 in all 4 detections.
   - **1624:** two LOW RESPONSIVITY columns cover every row of the 178-row cutout in every exposure. Away from the
     star they lie 400–3,080 DN below their neighbours in the NDs and 434–1,111 DN below in the detections.
   - **1625:** no bad bits.
3. **Comparison stars** (image-based, as planned). The 180″ cutouts held 0 (1624) and 1 (1625) stars in the planned
   range, below the required 3. **So every exposure of both sources is unknown, and the planned verdict is
   undetermined for both.**
4. **psfcat.** The same images and photometry feed DR24, so this is **not** an independent observation.
   - The target is missing (no source within 1.5″) in all 9 NDs of 1624 and all 3 NDs of 1625.
   - It is matched with flags 0 in every detection: PSF magnitude spread −0.014 to +0.005 mag for 1624 and
     −0.054 to +0.053 mag for 1625.
   - A missing row is not taken as zero flux.

## 3. Deviations (post hoc; `followup.json`, `blend_1625.json`, `verdict.json`; no new download)
- **D1 — comparison stars from the whole-quadrant psfcat**, replacing the failed set.
  - **Selection:** stars within 0.5 mag of the target's PSF magnitude, flags 0 in every bracket exposure, more than
    15″ away. Stars flagged in the tested exposure are excluded. This gives 164 stars (1624) and 243 (1625).
  - **1624:** all 164 are present in all 13 exposures; median offsets are −0.003 to +0.007 mag.
  - **1625:** median offsets are −0.020 to +0.026 mag.
  - **Caveat:** the D1 scale comes from PSF magnitudes and divides an aperture flux. It corrects zero point and
    transparency, but **not** seeing-dependent aperture loss.
  - **1625 under D1:** every exposure is unknown, because a brighter psfcat source sits 6.6″ away, inside the
    annulus. Without that contamination test, D1 would have classed two 1625 NDs "normal" on contaminated ratios
    (R 3.5 and 6.1). D3 replaces D1 for 1625.
- **D2 — 1624 column-interpolated ratio R_i** (the PLAN 2.2 rule with the D1 scale).
  - **Method:** the aperture's bad pixels are replaced by the row's nearest pixels outside the image's own bad mask.
    The same offsets are interpolated in the 4 detections.
  - **Review fix:** a first run used neighbours that were themselves bad; the stored values are corrected.
  - **Results per ND** (`verdict.json` → `nd_R_i_by_expid`): R_i 0.68–1.22.
    - **8 of 9 are ≥ 0.75** and lack a psfcat row, which is the post hoc "data problem supported" rule.
    - **206350717** is undetermined: R_i 0.68, σ ≈ 0.22 from the same-offset bracket spread.
    - **Depth:** the PLAN rule for this class has no depth condition. Measured against the frozen dim ratio 0.398,
      3 NDs are < 3σ above it (206250648, 206350476, 206350717). So 6 of 9 NDs exclude a 1-mag drop at 3σ, and 3 do
      not.
- **D2b — position against the columns** (psfcat x; 2nd-order fit rms ≤ 0.071 px).
  - **Distances:** the target lies 0.00–1.70 px from a LOW RESPONSIVITY column in all 9 NDs and 2.11–11.06 px away
    in all 4 detections. That split between 1.70 and 2.11 px was read off the data, post hoc.
  - **Star test:** the star list comes from detection 206550548, so that frame is not tested. In the other 12
    frames, no star other than the target lies within 1 px of the column. The target lies within 1 px in 7 NDs.
  - So the test has **no independent on-column stars** and shows only the target's own pattern.
- **D3 — 1625 neighbour test.**
  - A psfcat source 3.85–3.87 mag brighter (g 14.37–14.39) lies 6.55–6.71″ away in every exposure.
  - **Wing estimate:** the same aperture placed at the target's position rotated 90°, 180° and 270° about the
    neighbour, without a PSF fit. σ uses the random-aperture noise, the wing spread and the bracket scatter. The
    PLAN's control-MAD term is dropped (deviation).
  - **Wing-subtracted ratios:**
    - 1.46 ± 0.69 (214237995), 2.01 ± 0.31 (214543142), 0.76 ± 0.41 (214647008);
    - in the detections 0.88–1.11.
    - All 3 NDs are unknown (too shallow): a 1-mag drop would be under 3σ.
  - **Blend mechanism (same-exposure):** the one other faint star (g 18.55) with a 4.33-mag-brighter neighbour at
    5–8″ is missing in both ND catalogues that cover it. Isolated similar stars inside each frame are missing at
    4/157, 1/161 and 2/58. That is n = 1 pair.

## 4. Answers per source

### 1624 / 708205400009336: **undetermined**; a specific data problem is supported, not proven
- **Planned:** undetermined; every exposure is unknown (no comparison stars).
- **Brightness drop confirmed in the images?** No. Under D1 the 4 detections are normal (R 0.97–1.06, σ 0.04). No
  ND is judgeable: all 9 carry LOW RESPONSIVITY pixels in the aperture. The raw ND ratios (−0.51 to 0.60) include
  the depressed column and are not used as a brightness.
- **A specific data problem that explains the flag?** Supported, not proven:
  - the official mask flags a detector column in every exposure;
  - the target lies within 1.70 px of it in all 9 NDs, and ≥ 2.11 px away in all 4 detections (split post hoc);
  - psfcat lacks the target in those 9 NDs;
  - all 164 comparison stars are present (D1);
  - R_i ≥ 0.75 in 8 of 9 (D2), but only 6 of 9 exclude a 1-mag drop at 3σ.
- **Which exposures can be judged?**
  - judgeable: the 4 detections (D1);
  - unknown: all 9 NDs;
  - data problem supported (D2): 8;
  - undetermined: 206350717.
- **Missing to close it:** a flux measurement that does not use the bad columns, either a mask-aware PSF fit or
  difference-image forced photometry. Neither was built (owner rule); the ZTF forced-photometry service needs an
  account.

### 1625 / 709212400017577: **undetermined**; consistent with a wide-seeing blend, not tested to a decision
- **Planned:** undetermined; every exposure is unknown. **D1:** every exposure unknown (neighbour in the annulus).
- **Brightness drop confirmed?** No. Under D3 none of the 3 NDs is judgeable. The wing-subtracted ratios (0.76–2.01)
  are not low, but their errors (0.31–0.69) cannot separate a 1-mag drop.
- **A specific data problem?** Consistent with one, not shown:
  - seeing 4.05–4.64″ in all 3 NDs, against 1.86–2.82″ in the detections;
  - a 3.85–3.87-mag-brighter star 6.6″ away;
  - psfcat keeps the neighbour and lacks the target in exactly these 3 NDs;
  - one comparable pair behaves the same way.
- **Which exposures can be judged?** Judgeable: the 4 detections (D3, normal). Unknown: the 3 NDs.
- **Missing:** a good-seeing exposure inside HMJD 59896.38–59900.47 (the window was not searched further; no
  other-band epochs), or a neighbour-subtracted PSF / difference-image measurement.

- **Neither source** gives real-dimming evidence, so no variable-star, occultation or known-event check was opened.

## 5. #385 batch, for completeness (`../v16/summary.json`)
- **Queued and started:** 8 areas, each with 500 calibration groups (4,000).
- **Passed calibration:** 3 areas (1627, 1624, 1625).
- **Searched:** 5,995 groups, 11,919 series:
  - 1627: 2,064 groups, 9 runs, 0 survivors;
  - 1624: 1,887 groups, 16 runs, 1 survivor;
  - 1625: 2,044 groups, 5 runs, 1 survivor.
- **Failed areas:** 5 (187, 1626, 2274, 1358, 1357). Their 5,914 search groups were never read.
- **Cumulative searched with v13:** 17,773 groups (11,778 + 5,995 from the CHANGELOG; not stored as one number).
  None gives an image-confirmed target event. That is not the same as "no star dimmed".
- The 28 flagged runs removed by vetting (30 flagged, minus the 2 survivors' runs) were not image-checked.

## Deviations
- **D1, D2 / D2b, D3:** as above.
- **Review fixes:**
  - **D2 interpolation:** neighbours must lie outside the image's own mask. R_i moved by −0.09 to +0.14 (not
    stored); no class changed.
  - **D2b:** the star-list frame is not tested.
  - **D1:** flagged comparison stars are excluded; the bracket reference is scaled like R_i; a brighter source in
    the annulus makes the exposure unknown.
  - **Blend test:** a footprint proxy is used for isolated stars.
  - **`check_script.py.txt`:** fixes in branches that did not execute in the stored run, so `check.json` is
    unchanged.
- **Pair test:** first run with the target counted as a pair; fixed before recording.
