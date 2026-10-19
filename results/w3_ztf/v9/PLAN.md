# W3 ZTF v9 — official single-exposure PSF catalogues (psfcat) as a data entry (plan, fixed before download)

Owner authorization: the round after #352. Scope: no new photometry code, no relaxed quality rule, no new search.
The repository had no earlier psfcat or ZFPS work (grep for `psfcat`, `ZFPS`, `forcedphot` in scripts and results:
none).

## 1. Sample (fixed now, from `../v8/image_trial_selection.json`; metadata only)
- The same 3 calibration stars: objectids 672208200000073, 672208200000239 and 672208200002085. All lie in field 672,
  CCD 8, quadrant 2, so one psfcat file per exposure serves all three.
- **Overlap exposures:** 6 of the 8 v8 control exposures, at indices `linspace(0, 7, 6)` by MJD: 193335993,
  195727655, 206514034, 206813989, 309032905 and 309627976.
- **Gap exposures:** 6 of the 20 v8 gap exposures, at indices `linspace(0, 19, 6)` by MJD: 212252128, 215545043,
  228830788, 241116549, 248752283 and 249053293.
- Products are located from the official metadata (`filefracday`, field, CCD, quadrant, filter) through the IBE
  product path `…_psfcat.fits`. Each file is opened to confirm it is readable.

## 2. Matching (by sky position, not by objectid)
- The target is the nearest psfcat source within 1.5″ of the star position in `../v6/sample.json`.
- The match is **ambiguous** if a second psfcat source lies within 3″, and **missing** if no source lies within 1.5″.
  A missing source is never read as a disappearance or as an upper limit.
- Stored per match: separation, the second-nearest distance, exposure time, filter, the psfcat flag columns,
  magnitude and error with their header and column definitions, and the image infobits from TAP.

## 3. Overlap test (fixed before any comparison is looked at)
- **Pairs:** star × overlap exposure where the DR24 light-curve API has a catflags = 0 detection at that expid, and
  the psfcat match is unique with flags = 0.
- **Δ** = psfcat magnitude − DR24 magnitude. The psfcat magnitude is computed as the header and the IRSA
  documentation define it (zero point; any colour term stated). No coefficient is fitted to the data.
- **Pass if all hold:**
  - n ≥ 8 pairs;
  - |median Δ| ≤ 0.05 mag (the release notes quote light-curve bias corrections up to about 0.025 mag);
  - robust scatter (1.4826 × MAD) of Δ ≤ max(0.03 mag, 2 × median quoted combined error).
- Only overlap exposures (with DR24 detections) enter the test, not gap exposures.

## 4. If the overlap test passes
- **Gap fill:** psfcat for gap exposures of the quadrant (infobits 0, the star ≥ 10″ inside the footprint, no DR24
  catflags = 0 detection for at least one of the 3 stars).
  - At most 60 distinct exposures, evenly spaced by MJD if more; ≤ 200 MB in total; each file downloaded once and
    reused for all 3 stars.
- **Usable observation:** a unique match with flags = 0 in an infobits = 0 exposure.
- **Reported per star:** the usable gap observations added, and the count of usable observations in the gap before and
  after the fill, against the frozen 20-epoch baseline. This is indicative only: the baseline in the pre-screen is
  counted per event.

## 5. If the overlap test fails, or for exposures psfcat cannot fill
- Check whether difference images (`scimrefdiffimg.fits.fz`) exist for the gap exposures (HTTP HEAD only), which
  ZFPS needs.
- Prepare a ZFPS request package for the 3 stars:
  - coordinates, a JD range, the official request form and result-reading steps;
  - quality and control checks fixed in advance.
- Write a registration-request draft for the owner. Nothing is sent and no account is created.
