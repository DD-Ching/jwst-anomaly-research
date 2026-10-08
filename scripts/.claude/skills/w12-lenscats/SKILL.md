---
name: w12-lenscats
description: Recipe for the W1/W2 dark-deflector test on published strong-lens catalogues (negative-mass or Ellis-wormhole lens = multiple images with no visible deflector) - catalogue merge, coverage and depth, per-class deflector tests, vetting and limits - with its failed-approach rules. Use when changing src/jwst_anomaly/lenscats.py, scripts/w12_lenscats.py or scripts/w12_hsc_probe.py, or adding a lens catalogue or deeper imaging.
---

# W1/W2 in published lens catalogues

Rules: root CLAUDE.md, `scripts/CLAUDE.md`. Decisions: D-047, D-056 (catalogue test), D-060 (HST catalogue probe).
Results: docs/exotic_limits.md "W1/W2 in published lens catalogues".

## Commands
- `python scripts/w12_lenscats.py screen --out <dir>` then `vet --out <dir>` (Data Lab TAP; Tractor boxes cached in
  `<dir>/tractor_cache`, brick summary `bricks_dr10_south.csv`).
- `python scripts/w12_hsc_probe.py` (HST Hubble Source Catalog probe, D-060).

## Method
- Merge lenscat, Euclid Q1 Discovery Engine and SuGOHI by position (3″), pinned by sha256.
- Coverage and depth from the survey's footprint product (`ls_dr10.bricks_s`), never from "a source nearby".
- Only selections that can show a dark lens give a limit: lensed-quasar pair test (two PSF images ≥ 2″ apart;
  search circle about the images' centroid holding every image) and radio-interferometric systems. Galaxy-finder
  and sub-mm selections are insensitive.
- "faint galaxy", "blended", "too close" and position/mask problems are undecided and leave N; k and N per class.
- Limits s₉₅(k) / N, completeness assumed (state it); vet "none" systems against SIMBAD, published z_l, clusters,
  literature before anything else.

## Failed approaches (rules)
- A limit is valid only over systems where the test could have found the signal — prove it by injection.
- A recovery factor from deleting deflectors and re-running the same code is 1 by construction.
- Coverage from detected sources drops the dark configuration itself.
- Evaluate exclusion flags on every system, not only on flagged ones.
- Count k with the same selection mask as each class's N.
- A fold or cusp pair of a quad leaves the lens outside the pair circle.
- Data Lab TAP takes no table uploads or q3c: batch box ORs; split boxes at RA 0/360 and use full RA near the poles.
- lenscat types cluster-survey entries as "galaxy", has AGEL declination and SPT position errors and rounded
  positions, and keeps rejected candidates; "no lens redshift" ≠ "no lens".
- HST *catalogue* photometry cannot decide a dark deflector (a "none" is more likely a missed lens, 0.54): use
  PSF-subtracted image models or HSC PDR3 photometry.
- Known biases of the assumed completeness: unrelated faint sources remove dark-lens systems as "faint galaxy"
  (efficiency ~0.98); a colourless pair takes a compact PSF-typed lens as an image.
- LS DR10 finds no lens galaxy in real quasar lenses at 1.9–2.6″: 0/5 controls (lens light absorbed into the image
  PSFs). A "none" there is uninformative; always run a control sample of known lenses through the same chain.
- VizieR writes some coordinates as sexagesimal strings (`RA1` "h:m:s"); parse them, and fail loudly if nothing parses.
- Lemon et al. 2023 `z2`: a blank `n_z2` (or "zqso=") is a second quasar redshift; "z_lens="/"zgal=" are not.
- VizieR ASU-TSV headers carry the request time: pin the data lines, not the raw file.
