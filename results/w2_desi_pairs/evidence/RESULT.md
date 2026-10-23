# W2 DESI pairs, evidence audit: of the 13 rejected pairs, 3 have stronger ordinary explanations, 2 are disputed and 8 stay open; J2244−0550 is a likely lens in its own source paper and has public HST and Keck AO imaging with no published use found; the DESI spectral-classifier line stops

- **Plan:** `PLAN.md` (claim 0162310, 01:48 UTC).
- **Scope:** a bounded audit of identity, published evidence and follow-up value.
  - No spectrum was downloaded or read, and #472's statistic, threshold and failed verdict are unchanged.
  - No new classifier or PSF tool was built.
- **Tables:**
  - the evidence table is `EVIDENCE.md`, with its machine-readable form `evidence.json`;
  - the inputs are `identity.json`, `literature.json` and `archives.json`, written by `scripts/w2_desi_evidence.py`.
- **Labels:**
  - positions, TARGETIDs, magnitudes and archive rows: observed;
  - matches and flags: derived;
  - classes: catalogued or published, quoted with their basis in `literature.json`;
  - the readings in §5: this audit's judgement.
- **Not found:** "not found" means not found in the searches listed in `literature.json` → `searched`, not "does
  not exist".

## 1. Identity (`EVIDENCE.md` §A, `identity.json`)
- **12 of 13 match the catalogued pair:** the separation agrees within 0.5″ (#472's rule) and each DESI fibre lies
  within 0.25″ of its own Gaia DR3 source.
- **J2351+1343 is not.**
  - The SQLS companion has no SQLS redshift and no DESI spectrum. It is very red (g 22.72, z 18.74), has a Gaia
    parallax at 3.2σ, and SLED calls the system QUASAR+STAR.
  - DESI's second fibre is on a third, faint source (g 23.45) 2.05″ from the quasar, which is 3.87 mag brighter in
    Tractor g. That fibre's spectrum may be crosstalk.
  - #472 matched the pair only because its separation fell within the 0.5″ rule.
- **Fibre overlap:** 6 pairs have fibre separations below the 1.5″ fibre diameter (0.82–1.31″). Their DESI spectra
  mix the two components. Jing et al. class one of them, J0909−0749, as a "Contaminant".
- **TARGETIDs:** for J1428 and J1511, Jing et al.'s catalogue lists another TARGETID for one component, 0.06″ from
  #472's, so DR1 may hold two coadds of it.

## 2. What the labels rest on (the priority checks)
- **J2244−0550: the source paper contradicts itself.**
  - Lemon et al. 2023 (arXiv:2206.07714, only v1; LaTeX source read) discuss it in §4.2 "Likely Lensed Quasars"
    (§4.2.16): "almost identical profiles", no lens residuals in Pan-STARRS or LS.
  - Its own table 1 gives "NIQ"; VizieR, built from the journal's electronic version, gives "UQP".
  - The published MNRAS text was not read (OUP returned HTTP 403), so a change in the journal version is not
    excluded.
  - The table class does not separate the two groups: 13 of the 17 §4.2 likely lenses carry VizieR class UQP (§3).
- **J0130+0725: "dual favoured, lensing not excluded" is kept.**
  - Chen et al. (arXiv:2607.00093 §3.1.1) used the same DESI DR1 spectra, because their own Palomar slit blended the
    two images.
  - Independent of DESI:
    - Lemon et al. report absorption differences and read a projected pair, while saying one source "cannot
      necessarily be ruled out";
    - the D-064 HST F814W limit is ≈ 23 with no lens light.
- **J0728+2607: already settled as far as public data go; #472 added nothing.**
  - Lemon et al.'s NOT spectra show "large differences in their emission line profiles and continua, thus making it a
    likely physical quasar pair".
  - The D-064 HST check found no lens light (≥ 1.8 mag below any ordinary lens).
  - Listing it as a #472 target repeated finished work.
- **J0740+2926 and J1035+0752 (the Hennawi overlap): "binary" rests on not seeing a lens galaxy.**
  - Hennawi et al. 2006 §5.1 say the APO discovery spectra are mostly "too low a signal to noise ratio to make
    convincing arguments for or against the lensing hypothesis based on spectral dissimilarity". These two are not
    among the five with Keck spectra.
  - The binary call rests on deeper UH88 images without a lens galaxy, and "does not strictly speaking confirm the
    binary hypothesis".
  - The catalogue χ² (4.5 and 41.9) is a photometric colour-similarity statistic, not a spectral test.
  - SQLS DR3 copies the class from Hennawi.
  - Later papers dispute both:
    - Jing et al. class J1035 as a lensed-quasar candidate (LQC) from the same DESI spectra;
    - Ji et al. (arXiv:2412.08397; DESI EDR + SDSS) list J0740 as a lensed-quasar candidate.
  - #472 used both as negative controls (D = 1.04 and 1.35). Their "binary" truth label is therefore weaker than #472
    assumed, which adds to the gate's failure being uninformative.
- **J1511+3357 and J1552+3009, SQLS "Binary QSO": the spectra are lens-like and the label is again a non-detection.**
  - SQLS IV (Inada et al. 2010, §4.2) took Subaru spectra of both components:
    - J1511: "similar shapes" and "similar" continua;
    - J1552: the fainter component's MgII is "slightly broader", with similar continua.
  - Both are called binaries because no galaxy is seen to M_I ≈ −19.9 and −19.8 in deep Subaru I images (0.6″
    seeing).
  - SQLS §4.1 itself: "our selection process would reject 'dark lenses'".
- **J1613+1708 and J2351+1343, "no lensing object": imaging only, no SQLS spectra.**
  - For J1613, that both components are quasars at one redshift rests only on overlapping DESI fibres (0.82″).
  - For J2351, see §1.

## 3. Record corrections
- **#472 wording** (numbers, T and verdict unchanged) in `../RESULT.md`, `../PLAN.md` (Errata), CHANGELOG and the
  D-064 addendum:
  - D ≤ T means "no difference above threshold", not "one source";
  - 9 / 9 is a G2 pass count, not an accuracy for "one source" (6 of 7 binaries also met D ≤ T).
- **D-064 sample, recorded and not re-run here:**
  - Lemon's UQP class includes 13 of the 17 systems the paper discusses as likely lenses: J0138+4841, J0221+0555,
    J0326−3122, J0332−6608, J0346+2154, J0916−2848, J1008−2911, J1019−1322, J1033−8249, J1036−8544, J1103−1005,
    J2244−0550 and J2255+8009.
  - D-064 maps UQP to "rejected", so these 13 are counted as rejected pairs.
  - The SQLS and Hennawi "binary" labels of J0740, J1035, J1511 and J1552 rest on not seeing a lens galaxy. They are
    therefore not independent of the dark-lens hypothesis W2 asks about.
- **#472's target list:** J2351+1343's DESI pair is not the catalogued pair.

## 4. Up to 3 open pairs: is there public data not yet used that fills the gap? (`archives.json`)
- **Selection rule:** applied to the literature evidence only, never to #472's score.
  - Independent resolved spectra are reported as similar, and the only argument against a lens is a non-detection of
    the lens galaxy.
  - Ranked by how lens-like the spectra are and how weak the imaging is.
  - Chosen: J2244−0550 ("almost identical"; ground imaging only), J1511+3357 (similar Subaru spectra; "binary" from a
    non-detection) and J1202+0703 (similar except one absorber; ground imaging).
  - **Deviation:** the archive queries for J2244 and J1202 ran while the third pair was still being chosen. The rule
    above does not use archive results; J1202 has no new imaging and was still kept.
- **Queries (metadata only, 2026-10-23):**
  - MAST CAOM, all missions, footprint containing the position within 3″;
  - ESO ObsCore (same footprint test);
  - KOA NIRC2, OSIRIS, KCWI, LRIS, ESI and MOSFIRE (pointing within 60″).
  - **Not checked:** CADC (CFHT, Gemini, Subaru) returned HTTP 503 or a connection error twice; SMOKA; HSC (needs an
    account); Euclid DR1. Euclid Q1's fields do not cover these positions.

| pair | gap | public data found that no published paper is known to use | access | could it change the call? |
|---|---|---|---|---|
| J2244−0550 | lens light at 1.0″ separation | HST ACS/WFC F814W, 674 s, 2024-08-31, GO 17308 (Lemon, "An HST Gap Program for lensed Quasars"), `jf56c1010`; HST WFC3/IR F160W 422 s, WFC3/UVIS F475W and F625W 450 s each, 2024-07-05, GO 17455 (Shen, "A Snapshot Survey of Sub-arcsec Dual Quasars and Lenses at z>2"), `if8a5k010/020/030`; Keck OSIRIS imager (AO), Kp, 8 × 180 s, 2020-11-01, C276 (Djorgovski); Keck LRIS 1.5″ long slit, 2 × 300 s on each of 2020-10-18 and 2022-10-24 (GraL) | MAST `dataRights` PUBLIC; KOA proprietary interval 18 months, now past, download not tried | **yes**: a lens galaxy in F160W / F814W would make it an ordinary lens; a calibrated non-detection keeps it a lens-like pair without a visible deflector, the W2 case |
| J1511+3357 | lens light vs host light at 1.1″; a time delay | **none found** (MAST: SDSS, PS1, GALEX, TESS only; nothing in ESO or KOA) | — | only the same data could be re-analysed: no new observational information found (not confirmed in CADC or SMOKA) |
| J1202+0703 | lens light at 1.2″ | Keck LRIS 1.5″ long slit, 2 × 600 s on each of 2022-06-23 and 2023-03-28 (GraL, Djorgovski); **no high-resolution imaging found** | KOA proprietary interval 18 months, now past, download not tried | weakly: seeing-limited spectra of a 1.2″ pair give another epoch of flux ratios, not a lens-galaxy test |

- **HST resolution:** the ACS/WFC (0.05″ pixels) and WFC3/IR (0.13″ pixels) scales are from the STScI instrument
  handbooks. Both resolve a 1.0″ pair, which the ground-based imaging Lemon et al. used could not.
- **Cost of the J2244 test:**
  - It needs PSF-subtracted photometry of a 1.0″ pair. The D-064 HST chain was only run on 2.06″ and 2.15″ pairs.
    The close-pair chain failed its no-lens negative control (D-060 addendum, S/N 13).
  - So the test is not "run the existing tool". It needs a pre-registered PSF test with a known lens and a no-lens
    pair at ~1″ from the same programmes as controls, and it stops if those controls fail.
  - Both HST programmes are the proposers' own follow-up of such candidates, so a result may already be in
    preparation.

## 5. Readings (`EVIDENCE.md` §B)
- **Stronger ordinary explanation (3):**
  - J0728+2607: binary;
  - J0130+0725: projected or dual favoured, lens not excluded;
  - J2351+1343: the catalogued system is a quasar plus a star, and the DESI pair is a different pair.
- **Disputed (2):** J0740+2926 and J1035+0752. Catalogued as binaries from lens non-detection and low-S/N spectra;
  later lens candidates.
- **Open (8):**
  - J2244−0550: a likely lens in its source paper;
  - J1511+3357 and J1552+3009: "binary" from a non-detection, with similar spectra;
  - J1202+0703, J0133+0816 and J1428+0500: UQP;
  - J1613+1708: imaging-only rejection, DESI fibres overlap;
  - J0909−0749: DESI pair unusable.
- **Not shown by anyone:** that any of these is a dark lens. Similar spectra do not show one source, and no visible
  lens does not exclude an ordinary lens.

## 6. Cost and review
- **Wall clock:** 01:48 UTC to the end recorded in `cost.json`.
- **Downloads:** no spectrum. Archive and catalogue metadata (TAP, VizieR, MAST) are small. The literature
  subagents fetched arXiv PDFs and TeX sources, about 0.2 GB on scratch disk, deleted after the round (`cost.json`).
- **Deviation:** one of the three literature subagents (SQLS / Hennawi) stalled on a repository-access request and
  was stopped. Its downloaded papers were read directly instead.
- **Review:** see the PR.
