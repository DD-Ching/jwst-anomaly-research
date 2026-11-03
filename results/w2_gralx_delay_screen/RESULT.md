# Lensed quasars with published delays: does any lack a lens galaxy in deep high-resolution imaging? (2026-11-03)

Owner-requested screen, at most 60 minutes, published measurements and image analyses only; it stops at candidate
screening.

**Population:** the 73 systems of Gaia GraL X Table 2 (arXiv:2603.12810v1, App. B).

**Question:** does any system have both:
- a reliable measured delay;
- high-resolution imaging deep enough to show an ordinary lens galaxy, in which none is detected?

**Sources** (in this directory):
- `table2.json`: Table 2 transcribed from the pinned LaTeX source;
- `coords.json`, `dedup.json`;
- `deflectors.json`: one cited source per system, curated by hand;
- `strides_lenslight.json`: the STRIDES lens-light table, transcribed;
- `summary.json`;
- script `scripts/w2_gralx_delay_screen.py`.

**Labels:**
- Delays and classes: published.
- Matches, counts and colour differences: derived.

## Verdict

**Zero systems qualify, so no named object is delivered.**

**Main reason:**
- 72 of the 73 have a deflector seen in imaging or spectroscopy; 46 of them have a lens redshift.
- The only system without a detected deflector is APM08279+5255. It fails both conditions:
  - its Table 2 "delay" is a reverberation-mapping lag, not an inter-image delay;
  - its high-resolution data are dominated by a compact, very bright quasar, and no limit was found.
- No system carries a quantitative non-detection limit. So there is no case in which to bound the redshift, M/L and
  depth of an ordinary lens.

**Selection effect:** delay monitoring targets confirmed lenses. Most already have HST photometry or a lens redshift
before their delays are measured.

## 1. Population and dedup
- **Table 2:** 73 systems, with 183 delays that carry an error and 13 that do not.
  - The arXiv PDF prints only part of the table. The LaTeX source holds all rows, most of them commented out.
- **Dedup** (`dedup.json`): 14 of the 73 were already in the repo's W1/W2 records (`w12_niq`, `w2_desi_pairs`;
  D-064), all as **control lenses**:
  SDSS0246-0825, J0607-2152, SDSS0806+2006, SDSS J0832+0404, SDSS J0924+0219, FBQ0951+2635, SDSS1001+5027,
  PG1115+080, SDSS1226-0006, SDSSJ1335+0118, SDSSJ1455+1447, SDSSJ1515+1511, SDSSJ1620+1203 and J2205-3727.
  - None was tested there for a dark deflector. D-064's LS DR10 test found the lens in 0 of 5 real controls, so its
    "none" carries no information.
  - D-056's decided systems are not stored as a list; of those DECISIONS.md names, J1515+1511 is in this population.
  - No rejected pair, binary or NIQ system from D-064, J2244 or J0740 is among the 73.

## 2. Deflector classes (`deflectors.json`, `summary.json`)
| Class | n | Systems |
|---|---|---|
| Galaxy seen | 72 | all others; 46 with a lens redshift |
| Quantitative non-detection limit | 0 | — |
| Data too shallow or blended | 1 | APM08279+5255 |
| Catalogue field blank | 0 | — |

**APM08279+5255:**
- HST NICMOS's third image C lies between A and B and is "rather than being the lensing galaxy" (Ibata+1999).
- Keck high-resolution imaging only suggests a massive lens "at z~3" (Egami+2000).
- Feruglio+2017 call it "the, so far undetected, lens galaxy".
- No limit was found in these abstracts.

**Weakest detections**, still counted as "seen":
- **HST lens light from the automated STRIDES fit only** (Schmidt+2023), plus VHS Ks for 2M1134: DESJ0029-3814,
  WISEJ0259-1635, J0420-4037, GRAL J1131-4419 and 2MASSJ1134-2103.
- **Ground imaging only, no lens redshift:** DESJ0501-4118, J050129.5-073309 (a "candidate lensing galaxy" in LS DR10),
  J2124+1632, SDSS1650+4251 and Q2138-431 (a 5σ VISTA J detection, Hawkins+2021).

## 3. Delays checked against their originals
The original delay was read only where the deflector condition might hold or the Table 2 entry has no error. For
every other system the conjunction already fails on the deflector.

| System | Table 2 | Original | Reliable lens delay? |
|---|---|---|---|
| APM08279+5255 | AB 900 d | Si IV / C IV reverberation lag of ~900 rest-frame days (Saturni+2016) | no: not an inter-image delay |
| SDSS J1640+1932 | 23–25 d, no errors | SIE prediction from CFHT imaging (Wang+2017) | no: model |
| Q2138-431 | AB 410 d | "image B to lead image A by around a year, consistent with 410 days from the mass model" (Hawkins+2021) | no: model value, no error |
| GRAL J1131-4419 | BC 3.5 ± 0.7, BD -10.1 ± 0.8, CD -13.6 ± 2.6 | (A+B)–C 3.5, (A+B)–D -10.1, C–D -13.6 d; covariance diagonal 0.5, 0.6, 0.7 d² (TDCOSMO XVII); 7 % best single, 4 % combined | yes |

Table 2 relabels the blended pair A+B as B, and lists 2.6 d for C–D, where the original covariance gives √0.7 d.

## 4. Closest case, not delivered: GRAL J1131-4419
- **Delay:** reliable (§3). HST WFC3 imaging exists.
- **Spectra:** Keck/LRIS spectra "show no evidence of the central lens galaxy" (Wertz+2019). The model puts the lens at
  (-0.424″, -0.744″) from A.
- **Imaging:** the only detection is the automated STRIDES lens-light fit:
  - n = 4 (fixed), θ_e 0.29″, I_e(F160W) 18.31 ± 0.03 and I_e(F814W) 24.28 mag arcsec⁻², plus an n = 1 component;
  - F814W − F160W = 5.97 mag.
  - The colour is not a usable test. Across that table, F475X − F814W is 1.02–1.04 mag in all 44 rows (the blue band
    is tied), and F814W − F160W spans 1.37–6.96 mag (median 3.475) for clear lenses too.
- **Why it does not qualify:** a published image analysis reports lens light, so it is not a non-detection.
- **What would move it:** a published, person-checked PSF-subtracted HST F160W analysis that finds no lens light at the
  model position and states a limit. A spectroscopic lens redshift would instead confirm "seen".

## 5. Not applicable
The redshift, M/L and depth check of the ordinary lenses still allowed needs a quantitative non-detection limit. There
are 0 such systems.

## 6. Limits of this screen
- **How classes were read:**
  - One cited source per system, read from abstracts or full text by three research subagents.
  - The main session re-read the decisive cases:
    - APM08279+5255, Q2138-431 and SDSS J1640+1932 (abstracts);
    - GRAL J1131-4419 (the STRIDES and TDCOSMO XVII LaTeX sources, and the Wertz+2019 abstract).
- **Weaker sources:** several detections rest on CASTLES photometry pages, and some lens redshifts on Mg II absorbers
  (noted per system).
- **Delay originals** were read for 4 systems only.
- **Not accessed:**
  - the Cambridge lensed-quasar database (HTTP 500);
  - the SIMBAD CDS host (the Harvard mirror was used);
  - the COSMOGRAIL XIX PDF (too large);
  - Bekov+2024 and Rivera+2024 (not on arXiv).
