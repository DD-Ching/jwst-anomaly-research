# Lensed quasars with published delays: does any lack a lens galaxy in deep high-resolution imaging? (2026-11-03, revised 2026-11-04)

Owner-requested screen using published measurements and image analyses only; it stops at candidate screening.

**Population:** the 73 systems of Gaia GraL X Table 2 (arXiv:2603.12810v1, App. B).

**Question:** does any system have both:
- a reliable measured delay;
- high-resolution imaging deep enough to show an ordinary lens galaxy, in which none is detected?

**Revision (2026-11-04, owner review):**
- Keeley et al. 2026 (Phys. Rev. D 114, 043043, doi:10.1103/53tj-sjlg) Tables I and VI correct the evidence for
  GRAL J1131-4419, WISE J0259-1635 and 2M1134-2103 (§4).
- The earlier suggestion to redo J1131's HST analysis is withdrawn.

**Sources** (in this directory):
- `table2.json`: Table 2 transcribed from the pinned LaTeX source;
- `coords.json`, `dedup.json`;
- `deflectors.json`: cited sources per system, curated by hand;
- `strides_lenslight.json`: the STRIDES lens-light table, transcribed;
- `keeley2026_tables.json`: Keeley+2026 Tables I and VI, transcribed;
- `summary.json`;
- script `scripts/w2_gralx_delay_screen.py`.

**Labels:**
- Delays and classes: published.
- Matches, counts and colour differences: derived.

## Verdict

**Among these 73 systems, none was found with both a reliable inter-image delay and deep imaging that gives a
quantitative non-detection limit for an ordinary lens galaxy.** No named object is delivered.

**Basis:**
- 72 of the 73 have a deflector seen in imaging or spectroscopy.
  - 47 of them have a lens redshift as cited.
  - One of those, J1131's 0.472, is quoted from Keeley+2026 Table I; its original spectrum is unpublished.
- The only system without a detected deflector is APM08279+5255. It fails both conditions:
  - its Table 2 "delay" is a reverberation-mapping lag, not an inter-image delay;
  - its high-resolution data are dominated by a compact, very bright quasar, and no limit was found.
- No system carries a quantitative non-detection limit.
  - Seven systems have no deflector detected in Keeley+2026's mid-IR fits, but no limit is given, and each has a
    deflector detected in other data (§4).

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
| Galaxy seen | 72 | all others; 47 with a lens redshift as cited |
| Quantitative non-detection limit | 0 | — |
| Data too shallow or blended | 1 | APM08279+5255 |
| Catalogue field blank | 0 | — |

**APM08279+5255:**
- HST NICMOS's third image C lies between A and B and is "rather than being the lensing galaxy" (Ibata+1999).
- Keck high-resolution imaging only suggests a massive lens "at z~3" (Egami+2000).
- Feruglio+2017 call it "the, so far undetected, lens galaxy".
- No limit was found in these abstracts.

**Weakest detections**, still counted as "seen":
- **HST lens light from the automated STRIDES fit only** (Schmidt+2023): DESJ0029-3814 and J0420-4037.
  - WISEJ0259-1635, GRAL J1131-4419 and 2MASSJ1134-2103 were listed here before the revision. Keeley+2026 Table VI also
    fits their deflector light (§4).
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

## 4. Correction from Keeley+2026 Tables I and VI (`keeley2026_tables.json`)
In Table VI, Sérsic n is fixed to 4 and the deflector centroid is shared across filters.

| System (Table I name) | Table VI deflector-light fit (JWST/MIRI) | Other detections | Lens redshift |
|---|---|---|---|
| GRAL J1131-4419 (GRAL1131) | F560W R_s 0.43″, q 0.60, θ 173°, centroid (-0.073, 0.774)″; F1280W R_s 0.072″, q 0.40, θ 165° | HST STRIDES fit (I_e F160W 18.31 ± 0.03) | 0.472: Table I quoted value |
| WISE J0259-1635 (J0259) | F560W R_s 0.08″, q 0.39, θ 311°, centroid (0.028, 0.687)″ | HST STRIDES fit | Table I lists 0.904, not adopted |
| 2M1134-2103 (2M1134) | F560W R_s 0.94″, q 0.59, θ 237°, centroid (1.196, 1.555)″ | VHS Ks (Rusu+2019); HST STRIDES fit | Table I lists 0.66, not adopted |

- **J1131's lens redshift:** 0.472 is a Table I quoted value. Footnote e cites "Sluse et al. (2025, in prep)".
  - The original spectrum is unpublished, so it was not verified.
  - Wertz+2019's Keck/LRIS spectra "show no evidence of the central lens galaxy".
- **J0259 and 2M1134 redshifts:** Table I gives them without a separate reference (footnote a: the discovery paper).
  They were not checked and are not counted.
- **Classes:** the three systems stay "galaxy seen". Their evidence is now HST plus JWST/MIRI. This replaces the
  earlier "automated STRIDES fit only" reading.
- **Same table, other systems:** Table VI also fits the deflector light of J2205-3727 (F560W, F1280W) and of five
  lenses already seen in HST: HE0435-1223, SDSS J0924+0219, PG1115+080, RXJ1131-1231 and WFI2033-4723.
- **Seven systems absent from Table VI.** They are in Table I but not in Table VI, whose caption says "Lenses not
  included in this table did not have detected deflectors or satellites":
  J0607-2152, J0659+1629, H1413+117, J1537-3010, PS J1606-2333, WGDJ2026-4536 and GraLJ203802-400815.
  - The table gives no limit.
  - Each has a deflector detected in other data (`deflectors.json`), so its class is unchanged.
- **STRIDES colours:** they are not used as evidence. In that table F475X − F814W is 1.02–1.04 mag in all 44 rows,
  i.e. tied.

## 5. Not applicable
The redshift, M/L and depth check of the ordinary lenses still allowed needs a quantitative non-detection limit. There
are 0 such systems.

## 6. Limits of this screen
- **How classes were read:**
  - Cited sources per system, read from abstracts or full text by three research subagents.
  - The main session re-read the decisive cases:
    - APM08279+5255, Q2138-431 and SDSS J1640+1932 (abstracts);
    - GRAL J1131-4419 (the STRIDES and TDCOSMO XVII LaTeX sources, and the Wertz+2019 abstract);
    - Keeley+2026 Tables I and VI.
- **Keeley+2026:**
  - The tables come from the arXiv v3 LaTeX source (arXiv:2511.07765v3).
  - The DOI metadata were checked with Crossref.
  - The table numbers follow the source's counters: Table I is the lens sample and Table VI the deflector-light fits.
  - The published PRD text was not compared line by line.
- **Weaker sources:** several detections rest on CASTLES photometry pages, and some lens redshifts on Mg II absorbers
  (noted per system).
- **Delay originals** were read for 4 systems only.
- **Not accessed:**
  - the Cambridge lensed-quasar database (HTTP 500);
  - the SIMBAD CDS host (the Harvard mirror was used);
  - the COSMOGRAIL XIX PDF (too large);
  - Bekov+2024 and Rivera+2024 (not on arXiv).
