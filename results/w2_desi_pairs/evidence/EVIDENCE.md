# Evidence table: the 13 rejected lensed-quasar pairs with DESI DR1 spectra of both components

- **Sources:**
  - machine-readable `evidence.json`, built by `scripts/w2_desi_evidence.py build` from:
    - `identity.json`: DESI zpix rows, Gaia DR3 and LS DR10 Tractor matches;
    - the pinned D-064 catalogue rows;
    - `literature.json`: classes, quotes and references, with their keys;
  - `archives.json` holds the archive metadata for the 3 pairs of RESULT §4.
- **Labels:**
  - positions, TARGETIDs, magnitudes and archive rows: observed;
  - matches and flags: derived;
  - classes: catalogued or published;
  - readings: this audit's judgement.
- **Not used:** no spectrum was read and #472's statistic was not rerun.

## A. Identity (DESI DR1 `zpix`; separations in arcsec, catalogued / DESI fibres)
| pair | TARGETID a (RA, Dec) / TARGETID b (RA, Dec) | sep | DESI z a / b | programme | distinct Gaia DR3 match ≤ 0.5″ for both | Δg (Tractor) |
|---|---|---|---|---|---|---|
| J0130+0725 | 39627965927001229 (22.58586, +7.42091) / 39627965927001231 (22.58597, +7.42149) | 2.06 / 2.13 | 1.5530 / 1.5454 | dark/dark | yes | 0.7 |
| J0133+0816 | 39627983916368948 (23.26933, +8.27903) / 39627983916368957 (23.26976, +8.27880) | 1.55 / 1.74 | 1.2790 / 1.2793 | dark/dark | yes | 0.07 |
| J0728+2607 | 39628394240938051 (112.11564, +26.11704) / 39628394240938063 (112.11614, +26.11741) | 2.15 / 2.07 | 1.0259 / 1.0259 | dark/dark | yes | 0.81 |
| J0909−0749 | 39627600242411562 (137.49474, −7.81790) / 39627600242411558 (137.49451, −7.81783) | 0.81 / 0.84 | 1.0774 / 1.0790 | dark/dark | yes | 0.61 |
| J1202+0703 | 39627956565313598 (180.61815, +7.05703) / 39627956565313594 (180.61796, +7.05731) | 1.21 / 1.21 | 2.1912 / 2.1871 | dark/dark | yes | 0.77 |
| J1428+0500 | 2781075411566612 (217.23101, +5.00552) / 39627909081598677 (217.23077, +5.00609) | 2.23 / 2.20 | 1.3749 / 1.3764 | bright/dark | yes | 0.19 |
| J2244−0550 | 39627651668773421 (341.04736, −5.84732) / 39627651668773415 (341.04721, −5.84706) | 1.02 / 1.07 | 2.8438 / 2.8430 | bright/bright | yes | 0.46 |
| J1613+1708 | 39628199625230574 (243.33351, +17.14416) / 39628199625230565 (243.33332, +17.14430) | 0.78 / 0.82 | 1.5446 / 1.5475 | dark/dark | yes | 0.56 |
| J1511+3357 | 2305843018827050187 (227.79105, +33.95045) / 39632961548586204 (227.79099, +33.95081) | 1.10 / 1.31 | 0.7996 / 0.7997 | backup/dark | yes | 1.33 |
| J1552+3009 | 39628486805029505 (238.10671, +30.15088) / 39628486805029507 (238.10674, +30.15053) | 1.25 / 1.27 | 0.7519 / 0.7513 | dark/dark | yes | 0.35 |
| J0740+2926 | 39628469268646194 (115.05602, +29.44677) / 39628469268646191 (115.05596, +29.44603) | 2.64 / 2.68 | 0.9800 / 0.9779 | dark/dark | yes | 1.36 |
| J1035+0752 | 39627980183437464 (158.83006, +7.88238) / 39627980183437471 (158.83067, +7.88281) | 2.66 / 2.67 | 1.2202 / 1.2202 | dark/dark | yes | 1.17 |
| J2351+1343 | 39628119941842459 (357.78664, +13.72301) / 39628119941842445 (357.78605, +13.72296) | 1.86 / 2.05 | 1.2246 / 1.2227 | dark/dark | **no** | 3.87 |

Full names of the SQLS / Hennawi pairs:
- J1613 = SDSS J161320.01+170839.2;
- J1511 = J151109.85+335701.7;
- J1552 = J155225.62+300902.0;
- J0740 = J074013.42+292645.8;
- J1035 = J103519.22+075256.3;
- J2351 = J235108.70+134320.8.

**Identity findings:**
- **J2351+1343 is not the catalogued pair.**
  - The SQLS companion (g 22.72, z 18.74; Gaia DR3 parallax/error 3.2) has no DESI spectrum and no SQLS redshift. SLED calls the system QUASAR+STAR.
  - DESI's second fibre sits on a different, faint source (g 23.45) 2.05″ from the quasar, which is 3.9 mag brighter.
  - That fibre has no Gaia source within 0.5″; its nearest Gaia source is the quasar, at 1.94″.
  - The pair entered #472 only because its DESI separation fell within 0.5″ of the catalogued one.
- **Close pairs:** 6 pairs have DESI fibre separations below 1.5″, the fibre diameter: J1613 0.82″, J0909 0.84″,
  J2244 1.07″, J1202 1.21″, J1552 1.27″ and J1511 1.31″.
  - Their fibres overlap, so each spectrum carries light from both components.
  - J0133 (1.74″) completes the 7 at ≤ 2″ (`close_le_2arcsec`).
  - Pairs with large Δg (J1511 1.33, J0740 1.36, J1035 1.17) risk crosstalk mainly in the fainter spectrum.
- **Other DESI TARGETIDs** (`literature.json` → `jing_targetid_check`):
  - Jing et al.'s catalogue lists the same TARGETID pairs as #472 for 11 of the 13.
  - For J1428 and J1511 it lists another TARGETID (dark programme) for one component, 0.06″ from #472's. DR1 may therefore hold two coadds of that component.
- **Extra Tractor sources within 6″** (`identity.json`), not lens-galaxy detections:
  - J2244: a g-only point source (g 22.75) 1.8–2.2″ from both images;
  - J1035: a red point source (g 23.52, z 21.47) 2.1″ from both images, off the line joining them;
  - J0133 and J0728: faint REX sources 4–6″ away.

## B. Classes, evidence and what is missing
"same DESI" = the source used the DESI DR1 spectra #472 holds. Keys are defined in `literature.json` → `sources`.

| pair | original class (basis) | later classes (same DESI?) | independent evidence | reading | missing |
|---|---|---|---|---|---|
| J0130+0725 | Lemon2023 UQP. NOT spectra differ through absorption in the fainter image; text: "likely a projected pair"; one source not ruled out | Chen2026 dual quasar, lens "cannot be completely ruled out" (yes); Jing2025 QP (yes) | NOT and eBOSS spectra; HST F814W: no lens light to ≈ 23 (D-064) | ordinary favoured: projected / dual; lens not excluded | an independent velocity or profile difference that absorption cannot explain |
| J0133+0816 | Lemon2023 UQP. Small CIII] / MgII width differences; deep CFHT riz shows no lens | Jing2025 QPC (yes); Dawes A (no) | NOT spectra; CFHT imaging | open | high-resolution imaging; higher-S/N resolved spectra |
| J0728+2607 | Lemon2023 UQP, but the text says "large differences in ... profiles and continua ... likely a physical quasar pair" | Jing2025 QP (yes); Dawes C (no) | NOT spectra; HST F814W: no lens light (empirical PSF, D-064) | ordinary favoured: binary; nothing new to do | none for W2 |
| J0909−0749 | Lemon2023 UQP. Continua "significantly different"; PS1 shows no lens | Jing2025 **Contaminant** (yes; fibres overlap) | NTT spectra; VLASS 2.6 mJy | open; DESI pair unusable | resolved imaging and spectra at 0.8″ |
| J1202+0703 | Lemon2023 UQP. "Similar continua and profiles, except for absorption in the CIV line of the fainter trace"; no lens | Jing2025 QPC (yes); Dawes A (no) | NTT spectra; Keck LRIS 2022/2023 in KOA (unpublished) | open | high-resolution imaging (none found, RESULT §4) |
| J1428+0500 | Lemon2023 UQP. MgII differs through absorption in the fainter image; LS fits two PSFs; Lemon2019 NIQ | Jing2025 QP (yes, one component another TARGETID); Dawes A (no) | WHT spectra; no HST (D-064, 2026-10-09) | open | high-resolution imaging |
| J2244−0550 | **Lemon2023 conflict:** §4.2.16 "likely lensed" ("almost identical profiles"; no lens in PS1/LS); table NIQ / UQP | Jing2025 QPC (yes); Dawes A (no) | NTT spectra; **HST 2024 and Keck AO 2020 exist, no published result found** (RESULT §4) | open, lens-leaning in its source paper; misfiled as "rejected" (RESULT §3) | lens light in high-resolution imaging (the data exist) |
| J1613+1708 | SQLS DR7 "no lensing object": UH88 I imaging only, **no SQLS spectra** | Jing2025 QPC (yes; fibres overlap) | none beyond DESI | open; that both components are quasars rests on overlapping fibres | resolved spectra; high-resolution imaging |
| J1511+3357 | SQLS DR5 "binary QSO". Subaru spectra **similar** (line shapes and continua); "binary" because no galaxy is seen to M_I ≈ −19.9 | Jing2025 QPC (yes, one component another TARGETID); Kayo & Oguri 2012 binary (no) | Subaru spectra and I imaging | open: the label is a lens-galaxy non-detection, which a dark lens also gives | high-resolution imaging or a time delay (none found, RESULT §4) |
| J1552+3009 | SQLS DR5 "binary QSO". Subaru spectra similar, MgII of the fainter image "slightly broader"; no galaxy to M_I ≈ −19.8 | Jing2025 QPC (yes); Kayo & Oguri 2012 binary (no) | Subaru spectra and imaging; UKIRT K | open, as J1511 | high-resolution imaging or a time delay (not checked) |
| J0740+2926 | Hennawi2006 binary; SQLS "QSO pair" (cites Hennawi). APO spectra too low S/N for a spectral argument; no lens in deeper UH88 images; catalogue χ² 4.5 is a **colour** similarity | Jing2025 QP (yes); Ji2024 lens candidate (DESI EDR + SDSS; no); #472 control D = 1.04 (yes) | UH88 imaging; low-S/N APO spectra | **disputed** | higher-S/N resolved spectra not from DESI; high-resolution imaging |
| J1035+0752 | Hennawi2006 binary; SQLS "QSO pair, no lensing object" (UH88 I, H). Same basis as J0740; colour χ² 41.9 (dissimilar colours) | **Jing2025 LQC, lens candidate** (yes); #472 control D = 1.35 (yes) | UH88 imaging; low-S/N APO spectra | **disputed** | as J0740 |
| J2351+1343 | SQLS DR3 "no lensing object" (UH88 i imaging); the companion has no redshift | Jing2025 QPC; SLED QUASAR+STAR (yes) | LS: the companion is very red; Gaia parallax at 3.2σ | ordinary for the SQLS system (quasar + star); the DESI pair is a different pair (§A) | none for the SQLS system |
