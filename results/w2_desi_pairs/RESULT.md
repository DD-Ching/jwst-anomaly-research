# W2 DESI pairs: the frozen spectral test does not tell binaries from lenses at DESI S/N (control gate failed: binaries flagged 1 / 7); the 13 rejected pairs with DESI spectra of both components were not read

- **Plan:** `PLAN.md`, frozen at 2279b18 before any spectrum was retrieved. The tool and its unit tests were committed
  at 7dd53f3 / d95e0ff, before any spectrum was retrieved.
- **Labels:**
  - DESI spectra: observed;
  - D, s, Δv and the classes: derived;
  - catalogue roles: catalogued;
  - T = 3: ASSUMPTION.
- **Scope:**
  - A bounded trial of one route, not a search.
  - It never tests whether a deflector is dark.
  - Nothing here is a candidate.

## 1. Data
- **Sample** (`pairs.json`, PLAN §2): 106 catalogued lenses, 123 rejected pairs (D-064 pinned inputs) and 27 Hennawi et
  al. 2006 binaries with θ ≤ 10″ and |Δv| < 2000 km/s.
- **DESI DR1 pairs:**
  - 9 lenses, 7 binaries and 13 rejected pairs have QSO spectra of both components at the catalogued separation
    (within 0.5″) or positions (within 1″);
  - zpix via NOIRLab Data Lab TAP; spectra via SPARCL, anonymous, DESI-DR1.
  - The availability count in `availability.json` (20 rejected, 12 lenses, 17 binaries) was made before the pair
    rule. It counted any ≥ 2 QSOs in the box, and the Hennawi count was before the |Δv| cut.
- **Retrieved** (`fetch_log.json`): only the 16 control pairs.
  - 32 spectra; the line windows are kept in `spectra/` (610 kB in total, with SPARCL ids).
  - The rejected pairs' spectra were **not retrieved**.

## 2. Control gate (`gate.json`, `controls.json`, `controls_sheet.png`; PLAN §4)
| | usable | D ≤ 3 ("one source") | D > 3 ("two spectra") | rule | result |
|---|---|---|---|---|---|
| G1 | 9 lenses, 7 binaries | | | ≥ 6 each | pass |
| G2 lenses | 9 | 9 (1.00; Wilson 0.70–1.00) | 0 | ≥ 80 % one source | pass |
| G3 binaries | 7 | 6 | 1 (0.14; Wilson 0.03–0.51) | ≥ 70 % two spectra | **fail** |

- **Lenses:** D = 1.02–1.30.
- **Binaries:** D = 1.04–1.60, except SDSSJ1508+3328 at 3.65 (Hβ χ²_ν 5.42).
- **Verdict:** the frozen statistic **is not validated**. As PLAN §4 requires:
  - the rejected pairs' spectra were not retrieved or read;
  - no threshold, line list or statistic was changed;
  - no second statistic was tried.
- **The failure is not fibre crosstalk.** The binaries are 2.7–7.3″ apart (DESI fibres are 1.5″); the lenses are
  0.6–2.1″.
- **Two of the 13 rejected pairs are themselves Hennawi binaries.** SDSSJ0740+2926 and SDSSJ1035+0752 sit in both
  lists. As negative controls their D is 1.04 and 1.35, so the test would have called these confirmed binaries "one
  source".

## 3. What the sheet shows (eye reading of `controls_sheet.png`; not used for the verdict)
- **The tool works as intended:**
  - the lens profiles overlap after scaling;
  - narrow absorption lines coincide in J1037+0018, J1550+0221 and J121645.93+352941.6.
- **Binaries:**
  - 2QZJ1056-0059 (CIV), SDSSJ1432-0106 (CIV) and SDSSJ1508+3328 (Hβ) differ visibly; only the last exceeds T;
  - the others look alike within the noise.
- **Reading:** a mean per-pixel χ² dilutes a broad, moderate shape difference when one spectrum is noisy. At DESI S/N
  this statistic lacks power. A causal reading beyond this is not tested.

## 4. Post hoc observations (not used; any of them would be a new pre-registered method with fresh controls)
- **Profile cross-correlation Δv** (reported only, `controls.json`):
  - lenses: 0 to 86 km/s in absolute value;
  - binaries: 170–579 km/s for 6 of 7, and 0 km/s for 2QZJ1056-0059.
  - These are the same 16 controls, so this cannot validate a Δv rule.
- **Catalogue redshift differences** are not a reliable separator: the lens J121645.93+352941.6 has DESI
  z_a − z_b = −0.0054 (`controls.json`). The rejected pair J0130+0725 has DESI z 1.5530 and 1.5454 (`pairs.json`),
  a difference of the same order, so its redshift difference alone decides nothing.

## 5. What this round adds, and what is missing
- **Added:**
  - DESI DR1 holds QSO spectra of both components for 13 rejected pairs, including J0130+0725, J0728+2607 and
    J1428+0500 from D-064's untestable list;
  - the 16-pair control set (9 lenses, 7 binaries) and the tool to read them;
  - a measured failure of the simplest profile test.
- **Not added:** no rejected pair is classified, so the W2 binary-vs-lens question for them stays open.
- **Missing for a decision**, as a new pre-registered round:
  - a statistic with power on broad, moderate differences (e.g. binned profiles, velocity offsets, or narrow
    absorbers present in only one component);
  - fresh controls not seen here, e.g. more binaries from the DESI DR1 pair catalogue (arXiv:2505.03103). That
    catalogue's own classes rest partly on spectral similarity, so a control set must avoid circularity.
- **Known literature** (PLAN §0): J0728+2607 has strongly different NOT spectra (Lemon et al. 2023). J0130+0725's
  nearly identical spectra leave lensing open (arXiv:2607.00093).

## 6. Cost
- Wall clock 23:47–about 01:30 UTC including review: route choice 23:47–00:08; tool and pairs to 00:38 (a serial
  pairs run timed out after 25 min and was replaced by batched queries); gate 00:43.
- Downloads: about 20 MB, an estimate (not stored): SPARCL JSON is about 0.47 MB per spectrum (32 kept, plus a few in
  the aborted first fetch and a one-spectrum probe), and the VizieR tables and TAP rows are small.
- At most 4 concurrent requests in the availability count; otherwise serial.
- **Deviation:** the first `fetch` raised a SPARCL 502. The partial spectra were deleted, retries and per-system
  logging were added, and the 16 pairs were fetched again; no spectrum was analysed in between.
