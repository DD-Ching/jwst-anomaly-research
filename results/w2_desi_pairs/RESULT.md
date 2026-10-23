# W2 DESI pairs: the frozen spectral test flags 1 of 7 catalogued binaries (gate needs ≥ 70 %), so it is not validated; 11 of the 13 rejected pairs with DESI spectra of both components were not retrieved or read (2 read only in the control role), and none is classified

- **Plan:** `PLAN.md`, frozen at 2279b18. The tool and unit tests were committed at 7dd53f3 / d95e0ff.
- **Deviation (probe before freezing):** about 00:07 UTC, before the plan commit, one SPARCL `find` and `retrieve`
  probe read the record format. It used the two spectra of catalogued lens J0440-0905 (a positive control).
  - Only field names and array lengths were printed. No flux was examined or analysed.
  - The two spectra were retrieved again with the others at 00:40.
  - Apart from this probe, no spectrum was retrieved before the freeze.
- **Labels:**
  - DESI spectra: observed;
  - D, s, Δv and the classes: derived;
  - catalogue roles: catalogued;
  - T = 3: ASSUMPTION.
- **Scope:** a bounded trial of one route, not a search. It never tests whether a deflector is dark, and nothing here
  is a candidate.

## 1. Data
- **Sample** (`pairs.json`, PLAN §2): 106 catalogued lenses and 123 rejected pairs (D-064 pinned inputs), plus 27
  Hennawi et al. 2006 binaries with θ ≤ 10″ and |Δv| < 2000 km/s.
- **DESI DR1 pairs** (zpix via NOIRLab Data Lab TAP; spectra via SPARCL, anonymous, DESI-DR1):
  - 9 lenses, 7 binaries and 13 rejected pairs have QSO spectra of both components, at the catalogued separation
    (within 0.5″) or positions (within 1″).
  - The earlier availability count in `availability.json` (20 rejected, 12 lenses, 17 binaries) counted any ≥ 2
    QSOs in the box, and the Hennawi entries before the |Δv| cut.
- **Retrieved** (`fetch_log.json`): only the 16 control pairs, i.e. 32 spectra.
  - Line windows are kept as text in `spectra/<role>_<name>.csv.gz` (16 files, ≤ 66 kB each); SPARCL ids are in
    `fetch_log.json`.
- **Overlap not foreseen by the plan** (`gate.json` → `overlap_rejected_with_controls`): two rejected pairs,
  J074013.42+292645.8 and J103519.22+075256.3, have the same DESI TARGETIDs as the Hennawi controls SDSSJ0740+2926 and
  SDSSJ1035+0752.
  - Their spectra were therefore retrieved and read, in the control role. D = 1.04 and 1.35: "one source" for two
    catalogued binaries.
  - The other **11 rejected pairs were not retrieved or read.**
  - A future round must remove such overlaps before the gate.

## 2. Control gate (`gate.json`, `controls.json`, `controls_sheet.png`; PLAN §4)
| | usable | D ≤ 3 ("one source") | D > 3 ("two spectra") | rule | result |
|---|---|---|---|---|---|
| G1 | 9 lenses, 7 binaries | | | ≥ 6 each | pass |
| G2 lenses | 9 | 9 (1.00; Wilson 0.70–1.00) | 0 | ≥ 80 % one source | pass |
| G3 binaries | 7 | 6 | 1 (0.14; Wilson 0.03–0.51) | ≥ 70 % two spectra | **fail** |

- **Lenses:** D = 1.02–1.30.
- **Binaries:** D = 1.04–1.60, except SDSSJ1508+3328 at 3.65 (Hβ χ²_ν 5.42).
- **The miss is per line too:** the highest single-line χ²_ν among the 6 unflagged binaries is 2.41 (2QZJ1056-0059
  C IV), below T (`controls.json`).
- **Verdict:** the frozen statistic **is not validated**. As PLAN §4 requires:
  - 11 of the 13 rejected pairs were not retrieved or read; the 2 that share TARGETIDs with controls were read only
    in the control role (§1);
  - `apply` and a rejected-role `fetch` refuse to run;
  - no threshold, line list or statistic was changed, and no second statistic was tried.
- **Separations, by role** (PLAN §4 close-pair list, added at review: `gate.json` → `close_pairs_le_2arcsec`):
  - 8 of 9 lenses are at ≤ 2″ (0.6–1.7″; J1550+0221 at 2.05″); no binary is (2.7–7.3″).
  - Separation is therefore confounded with role. Fibre crosstalk at ≤ 2″ (fibres 1.5″) may contribute to the lenses'
    low D, so the G2 pass is not clean evidence of lens identity.
  - Crosstalk is unlikely to explain the binaries' low D, given their 2.7–7.3″ separations.

## 3. What the data show (descriptive; not used for the verdict)
- **Sheet** (eye reading of `controls_sheet.png`):
  - the lens profiles overlap after scaling;
  - narrow absorption lines coincide in J1037+0018, J1550+0221 and J121645.93+352941.6;
  - among the binaries, 2QZJ1056-0059 (C IV), SDSSJ1432-0106 (C IV) and SDSSJ1508+3328 (Hβ) differ visibly, and
    the others look alike within the noise.
- **Catalogued velocity differences:** the 7 binaries' Hennawi |dV| is 10–350 km/s (`pairs.json` comments). The unit
  test showed power at an 800 km/s offset (`tests/test_w2_desi_pairs.py`). PLAN §1's assumption that binaries differ
  in profile and velocity did not hold strongly enough here for this statistic.
- The line S/N on the weaker side overlaps between the roles, so S/N as the cause of the miss is not tested.
- **Profile cross-correlation Δv** (reported only, `controls.json`; post hoc, same 16 controls):
  - lenses 0–86 km/s in absolute value;
  - binaries 170–579 km/s for 6 of 7, and 0 km/s for 2QZJ1056-0059.
  - No rule is drawn from these controls.
- **DESI pipeline redshift differences:**
  - lens J121645.93+352941.6 has z_a − z_b = −0.0054 (`controls.json`);
  - rejected pair J0130+0725 has z 1.5530 and 1.5454 (`pairs.json`), a difference of the same order, so its redshift
    difference alone decides nothing.

## 4. Review notes (none changes the gate)
- **Line S/N:** computed on each spectrum's own good pixels rather than on pixels good in both. The borderline lines
  are J121645 C III] (10.1, a lens) and SDSSJ1432 C III] (8.5; counting it gives D = 1.40 per the reviewer, not
  stored).
- **Continuum fit:** requires ≥ 10 continuum pixels with both sides present, a condition PLAN does not list. It never
  triggered.
- **Redshift source:** z and zwarn are from the first `zpix` row per TARGETID. SPARCL then keeps the coadd with the
  most good pixels, and that coadd's own redshift and zwarn are not re-checked.
- **Hβ windows:** the profile window ends just short of [O III] 4960 and the red continuum window starts just past
  [O III] 5008. This touches only SDSSJ1508+3328, the one flagged binary.

## 5. What this round adds, and what is missing
- **Added:**
  - DESI DR1 holds QSO spectra of both components for 13 rejected pairs, among them J0130+0725, J0728+2607 and
    J1428+0500 from D-064's untestable list;
  - a 16-pair control set with a reusable tool;
  - a measured failure of the simplest profile test (binaries flagged 1 / 7).
- **Not added:** no rejected pair is classified, so the W2 binary-vs-lens question for them stays open.
- **Missing for a decision** (a new pre-registered round):
  - a statistic that flags binaries with small catalogued velocity differences, e.g. binned profiles, velocity offsets
    or narrow absorbers found in one component only;
  - fresh controls not seen here, with role and separation matched and no overlap with the targets;
  - more binaries, e.g. from the DESI DR1 pair catalogue (arXiv:2505.03103), avoiding circularity with its
    spectral-similarity classes.
- **Known literature** (PLAN §0):
  - J0728+2607: strongly different NOT spectra (Lemon et al. 2023).
  - J0130+0725: nearly identical spectra; a dual quasar, with lensing not excluded (Chen et al., arXiv:2607.00093,
    which used the same DESI DR1 spectra).

## 6. Cost
- **Wall clock** 23:47 to about 01:35 UTC including review (end estimate; `cost.json`):
  - route choice 23:47–00:08;
  - tool and pairs to 00:38 (a serial pairs run timed out after 25 min and was replaced by batched box queries);
  - gate 00:43.
- **Downloads:** about 20 MB, an estimate (not stored). SPARCL JSON is about 0.47 MB per spectrum (32 kept, plus a
  few from the aborted first fetch and the probe); the VizieR tables and TAP rows are small.
- **Deviations:**
  - The availability count ran 4 concurrent Data Lab queries; scripts/CLAUDE.md asks for batched queries.
  - The first `fetch` raised a SPARCL 502. Its partial spectra were deleted, retries and per-system logging were
    added, and the 16 pairs were fetched again; no spectrum was analysed in between.
