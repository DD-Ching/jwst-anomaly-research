# Screening of other directions for newly public observations (2026-10-24): no direction qualifies now; recommendation: none

- **Scope:** owner request 2026-10-24. Up to three questions from other project lines that newly public, or so far
  unused, observations could answer; recommend zero or one.
- **Method:** literature, product documentation and archive metadata only. No target analysis, no new tool.
- **Limits:** 45 min (00:10 UTC onward); 100 MB through `scripts/net_ledger.py`. Used: 2.6 MB (`net_ledger.jsonl`).
  Web searches ran through the search tool, outside the container, and are not metered.
- **Orientation:**
  - `main` at ee2bd8e. Open PRs: only needs-human ones, plus dependabot #442.
  - The active agent work is the D-034 reference-catalogue routes (#484–#490); none of the three items below
    overlaps with it.
  - Candidates were taken from the "Reopen only if" column of `docs/state/pending.md` and from the vetting records
    whose stop condition names new data.
  - Paused lines stay paused: J2244 image modelling (D-060 addenda), DESI classifier, X4, W3 / ATLAS. Only reopen
    conditions were checked.
- **Labels:** archive metadata: observed (`metadata/*.csv`). Depth figures: estimates (`estimates_a370_1565.json`).

## Comparison

| | A. Sunrise `n0153` (WHL0137−08) | B. Abell S1063 anti rows 1316 / 1557 | C. Abell 370 anti row 1565 |
|---|---|---|---|
| **Open question** | Is the 0.64 mag F150W brightening (2022→2025) a nuclear transient (fades back) or AGN variability (keeps varying)? Record: PR #61 `docs/candidates/n0153.md` | Do 1316 / 1557 have their own redshifts? A foreground or member redshift would remove a row from the anti count. Record: `results/arcs/abells1063-cats/vetting/RESULT.md` | Is 1565 a z ≈ 2.9 background source with an anti orientation, or an unlensed z 0.386 knot in the ionised nebula seen by MUSE? Lagattuta+2019 attribute similar z ≈ 0.38 emission (CL49) to a jellyfish tail; whether it is the same structure is not checked Record: PR #199 `docs/candidates/abell370-1565.md` |
| **Stated reopen condition** | A new public JWST epoch of the field | The full 116.28VK.001 data set, or a redshift catalogue from it | A redshift for 1565, or imaging that fixes F150W/F200W to ≤ 0.15 |
| **New data found** | None. MAST CAOM: every JWST product at the position is 2282 o010 (2022-07-30) or 6882 o052 (2025-07-08); no NIRSpec / NIRISS spectrum (`metadata/sunrise_n0153.csv`) | Nothing new since the 2026-10-16 check. ESO ObsCore: 20 / 39 products public, the latest released 2026-06-16; 19 more are released 2026-11-20 to 2027-03-12 (raw: 38 / 77 files, same dates) (`metadata/eso_116_28VK_*.csv`) | Public but unused: GO 5890 F150W `jw05890-o001_t001_nircam_clear-f150w`, 37,106 s (public 2025-12-25); F210M from 3538 o002/o013/o014/o017, 2883 o004 and 6882 o053 (30,257 s); 6882 o053 F150W (1,245 s); 3362 NIRISS F090W WFSS `jw03362-o002_t001_niriss_f090w-gr150{c,r}` (1,117 / 687 s), which covers Hα at z 0.386 (0.910 µm) (`metadata/a370_1565_*.csv`) |
| **Enough to separate ordinary from hypothesis?** | No data | Not yet | **Marginal.** Estimated ratio error now: 0.14 if r = 0.47, 0.21 if r = 1.0. After 2027-01-28 (9371 and 7345 F200W public): 0.13 / 0.16. These are optimistic white-noise scalings (`estimates_a370_1565.json`), and the threshold is met only for r = 0.47. A colour tests the photo-z, not the orientation, and gives no redshift. The F090W WFSS is short (1.8 ks over two grisms): only a compact Hα of high equivalent width would show; a non-detection says nothing (estimate, not stored) |
| **Existing method and validation** | `scripts/transient_forced.py` (on `main`) + `epoch_difference.py` (on #61's branch only); controls: 388 / 378 sources | ESO SODA cutouts + `scripts/arcs_vet_spectra.py` (on `main`, review-fixed 2026-10-16) | Forced photometry with random-aperture noise (#199). **Its scripts were scratch and not kept.** WFSS needs a grizli extraction with contamination modelling, which has never been run in this repo |
| **Done elsewhere?** | Not searched in the literature (an anonymous 27 mag source); MAST shows no spectrum | No published redshift for 1316 / 1557 at the 2026-10-16 check; not re-searched this round | No "abell370" entry in the DJA v7 imaging index (rechecked 2026-10-24, 2.5 MB; other naming not checked). CANUCS DR1 (arXiv:2506.21685) gives no spectroscopic z for 1565. Cycle 5 LUCID (GO 10341) has planned NIRSpec MSA in Abell 370 for 3 < z < 6 Lyα emitters; whether 1565 is a target is not public |
| **Minimal test, cost, stop** | — | After the releases: SODA cutouts of the new OBs, about 15 MB each, and the existing stack and limit code; about 1–2 h. Stop if the limits stay far above typical line strengths, e.g. rest EW ~100 Å (a suggested threshold, not frozen) | Cutouts (byte ranges) of 6–9 level-3 mosaics: about 0.3–0.5 GB, 2–3 h including re-implementing the forced photometry and cross-programme astrometry. Stop if the measured random-aperture error of the ratio is > 0.15, or the alignment residual is > 0.03″ |
| **Verdict** | **Not qualified:** no new data | **Not qualified yet:** data not public until 2026-11-20 | **Not qualified now:** new data exist, but the estimated precision is marginal and the outcome is likely "undetermined" again |

## Recommendation
- **None now.** No direction has new public data that is likely to settle its open question at a cost in line with
  the expected gain.
- **What is missing:**
  - A: any new JWST epoch or spectrum of WHL0137 at n0153.
  - B: the 116.28VK.001 OBs released from 2026-11-20 to 2027-03-12 (or a published catalogue).
  - C: red-side depth (the 9371 and 7345 F200W, public 2026-12-29 to 2027-01-28) or a spectrum.
    - Even then the estimated precision is 0.13–0.16, so C qualifies only if the owner accepts a colour test that
      cannot give a redshift.
- **Natural re-screen points (no schedule created):**
  - S1063 after 2027-03-12, or earlier once the OBs at 1316 / 1557 are public;
  - Abell 370 1565 after 2027-01-28;
  - Euclid DR1, staged: first part announced for November 2026, full release mid-2027 (ESA notice of 2026-06-15),
    for the W2 / W5 lines whose reopen conditions name it.

## Sources
- MAST CAOM TAP (`dbo.ObsPointing`) and ESO TAP (`ivoa.ObsCore`, `dbo.raw`), queried 2026-10-24 through the
  ledger: `metadata/`.
- JWST programme pages: GO 3362 "JWST in Technicolor" (PI Muzzin; NIRISS F090W WFSS of Abell 370, MACS0416,
  MACS1149), https://www.stsci.edu/jwst-program-info/download/jwst/pdf/3362/; GO 10341 LUCID (PI Morishita),
  https://www.stsci.edu/jwst-program-info/download/jwst/pdf/10341/.
- CANUCS / Technicolor DR1, arXiv:2506.21685. Lagattuta et al. 2019 (Abell 370 MUSE mosaic), arXiv:1904.02158.
- Euclid DR1 timeline (ESA), https://www.cosmos.esa.int/web/euclid/dr1-timeline (accessed 2026-10-24).
- DJA imaging index v7, https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html (no "abell370" string).
