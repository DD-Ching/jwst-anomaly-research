# W3 ZTF × ATLAS forced photometry: the documented shared Data Lab token works; 4 / 4 stars returned usable reduced-mode photometry; ATLAS has good visits in the decisive window of 11 of 19 ZTF r / g common gaps (6 seasonal gaps empty); effect on decisions not measured

- **Pre-registration:** the plan, the 4-star selection and the per-gap measures were committed at 07:35:16 (e84e4dc).
  That was before the first query (07:36:22).
- **Labels:**
  - **observed:** the ATLAS photometry and the service responses;
  - **derived:** counts and windows;
  - **simulated:** the gaps come from injected events of #459 / #460, with real ZTF sampling.
  - Nothing here is a real dimming or a candidate. This is an access and coverage trial on 4 stars, not a validation.
- **Scope:**
  - Not done: a ZTF re-read, a search group, model fitting, or a change to any ZTF rule.
  - The token was used only in memory and was never printed, stored or committed.
  - Only our own 4 tasks were polled and then deleted; no other task or queue list was viewed.

## 1. Four separate answers
| question | answer | evidence |
|---|---|---|
| Does the documentation allow shared access? | **Yes.** NOIRLab's ANTARES notebook (Chapter 3) says "You may use the API token for the Astro DataLab account". It encourages a personal account, because the queue is shared. The ATLAS server pages state no prohibition. | PLAN §1 |
| Did the queries succeed? | **Yes.** 4 / 4 POSTs gave HTTP 201 and no 429. Tasks started 0–11 min after submission and ran 2.5–10 min. Our tasks were deleted after download (HTTP 204). | `tasks.json`, `cost.json` |
| Did we get valid data? | **Yes, with a caveat.** 970–2,105 rows per star; 91–94 % pass the FAQ quality cuts. Reduced-mode (total-flux) baselines are c 14.68–15.36 and o 14.45–14.67 AB (ZTF r 14.53–14.70). Quoted errors are 0.5–0.8 %, but the robust scatter of passing rows is 1.0–2.3 % (chi/N 6–25), so the errors are about 2–4× too small. That is still ≪ the 10 % precision bar. 2–14 negative-flux rows per star are kept. | `analysis.json` `stars` |
| Does ATLAS add discriminating information? | **Not measured.** Coverage only: 11 of the 19 common gaps have a quality-passing, precise visit inside the decisive window, where a normal point would decide the edge. Whether that would decide the simulated events needs a paired injection, which needs a ZTF replay; the brief excluded that. | §2 |

## 2. Per gap (`analysis.json` `gaps`, `summary`; `lightcurves.png`)
- **Classes:**
  - 11 **decisive**: a passing, precise visit in the decisive window;
  - 2 **gap only**: such visits elsewhere in the gap;
  - 6 **none**: no ATLAS observation in the gap.
- **By star:**
  - 601209200002983 (4991): 1 / 7 decisive;
  - 601213300003718 (4991): 4 / 6;
  - 632215200005622 (2344): 3 / 3;
  - 632211200001819 (2344): 3 / 3.
- **Why the 4991 gaps stay empty:** 6 of 601209200002983's 7 gaps run from late February to May, when the field
  (RA ≈ 11°) approaches solar conjunction. ATLAS, also ground-based, has no data in 5 of them and one visit in the
  sixth. The contact sheet shows the two surveys' seasons ending together.
- **Where ATLAS fills:** the non-seasonal ZTF gaps. Decisive gaps by the r-side class of v24:
  - A (no ZTF r exposure): 9 / 13;
  - B (only flagged r exposures): 1 / 5;
  - C1 (DR24 product gap): 1 / 1.
- **Baseline:** every gap has passing precise ATLAS visits within 60 d before it and within 60 d after it.
- **What the in-gap data show:** 464 passing, precise in-gap points; median relative flux 1.002, 1st–99th percentile
  0.93–1.10.
- **Outliers, observed but not vetted or interpreted:**
  - 632215200005622, MJD 60375.64 (2024-03-11): three exposures of one visit at 1.46–1.52 × baseline, and a fourth
    with uJy = 9999 at another detector position.
  - 632211200001819, MJD 59589.67: one exposure at 0.67, with mag5sig 17.95 and Sky 18.25. That night is brighter
    and shallower, but still inside the cuts.
  - A future rule would need visit-level consistency (the 4 exposures per night) before it trusts single points.

## 3. Time, band and flux checks
- **Time:**
  - ATLAS MJD is the exposure start (UTC). Plus 15 s and the heliocentric correction, it matches the ZTF series' mid-exposure HMJD.
  - The correction is at most a few minutes against gaps of days.
- **Band:**
  - ATLAS c (cyan) and o (orange) are compared, each against its own median; neither is converted to ZTF r.
  - The W3 model is achromatic; chromatic blending is not modelled.
- **Flux:**
  - Reduced mode is total flux on the target images (tphot single-PSF fit).
  - Difference mode was not used: its flux is measured against templates that changed (FAQ), so it is not a brightness
    ratio.
- **Saturation and blending:**
  - Not expected at m ≈ 14.5. All four stars pass #459's 16″ clean-aperture cut, and ATLAS's FWHM is about 3.7″.
  - The high chi/N suggests PSF mismatch for bright stars, so treat the quoted errors with care.

## 4. Is a 40-star test worth it?
- **For:**
  - Access works.
  - The data are precise.
  - ATLAS covers 11 / 19 decisive windows here, including every non-seasonal gap of the 2344 stars.
- **Against, or not yet known:**
  - Seasonal gaps are not fillable.
  - The effect on the #459 criterion (≥ 50 % of 63 restored) is unmeasured.
  - A 4-star sample does not generalise.
- **A useful next step needs two things the brief excluded:**
  1. ATLAS for all 40 stars: about 36 more tasks on the shared queue, at about 5–10 min each if run one at a time, and
     roughly 15 MB.
  2. A paired injection, with the ZTF replay (about 1.25 GB) and an ATLAS noise rule fixed beforehand.
  This is the owner's decision.
- **At search scale (thousands of stars) the shared Data Lab queue is not appropriate.** NOIRLab asks for personal
  accounts, so a personal ATLAS account (the owner's) would be needed.

## 5. Re-running with a personal account (not done; no registration was made)
- The script reads the token from the environment variable `ATLAS_FP_TOKEN`. A personal token comes from
  `POST /api-token-auth/` with the account's username and password (API guide).
- Then run `fp_script.py.txt submit <key>` and `poll <key> <dir>`, then `analysis_script.py.txt <dir>`. The targets
  and time ranges are in `selection.json` (`requests`).

## 6. Deviations
- **Plan text:** the §1 check time in PLAN read 07:33–07:43 at the plan commit. It was corrected to 07:32–07:35, the
  actual time; the content is unchanged.
- **Added after the data were seen:**
  - `analysis_script.py.txt` gained two outputs after the test star was analysed: the robust scatter and the r-side
    class / in-gap summary. These are additions, not changed definitions.
  - The contact sheet is an addition.
- **#460 wording clarified in this PR:** its "64 / 64 truly dim" are dims in the injected model (simulated-injection
  measurements), not real astronomical events.

## 7. Cost (`cost.json`)
- 2.1 MB in all (0.85 MB ATLAS results, 1.2 MB documentation pages), against a 500 MB budget.
- 07:31–08:10 UTC.
- The raw ATLAS files are not committed; their sha256 are in `tasks.json`.
