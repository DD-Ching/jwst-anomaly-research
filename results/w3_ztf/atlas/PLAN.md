# W3 ZTF × ATLAS forced photometry — does ATLAS observe the ZTF r / g common gaps? (≤ 4 stars, access trial)

Owner brief 2026-10-22 after #460: 2 h, ≤ 500 MB, ≤ 4 targets.
- **Not in scope:** changes to the ZTF edge rules, model fitting, ZTF partition re-reads, ASAS-SN, X4 and the other
  paused units.
- **Never done:** e-mail, expert inquiries, schedule or payment changes, new accounts.
- **Credentials:** never printed, stored or committed.
- **Other users:** never look at, cancel or modify their tasks.

Session start 07:31 UTC.

## 0. Prior state (checked 07:31–07:33)
- `main` at 971b613 (#460 merged). Open PRs: #458 (MACS0416 claim), #442 (dependabot), and the `needs-human` backlog.
  None concerns ATLAS or W3 / ZTF.
- Earlier ATLAS records:
  - `../altdata/PLAN.md` / `RESULT.md`: "needs a registered account", not queried.
  - SOURCES.md: the forced-photometry server.
  - No ATLAS data or code exists in the repository.
- #460 left 74 W3 sides insufficient under the joint r + g rule (`../joint/worth_*.json`). These are the ZTF r / g
  common gaps of simulated injected events. The gap structure is real sampling; the events are simulations.

## 1. Official access (checked 07:33–07:43; texts saved in the session, token masked)
- **Shared access, documented:** NOIRLab's ANTARES devkit notebook, "Lightcurves and Photometry Search", Chapter 3
  (https://nsf-noirlab.gitlab.io/csdc/antares/devkit/notebooks/LightcurvesPhotometrySearch/), states as fetched
  2026-10-22: "You may use the API token for the Astro DataLab account on the ATLAS Forced Photometry server which is
  provided below, but we encourage you to create your own ATLAS FP account ... since users of the DataLab account will
  share the same queue".
- **No ban found:** the ATLAS server's API guide, FAQ, home and Output pages say nothing against using it. Their
  banner, dated 2026-09-09 07:45 UTC, reads "Normal operations have resumed".
- **Queue and limits** (API guide):
  - tasks are POSTed to `/queue/`;
  - throttling returns 429 with a `Retry-After` header;
  - the task URL is polled until `finishtimestamp` is set;
  - `result_url` is short-lived;
  - the example deletes the task after download.
  `queuepositions.json` and the task list are **not** called, because on a shared account they show other users'
  tasks.
- **Product** (Output page):
  - MJD is the start of the 30-s exposure (UTC, not barycentric);
  - uJy and duJy, whose sign is "always statistically meaningful", so negative fluxes are kept;
  - F (c / o), err (tphot code), chi/N, x, y, maj, min, apfit, mag5sig, Sky, Obs.
- **Modes:** difference images are the default and are measured against a template that changed at MJD ~58417 and
  ~58882 (FAQ). Reduced mode (`use_reduced`) applies tphot to the target images and gives total flux. tphot fits one
  PSF plus sky and needs care in crowded fields; ATLAS pixels are 1.86″ with a FWHM of about 2 pixels.

## 2. Targets and requests (`select_script.py.txt` → `selection.json`; fixed before any query)
- **Gaps:** the 74 joint-insufficient W3 sides, each with its interval:
  - **common gap:** from the joint r + g edge to the next r / g normal point (`U_J`);
  - **decisive window:** from the joint edge to 0.25 × the joint span. A normal point inside it would decide the
    side.
- **Rule:** per area, rank by the most such sides, then brighter ZTF refmag_r, then key. Take 2 per area.
  - **Chosen:** 601209200002983 (4991, 7 sides), 601213300003718 (4991, 6), 632215200005622 (2344, 3) and
    632211200001819 (2344, 3). That is 19 gaps in all.
  - All four pass #459's clean-aperture cut (neighbours within 16″ add < 10 % of G).
- **Requests:**
  - One task per star, **reduced mode** (total flux), `send_email` False.
  - `mjd_min` / `mjd_max` cover the star's gaps ± 60 d, so baseline on both sides is included (`selection.json`).
  - The test star is submitted first. The other three are submitted only if the test product is usable: it parses
    with the documented columns, and at least one c or o row passes the quality cuts.
- **Token:** read at run time from the official notebook page into an environment variable. It is never printed,
  logged, stored or committed.
- **Tracking:** task ids, HTTP codes, time stamps, row counts, byte counts and the result sha256 go to `tasks.json`.
  - Only tasks submitted here are polled, at ≥ 30-s intervals; 429 waits for `Retry-After`.
  - Nothing is resubmitted. A finished task of ours is deleted after its result is saved, as in the API guide.
  - Tasks still queued at 09:05 UTC are recorded for continuation; that is not "no data".
- **Stop:** an authentication refusal (401 / 403), a withdrawn token or no connection. The response code and detail are
  recorded; no other token and no workaround.

## 3. What is measured (fixed now; no injection, no truth)
- **Time:** ATLAS MJD + 15 s, converted to heliocentric MJD at the star (as `w3_ztf.exposure_times`), to compare with
  the ZTF series times.
- **Quality:** the FAQ's filters (Rest et al.):
  - duJy < 10000, err == 0;
  - 100 < x, y < 10460;
  - 1.6 < maj, min < 5;
  - −1 < apfit < −0.1;
  - mag5sig > 17, Sky > 17;
  - filter c or o only.
  Saturation is not expected at m ≈ 14.5. Blending is limited by the clean-aperture cut.
- **Brightness:**
  - relative flux = uJy / median(uJy) of the passing rows, per filter;
  - σ = duJy / that median;
  - low-S/N and negative values are kept with their errors, and a missing epoch is never zero;
  - precise = σ ≤ 0.10.
- **Per gap:**
  - ATLAS exposures in the common gap, overall and passing quality, and the passing ones that are also precise;
  - visits, where exposures closer than 0.3 d count as one visit;
  - the same counts inside the decisive window;
  - whether passing precise visits exist in the 60 d before and after the common gap (baseline).
- **Gap class:**
  - **decisive:** ≥ 1 passing precise visit in the decisive window;
  - **gap only:** such a visit in the common gap but not in the decisive window;
  - **low quality only:** observations in the gap, none passing and precise;
  - **none:** no ATLAS observation in the gap.
- **Not done:**
  - A paired injection demo needs the injected ZTF series, i.e. a ZTF re-read, which the brief excludes.
  - So whether a decisive ATLAS point would read normal or dim in a given simulated event is not evaluated.
  - No calibration to ZTF r, and no PSF or photometry code of our own.
- **Labels:**
  - observed: ATLAS photometry;
  - derived: counts;
  - simulated: the gaps come from injected events.
  Four stars are not a validation of the 40.
