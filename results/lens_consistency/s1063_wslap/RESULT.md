# Abell S1063 WSLAP+ JWST 2025: input check and spec-z gate K (2026-10-24)

Next field of the fitted-model route (SPT-CL J0546 handoff). Model: WSLAP+ JWST 2025 "Case3" (Zenodo
10.5281/zenodo.17726321, CC-BY-SA-4.0), α maps at D_ds / D_s = 1. Image list: the 36 Chan+2020 images of the 12
systems with a MUSE z_spec (`abells1063_images/chan20_images.ecsv`). Plan frozen in be21669 together with the
functions and tests, before any map value was read. Script `scripts/s1063_wslap.py` (tests
`tests/test_s1063_wslap.py`); numbers from `summary.json` (frozen), `posthoc.json` (post hoc,
`scripts/s1063_wslap_posthoc.py`); per-image residuals `imageplane_zspec.ecsv` (mean source) and
`imageplane_zspec_wls.ecsv` (WLS source), provenance model_prediction.

| Check (PLAN rule) | Found | Status |
|---|---|---|
| Map header | 4200 × 4200, 0.06″ / px, RA---TAN / DEC--TAN, no BUNIT | read |
| C: coverage | 36 / 36 images inside the map | **passes** |
| U: unit at z_spec | median per-system backtrace rms: A (arcsec) 0.369″, P (pixels × 0.06″) 21.14″; threshold ≤ 1.0″ | **passes** (A) |
| K: spec-z gate, mean source | rms 4.197″, 94.4 % of 36 within 3″; rule rms ≤ 1.0″ and ≥ 80 % | **fails** (rms) |
| K-WLS (declared, not a gate) | rms 0.581″, 100 % within 3″ | reported |
| R: DJA v7 field | `abells1063` listed (35 fields) | **passes** |

**Reading.** As frozen, gate K fails on its rms condition, so the extra-image stages do not start from this stage.
The unit and normalisation are decided (reading A, the file names' arcsec at D_ds / D_s = 1). The failure is one
system: system 10 (z 0.73, all three images Chan class 2) gets a single predicted image from its mean back-traced
source (`n_predicted` 1, all three images `shared_match`; residuals 19.39″, 15.74″, 1.28″). The other 11 systems have
mean-source rms 0.10–1.10″ (`reported.per_system`). Under the pre-declared WLS source, system 10 gets 3 predicted
images (rms 1.46″) and all 12 systems lie at 0.07–1.46″.

**Comparison (same images, same z, same mean-source code):** the CATS v4.1 map gave rms 18.5″ (D-035 addendum b);
the WSLAP+ JWST 2025 map gives 4.20″, all of it from system 10 (below). This is a reproduction check of a fitted model:
WSLAP+ used BUFFALO identifications and its magnification maps include z 0.73 and 1.261 (systems 10 and 2), so the
images are not held out.

**Post hoc (added after K failed; no verdict depends on it; `posthoc.json`).**
- Without system 10, mean source: rms 0.506″, 33 / 33 within 3″ (K's thresholds would pass); WLS: 0.416″.
- System 10's single predicted image matches the near-critical failure mode seen on the CATS map (a mean of traced
  positions outside the caustic, D-035 addendum b). That it is the estimator, not the model, is a hypothesis: the
  backtrace rms of system 10 is also the largest under reading A (3.87″, `check_U.per_system`), and a different
  WSLAP+ identification of system 10 is not excluded (the BUFFALO list was not read).

**Data:** only the two α members were fetched (one range request each over 131,270,062 / 131,129,129 deflated bytes,
CRC-32 matched the zip directory, sha256 in `summary.json`), held in the session scratchpad and deleted after the run;
the other 12 members (≈ 750 MB) were not fetched. Transfer ≈ 262 MB, reason stated in the plan.

## Next
- A frozen plan that uses the WLS source as its gate, citing this K-WLS number (the Bullet precedent, D-034 addendum
  k), then gate C2 on the DJA v7 `abells1063` catalogue and the extra-image predictions. The WLS number was declared
  before the run, but making it a gate after the mean-source gate failed is a deviation and must be labelled so.
