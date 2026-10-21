# Line C (2026-10-21): W2 at HST resolution for blended or too-close lensed quasars — result

Plan: `PLAN.md` (frozen, commit 2a44d9b, before any target residual). Code: `scripts/w2_hst_blend.py` (classes,
coverage, select), `scripts/w2_hst_blend_chain.py` (rules), `scripts/w2_hst_blend_run.py` (download, measure,
psf, posthoc, summary, sheet); tests `tests/test_w2_hst_blend.py`. Labels: HST pixels and catalogue values
`observed`; fits, residuals, limits `derived`; injections and PSF-mismatch stamps `simulated`; lens predictions
`model_prediction`; thresholds ASSUMPTIONs (`params.json`). Numbers below are in `summary.json`,
`measurements.ecsv`, `gate.json`, `injections.ecsv`, `psf_mismatch.{ecsv,json}`, `posthoc.json`,
`classes_summary.json` or `selection.ecsv` unless marked "not stored".

**Bottom line.** The control gate passes as frozen (9 of 10 known lenses show light between the images,
efficiency 0.9, one-sided 95 % lower bound 0.61). The PSF-mismatch simulation shows the Moffat statistic can
report "light between images" with no lens (S/N 13.0 in the one control with an empirical PSF). So, by the frozen
rule, a target's "light" is "PSF mismatch not excluded", not "deflector seen". Of the 10 targets, 6 are "light
between images (PSF mismatch not excluded)", 1 is "none between images" (DESJ0405-3308), and 3 are undecided.
No system is a dark-lens candidate. The one "none" has a compact residual at the images' centroid in a post hoc
check (S/N 30.7), which an ordinary faint lens explains.

## 1. Class membership and selection (`derived`)
- The D-056 chain was rebuilt for the 444 galaxy-scale quasar and radio systems, 343 of them covered. It
  gives the same class counts as D-056: quasar blended 242, too close 65 (the 307), radio 1 and 1
  (`classes_summary.json`). Cost: about 1 min and 28 MB of Data Lab and catalogue downloads (`downloads.log`;
  not stored in the repo). The Faber–Jackson calibration did not need rebuilding.
- HST coverage comes from one MAST CAOM TAP query per filter (the TAP rejects OR-ed geometry, and the position
  services take about 6 s per object; not stored). 89 of the 307 lie in F814W or F160W imaging; 81 lie in F814W
  ACS/WFC or WFC3/UVIS imaging (`hst_coverage.ecsv`).
- Selection, frozen in `PLAN.md` §2:
  - 10 targets with no catalogued lens z, the deepest F814W first;
  - 10 same-programme controls with a catalogued lens z.
  None of the 20 has a catalogued θ_E, so the pair check is "unchecked" for all 20.

## 2. Controls and gate (PLAN §4.1)
| control | N img | max sep (″) | S/N | status | m_lim r_eff 0.3″ |
|---|---|---|---|---|---|
| SDSSJ0246-0825 | 3 | 2.24 | −3.5 | none between images | 21 |
| HE0047-1756 | 2 | 1.43 | 134.3 | light | 22 |
| SDSSJ0806+2006 | 2 | 1.49 | 95.4 | light | 24 |
| WG0214-2105 | 4 | 1.85 | 42.5 | light | 23 |
| SDSSJ1620+1203 | 2 | 2.82 | 27.2 | light | 21 |
| WGD2038-4008 | 4 | 2.87 | 57.7 | light | 20 |
| J1042+1641 | 4 | 1.68 | 35.4 | light | 23 |
| WISEJ0259-1635 | 4 | 1.58 | 28.6 | light | 23 |
| DESJ0408-5354 | 4 | 4.19 | 92.7 | light | 23 |
| SDSSJ1251+2935 | 3 | 1.88 | 52.9 | light | 21 |

- **Gate: passed**, n_light = 9 of n_decided = 10. The efficiency is 0.9, exactly at the threshold; the one-sided
  95 % Clopper–Pearson lower bound is 0.61.
- The miss, SDSSJ0246-0825: the image set holds three point sources (F814W 16.94, 18.09 and 20.55). The largest
  separation is 2.24″, against the LS pair's 1.10″. On the contact sheet the third source lies away from the bright
  pair, so it is probably an unrelated point source taken as an image. That enlarges the aperture, and the residual
  is dominated by over-subtracted wings (raw S/N −25.3). The fitted Moffat α sits on its lower bound (1.2). This is
  an eye reading of `contact_sheet_controls.jpg`; the rule was not changed.
- In the other nine the contact sheet shows an extended residual between the images.

## 3. Injections through the whole chain (PLAN §4.2, `simulated`)
- **Grid.** Sérsic n = 4 lens, convolved with the fitted Moffat. Recovered over the 10 controls (`summary.json`):

  | r_eff | F814W 19 | 20 | 21 | 22 | 23 | 24 | 25 |
  |---|---|---|---|---|---|---|---|
  | 0.3″ | 10 | 10 | 9 | 6 | 5 | 1 | 0 |
  | 0.6″ | 10 | 10 | 9 | 6 | 4 | 0 | 0 |

  Per-control m_lim(0.3″) runs from 20 to 24, median 22.5.
- **Blending loss.** The recovered fraction of the injected flux has a median of 0.78 (p10–p90 0.58–0.98) at
  0.3″ and 0.58 (0.43–0.85) at 0.6″. The rest falls outside the aperture or goes into the point-source fits and the
  halo correction. In the 7 decided targets the medians are 0.75 and 0.54.
- **PSF mismatch** (`psf_mismatch.{ecsv,json}`):
  - Of the 5 controls tried (2200-px cutouts), only SDSSJ0806+2006 has ≥ 3 usable Gaia stars (5). The others had
    0, 1, 1 and 2.
  - Its images were rebuilt with the empirical PSF, with no lens, and measured by the Moffat chain. The result is
    "light between images" at S/N 13.0. The same stamp made with Moffat images gives −0.3.
  - With the empirical PSF, the injected lens is still recovered to F814W 24 at both r_eff, as in the Moffat
    reference. Its S/N is lower, though: at F814W 19 and r_eff 0.3″ it is 688, against 1282 for the reference.
  - **Frozen rule triggered:** `psf_mismatch_not_excluded` = true. A false positive from PSF mismatch alone
    means the control gate cannot by itself show that lens light is what the chain measures.
- **Dimming or extinction** (`measurements.ecsv`, Δ = m_lim(0.3″) − predicted lens F814W):
  - controls: Δ_typ from 0.31 to 3.96 mag, Δ_faint from −1.45 to 2.20 mag;
  - Δ_faint < 0 means a 2-rms under-luminous ordinary lens (F814W − z = 0.6) is fainter than the injection limit.
    That holds for 3 of 10 controls (SDSSJ0246-0825, HE0047-1756, SDSSJ1251+2935).
  - No second band: F160W is outside this unit.

## 4. Targets (PLAN §4.3; the gate passed)
| target | N img | max sep (″) | S/N | frozen outcome | m_lim 0.3″ | Δ_typ / Δ_faint | post hoc S/N ÷ scaled FP |
|---|---|---|---|---|---|---|---|
| HS2209+1914 | 2 | 1.04 | 32.8 | light between images (PSF mismatch not excluded) | 22 | 0.25 / −1.51 | 0.15 |
| J1817+2729 | 4 | 1.80 | 83.9 | light between images (PSF mismatch not excluded) | 23 | 2.44 / 0.68 | 4.1 |
| J0203+1612 | 3 | 2.88 | 50.3 | light between images (PSF mismatch not excluded) | 22 | 2.46 / 0.70 | 32.9 |
| DESJ0405-3308 | 4 | 1.41 | 0.5 | none between images | 23 | 1.91 / 0.15 | — |
| J0530-3730 | 3 | 1.04 | −59.0 | undecided: fit inconclusive (over-subtraction) | — | — | — |
| J1537-3010 | 4 | 3.29 | 299.9 | light between images (PSF mismatch not excluded) | 24 | 4.75 / 2.99 | 13.9 |
| PSJ1606-2333 | 4 | 1.73 | 41.3 | light between images (PSF mismatch not excluded) | 23 | 2.35 / 0.59 | 0.76 |
| J2257+2349 | 2 | 1.69 | 171.5 | light between images (PSF mismatch not excluded) | 21 | 0.30 / −1.46 | 30.0 |
| J0146-1133 | 1 | — | — | undecided: unresolved (< 2 point sources within 3″) | — | — | — |
| DESJ0407-5006 | — | — | — | undecided: no HAP cutout ("Missing HAP files") | — | — | — |

- **Post hoc diagnostic 1** (`posthoc.json`; not a frozen rule).
  - The S/N 13.0 false positive is scaled to each system by total image flux and aperture error. This assumes the
    same fractional mismatch, which is untested across instruments and visits.
  - Four "light" targets exceed the scaled false positive: J0203+1612 ×32.9, J2257+2349 ×30.0, J1537-3010 ×13.9 and
    J1817+2729 ×4.1. On the contact sheet these four show an extended residual galaxy between the images.
  - HS2209+1914 (×0.15) and PSJ1606-2333 (×0.76) do not. HS2209+1914's "light" is the halo correction turning a
    raw S/N of −81.1 positive, with α on its bound.
  - Among controls, 8 of 9 detections exceed the scaled level. The exception is WISEJ0259-1635 (×0.56), a known
    lens. So this scaling cannot decide a system either.
- **Post hoc diagnostic 2** (`posthoc.json`). DESJ0405-3308, the only "none", has a residual of S/N 30.7 in a 0.2″
  circle at the images' centroid, with F814W 24.0 inside the circle. The contact sheet shows a compact source
  there. In the frozen aperture, the negative halo residuals of the four bright images cancel it: α is at its bound
  and the raw S/N is −3.4.
  - The ordinary reading is a faint, compact lens galaxy that the frozen statistic misses: model limits (PSF wings),
    depth (Δ_faint only 0.15 mag) and blending.
  - It is not a dark-lens candidate. The frozen route is still `/vet-candidate` (literature lens light, lens z),
    which was not run in this unit.
- **Undecided.** J0530-3730: the halos are over-subtracted; the brightest pair is 16.7 and 17.3 mag. J0146-1133: one
  point source within 3″. DESJ0407-5006: has CAOM F814W, but Hapcut has no HAP product.

## 5. What this answers, and what it does not
- **Answered: can PSF-subtracted HST models of the D-056 blended or too-close quasars reach the D-060 efficiency on
  known lenses?** Yes, but only just: 0.9 exactly, with a 95 % lower bound of 0.61, on 10 same-programme
  controls. The run also shows that the Moffat N-image statistic is not specific: PSF mismatch alone gave S/N 13 in
  the one simulated control with an empirical PSF.
  - So "light between images" does not decide a deflector, and none of the 6 "light" targets is counted as
    "deflector seen".
  - None of the 10 targets is a dark-lens candidate. The one "none" has a post hoc compact central residual.
- **Not answered.** Whether any of the 307 lacks a lens galaxy. This needs a statistic that separates lens light
  from PSF mismatch, which in turn needs the following.
  - **PSF:** a PSF per image matched to the drizzled HAP-MVM product. Reuse candidates are STScI's empirical PSF
    libraries for ACS/WFC and WFC3/UVIS (via photutils) and the `w12_niq_epsf` star stack, which here had stars
    in 1 of 5 fields.
  - **Joint model:** a joint lens-light plus point-source model, e.g. lenstronomy with PSF iteration.
  - **Cost:** both need new packages (photutils and lenstronomy are not installed). The lenstronomy cost is about
    minutes of CPU per system (not measured). Per the authorization, they were not built here.

## 6. Deviations and implementation notes
- **Deviation (noise).** The "MAD × correlation" term of the per-pixel noise uses the pixel MAD of the 250-px
  cutout, not the MAD of the fit residual (`noise_function`). The other term of the max, the empirical
  sky-aperture scatter, is as in D-064.
- **Deviation (outputs).** There are two contact sheets, `contact_sheet_controls.jpg` and
  `contact_sheet_targets.jpg`, not one `contact_sheet.jpg`. `summary.json`, `cpu_targets.json` and `posthoc.json`
  were added.
- **Post hoc.** Everything in `posthoc.json` (scaled false positive, 0.2″ centre aperture) is post hoc.
- **Order of runs.** Controls first ran without injections (gate only), then again with injections; the statuses
  were identical. The injection Sérsic is 5 × 5 subsampled and PSF-convolved, whereas D-064 used an unconvolved
  pixel-centre Sérsic. This was found during the offline tests, before any data were run: an unconvolved n = 4 cusp
  is taken as a point source.
- **Fit bounds.** The fitted Moffat α hits its lower bound (1.2) in 5 of 18 fits.
- **Downloads:** 1,935,377,411 bytes in 176 files (`summary.json`). The 2200-px wide cutouts for the PSF test took
  1,745,375,040 bytes in 36 files, because Hapcut returns every filter and the combined-skycell products. Raw
  files were deleted after use; the log is in the scratchpad (`w2-hst/downloads.log`).
- **CPU:** controls with injections 43.8 s (`gate.json`); targets with injections 27.4 s (`cpu_targets.json`). The
  PSF stage was not stored (about 6 min wall clock, mostly download).
