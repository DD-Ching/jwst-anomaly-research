# Line C (2026-10-21): W2 at HST resolution for blended or too-close lensed quasars — result (v2)

Plan: `PLAN.md` (frozen in commit 2a44d9b, before any target residual). Code: `scripts/w2_hst_blend.py` (classes,
coverage, select), `scripts/w2_hst_blend_chain.py` (rules), `scripts/w2_hst_blend_run.py` (download, measure,
psf, verdict, posthoc, summary, sheet). Tests: `tests/test_w2_hst_blend.py`.

**Labels.**
- Data: HST pixels and catalogue values are `observed`. Fits, residuals and limits are `derived`. Injections and
  PSF-mismatch stamps are `simulated`. Lens magnitudes are `model_prediction`. Thresholds are ASSUMPTIONs
  (`params.json`).
- **[pre-registered]** marks tests frozen in `PLAN.md`.
- **[owner decision, after freezing]** marks the owner's rules of 2026-10-21, applied after the plan was frozen.
- **[post hoc]** marks diagnostics that are not in the plan. They are never used as rules.

**Sources of numbers.** Every number comes from `summary.json`, `gate.json`, `measurements.ecsv`, `injections.ecsv`,
`psf_mismatch.{ecsv,json}`, `verdict.json`, `selection_diff.ecsv`, `classes_summary.json`, `download_status.json`
or `posthoc.json`. The v1 versions carry a `_v1` suffix. Anything else is marked "not stored".

## Verdict

**Method validation FAILED** (`verdict.json`).

- (a) **Recovery on known ordinary lenses [pre-registered gate].** After the fix, 9 of 10 controls show light
  between the images. The efficiency is 0.9 (one-sided 95 % lower bound 0.61), so the gate passes as frozen.
- (b) **False positive on negative controls [pre-registered PSF-mismatch test].** Only one negative control could
  be built: a stamp with no lens light, whose images were made with an empirical PSF, for SDSSJ0806+2006. Through
  the Moffat chain it gives **S/N 13.0**, which is ≥ the frozen 5σ threshold. The same stamp made with Moffat images
  gives −0.3. Only 1 of the 5 controls tried had the ≥ 3 Gaia stars that an empirical PSF needs.
- **[owner decision, after freezing].** A negative control that reaches the detection threshold through PSF
  mismatch alone fails the validation. This is stricter than the frozen wording "PSF mismatch not excluded". As a
  result:
  - the target residuals are **not interpreted**;
  - they are reported only as stored measurements, with the status "not interpreted (method validation failed)";
  - no dark-lens candidate list is made;
  - no second PSF method was added and no threshold was changed.
- **Gap (named; the line stops here).** A per-image PSF matched to the drizzled HAP-MVM products, or a joint
  lens-light plus point-source model, is needed before this test can decide anything. See §6.

## 1. Implementation fix v1 → v2: F814W depth counted each exposure once (coordinator review)
- **The error.** v1 summed `t_exptime` over every contained level-3 CAOM row. HST level 3 holds overlapping
  products of the same exposures: HAP-MVM skycells (every visit, repeated in each overlapping skycell), HAP-SVM
  visits and pipeline association products. In all 23 systems involved, v1 over-counted the depth by factors of 1.9 to
  7.0 (`selection_diff.ecsv`).
- **The fix.** `visit_exposure` counts each distinct visit (ipppss, instrument) once, in this order of preference:
  - the visit's SVM exposure;
  - else the sum of its association products;
  - else the sum of its individual exposures.

  MVM products are never counted. The association-or-exposure step was made explicit in review round 2. Recomputed
  under the explicit rule, no depth changes for any of the 82 systems in F814W coverage (recomputed in session; not
  stored), and `selection.ecsv` is byte-identical. The reason is that every visit mixing the two kinds also has an
  SVM product, which takes precedence. The rule has unit tests. This is an
  implementation fix, not a re-tuning: the 1000 s threshold, the order and the caps are as frozen.
- **What changed** (`selection_v1.ecsv` → `selection.ecsv`, `selection_diff.ecsv`):

| | v1 | v2 |
|---|---|---|
| targets ≥ 1000 s | 10 | 7. J0203+1612, J2257+2349 and J0146-1133 drop out: each has one 674 s visit (17308), counted as 4044–4718 s in v1 |
| controls | 10 | 10. SDSSJ1620+1203, WGD2038-4008 and J1042+1641 drop out: they still pass 1000 s but fall behind the cap of 10 in the new order. SDSSJ0248+1913, 2M1310-1714 and SDSSJ1330+1810 come in |
| example depths | SDSSJ0246-0825 67232 s; HS2209+1914 28476 s | 21648 s; 14238 s |
| gate (pre-registered) | 9/10 = 0.9, passed | 9/10 = 0.9, passed |
| negative control (pre-registered) | S/N 13.0 (SDSSJ0806+2006) | S/N 13.0 (the same stamp; see below) |
| target interpretation | frozen wording: 6 "light (PSF mismatch not excluded)", 1 "none", 3 undecided | not interpreted (owner decision) |

- **Re-run.** The download, the controls with injections, the PSF test in the v2 order, the targets, the verdict,
  the summary and the sheets were all re-run.
- **Unchanged inputs.** The 7 controls and 6 targets that are in both versions give identical measurements:
  positions, amplitudes, S/N and limits (checked in session; values in `measurements_v1.ecsv` and
  `measurements.ecsv`).
- **PSF test reuse.** The v2 order tries the same first four controls as v1: SDSSJ0246-0825, HE0047-1756,
  WG0214-2105 and SDSSJ0806+2006, with star counts 0, 1, 1 and 5. The fifth try is new: SDSSJ0248+1913, with 1 star,
  downloaded in d9dab07.
  - Since review round 2, a reused row is accepted only if the control measurement it was built from equals the
    current one (`reuse_ok`). The compared inputs are the cutout file, image positions and amplitudes, Moffat γ and
    α, sigma_pix and the zero point.
  - All 5 reused systems pass (`psf_mismatch.json`: `rows_reused`, with `reuse_refused_input_mismatch` empty). The
    first four come from `psf_mismatch_v1.ecsv` with `measurements_v1.ecsv`; SDSSJ0248+1913 comes from d9dab07's
    table with the current `measurements.ecsv`.
  - The old labels were normalised: 3 rows read "none" and now read "not detected between images".
  - No number changed and no wide cutout was downloaded again.

## 2. Class membership and coverage (`derived`)
- **Class rebuild.** The D-056 chain was rebuilt for the 444 galaxy-scale quasar and radio systems (343 of them
  covered). It gives the same class counts as D-056: quasar blended 242, too close 65 (the 307), radio 1 and 1
  (`classes_summary.json`).
- **Coverage.** One MAST CAOM TAP query was made per filter. 81 of the 307 lie in F814W ACS/WFC or WFC3/UVIS
  imaging (`hst_coverage.ecsv`). None of the selected systems has a catalogued θ_E, so the pair check is
  "unchecked" for all of them.
- **Lens redshift is not lens presence.** "No catalogued lens z", the target criterion, does not mean "no lens".
  DESJ0405-3308, a target, has a published deflector (see §5).

## 3. Controls [pre-registered] (`gate.json`, `measurements.ecsv`, `contact_sheet_controls.jpg`)
| control | N img | max sep (″) | S/N | status | m_lim r_eff 0.3″ |
|---|---|---|---|---|---|
| SDSSJ0246-0825 | 3 | 2.24 | −3.5 | not detected between images | 21 |
| HE0047-1756 | 2 | 1.43 | 134.3 | light between images | 22 |
| WG0214-2105 | 4 | 1.85 | 42.5 | light between images | 23 |
| SDSSJ0806+2006 | 2 | 1.49 | 95.4 | light between images | 24 |
| SDSSJ0248+1913 | 4 | 1.68 | 8.2 | light between images | 21 |
| WISEJ0259-1635 | 4 | 1.58 | 28.6 | light between images | 23 |
| DESJ0408-5354 | 4 | 4.19 | 92.7 | light between images | 23 |
| SDSSJ1251+2935 | 3 | 1.88 | 52.9 | light between images | 21 |
| 2M1310-1714 | 2 | 5.02 | 208.5 | light between images | 22 |
| SDSSJ1330+1810 | 4 | 1.72 | 112.3 | light between images | 22 |

- **SDSSJ0246-0825.** On the sheet, the third "image" lies away from the bright pair. It is probably an
  unrelated point source (eye reading); the residual is dominated by over-subtracted wings (raw S/N −25.3).
- **SDSSJ0248+1913.** It passes at S/N 8.2 on a noisy residual.
- **2M1310-1714.** Its two point sources are 5.02″ apart, and a bright galaxy sits at the edge of the fit stamp
  (eye reading).
- **Fit bounds.** The fitted Moffat α is at its lower bound in 4 of 16 fits.
- **Gate caveat.** Since (b) shows that PSF mismatch alone can reach the threshold, the gate's 9/10 does not show
  that lens light is what the chain measures.

## 4. Injections and losses [pre-registered] (`injections.ecsv`, `summary.json`; `simulated`)

Whole-chain injections into the 10 controls (Sérsic n = 4, convolved with the fitted Moffat). Recovered of 10:

| r_eff | F814W 19 | 20 | 21 | 22 | 23 | 24 | 25 |
|---|---|---|---|---|---|---|---|
| 0.3″ | 10 | 10 | 10 | 7 | 4 | 1 | 0 |
| 0.6″ | 10 | 10 | 10 | 6 | 3 | 0 | 0 |

- **Limits.** Per-control m_lim(0.3″) runs from 21 to 24, median 22.0. In v1 it ran from 20 to 24, median 22.5.
- **Blending loss.** The recovered fraction of the injected flux has a median of 0.81 (p10–p90 0.62–0.93) at 0.3″
  and 0.63 (0.47–0.71) at 0.6″.
- **PSF mismatch** (`psf_mismatch.{ecsv,json}`).
  - The no-lens stamp gives S/N 13.0 with the empirical PSF and −0.3 with Moffat images. These are the negative
    controls of the verdict.
  - Injected lenses are recovered to F814W 24 with both PSFs. At F814W 19 and r_eff 0.3″ the S/N is 688 with the
    empirical PSF, against 1282 with Moffat images.
- **Dimming / extinction.** Δ = m_lim(0.3″) minus the predicted ordinary lens magnitude (`measurements.ecsv`).
  - Δ_typ runs from 0.29 to 4.74 mag; Δ_faint from −1.47 to 2.98 mag.
  - Δ_faint < 0 means a 2-rms under-luminous lens would be fainter than the limit. That holds for 4 of 10
    controls: SDSSJ0246-0825, HE0047-1756, SDSSJ0248+1913, SDSSJ1251+2935.
- **Second band.** None; F160W is outside this unit.

## 5. Targets: stored measurements, not interpreted [owner decision, after freezing]
| target | status | N img | max sep (″) | stored S/N | stored m_lim 0.3″ |
|---|---|---|---|---|---|
| HS2209+1914 | not interpreted (method validation failed) | 2 | 1.04 | 32.8 | 22 |
| DESJ0405-3308 | not interpreted (method validation failed) | 4 | 1.41 | 0.5 | 23 |
| DESJ0407-5006 | not interpreted (method validation failed); no HAP cutout | — | — | — | — |
| J0530-3730 | not interpreted (method validation failed) | 3 | 1.04 | −59.0 | — |
| J1537-3010 | not interpreted (method validation failed) | 4 | 3.29 | 299.9 | 24 |
| PSJ1606-2333 | not interpreted (method validation failed) | 4 | 1.73 | 41.3 | 23 |
| J1817+2729 | not interpreted (method validation failed) | 4 | 1.80 | 83.9 | 23 |

- **What is stored.** The chain's own status strings are kept, uninterpreted, in the column
  `chain_output_uninterpreted`. No target is read as having, or lacking, lens light. No dark-lens candidate list is
  made.
- **DESJ0405-3308 [post hoc literature check by the coordinator].** It has a published deflector: a photometric
  z_d ∼ 0.3 from Gilman et al. 2020a (MNRAS 491, 6077), as quoted in Nierenberg et al. 2023 (arXiv:2309.10101, §3),
  accessed 2026-10-21.
  - The chain output for it is "**not detected** between images" (v1: "none between images"). The chain did not
    detect the light of a known lens.
  - This is evidence about the tool. It says nothing about whether the lens galaxy exists.
  - It also shows that "no catalogued lens z" ≠ "no lens".
- **DESJ0407-5006.** Both requests to Hapcut, in v1 and v2, got an explicit no-data answer: astroquery
  `NoResultsWarning: Missing HAP files for input target. Cutout not performed.`.
  - v2: recorded in `download_status.json`, status "no HAP product".
  - v1: the same warning, followed by "0 files", is in the session's stdout log (scratchpad `w2-hst/download.log`;
    not stored in the repo).
  - Neither was a failed request.

## 6. Post hoc diagnostics [post hoc] (not used to interpret targets)
- `posthoc.json` (v2, controls only) and `posthoc_v1.json` (v1, controls and targets) hold two diagnostics: the
  flux-scaled PSF-mismatch false positive and a 0.2″ centre-aperture residual.
- The scaling uses one reference system for all three quantities (the system with the largest no-lens S/N). With
  one ePSF system, its numbers are unchanged from v1. For example, SDSSJ0806+2006 has a reference S/N of 13.0, and
  HE0047-1756 a scaled value of 57.2 in both versions.
- These diagnostics are stored only. **They are not used to interpret any target.** The v1 RESULT's readings of
  targets as "light" or "extended lens" are withdrawn.

## 7. Review fixes (coordinator `/code-review`, 2026-10-21)
1. **Depth double counting:** fixed (§1).
2. **Download errors:** `hapcut_with_retry` now separates "download failed" (retried, then fails loudly; status
   "download failed") from an explicit no-data answer ("no HAP product"). The status is in `download_status.json`.
3. **`measure --only`:** it now replaces only the systems it re-measures and recomputes the gate on the full
   control set. `--only` was not used in any committed run, so no result is affected.
4. **Post hoc reference:** the false-positive scaling uses one reference system for all three quantities. With one
   ePSF system the numbers are unchanged.
5. **`cmd_summary`:** an empty role (the gate-failed path) no longer crashes it. Tested.
6. **Retries:** the COUNT(*) TAP check, the 2200-px Hapcut call and the Gaia query in `cmd_psf` retry with back-off
   and fail loudly. `cmd_psf` writes partial results before any error is re-raised.

**Round 2 (coordinator `/code-review` of d9dab07).**
1. **All-NaN no-lens S/N:** it now gives null ("not measurable") instead of raising. All JSON is written with
   `allow_nan=False`. A failure while writing partial PSF results no longer masks the original error.
2. **`verdict` is idempotent:** `status` is always re-derived from `chain_output_uninterpreted`. The column is
   created or filled for rows that lack it. A passing verdict restores the chain status.
3. **`posthoc` scope:** it filters on the chain status, not the verdict relabel. It is restricted to controls on
   purpose while the method is not validated, and the scope is stated in the provenance field. Its numbers are
   unchanged; only the provenance text changed.
4. **Download status, overwrite and partial sets:** a later success overwrites an earlier "download failed". A
   cutout set counts as complete only after a recorded "ok"; a partial set is deleted and redone.
5. **Download status, timing and "not attempted":** the status is written after every system, and an attempt in
   progress is recorded as such. A system with no record, or with an interrupted attempt, is labelled "not
   attempted".
6. **Checked PSF reuse:** see §1. Labels were normalised; no number changed.
7. **`measure --no-inject`:** it now drops the injection rows of re-measured systems. It was used only once, for the
   first v1 control pass. That pass ran before any injection file existed, and the full run with injections then
   replaced its output. No committed file is affected.
8. **Labels from `verdict.json`:** the target-sheet title and the `posthoc` provenance are read from
   `verdict.json`. With the failed verdict the sheet text is the same as committed, so the sheets were not
   regenerated.
9. **Association or individual exposure:** the rule is now explicit; the selection is unchanged (§1).

- **Re-run in round 2:** selection (identical), the 250-px download (needed by `posthoc`; DESJ0407-5006 again gave an
  explicit "no HAP product"), PSF with checked reuse, verdict, posthoc and summary.
- **What changed in committed files:** labels and provenance text only.
  - `psf_mismatch.ecsv`: 3 status labels.
  - `psf_mismatch.json` and `summary.json`: `rows_reused` (was `rows_reused_from_v1`, now 5 systems) and
    `reuse_refused_input_mismatch` ([]).
  - `posthoc.json`: provenance text.
  - `measurements.ecsv`: string columns are written as plain strings, with identical values.
  - `download_status.json`: UTC times.
  - `summary.json`: the download totals.

## 8. Deviations, gaps and costs
- **Deviation (noise).** The "MAD × correlation" noise term uses the pixel MAD of the cutout, not the MAD of the fit
  residual. The empirical sky-aperture term is as in D-064.
- **Deviation (label) [owner decision].** A sub-threshold residual reads "not detected between images" in v2. The v1
  files say "none between images".
- **Deviation (rule) [owner decision, after freezing].** The method-validation rule (`method_validation`,
  `verdict.json`).
- **Deviation (outputs).**
  - Two contact sheets instead of one. In the target sheet only the stored S/N appears, marked "not interpreted".
  - Files added: `summary.json`, `cpu_targets.json`, `verdict.json`, `download_status.json`,
    `selection_diff.ecsv`, and the `_v1` copies.
  - PSF-test rows for systems with unchanged inputs were reused after the `reuse_ok` check (§1).
- **Injection Sérsic.** It is 5 × 5 subsampled and PSF-convolved; D-064 used an unconvolved pixel-centre Sérsic. This
  was fixed in the offline tests, before any data were run.
- **Gap.**
  - **PSF:** a per-image PSF matched to the HAP-MVM drizzle, such as STScI's empirical PSF grids for ACS/WFC and
    WFC3/UVIS (via photutils).
  - **Joint model:** a joint lens-light plus point-source model (lenstronomy).
  - **Success test:** either must keep no-lens negative controls below 5σ and controls at ≥ 0.9.
  - **Cost:** both are new dependencies (neither is installed), with minutes of CPU per system (not measured).
    Not built, per the authorization and the owner's decision.
- **Negative controls are scarce.** Gaia-star ePSFs found ≥ 3 stars in 1 of 5 HAP fields (88″). A library PSF or
  an off-target stellar field would be needed for more than one negative control.
- **Downloads:** 2,278,679,171 bytes in 356 files (`summary.json`), under the 4.5 GB cap.
  - The 2200-px cutouts for the PSF test took 1,978,194,240 bytes (42 files).
  - The 250-px cutouts took 197,064,000 bytes (304 files: v1, v2 and the round-2 re-download for `posthoc`).
  - Raw files were deleted after use; the log is in the scratchpad (`w2-hst/downloads.log`).
- **CPU:** v2 controls with injections 37.0 s (`gate.json`); v2 targets 22.1 s (`cpu_targets.json`). v1: 43.8 s and
  27.4 s (`gate_v1.json`, `cpu_targets_v1.json`). The PSF stage is not stored.
