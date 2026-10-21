# Line C (2026-10-21): W2 at HST resolution for blended or too-close lensed quasars — frozen plan

Frozen and committed before any target residual is computed or looked at. Owner authorization: interactive,
2026-10-21 (coordinator plan `docs/plans/2026-10-21-parallel-round.md`, line C). Code: `scripts/w2_hst_blend.py`
(CLI) and `scripts/w2_hst_blend_chain.py` (rules as functions, `Params`); tests `tests/test_w2_hst_blend.py`.
Every threshold below is an ASSUMPTION named in `w2_hst_blend_chain.Params`.

Question: in the D-056 lensed-quasar systems that Legacy Surveys DR10 could not decide ("blended": < 2 PSF
images within 3″; "too close": pair < 2″), does archival HST F814W show light between the quasar images (an
ordinary lens galaxy)? A "none" is an anomaly to vet, never evidence of a dark or exotic lens.

## 1. Class membership (rebuilt, `derived`)
- `w2_hst_blend.py classes`: the unchanged D-056 code path (`w12_lenscats.load_systems`, galaxy-scale cut,
  `selection_class`, `lenscats.brick_coverage` on `ls_dr10.bricks_s`, Tractor boxes, `deflector_test`), with
  Tractor queried only for the 444 galaxy-scale quasar and radio systems (343 covered). "blended" and "too close"
  depend only on each system's own Tractor sources, not on the Faber–Jackson calibration, so the calibration is
  not rebuilt and the deflector / faint galaxy / none labels of this rebuild are not used.
- Reproduced exactly: quasar blended 242, too close 65 (= the 307), radio 1 and 1
  (`classes_summary.json`). The 307 are in `classes.ecsv`.

## 2. Imaging coverage and selection (computable: `select_systems`)
- Coverage (`w2_hst_blend.py coverage`, `observed`): one MAST CAOM TAP query per filter (F814W, F160W) for all
  calibrated (level 3) HST images; matched locally (observation centre within 0.03°, then the position inside
  `s_region`). 89 of the 307 lie in an F814W or F160W image; 81 in F814W from ACS/WFC or WFC3/UVIS
  (`hst_coverage.ecsv`). Not used: the HSC v3 catalogue and LS DR10 (shown non-decisive, D-060, D-064).
- Rule, applied per system, all conditions joined by "and":
  1. D-056 quasar class, status "blended" or "too close";
  2. inside at least one calibrated F814W image from ACS/WFC or WFC3/UVIS (0.04″/px HAP-MVM drizzle);
  3. summed exposure of those F814W images ≥ 1000 s (depth).
- Role from the existing model information only:
  - **target**: no catalogued lens redshift (`lens_z_known` False);
  - **control**: a catalogued lens redshift (a published deflector) **and** F814W from a programme that also
    supplied a selected target's F814W (same imaging type: instrument, depth, observing strategy).
- Order: summed F814W exposure, deepest first; ties by `system_id`. First 10 targets, then the first 10 controls.
- Nothing about how any image looks enters the rule. Result (`selection.ecsv`):

| role | systems (system_id) |
|---|---|
| control | SDSSJ0246-0825 (L00031), HE0047-1756 (L00005), SDSSJ0806+2006 (L00058), WG0214-2105 (L00023), SDSSJ1620+1203 (L00165), WGD2038-4008 (L00190), J1042+1641 (L00093), WISEJ0259-1635 (L00036), DESJ0408-5354 (L00041), SDSSJ1251+2935 (L00114) |
| target | HS2209+1914 (L00202), J1817+2729 (L00178), J0203+1612 (L00022), DESJ0405-3308 (L00038), J0530-3730 (L00049), J1537-3010 (L00155), PSJ1606-2333 (L00161), J2257+2349 (L00207), J0146-1133 (L00018), DESJ0407-5006 (L00039) |

- None of the 20 has a catalogued θ_E, so the pair-separation check (D-064, 0.5″) cannot run: all 20 are
  "unchecked" and counted as such. LS "too close" separations are not used as the expected separation
  (ground-based blends).
- Not redone: the D-064 NIQ pairs (none is in this class).

## 3. Chain (`measure`, per system, identical for controls, targets and injections)
1. Data (`observed`): MAST HAP-MVM cutout (`astroquery.mast.Hapcut`, 250 × 250 px = 10″) at the catalogue
   position; the single-skycell F814W file with most valid pixels within 3″ (`w12_niq_hst.pick_image`).
2. Point sources: local maxima (1-px Gaussian smoothing, 0.15″ minimum separation) within 3″ (D-056
   image_radius) that are ≥ 20 σ_pix above the background, ≥ 1 % of the brightest peak and point-like:
   concentration flux(r ≤ 2 px) / flux(r ≤ 5 px) ≥ 0.7 × that of the brightest peak (ASSUMPTION: the brightest
   peak of a lensed quasar is an image). Extended peaks are never images. At most 4 (brightest first).
3. Fit: N Moffat point sources with one shared shape plus a constant background on a stamp (images' extent
   + 1″) (generalises `w12_niq_hst.two_psf_fit`).
4. Halo correction: per image, the radial residual profile on the side facing away from the other images'
   centroid, out to the nearest other image, is subtracted (equals `w12_niq_hst.remove_halos` for N = 2; unit test).
5. Statistic: residual flux in the smallest circle about the images' centroid holding every image (D-056
   search circle; the pair-diameter circle for a double), outside 0.3″ image cores. Error per pixel: the larger
   of the residual MAD × drizzle-correlation factor and the scatter of same-area empty-sky apertures in the
   cutout (`w12_niq_hst.noise_correlation`, `empirical_aperture_noise`; D-064 rule). S/N = flux / (σ √n_pix).
6. Status (`classify`): < 2 point sources → "unresolved"; largest image separation < 0.8″ → "too compact for the
   core mask"; (N = 2 with an expected separation off by > 0.5″ → "pair mismatch"; never applies here);
   noise not finite → "noise not estimable"; S/N ≤ −5 → "fit inconclusive"; S/N ≥ 5 → **"light between
   images"**; else **"none between images"**. Only the last two are decided.

## 4. Order of work and gates
1. **Controls first.** Gate (`control_gate`, D-060's revisit condition): efficiency = n(light) / n(decided)
   ≥ 0.9 with ≥ 5 decided controls. If the gate fails, that is the result: target residuals are not computed.
2. **Injections (`simulated`), whole chain** (point-source finding, fit, status) on every decided system's real
   cutout: Sérsic n = 4, ellipticity 0.2, r_eff 0.3″ and 0.6″, 5 × 5 subsampled and convolved with the system's
   fitted Moffat PSF, F814W = 19, 20, …, 25, at the SIS position on the line between the images (double, D-064
   formula) or the images' centroid (N ≥ 3). Recovered (`recovered`): status "light between images", same N as
   without injection, S/N ≥ 5 above the no-injection value. Limit m_lim(r_eff) = faintest magnitude of the
   unbroken run of recoveries from the bright end (`contiguous_limit`).
   - **Blending loss:** recovered flux fraction (injected-run residual minus baseline, over injected flux) per
     system and magnitude, and m_lim against the largest image separation. Reported, no threshold.
   - **PSF mismatch (Moffat against an empirical PSF):** for up to 3 controls (selection order) whose 2200-px
     HAP cutout holds ≥ 3 usable Gaia DR3 stars (`w12_niq_epsf` star rules and `build_epsf`): simulated stamps
     with white noise at the system's σ_pix and its N images placed with the empirical PSF at the fitted
     positions and fluxes, with no lens and with the same lens grid (lens convolved with the ePSF), run through
     the Moffat chain; the same with Moffat-made images as the no-mismatch reference. Reported: false-positive
     S/N with no lens and the m_lim shift. **Rule:** if any no-lens ePSF simulation gives S/N ≥ 5, a target's
     "light between images" is worded "light between images (PSF mismatch not excluded)", not "deflector seen".
   - **Dimming / extinction:** m_lim(0.3″) against the predicted ordinary lens F814W (`predicted_f814w`: D-056
     Faber–Jackson in LS z, a = 20.41, k = 0.67, slope −10, rms 0.88, n = 605, `model_prediction`; SIS θ_E =
     half the largest HST image separation; z_l catalogued else 0.5, z_s catalogued else 2; F814W − z = 0.6,
     the faint end of D-064's 0–0.6). Δ_typ = m_lim − m_typ and Δ_faint = m_lim − (m_typ + 2 rms) state how much
     dimming or extinction beyond an ordinary lens the test still excludes. A second band (F160W) is not part of
     this unit: the chain is not validated at WFC3/IR sampling.
3. **Targets last**, only if the gate passed, through the same chain and injections.

## 5. Outcome wording (frozen)
- Gate failed: "the chain is not validated on known lenses (efficiency k/n); no target is interpreted".
- Per target (gate passed):
  - "light between images" → "light between the images in HST F814W (S/N x): decided, deflector seen (HST)"
    (or "PSF mismatch not excluded" under the rule above);
  - "none between images" → "no light between the images to F814W = m_lim (injection-calibrated, r_eff 0.3″),
    Δ_typ / Δ_faint mag below an ordinary lens". Ordinary explanations stay open: insufficient depth, a compact
    lens taken as an image, a binary or dual quasar, model limits, blending, extinction. Goes to `/vet-candidate`
    before any interpretation; never "dark lens";
  - any other status → undecided, with its reason.
- Decided systems are reported; folding them into the D-056 limits is the coordinator's decision.

## 6. Budgets and outputs
- Downloads ≤ 2.5 GB, logged (URL, bytes, sha256) in the scratchpad `w2-hst/downloads.log`; raw cutouts deleted
  after use. One CPU-bound process at a time under `nice -n 10`, `OMP_NUM_THREADS=1`. Hard stop 03:30 UTC.
- Outputs here: `measurements.ecsv` (controls and targets), `injections.ecsv`, `psf_mismatch.ecsv`,
  `gate.json`, `params.json`, `contact_sheet.jpg`, `RESULT.md`. Deviations are labelled "deviation" or "post hoc".
