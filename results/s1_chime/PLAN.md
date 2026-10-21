# S1 burst twins on CHIME/FRB Catalog 2: frozen plan for a ≤ 30-event validation (2026-10-21)

Line B of `docs/plans/2026-10-21-parallel-round.md`, owner authorization 2026-10-21 (interactive). This is a small
validation of cross-position burst pairing on CHIME dynamic spectra, **not a search**: no all-catalogue pairing, no
generalisation of the sample's rates to the sky or the catalogue, no Fermi re-run, no energy study.

Frozen at commit time, **before** any dynamic spectrum of the sample below is downloaded or looked at. Only the
catalogue table, the VOSpace listing and one format-check file outside the sample (FRB20180725A, `excluded_flag = 1`,
single component) have been read.

## Prior work (not redone)
- Issue #153 (owner account, 2026-10-11, research note): the data route, the file format and a 307-burst pilot.
  Format facts below come from #153 and arXiv:2601.09399 §3.3. The pilot's failed approaches are rules here:
  compare over the union of both emission regions; the χ² gain per band must be > 0 and detected; time reversal is
  not a fair null (pulse-shuffled surrogates, as D-071).
- D-071 (Fermi GBM S1): the statistic `burst_twins.xcorr_max`, the pulse rule `burst_twins.find_pulses`, the
  surrogate `burst_twins.pulse_shuffle`, the χ² `burst_twins.twin_chi2` (secondary here) and the ρ* rule are reused.
- #167 (M2 re-screen) and #169 (TTE, E_peak): not touched.

## Access (this round; ≤ 10 min confirmation of the #153 route)
- Listing `https://cadc-west-01.canfar.net/vault/nodes/AstroDataCitationDOI/CISTI.CANFAR/25.0066/data/dynamic_spectra/data?limit=10000`:
  HTTP 200, 3,651,493 bytes, 4,536 `<TNS>_stokesi_dynamic_spectrum.h5` entries (#153 counted 4,535), total
  56.4 GB, median 11.8 MB (`length` property). 4,536 of the table's 4,539 events have a file.
- Catalogue CSV via `.../vault/files/.../data/table/chimefrbcat2.csv`: HTTP 200, 4,057,396 bytes, sha256
  5108ada7… (matches the E1 pin `data/manifests/e1_events.ecsv`).
- Format check `.../data/dynamic_spectra/data/FRB20180725A_stokesi_dynamic_spectrum.h5`: HTTP 200 after one
  redirect, 10,200,204 bytes, sha256 18d3e40f…. The redirect went to `ws-uv.canfar.net/minoc` (503 in #153).
  `data`, `model`, `flag`, `times`: 16384 × 162 (float32, gzip); `good_freq`, `calibration/{good_freq,
  flux_conversion_factor (Jy), spectrum}`, `statistics/*`, `index_map/{freqs, times}`; attrs `res_time` 0.98304 ms,
  `res_freq` 24.414 kHz, `dm_incoherent`, `pulse_emission_region`, `burst_parameters_json`, `beam_number`. The window
  is `num_time` bins (162 here = 159 ms; `window_data.window` = 0.08 s); `num_time` may differ per file.
- All downloads are logged (URL, bytes, sha256) in the session scratchpad `s1-chime/downloads.log` and pinned in
  `data/manifests/s1_chime_dynspec.ecsv`. Files are streamed to a scratch directory, reduced and deleted.

## Sample (frozen; non-outcome criteria only) — `results/s1_chime/sample.ecsv`
Rule (`s1_chime_lib.select_sample`): per event n_comp = number of table rows (fitburst sub-bursts), snr = max
`snr_fitb`. Clean = `excluded_flag = 0` and `sidelobe_flag = 0` and n_comp ≥ 2 and a file in the listing.
- Repeaters: the 3 repeaters with the most clean events, 5 highest-snr clean events each:
  FRB 20220912A (FRB20221116F, 20230314D, 20221216B, 20221114B, 20230826A);
  FRB 20180916B (FRB20181222A, 20210712A, 20190825A, 20190604F, 20190605B);
  FRB 20190208A (FRB20200314A, 20211229A, 20211228A, 20190208A, 20211009B).
- Non-repeaters: the 15 highest-snr clean events: FRB20190423A, 20221215H, 20200116B, 20220316B, 20201010B,
  20230911B, 20200315A, 20200615C, 20220426A, 20181226B, 20190111A, 20190124C, 20200207A, 20200318B, 20230205A.
- 30 events, 373.6 MB in total (listing sizes). Ties on snr break by TNS name.

## Processing rules (frozen; code `scripts/s1_chime_lib.py`, `Params`; tests `tests/test_s1_chime.py`)
1. **Dedispersion:** none of our own; the files are incoherently dedispersed at `dm_incoherent` (fitburst DM).
2. **Downsampling:** none; native 0.98304 ms bins.
3. **Sub-bands:** 400–500, 500–600, 600–700, 700–800.2 MHz. Per band and bin, S_b(t) = mean over valid channels of
   `flux_conversion_factor` × `data` (Jy); valid = `good_freq` & `calibration/good_freq` & finite positive factor &
   `flag[f, t]`. The fitburst `model` is reduced the same way.
4. **Masking:** a band with < 10 % valid channels (or < 20 valid off-pulse bins) is missing (zero, not used).
5. **Crop and noise:** crop = union of `pulse_emission_region` ± 10 bins. Off-pulse = outside the crop; its median
   is subtracted from the data; σ_b = 1.4826 × MAD of the off-pulse data (radiometer noise; ASSUMPTION: no
   self-noise).
6. **Eligibility (D-071):** `burst_twins.find_pulses` on the cropped 4-band curve finds ≥ 2 pulses (prominence ≥ 5σ
   of sqrt(σ_peak² + σ_base²), ≥ 2 bins apart).
7. **Statistic:** ρ = `burst_twins.xcorr_max` of the two crops (global norms, zero padding, s = 1 only; band ratios
   must agree in Jy).
8. **Prescreen ρ\*:** 0.95 quantile of the per-catalogue maximum ρ over 100 pulse-shuffled surrogate catalogues
   (`burst_twins.pulse_shuffle` on every eligible crop; maximum over the different-source eligible pairs).
   Seed 20261021.
9. **χ² twin test (decides; issue-#153 fixed rule, `chi2_union`):** B = a_b × A shifted by the ρ lag + δ,
   δ ∈ {−0.5, −0.4, …, 0.5} bins (A Fourier-shifted; best p kept; 1 dof for δ). Window = union of both crops in B's
   frame, inside both full windows (else untestable = fail). One gain a_b ≥ 0 per band used in both (iterated
   weights 1/(σ_y² + a²σ_x² + (0.10 a x)²)). A band has signal if either member peaks ≥ 5σ in the window; each
   such band needs a_b/σ(a_b) ≥ 3. **Pass** = testable **and** ≥ 1 signal band **and** every signal band's gain
   detected **and** p ≥ 10⁻³.
   - Secondary, reported only: `burst_twins.twin_chi2` unchanged (one gain, 2σ bins) on the full windows.
   - Note on δ: added before any data were read, because the synthetic test of an exact copy at a 0.3-bin offset
     failed with an integer lag (FRB pulses are ~1 bin wide).
10. **Chain:** both eligible **and** ρ > ρ\* **and** χ² pass. No position cut (validation; repeaters share a position
    by construction) and no re-trigger veto (one file per event).

## Tests on the sample
- **(i) Recovery (simulated)** through the whole chain, `make_copy`: B' = r × g_b × shift(model_A) + (data_B −
  model_B), i.e. A's fitburst model copied into another sample burst's real residual noise (B's own burst removed
  by its model); A' = the real data of A. Every ordered (A eligible, host B ≠ A among the 30 reduced) once per set
  and ratio; shift = (B's region centre − A's) + U(−10, 10) bins (fractional, Fourier). Ratios r = 1, 0.5, 0.3.
  Seed 20261022. Gain sets (geometric mean of g over bands = 1):
  - **G0**: g = 1 (control).
  - **G_cal** (main; the per-sub-band mismatch from CHIME's real response): g_b = exp(l_X,b − l_Y,b) with l = ln(mean
    `calibration/spectrum` over the band) minus its band mean, for two distinct sample files X, Y. ASSUMPTION: the
    spread of measured calibrator response shapes between beams, epochs and calibrators stands in for the
    chromatic-beam imprint difference between two sky positions; it does not include the FRB's offset within its
    formed beam.
  - **G_beam**: Gaussian formed beam with FWHM ∝ 1/ν; ln g_b = −4 ln2 (u_B² − u_A²)(ν_b/600 MHz)², u ~ U(0, 0.5) per
    member (offset in units of the 600-MHz FWHM). ASSUMPTION (analytic; no CHIME beam model is used).
  - **P_path**: G_cal gains plus line-of-sight broadening of the copy: one-sided exponential τ = 1 ms (ν/600 MHz)⁻⁴
    and a boxcar for the intra-channel smearing of a DM difference of 300 pc cm⁻³ (8.3 µs × ΔDM × 0.0244 MHz ×
    ν_GHz⁻³). ASSUMPTION values.
  - Outcome per set and ratio: fraction of copies eligible, flagged (ρ > ρ\*), passed (whole chain); median ρ;
    fraction with |lag − true lag| ≤ 2 bins.
- **(ii) Chance similarity (derived):** all pairs of different sources among the eligible bursts: ρ, flag, χ² (pass
  regardless of ρ, and the whole chain).
- **(iii) Ordinary repeaters (derived):** all same-source pairs among the eligible bursts (same object, different
  events), with the same numbers, per repeater. Never read as "the same event replayed".

## Outcome rules (computable; `s1_chime_lib.verdict`; ASSUMPTION thresholds)
- Not testable if < 10 eligible bursts or no different-source pair.
- T1 recovery: passed fraction, G_cal, r = 1 ≥ 0.5.
- T2 chance: passed different-source pairs / all different-source pairs ≤ 0.01.
- T3 separation: median ρ of the G_cal r = 1 injections > the 0.95 quantile of the different-source ρ.
- **The S1 criteria carry over to CHIME on this sample** iff T1 and T2 and T3. Otherwise the failing test is named.
- T4 (reported, not part of the verdict): same-source pass fraction > 0.05 means the chain cannot tell same-object
  similarity from a twin.
- A pair that passes the whole chain is a screen statistic, not a candidate; at most it is listed for a contact-sheet
  look. Nothing is generalised beyond these 30 events.

## Outputs (`results/s1_chime/`)
`sample.ecsv`, `profiles.json` (reduced sub-band profiles, derived), `null.json`, `pairs.ecsv`,
`injections.ecsv.gz`, `summary.json`, `contact_sheet.png`, `RESULT.md`. Driver: `scripts/s1_chime.py`
(`select`, `fetch`, `run`).

## Budget
≤ 3 h for the line, hard stop 03:30 UTC; ≤ 5 GB downloads (sample 0.37 GB); one CPU-bound process under
`nice -n 10`, `OMP_NUM_THREADS=1`, ≤ 3 GB RAM.
