# S1 burst twins on CHIME/FRB Catalog 2: result of the ≤ 30-event validation (2026-10-21)

Frozen plan: `PLAN.md` (commit 2861381, before any sample dynamic spectrum was read). This is tool validation on
30 events, **not a search**. Nothing is generalised to the catalogue or the sky, and no limit is quoted.

**Verdict (frozen rule `s1_chime_lib.verdict`): the S1 criteria do not carry over to CHIME as frozen.** T1 fails:
recovery is 6.9 % at flux ratio 1 for the G_cal set, and ≥ 50 % was required. T2 holds (0 of 133 different-source
pairs pass) and T3 holds (median injected ρ 0.934 > the 0.95 chance quantile, 0.768). T4: 0 of 3 same-source pairs
pass. Source: `summary.json` (`verdict`, `verdict_inputs`).

**Re-run after /code-review (2026-10-21, implementation fix, no threshold changed).** `chi2_union` kept the sub-bin
δ with the highest p even when that δ's gains failed, so a δ that passed could be discarded. It now takes the best
p among the δ with detected gains (`select_delta`, unit test of the failure case). The v1 outputs are kept as
`summary_v1.json`, `pairs_v1.ecsv`, `null_v1.json`, `injections_v1.ecsv.gz`, `posthoc_v1.json` and
`posthoc_injections_v1.ecsv.gz`.
- Frozen chain: no pass changed among the 4,680 injections or the 136 real pairs, so every verdict input and the
  verdict are unchanged.
- Post hoc I3: 3 injections now pass. G0 at ratio 1 goes from 0.695 to 0.697, G_cal at ratio 0.5 from 0.513 to
  0.515, and G_beam at ratio 0.3 from 0.303 to 0.305.
- One real pair (FRB20190111A/FRB20200207A, ρ 0.703, not flagged) now reports a different p; it still fails.
- The numbers below are from the re-run.

New relative to issue #153 (prior work: route, format, 307-burst pilot): every number below. #153 did not run
D-071 eligibility, pulse-shuffled ρ\*, the union-window χ², injections through the whole chain, or same-repeater
pairs.

## Access (observed; `access.json`)
- The #153 route still works: listing HTTP 200 (4,536 files, 56.4 GB in total, median 11.8 MB), catalogue CSV
  HTTP 200 (sha256 matches the E1 pin), and one format-check file HTTP 200 after one redirect to
  `ws-uv.canfar.net/minoc`. No other entry was tried.
- Format: 16384 channels × `num_time` bins, at 0.98304 ms and 24.4 kHz, incoherently dedispersed and mean-normalised.
  It comes with the fitburst `model`, `flag`, `good_freq` and `calibration/*`. In the sample, `num_time` is 162 in
  24 files, and 62, 82, 102 or 326 in the others (`profiles.json`), i.e. windows of 61–320 ms.
- Sample: 30 files, 70 transfers, 881,197,121 bytes; with the CSV, the format check and the listing, 899,106,214
  bytes. Every repeat download had the same sha256. Pins: `data/manifests/s1_chime_dynspec.ecsv`. Raw files deleted.

## Reduction (derived; `profiles.json`)
- 29 of 30 files reduce. FRB20200314A has no `calibration` group, so all its bands are missing (ineligible).
- **Deviation (implementation fixes before any outcome; both earlier versions kept):**
  - v1 (`profiles_v1_unweighted.json.gz`) took the unweighted mean of fcf × data. A few channels with a huge fcf
    dominated the noise.
  - v2 (`profiles_v2_ivar.json.gz`) used 1/(fcf σ)² weights. A few channels with a tiny fcf then dominated, and the
    band S/N of FRB20181222A fell to ~2.
  - v3 (used) averages channels in S/N units and sets the band Jy scale to the median of fcf × σ.
- Two Dec-2018 files (FRB20181222A, FRB20181226B) have band σ of 362–1,064 Jy, against a sample median of 0.176 Jy
  (`summary.json` `band_sigma_median_jy`). **Deviation:** they, and FRB20200314A, are neither template sources nor
  hosts in the injections. They stay in the real-pair tests, which are scale-free.

## Frozen chain on the sample (`summary.json`, `null.json`, `pairs.ecsv`)
- Eligible (≥ 2 pulses at 5σ): 17 of 30. Repeaters: 5 of 15 (FRB 20220912A 3, FRB 20180916B 1, FRB 20190208A 1).
  Non-repeaters: 12 of 15. Per-burst flags are in `summary.json` `per_burst`.
- ρ\* = 0.870 (0.95 quantile of the per-catalogue maximum over 100 pulse-shuffled surrogate catalogues).
- **(ii) Chance, different sources:** 133 pairs. ρ median 0.485, max 0.882; 1 flagged; 0 pass χ² at any ρ; 0 pass the
  chain; 12 untestable (the union window leaves one file's window).
  - The flagged pair, FRB20221215H/FRB20200207A (ρ 0.882), fails χ² with p = 2.6 × 10⁻⁸³. On the contact sheet both
    are one dominant narrow pulse (look-alikes, as in #153).
- **(iii) Same source, different events:** 3 pairs, all FRB 20220912A. ρ median 0.454, max 0.485; 0 flagged; 0 pass
  χ²; 1 untestable.
- **(i) Recovery (simulated), frozen design I1:** A' is the real burst A, and B' is A's fitburst model copied into a
  host's real residual. 390 injections per set and ratio (`injections.ecsv.gz`).

  | set | ratio 1 passed | ratio 0.5 | ratio 0.3 | ratio 1: copy eligible / flagged / χ² pass given flag / median ρ |
  |---|---|---|---|---|
  | G0 | 0.069 | 0.064 | 0.021 | 0.872 / 0.738 / 0.094 / 0.940 |
  | G_cal | 0.069 | 0.059 | 0.031 | 0.877 / 0.744 / 0.093 / 0.934 |
  | G_beam | 0.067 | 0.062 | 0.031 | 0.877 / 0.751 / 0.089 / 0.944 |
  | P_path | 0.003 | 0.003 | 0.005 | 0.628 / 0.290 / 0.009 / 0.880 |

  - The lag is recovered within 2 bins in 99.7–100 % of ratio-1 copies for G0, G_cal and G_beam, and in 96.4 % for
    P_path.
  - G_cal gain size: measured calibrator-shape differences, RMS 0.179 per band between file pairs
    (`cal_shape_pair_diff_rms_per_band`).

## Post hoc diagnosis (labelled; does not change the verdict; `posthoc.json`, `posthoc_injections.ecsv.gz`)
- **D1:** the fitburst model does not reproduce the sub-band profiles of bright bursts. Over the crop, without a
  systematic, data-vs-model χ²/dof ranges from 1.02 to 135.0 across the 17 eligible bursts. See
  `posthoc_misfit.png`.
- **I2** (D-071 template design: both members carry model_A; A in Gaussian noise, the copy in the host's real
  residual): G_cal passes 0.236 / 0.136 / 0.064 at ratio 1 / 0.5 / 0.3, with χ² pass given flag 0.345 at ratio 1.
  P_path passes 0.036 at ratio 1.
- **I3** (as I2, but the host noise is Gaussian at the host's σ): G_cal passes 0.705 / 0.515 / 0.305, with χ² pass
  given flag 0.948. G0 passes 0.697 and G_beam 0.700. **P_path passes 0.100** (χ² pass given flag 0.289).
- **P_path kernels:** τ = 3.2 / 1.4 / 0.74 / 0.42 bins at 450 / 550 / 650 / 750 MHz. The ΔDM = 300 smearing boxcar is
  0.68 bins or less in every band, so it was not applied: the P_path loss comes from scattering alone
  (`P_path_kernels`).
- **R2** (all 29 usable files, no eligibility cut): same source 26 pairs, ρ median 0.411 (q95 0.622, max 0.671),
  0 above ρ\*, 0 pass χ². Different source 380 pairs, ρ median 0.513 (q95 0.806, max 0.953), 4 above ρ\*, 0 pass
  χ². Per repeater, the ρ median is 0.224 (FRB 20180916B), 0.438 (FRB 20190208A) and 0.469 (FRB 20220912A).

## Transferability (what the numbers support)
- **Chain mechanics:**
  - On exact copies in radiometer-like noise (I3), the chain recovers 70.5 % at ratio 1.
  - Losses come from copy eligibility and ρ\* rather than χ² (0.948 pass given flag).
  - The per-sub-band beam or calibration imprint at the tested size changes I3 recovery by ≤ 1 pp (G0 0.697,
    G_cal 0.705, G_beam 0.700). The free per-band χ² gain and ρ absorb it.
- **χ² against real structure:**
  - The D-071 χ² (scaled, shifted copy, 10 % systematic) rejects copies once either member carries non-radiometer
    structure at the sub-band level (I2 0.236; I1 0.069).
  - The ordinary source of such structure is fitburst misfit (D1). For a real twin it would be line-of-sight
    scintillation and drift differences. Their size for a twin is not measured here (no real same-signal pair
    exists in the release).
- **Scattering (line of sight):** 1 ms extra at 600 MHz (ν⁻⁴) cuts even I3 recovery to 0.100. A twin seen along
  another sightline with different scattering is therefore not recoverable by an s = 1, no-re-convolution χ².
  Smaller scattering differences were not tested.
- **DM smearing:** at 24.4 kHz channels, ΔDM = 300 pc cm⁻³ smears ≤ 0.68 bins, which is negligible at 0.98 ms.
- **Time resolution:** pulses are ~1 bin wide. The χ² needed a sub-bin lag (δ grid); without it, an exact copy at a
  0.3-bin offset fails (`tests/test_s1_chime.py`).
- **Window:** 61–320 ms; 12 of 133 chance pairs are untestable. The window holds the burst only, so the test
  compares profiles; it does not measure a delay.
- **Eligibility:** 9 of 15 repeater events are single-pulse at 5σ in four 100-MHz bands, and 1 has no calibration.
  Most same-object pairs are therefore not testable by the D-071 rule.
- **Same object vs same event:** no same-source pair passes the chain (0 of 3 frozen, 0 of 26 post hoc). Same-source
  ρ is not higher than different-source ρ (R2 medians 0.411 vs 0.513). On this sample, the chain does not mistake
  repeater similarity for a twin.
- **Bottom line:**
  - The tool is **not validated** for CHIME twins that differ by propagation.
  - It cannot search the catalogue.
  - A false-pass rate cannot be quoted beyond these pairs (0 of 133 frozen, 0 of 380 post hoc).

## Not done
No all-catalogue pairing, no Fermi re-run, no energy study, no s ≠ 1 stretch, no scattering re-convolution in the
twin model (each would be a design change for the owner).
