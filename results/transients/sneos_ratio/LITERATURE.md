# SN Eos 101.2 / 101.1 flux ratio: literature check (2026-10-25)

The check named in `PLAN.md` ("whether Allingham et al. 2026 or the SN Eos discovery paper already report this
ratio"), done after `RESULT.md` was frozen and merged (#541). Tool `scripts/sneos_literature.py`, tests
`tests/test_sneos_literature.py`; every number below is in `literature.json`. Published photometry is observed
(published); published |μ| are model_predictions; ratios are derived.

## What the papers report
- Neither paper states a 101.2 / 101.1 flux ratio or discusses microlensing of SN Eos (HTML full texts searched for
  "flux ratio", "microlens"). Coulter et al. 2026 (arXiv:2601.04156v2) Table S1 gives per-image AB magnitudes at 17
  band × epoch entries: HST PID 17504 (WFC3 F814W / F110W, 2024-03-11 to 03-18, difference against 2012 templates),
  JWST 6882 (2025-09-01) and 9493 (2025-10-08). Ratios here are computed from that table.
- Published |μ| (both quoted with a ~20 % systematic per image by the authors): Allingham et al. 2026
  (arXiv:2602.14074) Table 3, 25.79 / 27.49 → R_model 1.066; Coulter et al. Table S2, 25.46 / 29.88 → 1.174. Our
  pinned v2.1 maps give 1.162. Published model range used below: 1.066–1.174.
- Delay: our pinned map has 101.2 arriving 0.32 d **before** 101.1 (`sneos/summary.json`, Δt −0.324 d), matching
  Allingham et al. (−0.37 d; "101.2 must arrive before 101.1"). Coulter et al. Table S2 lists 101.2 0.99 ± 0.33 d
  after 101.1. The 2026-10-25 CHANGELOG line "101.2 trails by 0.32 d" was a wording error (corrected).

## Ratios from the published photometry (derived)
| Group | R = F(101.2)/F(101.1) | σ_ln (stat) | χ² / dof |
|---|---|---|---|
| JWST 2025, 8 entries | 1.103 | 0.007 | 95.0 / 7 |
| HST F814W 2024, 4 entries | 1.243 | 0.137 | 2.4 / 3 |
| HST F110W 2024, 5 entries | 1.805 | 0.039 | 5.7 / 4 |

- JWST 2025 sits inside the model range. Its entries span 1.009 (F410M) to 1.191 (F444W); χ² 95 / 7 means the
  scatter exceeds the quoted 0.01–0.05 mag errors (cause not tested).
- **HST F110W 2024 is 0.43 in ln above the top of the model range** (`offset_vs_model_range`), beyond the plan's
  0.30 floor (an ASSUMPTION), and differs from the JWST 2025 ratio by a factor 1.64 (not stored). F814W, with 3.5×
  larger errors, is within the floor.
- Flux-drop ratio 6882 → 9493 from the same table (comparable to our difference images): F277W 1.33, F356W 1.13,
  F444W 1.58 (σ_rel 0.12, 0.14, 0.11). Our 9493 pairs: F356W 1.37 ± 0.21, F444W 1.93 ± 0.41 (`summary.json`), each
  within about 1σ of these (not stored). Our combined 1.382 is above the published direct-image ratio 1.103.

## Reading (post hoc; not a candidate)
- The published model range is wider than the single R_model in `PLAN.md`; with it, our "agrees" stands.
- The 2024 F110W ratio is a flux-ratio offset in published photometry at an epoch our plan never covered.
  With a 101.1–101.2 delay of order 1 d, intrinsic SN variability cannot change the ratio between images; a ratio
  that differs between 2024 (rest-frame FUV, days after explosion) and 2025 (end of plateau) points to an effect
  that varies in time or wavelength. Ordinary candidates, untested: microlensing by intracluster stars at |μ| ≈ 25–30
  (stronger for the compact early photosphere), HST difference-photometry systematics (template, PSF, host
  residuals), and macro-model error (would not change with time). A difference between direct and flux-drop ratios
  in 2025 also hints at a time-varying ratio or host light in the direct photometry (not tested).
- Not run through `/vet-candidate`: it is a published measurement, not a screen hit of ours, and microlensing is
  the expected first explanation. An own measurement on the HST 17504 / 2012 images (frozen plan, injections) would
  test it.

## Limits
Table values transcribed by hand from the arXiv HTML (v2 of 2601.04156; 2602.14074 as served 2026-10-25); quoted
errors taken at face value; the two μ tables are model averages from different runs of the same group.
