# `arcs`: anti null with an empirical CEERS shape prior (exploratory)

Plan: `PLAN.md` (committed before the prior was read). Script: `scripts/arcs_shape_prior.py`; rule functions
unit-tested before the run (`tests/test_arcs_shape_prior.py`, including a check that a Gaussian-drawn prior
reproduces `lens_consistency.anti_probability` to 0.01 at g 0.2, 0.5, 1.5). Every number below is in
`summary.json`. Input: CEERS `jw01345-o001_t021_nircam_clear-f200w_cat.ecsv` (jwst 3.0.0) from `stpubdata` S3,
pinned by sha256 in `data/manifests/reference_sample.ecsv`, deleted after reading. Shapes `observed`; prior,
expectations and labels `derived`; thresholds ASSUMPTIONs (`summary.json` `params`). Figure `prior.png`
(`scripts/arcs_shape_prior_sheet.py`; looked at).

## Frozen outcome: **excess persists under the empirical prior** (gate passed)
| Prior | n | rms per component | Weak B: 9 anti vs | Strong B: 6 anti vs | Strong A: 10 vs | S1063: 4 vs |
|---|---|---|---|---|---|---|
| Gaussian σ_ε 0.26 (as run; `../pooled-tally/summary.json`, weak from `../null-calibration/summary.json`) | — | 0.26 | 9.14 | 2.273, P 0.029 | 3.222, P 0.0019 | 0.52, P 0.0021 |
| CEERS, `arcs` cuts (frozen) | 1108 | 0.233 | 7.42, P 0.65 (gate) | **1.392, P 0.0031** | 1.934, P 3.5 × 10⁻⁵ | 0.341, P 4.3 × 10⁻⁴ |
| CEERS, semimajor ≥ 4 px (post hoc) | 661 | 0.247 | 8.26, P 0.89 | 1.884, P 0.013 | — | — |
| CEERS, semimajor ≥ 6 px (post hoc) | 359 | 0.266 | 8.87, P 1.0 | 2.352, P 0.033 | — | — |

Gate (rule 3): weak-shear pool B two-sided P 0.654 ≥ 0.05 and P(anti | g = 0) 0.3317 (1/3 ± 0.01): passes.
Rule 4: pooled B P 0.0031 < 0.05.

## The 3 open single-component rows (descriptive, not a label)
| Row | g_min | e (1 − b/a) | P(anti) CEERS / Gaussian | Prior fraction with |ε| ≥ the row's |
|---|---|---|---|---|
| WHL0137 2194 | 0.24 | 0.795 | 0.027 / 0.048 | 0.016 |
| RXC J0600 581 | 0.39 | 0.688 | 0.0043 / 0.011 | 0.064 |
| MACS0308 2157 | 0.64 | 0.531 | 0.00057 / 0.00084 | 0.301 |

## Reading (narrowest statement)
- The CEERS F200W pipeline shapes under the `arcs` cuts have a thinner high-|ε| tail than the Gaussian σ_ε 0.26 null
  (fraction with e ≥ 0.5: 0.371 vs 0.437), so the empirical prior lowers the expected anti counts. A thin Gaussian
  tail does not produce the pooled anti excess or the S1063 count; the Gaussian null is, if anything, conservative
  against them for this sample.
- Post hoc, larger prior sources are more elongated (rms 0.247 at ≥ 4 px, 0.266 at ≥ 6 px; the strong rows' median
  semimajor is 6.1 px, not stored); at ≥ 6 px the expectation comes back to the Gaussian one (2.352 vs 2.273) and
  the excess stays at P 0.033. Size matching is approximate: the strong rows are magnified.
- This tests shape noise only. Addendum 9's clump-aware expectation (5.538, P 0.484; `../pooled-census/summary.json`) still explains the pooled count;
  the 3 single-component rows stay open. No candidate; not a lensing anomaly claim.

## Limitations
- One blank pointing (CEERS t021, 1108 sources); field-to-field variance of the prior not measured. Observed
  (PSF-convolved, noisy) shapes are used as the intrinsic prior and lensed again (ASSUMPTION); no redshift or
  magnitude matching to the lensed rows. Uniform phases (blank-field shear and PSF anisotropy neglected).

## Next
The deciding models' local shear error at the 3 rows (model-to-model spread of the predicted axis and |g|);
membership only with a spectrum.
