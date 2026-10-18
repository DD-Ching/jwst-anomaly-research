# What others have done, what this repository has done, and what is still unverified (2026-10-18)

Primary sources were checked on 2026-10-18 (abstract and HTML pages; quotes were extracted by a fetch model and not
checked character by character). Nothing here assumes that others missed what we look for.

| Work | What it does | Brightening and fading? | Same source at different positions? | Public outputs |
|---|---|---|---|---|
| DeCoursey+ 2024, JADES Transient Survey (arXiv:2406.05060) | GOODS-S, ~25 arcmin², epochs 2022–2023, 7 NIRCam bands, 5σ ~29.5–30 mag; difference images and inverted differences; 79 SN candidates, 9 non-SN removed, 3 spectroscopic types | **yes** ("both brightening and fading transients") | no | positions in Table 4; partial injection limits only |
| Fudamoto+ 2024 (arXiv:2404.08045) | Abell 370 Dragon arc, 2022 vs 2023; >40 microlensed stars | yes (positive and negative residuals) | no | Table 1 |
| Palencia+ 2026 (arXiv:2604.22702) | Same arc, ~10 epochs 2022–2024, each epoch minus the median of the others; 104 events | yes (inverted residuals) | no | Zenodo 10.5281/zenodo.19439926 |
| Stein+ 2021, similarity search (**arXiv:2110.13151**) | self-supervised representations of 42M Legacy Survey images; similar-looking objects | n/a | **no** (look-alikes) | code via github.com/georgestein/ssl-legacysurvey (MIT) |
| Stein+ 2021, lens mining (**arXiv:2110.00023**; the ID the owner gave is this paper) | the same representations used to find strong-lens *systems* (similarity to known lenses + linear classifier); 1,192 new candidates | n/a | **no** (lens-system classification) | predictions in the same repository |
| This repository | per-box injection-calibrated two-epoch search with a dimming class (D-027, batches 1–4: 47 calibrated boxes (12 + 12 + 10 + 13)); matched-control checks; published-model counter-image prediction; exotic-lens simulator (D-047) | yes | no (model-free proposal failed, D-100) | all results and rules in `results/` |

## Added value: measured vs unknown
- **Measured:** per-box completeness of our own search, including vanish (no flux left) at 26–27.5 AB (e.g. batch 4
  pooled 0.85–0.88), and its false-hit fraction on unchanged controls; how often rule-level changes in bright compact
  sources are within matched-control scatter (A370 h003, A370 h000, SMACS h005: all within).
- **Not shown:** that our search finds anything the JADES, Dragon-arc or VENUS teams' searches miss. Their fields
  overlap ours (A370 Dragon arc; VENUS fields), and MACS0417 h002 was already SN Athena. No comparison of selection
  functions has been run, so "complementary coverage" is a hypothesis.
- **Not available anywhere checked:** identification of the same physical source at different positions without a
  lens model (task 1); our own test of it failed (D-100).
- **New this round:** a calibrated achromaticity test that separates an intrinsic transient from magnification by any
  compact deflector (point mass, Ellis wormhole, negative mass) for a two-epoch, three-band fade
  (`results/transients/confirm_r2/RESULT.md`). On the one usable object its power is 0.18–0.59, so it could not decide.
