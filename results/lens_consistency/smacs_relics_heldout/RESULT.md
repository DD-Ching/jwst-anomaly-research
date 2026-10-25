# RELICS HST-era SMACS 0723 models vs the JWST-era Mahler+22 images: all three fail gate H (D-034 addendum aw)

Plan frozen in `PLAN.md` (commit 33067681) before any map value was read; script and tests committed before the run.
Numbers from `summary.json` (per-image rows in `heldout_<model>.ecsv`) unless marked post hoc (`posthoc.json`).
Labels: Mahler+22 and Golubchik+22 positions are observed; spectroscopic z observed, other z are ICLv2 `best.par`
model values (ASSUMPTION); residuals are model_prediction; thresholds are ASSUMPTIONs (`../mars_hff_heldout/PLAN.md`).

**Implementation fix (before any reading, run repeated):** the first run overwrote the arcs.dat spectroscopic z of
systems 1, 2, 3, 5 and 7 with 0 (they are not in `best.par`'s `z_m_limit`), so they dropped out (3 seen, 41 held-out).
`heldout_images` now keeps a listed z > 0, as PLAN.md says; a unit test covers it. Only the second run is reported.

## Result
Frame offset (Mahler → Golubchik / RELICS frame): −0.186″ E, +0.260″ N from 19 pairs within 1.0″. Seen: the 16
images of systems 1–5; held-out: 44 images of 16 JWST-era families, all type N (no new family has a seen image, so
each source is the mean of the family's own images); none outside a map.

| Model | Gate H rms | Within 3″ | H | Images > 3″ | Seen-family rms (reported, not a gate) |
|---|---|---|---|---|---|
| glafic v2 | 2.23″ | 84 % (of 44) | **fails** (rms) | 7 (6.1, 7.1, 10.2, 13.1, 17.1, 17.3, 21.1) | 2.43″, 69 % within 3″ (16) |
| zitrin-ltm-gauss v2 (Golubchik+22) | 2.73″ | 73 % | **fails** | 12 | 2.68″, 69 % |
| lenstool v1 | 7.28″ | 34 % | **fails** (not interpretable, below) | 29 | 9.27″, 19 % |

- glafic and zitrin miss their **held-out** images by about as much as their **seen** HST families (held-out rms
  2.23″ / 2.73″ vs seen 2.43″ / 2.68″). Under this estimator neither model is at the 1″ level even on its own
  families; for zitrin, Mahler+22 (§5.5; PDF text read 2026-10-25) quote Golubchik+22's own rms as about 2″. So the gate H fail reflects the models'
  precision, not a specific JWST-era family they cannot make.
- Images > 3″ in both glafic and zitrin: 6.1, 7.1, 17.1 (from `over_3arcsec`). They are not flags (frozen: no flag
  reading); system 7 is the z 5.17 family in the west, where Mahler+22 add a diffuse component.

## Post hoc (not in the plan): map readings on the seen families (`posthoc.json`)
- Sign readings (x, y, both flipped) make the seen rms of every model worse; the frozen reading (A) is the best for glafic
  and zitrin.
- lenstool v1: every sign reading gives ≥ 8.88″. Its |α| at the 17 seen positions is 1.14–1.56 × glafic's (median
  1.25; session check, not stored). Scaled by 0.8 its seen rms is 2.29″, 75 % within 3″ (reading `A_x0.8`), like
  the other two. So the lenstool v1 map as read (arcsec, D_LS/D_S = 1, the same z_lens and cosmology) is probably
  on a different normalisation. Its held-out numbers above are therefore not a model test. Its gate H was not
  re-run at 0.8: that scale was chosen after the run.

## Reading (narrow)
- No SMACS 0723 model tested so far is validated on images it was not fitted to: the three RELICS models fail
  gate H here, Chow+2023 fails even gate K (`../smacs_chow23/`), and ICLv2 is checked only on its own constraints
  (`../smacs_djaref/`). So there is no held-out-validated
  second model on SMACS 0723 and no two-model extra-image stage there. No flag.
- Caveat (PLAN.md): the glafic and lenstool constraint lists are not public (RELICS-Lenstool used seven systems,
  Mahler+22 §5.6), so some "held-out" images may have been their constraints. That would make a pass optimistic;
  it cannot explain these fails.
