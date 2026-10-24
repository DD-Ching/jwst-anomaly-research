# MARS HFF maps vs JWST-era images absent from their catalogue: MACS0416 reproduces, Abell 370 fails on rms (D-034 addendum ap)

Plan frozen in e5199c1c; functions and tests in 412de9dd; both before the run. Numbers from `summary.json` (per-image
rows in `heldout_<tag>.ecsv`) unless marked post hoc (`posthoc.json`). Labels: CANUCS and MARS positions and listed z
are observed (as the authors list them); residuals are model predictions; thresholds are ASSUMPTIONs (PLAN.md).

## Result
| Field | CANUCS images | Frame offset O (E, N) | Seen / held-out | Types E / N / M | Gate H rms | Within 3″ | H | Training K rms |
|---|---|---|---|---|---|---|---|---|
| Abell 370 | 115 | −0.033″, +0.006″ (94 pairs) | 75 / 40 | 2 / 25 / 13 | 14.81″ | 90 % | **fail** | 0.850″ |
| MACS0416 | 303 | −0.218″, +0.082″ (242 pairs) | 226 / 77 | 6 / 68 / 3 | 0.514″ | 100 % | pass | 0.513″ |

- **MACS0416: the HST-era MARS model reproduces the 77 held-out images** (34 families; rms 0.514″, every image
  within 3″; type E alone 6 images, rms 0.488″). The held-out rms equals the training rms of addendum ao.
- **Abell 370: gate H fails on rms**, with 36 / 40 within 3″. Four images are over 3″: 3b.1 (M, 3.85″), 5.3 (M,
  5.76″), 7.6 (N, 20.91″, 7 predicted images) and 24e.2 (M, z 4.916, 90.85″: the source predicts one image, so this
  is a missing-multiplicity miss, as for S1063 system 4 in addendum ao). Type E has only 2 images (rms 1.83″).
- No image fell outside the maps. Post hoc (`posthoc.json`): the 36 Abell 370 images within 3″ have rms 0.949″,
  median 0.138″; MACS0416 median 0.261″.

## Reading (narrow)
- This is the first test in this repository of a model on image positions it was not fitted to (a true held-out
  set relative to MARS, though not necessarily JWST-discovered images: CANUCS lists also carry MUSE-era systems MARS
  did not include).
- MACS0416: the free-form MARS map predicts new image positions as well as it fits its own; a validated second
  model for later extra-image / flux stages there.
- Abell 370: the four large residuals are not flags (frozen: no flag reading). Ordinary explanations not separated
  here: coarse 1.515″ map pixels near critical curves (24e.2 at z 4.9 gets one image), M-type anchors that include
  the held-out image itself, CANUCS z or family assignments, and model incompleteness away from the HST-era
  constraints. The CANUCS Abell 370 model itself failed gate K (rms 2.34″, addendum in `../abell370_canucsref/`),
  so it cannot arbitrate these four images.
- No cutouts were looked at: this stage reads catalogue positions only (no image pixels).

## Not done / limits
- MACS0717 and MACS1149: no JWST-era image list stored here; not tested.
- Mean-source and bilinear interpolation as in addendum ao; no bicubic or WLS variant run (not in the plan).

## Next
- MACS0416: MARS is a second model for a two-model extra-image comparison against the CANUCS reference route
  (`../macs0416_canucsref/`).
- Abell 370: a frozen check of 7.6 and 24e.2 against the BUFFALO model A maps (pinned, issue #465 route), which pass
  gate K2 0.910″, before any reading.
