# MACS0417: extra-image test with an arc-aware (catalogue-segment) reference flux (result)

Plan: `PLAN.md` (frozen and committed before the run). Numbers: `summary.json`, `control_macs0417.ecsv`,
`predicted_macs0417.ecsv`. Script: `scripts/m0417_arcref.py` (frozen rules unit-tested in
`tests/test_m0417_arcref.py`).

## Frozen outcome: gate C2 fails; the segment-total reference is not valid for MACS0417; 0 flags, nothing tested
| Step | Result | Rule |
|---|---|---|
| Gate C2, all controls | 12 / 30 testable recovered or confused (0.40) | n ≥ 8 and ≥ 0.8 |
| Gate C2, arc-referenced controls | 7 / 14 (0.50) | n ≥ 3 and ≥ 0.8 |
| Control classes (50 units) | 11 recovered, 1 confused, 18 absent, 10 inconsistent_reference, 8 no_reference, 2 undetectable | — |
| Predicted extra images (6) | family 8 `inconsistent_reference`; 9 and 17 `no_reference` (\|μ\| > 50, known before freezing); family 10: 2 `ambiguous`, 1 `undetectable` | flag = `absent` with C2 passing |

- The plan's first branch applies: segment totals do not scale with the model's |μ| here, so the route is dropped
  and none of the six predictions is read as present or absent. MACS0417's extra-image test stays untested.
- `depth_mag` 27.15; pixel scale 0.0312″ from the catalogue's own centroids; frame offset as `../macs0417_frame`.

## Why the controls fail (post hoc; not frozen rules)
- Split of the 18 `absent` controls (`summary.json` `posthoc_absent_controls`): 8 have their own segment's centroid
  outside the search radius (the 1a/1b/1c/1m knots share segment 3907, 3σ semi-major 9.2″, F200W 17.62 mag; 2a.2 and
  2b.2 sit in segment 3903, 9.9″; 5.2 and 6.2 in segment 3784, 65.1″); 6 have no containing segment; 4 have
  sources in the radius but too faint (6.3 and 8.1 are predicted from references 3784, 16.35 mag, and 2656, 21.24
  mag). Sizes and magnitudes: `own_a3_arcsec`, `ref_a3_arcsec`, `ref_mag` in `control_macs0417.ecsv`.
- Looked at (contact sheet of 1a.1, 1a.3, 2a.2, 5.2, 6.3, 8.1, 8.2, 10.1, 10.2, 15.1, 16.1, 16.2 with segment
  ellipses; not stored): 1a.1 lies on the bright lensed spiral, whose knots share one segment; 8.2 sits on a compact
  galaxy; 15.1 and 16.2 fall off the F200W footprint (16.2 is `no_reference`, not `absent`); 10.2 has no
  containing segment, and 5.2 has no compact segment of its own (it lies in the 65.1″ segment 3784).
- Reading: in this cluster core the pipeline segmentation blends lensed images with cluster members and with other
  images of the same galaxy, so a segment total is not the image's flux. This is a catalogue (deblending) effect,
  not a lensing result.

## Limitations and deviations
- Deviation (plan defect, noticed after the run): rule 6 has no footprint check, so the off-footprint control 15.1
  (`off_image` in `../macs0417_frame`) is classed `absent`. Without it the all-control fraction is 12 / 29, still far below 0.8; the gate outcome does not depend on it.
- Deviation: the script was meant to be committed before the run, but that commit was stopped by a lint error, so it
  lands with the result. After the run only lint fixes and the post hoc diagnostic columns and counts were added;
  the frozen-rule functions and their tests are unchanged.
- One model (CANUCS v1 maps via the committed `../macs0417_frame` tables); no injections; no limit.

## Next
- CANUCS v2 deblended image catalogues when public, or an arc-aperture sum on member-subtracted images, frozen and
  checked on these same controls first. Not pipeline segment totals (D-034 addendum 2026-10-22).
