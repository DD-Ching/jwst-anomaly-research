# MACS0417: extra-image test with an arc-aware (catalogue-segment) reference flux (result)

Plan: `PLAN.md` (frozen and committed before the run). Numbers: `summary.json`, `control_macs0417.ecsv`,
`predicted_macs0417.ecsv`. Script: `scripts/m0417_arcref.py` (frozen rules unit-tested in
`tests/test_m0417_arcref.py`).

## Frozen outcome: gate C2 fails; the segment-total reference is not valid for MACS0417; 0 flags, nothing tested
| Step | Result | Rule |
|---|---|---|
| Gate C2, all controls | 10 / 27 testable recovered or confused (0.37) | n ≥ 8 and ≥ 0.8 |
| Gate C2, arc-referenced controls | 6 / 12 (0.50) | n ≥ 3 and ≥ 0.8 |
| Control classes (50 units) | 9 recovered, 1 confused, 17 absent, 15 inconsistent_reference, 6 no_reference, 2 undetectable | — |
| Predicted extra images (6) | family 8 `inconsistent_reference`; 9 and 17 `no_reference` (\|μ\| > 50, known before freezing); family 10: 2 `ambiguous`, 1 `undetectable` | flag = `absent` with C2 passing |

- The plan's first branch applies: segment totals do not scale with the model's |μ| here, so the route is dropped
  and none of the six predictions is read as present or absent. MACS0417's extra-image test stays untested.
- `depth_mag` 27.15; pixel scale 0.0312″ from the catalogue's own centroids; frame offset as `../macs0417_frame`.
- Code review fix before this run: the first run used an undeclared 10″ centroid pre-cut in the counterpart search,
  which rule 1 does not have. That run also failed C2 (12 / 30, arc 7 / 14; not stored). The numbers here are from
  rule 1 as written (no distance cut).

## Why the controls fail (post hoc; not frozen rules)
- Split of the 17 `absent` controls (`summary.json` `posthoc_absent_controls`): 10 have their own segment's centroid
  outside the search radius, 1 has no containing segment (10.2) and 6 have sources in the radius that are too faint.
  The large segments involved: 3907 (3σ semi-major 9.2″, F200W 17.62 mag) holds the 1a/1b/1c/1m knot pairs; 3903
  (9.9″) holds 2a.2 and 2b.2; 3784 (65.1″, 16.35 mag) holds 3.1, 5.2, 6.2 and 16.3, and is the reference for 6.3, 15.3
  and 16.2; 3785 (16.0″, 16.20 mag) holds 16.2 and is the reference for 16.3 (`own_a3_arcsec`, `ref_a3_arcsec`,
  `ref_mag` in `control_macs0417.ecsv`).
- Looked at (contact sheet of 1a.1, 1a.3, 2a.2, 5.2, 6.3, 8.1, 8.2, 10.1, 10.2, 15.1, 16.1, 16.2 with segment
  ellipses; not stored): 1a.1 lies on the bright lensed spiral, whose knots share one segment; 8.2 sits on a compact
  galaxy, the reference (2656, 21.24 mag) for 8.1; 15.1 and 16.2 fall off the F200W footprint.
- Reading: in this cluster core single pipeline segments of 9.2–65.1″ (3σ; three of them 16.2–17.6 mag) hold several lensed knots
  or the image together with other light, so a segment total is not one image's flux. This is a catalogue
  (deblending) effect, not a lensing result. Cluster membership of these segments is not established here.

## Limitations and deviations
- Deviation (plan defect, noticed after the run): rule 6 has no footprint check, so the off-footprint control 16.2
  is classed `absent` rather than `off_image`. Without it the all-control fraction is 10 / 26 (computed in the session,
  not stored); the gate outcome does not depend on it.
- Rule 2 merges shared segments only within one family. Segment 3907 holds four families (1a, 1b, 1c, 1m), so the
  1a.3 and 1c.3 controls are predicted about 4× too bright and miss `recovered` narrowly (best ratio 0.333 and 0.317).
  Counting both as recovered gives 12 / 27, still a failure (computed in the session, not stored).
- Rule 6 excludes only the reference's segment from the search, so a catalogued sibling's segment could mask a
  missing image (it can hide an `absent`, never create one).
- The script and its tests were committed before the run (c4932e9, with two over-long lines that failed lint). After
  the run: lint fixes, the code-review fixes above, a positive-flux guard on references (no segment with S/N > 5 has
  flux ≤ 0, so no change), a single reference function, and the post hoc diagnostic columns and counts.
- One model (CANUCS v1 maps via the committed `../macs0417_frame` tables); no injections; no limit.

## Next
- CANUCS v2 deblended image catalogues when public, or an arc-aperture sum on member-subtracted images, frozen and
  checked on these same controls first. Not pipeline segment totals (D-034 addendum 2026-10-22).
