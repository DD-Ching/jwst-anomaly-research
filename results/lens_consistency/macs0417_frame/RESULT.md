# MACS0417: Jauzac et al. 2019 list → JWST frame offset; extra-image forced photometry again (result)

Plan: `PLAN.md` (frozen and committed before the run). Numbers: `summary.json`, `predicted_macs0417.ecsv`,
`control_macs0417.ecsv`. Script: `scripts/m0417_frame.py` (frozen rules unit-tested in `tests/test_m0417_frame.py`).

## Frozen outcome: frame offset established; gate C passes; 0 flags, but 0 of 6 predicted targets testable
| Step | Result | Rule |
|---|---|---|
| Offset O (JWST − RELICS, ΔRA cos δ, ΔDec) | (+0.101″, −0.292″) ± (0.001″, 0.001″), N = 471 of 728 selected RELICS sources | — |
| Gate F1 | passes: halves (even / odd id) differ by ≤ 0.002″ | N ≥ 30, error ≤ 0.02″, halves ≤ 0.03″ |
| Gate F2 | passes: median \|F200W peak − position\| of the 56 shifted images 0.038″ (unshifted 0.300″, `../macs0417_extra/offsets_posthoc.json`) | ≤ 0.18″ |
| Gate C (catalogued images as forced targets) | passes: 12 / 12 testable recovered (9) or confused (3); 41 `no_reference`, 1 `off_image` | ≥ 0.8 |
| Predicted non-catalogued images | 6: families 8, 9, 17 `faint` (μ 4.1, 3.4, 4.1), family 10 two `missing` (μ 6.1, −4.9) and one `candidate` (μ 2.8); all `no_reference` | flag = targeted image `absent` |

- F2 is consistent with H (the Jauzac list on the RELICS frame), but only at the level this diagnostic can show:
  `lc.peak_position` searches within 0.4″ of each position, so on thin arcs the brightest pixel within that radius
  need not be the catalogued knot, and a peak more than 0.4″ away cannot show up as a large residual.
- Post hoc: with one-to-one matches only (JWST sources matched once), N 465 and the same median O
  (`posthoc_unique_matches`); the quoted ±0.001″ is a formal error and ignores any position-dependent distortion.
- The testable controls come from families 5, 6, 11, 13 and 15, which have no extra predicted image. The families
  with extra predictions (8, 9, 10, 17) have no reference that passes the D-034 compact rule, so the extra-image test
  on MACS0417 stays untested. This is the plan's second branch: the next step is an arc-aware reference, not more
  frame work.
- With the shifted list the unpredicted catalogued images are 1b.4 and 11.1 (`summary.json`; unshifted only 1b.4), so
  the shift costs one model match (11.1) and the earlier run's family-11 prediction (μ 83) is no longer produced;
  family 11's control 11.2 is `recovered`. This one image points against the shift and is not explained here.

## Looked at (contact sheet of the 6 predictions and the 10 catalogued images of families 8, 9, 10, 17; not stored)
- The references 8.1, 9.1, 9.2, 10.1, 10.2 and 17.1 are thin arcs; 8.2 sits on a compact galaxy; 17.2 is faint.
- Predictions 8, 9, 17 and the family-10 `candidate` fall on empty sky; the two family-10 `missing` ones lie on and
  beside a bright clumpy galaxy. None was forced-tested, so none is read as absent or present.

## Limitations
- One model (CANUCS v1 maps, JWST frame assumed as before); no injections, so no limit.
- O is measured on RELICS sources of all morphologies with fwhm ≤ 0.5″ (pipeline `is_extended` not used as a cut).
- The frame offset is applied in this script only (a wrapped `m0417_extra.load_images`); `m0417_extra.py` and
  `../macs0417_images` are unchanged. D-034's `frame_offset_arcsec` moves the model and the list together, so an
  image-list-only offset field would be needed to reuse O there (not added).
- The F200W `_i2d.fits` is read twice (F2 diagnostic, then `run_field`); a cost, not a result issue.

## Next (optional)
- An arc-aware reference flux for families 8, 9, 10 and 17 (e.g. `halo_phot` totals or an arc-aperture sum),
  frozen and checked on the 12 testable controls first; or wait for the CANUCS v2 image catalogues.
