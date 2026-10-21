# MACS0417 / MACS1423: forced photometry at the images CANUCS v1 predicts beyond the catalogued ones (result)

Plan: `PLAN.md` (frozen and committed before the run). Numbers: `summary.json`, `predicted_*.ecsv`, `control_*.ecsv`;
post hoc diagnostic `offsets_posthoc.json` (`scripts/m0417_extra_offsets.py`). Script: `scripts/m0417_extra.py`.

## Frozen outcome: no flag; MACS1423 tested (no injections, so no limit), MACS0417 not testable by this route
| Field | Gate C (catalogued images as forced targets) | Predicted non-catalogued images | Flags |
|---|---|---|---|
| MACS1423 | passes: 11 / 11 testable recovered or confused (3 `no_reference`) | 1 `missing` (system 1, μ 17.3) `recovered`; 1 `candidate` (system 1, μ −8.8) `confused`; 2 `demagnified` (systems 1 and 2, μ 0.39 / 0.43, not tested) | 0 |
| MACS0417 | fails: 0 testable (55 / 55 `no_reference`) | 3 `missing` (family 10 twice, 11), 4 `faint` (8, 9, 10, 17), all `no_reference` | none read |

- MACS1423 system 2's extra image (the `extra_predicted` of the earlier run) is the demagnified central one
  (μ 0.43 < 0.5): expected undetectable, not tested by rule.
- MACS1423 system 1: the catalogued 1.3–1.5 are `unpredicted` by the joint source (`summary.json`), so its two
  non-catalogued predictions stand in for them; both have F200W flux there (recovered / confused).

## Post hoc (deviations from PLAN.md, labelled)
- **Why MACS0417 has no reference:** every family fails at least one D-034 reference rule. For the
  `extra_predicted` families (8, 9, 10, 11, 17) the 0.4″-recentred peak is 2.1–12.2 × the flux at the
  catalogued position (rule ≤ 2) in all 10 of their images (`offsets_posthoc.json` rows). The contact sheet (looked at,
  not stored) shows the 8, 11 and 17 images as thin extended arcs (f(0.2″)/f(0.4″) 0.29–0.68 at the peak; 9.1 / 9.2
  are compact by that rule, 0.71 / 0.77, and fail the peak / position rule only). Over all 56 published images, only 0.20 have
  peak / position ≤ 2, against 1.0 of the 17 MACS1423 images.
- **Frame:** the F200W peaks lie a median 0.30″ from the Jauzac et al. 2019 positions (median dx −0.09″,
  dy −0.27″), against 0.06″ in MACS1423 (`offsets_posthoc.json`). The coherent dy suggests a frame offset between
  the published (HST-frame) list and the JWST mosaic (not measured with stars; the median also includes extended
  arcs whose peak need not sit at the catalogued knot).
- **Contact sheet** of the 7 MACS0417 predictions and 10 catalogued images of families 8–17 (looked at): the 8,
  9 and 17 predictions (|μ| 3.4–4.2, `faint`) fall on empty sky; the two `missing` family-10 predictions (μ 10.7,
  −4.8) lie beside a bright clumpy galaxy, the third (μ 2.8, `faint`) on empty sky with an elongated galaxy in the stamp; the family-11 one (μ 83) lies on noisy sky near the arc itself. None was
  forced-tested, so none is read as absent or present.

## Limitations
- One model per field (CANUCS v1); its image catalogue is not public (ASSUMPTION as in `../macs0417_images`).
- Gate C uses the catalogued images' own predicted positions; the reference image of a family is among its
  targets, so the gate is partly circular: a loading / photometry check rather than an independent test of the flux
  prediction. Post hoc: all 11 testable MACS1423 controls passed, so leaving out one image per testable family
  (3 families) still leaves 8 / 8 (`control_macs1423.ecsv`).
- Deviation: relabelling the controls `control` means `system_search_radius` (which reads `observed` rows) gives
  systems 2 and 4 the 1.0″ base radius, while system 1 got 2.25″ (`search_arcsec` in `control_macs1423.ecsv`);
  the targets of system 1 used its residual-scaled radius. The plan said "the same `forced_check`".
- No injection-recovery was run: the MACS1423 "0 flags" covers 2 tested predictions and is not a limit.
- Family z for non-spectroscopic systems is the earlier run's CANUCS-fitted z.

## Next (optional)
- MACS0417: measure the Jauzac list → JWST frame offset on stars or compact knots, apply it as a
  `frame_offset_arcsec` of the image list, and re-run with an arc-aware reference (`halo_phot` totals), or wait for
  the CANUCS v2 image catalogues.
