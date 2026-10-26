# D115b Lenstool re-fit, stage V: Lenstool's χ² of the pinned model (result, 2026-10-26)

Plan `PLAN.md` (committed 10:15:35Z, before any run; its own "10:20Z" stamp is a typo). Tool `scripts/d115b_lenstool_refit.py`, rules tested in
`tests/test_d115b_lenstool_refit.py`. Numbers: `result.json` (gate), `posthoc.json` (diagnosis), Lenstool output
`chires_v0.dat.gz` and `lenstool_v0.log.gz`. All χ² values are `model_prediction` of the pinned `macs0416-canucs`
`best.par` (D-044) against the published image list (`observed`).

## Result
- **Gate V fails** (`gate_V` false). Lenstool 8.15.5 (`inverse 0`, `chi2 1`, image plane, `sigposArcsec 0.49`)
  gives χ²_pos 356.674 for the pinned `best.par`, against the header value 344.298 written by the CANUCS optimisation:
  +3.59 % (`rel_diff`), outside the 1 % tolerance (ASSUMPTION). All 303 images and 111 families are in
  `chires.dat`, with 0 warnings (post hoc count, `posthoc.json`).
- D115a contributes 4.60 and D115b 0.50 (post hoc reading, `posthoc.json`).

## Reading
Lenstool's own χ² evaluation of `best.par` does not reproduce the published χ² within the frozen tolerance, so a
Lenstool Δχ² for a freed 8515 σ has no validated reference yet. Per the plan, no stage F is run.

## Post hoc (after the gate failed)
- Post hoc: not the solver grid; `grille nombre` 256 and 512 give the same χ² 356.674111 as 128 (`posthoc.json`).
- Post hoc: the largest family terms are D52a 20.99, D52b 19.84, D34 16.37, D28 13.10 (`posthoc.json`). Not diagnosed further.
- Untested candidate causes: a different Lenstool version in 2024 (the header states none), the 6-decimal rounding of
  `best.par`, or a χ² setting of the optimisation run that `input.par` does not record.
- The port's fiducial χ²_pos for the same model is 330.81 (D115b record, "Grid step" row), so the three values
  are 330.81 (port), 344.30 (header) and 356.67 (Lenstool 8.15.5).

## Deviations
- **The tool and its tests were not committed before either run.** A lint failure (line length) stopped the commit
  in the same command, and this was noticed only afterwards. Only the plan was committed before the runs. The unit
  tests of the gate passed before the reported run (not stored); the gate function, not a reading, decided.
- First run (10:16Z): the builder kept `input.par`'s `multfile` name `gold_canucs_v2.fk5.dat`, which is not
  published. The plan specifies `multfile 1 arcs.dat`. Lenstool read no image and printed no χ² under `inverse 0`.
- A probe found Lenstool's runmode keyword `chi2`. `chi2 1` makes it print `INFO: chi2` and write `chires.dat`
  (a probe with the wrong `multfile` returned χ² 0). The builder now points `multfile` at `arcs.dat` and adds
  `chi2 1`; the reported run is the second (10:18Z). The gate and the tolerance are unchanged.
- A probe with `inverse 1` and nothing freed did not finish within 300 s and was stopped.

## Next
Find which definition differs before a re-fit: compare Lenstool's per-image χ² (`chires_v0.dat.gz`) with the port's
per-image χ² on the same back-traced sources, family by family. Then freeze a re-stated gate V (e.g. a Lenstool-internal
reference: Δχ² between Lenstool runs only) as a labelled deviation, or ask the CANUCS team which Lenstool version and
settings produced the header χ² (owner's call).
