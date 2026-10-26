# D115b Lenstool re-fit, stage V: Lenstool's χ² of the pinned model (result, 2026-10-26)

Plan `PLAN.md` (committed 10:15:35Z, before any run; its own "10:20Z" stamp is a typo). Tool
`scripts/d115b_lenstool_refit.py`, rules tested in `tests/test_d115b_lenstool_refit.py`. Numbers: `result.json` (gate),
`posthoc.json` (diagnosis after the gate failed), Lenstool output `chires_v0.dat.gz` and `lenstool_v0.log.gz`. All χ²
values are `model_prediction` of the pinned `macs0416-canucs` model (D-044) against the published image list (`observed`).

## Result
- **Gate V fails as frozen** (`gate_V` false). With the run file the plan specifies (`best.par` with its own
  `reference` line), Lenstool 8.15.5 (`inverse 0`, `chi2 1`, image plane, `sigposArcsec 0.49`) gives χ²_pos (`chimul`)
  356.67, against the header value 344.298308 written by the CANUCS optimisation. That is +3.59 % (`rel_diff`),
  outside the 1 % tolerance (ASSUMPTION). 111 families are read, with 0 warnings.

## Post hoc (after the gate failed; `posthoc.json`)
- **The difference is the rounded reference point.** `best.par` writes the WCS reference to 6 decimals
  (64.038217, −24.067501). `input.par` carries 64.03821668520617, −24.06750122335985, so the catalogue shifts by 1.3 mas
  against the potentials (not stored). With the `input.par` reference and nothing else changed, Lenstool
  gives 344.42, within 0.04 % of the header (not stored: (344.42 − 344.298308) / 344.298308).
- Not the solver grid: `grille nombre` 128, 256 and 512 give 356.67 each.
- With the rounded reference, the largest family terms are D52a 20.99, D52b 19.84, D34 16.37 and D28 13.10. D115a
  contributes 4.60 and D115b 0.50.
- So a 1.3 mas reference shift moves this model's χ² by 12.25 (356.67 − 344.42, not stored). Its high-magnification images make the
  image-plane χ² sensitive at the mas level. Any re-fit must keep the full-precision reference.

## Reading
Lenstool 8.15.5 reproduces the published χ² when it is given the reference the optimisation used. The frozen gate
failed on the rounded reference that the plan specified (356.67, vs 344.42 with the full-precision one). Per the plan, no stage F was run in this round.

## Deviations
- **The tool and its tests were not committed before the first two runs.** A lint failure in the same command stopped
  the commits, and this was noticed only afterwards. Only the plan was committed before them. The gate's unit tests
  passed before the second run (not stored), and the gate function decided the outcome.
- The reported numbers come from a re-run of the committed tool (commit `1c7cc76d`, after the review fixes). It reads
  `chimul` from `chires.dat` instead of Lenstool's `INFO: chi2` total, which is the same value here: 356.674111 vs
  356.67 (`chi2_info_line`).
- First run (10:16Z, not stored): the builder kept `input.par`'s `multfile` name `gold_canucs_v2.fk5.dat`, which is not
  published. The plan specifies `multfile 1 arcs.dat`. Lenstool read no image and printed no χ² under `inverse 0`.
- A probe found Lenstool's runmode keyword `chi2`. `chi2 1` makes it print a χ² and write `chires.dat`. A probe with the
  wrong `multfile` returned χ² 0 without a warning. The builder points `multfile` at `arcs.dat` and adds `chi2 1`. The
  gate and the tolerance are unchanged.
- A probe with `inverse 1` and nothing freed was stopped by a timeout (300 s, not stored).
- Lenstool's dof (plan: "reported, not gated") is not reported: `chi2 1` output carries none.

## Next
Stage V′ (a labelled deviation, frozen before its run): the same gate on the run file with `input.par`'s
full-precision reference. If it passes, stage F gets a frozen plan: a restricted `runmode inverse` with 8515's σ free
(the record's next step), Δχ² read against the stage V′ χ².
