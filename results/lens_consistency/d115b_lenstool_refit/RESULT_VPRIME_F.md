# D115b Lenstool re-fit, stages V′ and F: σ(8515) profile in Lenstool's own χ² (result, 2026-10-26)

Plan `PLAN_VPRIME_F.md` (committed 11:18:58Z, before any run of this round; its first "11:27Z" stamp was a typo,
corrected in the tool commit). Tool `scripts/d115b_lenstool_refit.py --stage vprime-f`, rules (`gate_v`, `rule_f`)
tested in `tests/test_d115b_lenstool_refit.py`, committed before the run. Numbers: `result_vprime_f.json` (gates),
`posthoc_vprime_f.json` (post hoc), Lenstool `chires_f1.00.dat.gz` and `chires_f0.90.dat.gz`. All χ² values are
`model_prediction` of the pinned `macs0416-canucs` model (D-044) against the published image list (`observed`).

## Result
- **Gate V′ passes** (`gate_Vprime` true): with `input.par`'s full-precision reference, Lenstool 8.15.5 gives
  χ²_pos 344.42 vs the header's 344.298308 (`rel_diff_vprime` 3.5 × 10⁻⁴; 111 families, 0 warnings). It confirms
  the post hoc stage V value with the committed tool; it is not a blind test.
- **Rule F as frozen: `costly`** (`rule_F`). The flag image (D115b.1 or D115b.2 within 0.7″) is present at f 1.00,
  0.90 and 0.85 and absent at 0.80, 0.75 and 0.70, so f_L* = 0.80 (96.48 km/s). Δχ²(0.80) = +3758.77 > 6
  (ASSUMPTION). The port's scan with all sources fitted found f* 0.75 (D115b record); Lenstool's f_L* is one grid
  step higher.

## Post hoc (after rule F was applied; `posthoc_vprime_f.json`)
- **The cost is one family row.** Δχ² is +3759.72 already at f 0.90, where the flag image is still present, and stays
  within 3755.72–3769.59 down to f 0.70. D14a's family row gives 3759.12 at f 0.90 (2.01 at f 1.00), while its four
  image rows sum to 10.68 (2.00 at f 1.00). So about 3748 of D14a's row is not in its image terms (not stored);
  the cause is not diagnosed. D14a.1 lies about 54″ from member 8515 (not stored).
- **Lost images carry no χ².** From f 0.90 down, D115a reports 2 warnings, and D115a.1 and D115a.2 have χ² 0.00. So
  D115a's χ² falls (−4.18 at f 0.90) where the port's "nothing refitted" run priced the lost image at thousands
  (5261 with σ × 0.7, D115b record).
- Without D14a's row, Δχ² is +1.43 at f_L* 0.80, +6.03 at 0.75 and +12.04 at 0.70 (`dchi2_without_D14a`). This
  does not price the lost D115a images either, so it is not a cost estimate.
- **F2 cost probe did not start.** Lenstool stopped at once: `ERROR: number of clumps with optimized parameters = 5
  but nlens_opt = 0` (`result_vprime_f.json` `stage_F2_probe`: `lenstool_error_lines`, status "lenstool stopped
  (signal 15)", 0.0 s). An optimised run needs
  `grille nlens_opt` and the freed potentials among the first `nlens_opt` (8515 is the 31st `potentiel` of
  `best.par`).

## Reading
The frozen outcome is `costly`. Post hoc, Lenstool's fixed-parameter barycentre χ² is not a usable cost measure for
this question: a family row jumps by about 3748 outside its image terms, and lost images score 0. The profile
neither supports nor removes the port's +12.28 (all sources fitted). The D115b verdict stays `inconclusive`;
Lenstool's f_L* (0.80) and the port's f* (0.75) agree within one grid step on where the flag image disappears.

## Deviations
- The first invocation (11:20Z) stopped at the F2 probe (Lenstool exit −15) before any output was written, so no
  number was read. The tool was changed to store V′ and F before the probe and to record a Lenstool exit (commit
  after the run of the plan, before the re-run), then the whole stage was re-run; only the re-run is reported.
- `/code-review` fixes after the run (stale `image.all` removed before each run, Lenstool's error lines kept, the
  post hoc writer and the `chires` copies made part of the tool, the caller's `OMP_NUM_THREADS` kept); the stage was
  re-run from that commit. V′, F and the post hoc numbers are unchanged, and the `chires` copies are byte-identical
  (not stored); only the probe's status text changed.
- The plan freed potentials with `limit` blocks only; Lenstool also needs `nlens_opt` (post hoc, above).

## Next
A stage F2 plan: Lenstool `inverse` with `nlens_opt` set and 8515 moved among the first `nlens_opt` potentials,
validated first by the same run with 8515's σ fixed (it should return near 344.42); its χ² must count lost images
(Lenstool's `forme` / penalty settings to be checked first). D14a's family-row term is diagnosed before any
Lenstool Δχ² is read.
