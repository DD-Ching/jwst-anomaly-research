# D115b Lenstool stage F2-B: budget probe of Lenstool's optimiser (result, 2026-10-26)

Plan `PLAN_F2.md` (committed 13:16Z before any run of this stage; its first "13:27Z" stamp was a typo, corrected in
the tool commit). Tool `scripts/d115b_lenstool_refit.py --stage f2b`, gate `gate_v2` and readers tested in
`tests/test_d115b_lenstool_refit.py`, committed (0ead6872) before the run. Numbers: `result_f2b.json` (gate),
`posthoc_f2b.json` and `burnin_f2b.dat.gz` (post hoc), `bayes_f2b.dat.gz` (header only). All values are
`model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044); member 8515 stays at its
`best.par` σ, so this stage gives no D115b number.

## Result
- **Gate V2: `undetermined`** (`gate_V2`). `inverse 3 0.5 20` with O1–O6 free (`nlens_opt 6`, `input.par`'s
  `limit 1`–`6`, 28 sampled parameters, counted from the `bayes.dat` header lines; count not stored) reached the 900 s wall cap (`status`, `runtime_s`
  900.1) on 4 OpenMP threads and wrote **0** samples to `bayes.dat` (`n_samples_bayes`). No Lenstool error line.
- So even this smallest BayeSys budget (ASSUMPTION) does not finish within one routine slot; Lenstool's own sampler
  at the CANUCS budget (`inverse 3 0.015 10000`) is out of reach for cloud routine runs.

## Post hoc (after the gate; `posthoc_f2b.json`)
- BayeSys wrote 40 burn-in lines in the 900 s (22.5 s per line; how many χ² evaluations a line stands for is not
  known). Their best ln(Lhood) is −89,515.45, against −302.92 for the fixed `best.par` model (stage V INFO line,
  χ² 356.67): the burn-in still samples the priors, as the plan expected.

## Reading
A Lenstool-sampler re-fit cannot start from `best.par`, and from the priors it needs far more than one routine slot
(0 samples in 900 s). Per the plan, F2 therefore needs a local optimiser started from `best.par`, planned
separately. Its cost per step is known: one fixed-parameter Lenstool χ² takes about 0.5 s wall (`result.json`
`runtime_s`, stage V). D115b stays `inconclusive`.

## Deviations
- The post hoc writer (`posthoc_f2b`, `--stage f2b-posthoc`) and its test were added after the gate was read; it
  only reads the stored burn-in file. `/code-review` fixes (stale `burnin.dat` removed before a run; seconds per
  line from the measured run time, not the cap) changed no stored number (re-run of the writer).
- Lenstool's `inverse 1` / `2` (non-BayeSys) modes were not checked in this stage.
