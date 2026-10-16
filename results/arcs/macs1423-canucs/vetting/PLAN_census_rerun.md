# MACS1423 clump-pair census: re-run with the fixed census code (plan, 2026-10-16, before the run)

**Why.** `PLAN_census.md` froze #273's census code at `43fbbbd`. #273's review fixes, merged later, changed
`arcs_vet_deblend.fit_models` (start flux and upper bound from `flux_scale`: faint sources whose filtered aperture sum
is ≤ 0 or near 0 were pinned) and the census null/count bookkeeping. The 2026-10-16 handoff says not to cite the
MACS1423 census numbers (P 0.259) until re-run.

**Change (only this).** `--s1063-scripts scripts` (main's `arcs_vet_census.py`, `arcs_vet_deblend.py` and their
imports) instead of the `43fbbbd` extract. Image, rows, anti sets, seeds, gates, σ_ε, fit radius and the outcome rule
are those of `PLAN_census.md`, unchanged. The MACS1423 script keeps its own null (no row is dropped from its k
counts; #273's 2184 drop is S1063-specific).

**Expectation (hypothesis).** The fix touches only rows whose filtered sum is ≤ 10 % of the positive-pixel sum, so most
classes should not change; the outcome changes only if the gate-passed `two-knot` set changes among the four anti
rows or f̂ moves by ≳ 0.1.

**What changes the plan.** Outcome unchanged ("explained") → the #282 conclusion stands and the numbers may be cited.
Outcome flips to "excess remains" → the PLAN_census.md "excess remains" branch (second band, photo-z of 2197/2322)
becomes the next step.

**Labels.** As `PLAN_census.md` (pixels observed; fits derived; injections simulated; P model_prediction;
exploratory).
