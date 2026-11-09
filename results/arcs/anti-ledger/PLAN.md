# `arcs`: anti-row ledger after vetting (plan, 2026-11-09, before the computation; post hoc)

**Question.** Addendum 9's post hoc split left the single-component, gate-passed class with 3 anti rows vs 0.631
expected (P 0.024; `../pooled-census/posthoc.json`). Addenda 13 and 19 have since classed WHL0137 2194
`model_sensitive` and RXC J0600 581 `member (model input)` under GLAFIC v3. Does any single-component excess remain
once the vetted rows are taken out?

**Hypothesis (ordinary).** The class count drops to 1 (MACS0308 2157) and P(≥ 1) is ≥ 0.05. **Alternative.** P stays
< 0.05 and 2157 alone is improbable under the shape-noise null at the class level.

**Data.** `../pooled-census/census.ecsv` only (stored; no new data, no network). Row probabilities `p_anti_row` as
stored. Poisson-binomial tail by `arcs_vet_noise.excess_p` (main, unchanged).

**Variants (frozen).**
- `reproduce`: all 54 non-two-knot gate-passed rows, k 3; must give 0.631 and P 0.0242 (stored) or the script stops.
- `primary`: drop 2194 and 581 (row and its probability), k 1.
- `drop_2194_only`: 581 kept as an arc (the Lenstool reading, addendum 19), k 2.
- `obs_only`: drop the two rows from the count, keep every probability (bias check), k 1.

**Outcome (frozen; ASSUMPTION threshold 0.05).** `primary` P ≥ 0.05 → `single-component excess not supported after
vetting (post hoc)`; P < 0.05 → `single-component excess remains (post hoc)`. The other variants are reported only.

**What changes the plan.** Not supported → the `arcs` anti line has no pooled excess left in the replication fields;
2157 stays an open single row (spectroscopy or a JWST grism only). Remains → state the class-level P for 2157 and
keep it as the one open row.

**Labels.** `p_anti_row` model_prediction (shape-noise null); statuses derived from the cited addenda; whole test post
hoc and exploratory (the split was defined after the counts were seen; aligned rows were never vetted for membership,
so dropping vetted anti rows only is a one-sided selection).
