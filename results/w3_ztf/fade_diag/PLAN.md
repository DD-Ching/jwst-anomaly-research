# W3 ZTF fade-control rate across areas: diagnostic plan (frozen before computing)

Written 2026-10-20 10:20 UTC, before any number below `RESULT.md` was computed. Post hoc with respect to the
v13–v22 rounds (their outcomes are known); offline on committed tables only. **No method change:** the frozen
v13 rule, the v6 pass criteria and the 10 % fade limit stay as they are; anything found here is input for the
owner's method decision (v22 RESULT section F).

## Question
In v22, 7 of 8 queued areas failed the fade control (10.2–13.8 % passing, limit 10 %), against ≤ 9.7 % in v16.
Is the fade-control rate (a) binomial scatter around one common rate, (b) tied to a per-area property known
**before** fetching light curves (inventory metadata), or (c) neither?

## Ordinary explanations to test first
- Binomial scatter at n = 1,000 trials per area (σ ≈ 0.9 points at p ≈ 0.09).
- Sampling: the fade model (3 mag drop over 5–20 d, exponential recovery 30–200 d; `control_model` in
  `scripts/w3_ztf.py`) looks more W3-like when sampling is dense or gaps fall at the drop.
- Sky position (ecliptic / Galactic latitude), which sets ZTF cadence and crowding.

## Data (observed metadata / derived decisions, all committed)
- Areas: every `results/w3_ztf/v1[6-9]/*/decision.json`, `v2*/*/decision.json` with `controls.fade.n` = 1,000 and
  an entry in `results/w3_ztf/v16/inventory.json` or `v22/inventory.json` (v14 8170, v15 183 / 3837 and v13 have
  no inventory or n = 200; excluded).
- Predictors (from the inventory, i.e. known before the fetch): `median_exposures`,
  `infobits_nonzero_in_footprint_median`, `product_ratio_median`, `score_min_over_cells`, mean
  `necessary_fraction` over the 10 cells; plus ecliptic and Galactic latitude of the HEALPix order-5 nested pixel
  centre; plus the round's inventory (v16 vs v22) as a grouping variable.

## Tests and decision rule (thresholds are ASSUMPTIONs)
1. Homogeneity: χ² of the 2 × k table (passed / not passed per area). H1 (heterogeneous) if p < 0.01.
2. For each predictor, Spearman ρ with the fade rate. A predictor is **reported as a candidate queue criterion**
   only if |ρ| ≥ 0.6 and p < 0.01 / (number of predictors) (Bonferroni).
3. Outcomes:
   - **A:** test 1 not heterogeneous → the v22 failures are scatter near the limit; the limit's margin is the issue.
   - **B:** heterogeneous and ≥ 1 predictor passes test 2 → report it, with how many v22 failures a cut on it
     would have skipped before fetching (no cut adopted).
   - **C:** heterogeneous, no predictor passes → the cause is not in the inventory metadata; say so.
- Any pattern seen after computing is labelled post hoc.
