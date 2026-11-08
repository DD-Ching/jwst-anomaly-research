# D-037 forced-photometry rules as functions: gate R passes on MACS1149 (defining field)

Plan: `PLAN.md` (frozen before the run). Script: `scripts/d037_rules_gate.py`; functions
`lens_consistency.repeated_pair_topology` and `lens_consistency.mu_model_dependent`, unit-tested before the run.
Numbers: `summary.json`; re-run `validate` output: `imageplane_macs1149-cats.ecsv`.

**Route check (passes).** `validate` re-run on `macs1149-cats`: rms 0.6731″ over 145 images (recorded 0.673″);
`shared_match` images 2.2/2.3, 16.1/16.2, 21.1/21.2, 1015.3/1015.4, 1050.2/1050.3, 1052.4/1052.5, as recorded in
docs/fields/macs1149.md.

**Gate R (passes).** On the three D-037 flags of the stored CATS forced table:

| Row | CATS μ (model_prediction) | Sharon nearest μ, sep | Shared pair | Labels (derived) | Hand (D-037) |
|---|---|---|---|---|---|
| 16 `absent` | 15.82 | 30.18, 0.9″ (ratio 1.91) | 16.1/16.2 | `model_topology` | topology / fitted z |
| 2 `confused` | −16.53 | none (CATS-only) | 2.2/2.3 | `model_topology` | topology |
| 1192 `confused` | 7.52 | 3.63, 0.44″ (ratio 2.07) | none | `mu_model_dependent` | model-dependent μ |

No other `absent` / `confused` row exists, so there are no rule-only rows. Diagnostic (not a flag): the `recovered`
row of system 21 lies within 8″ of the shared pair 21.1/21.2.

**Deviations (post hoc, after the first gate run; outcome unchanged).**
- The first run keyed the cross-model rows by rounded RA/Dec and mis-keyed system 16 (Sharon μ empty); the recorded
  run matches them within 0.05″ (`xmodel_row_match_arcsec`).
- After the review: the repeated-pair rule now needs the row within 8″ of **at least two** shared-match images of its
  system, not of every one, because `shared_match` does not say which predicted image is shared (a system with two
  separate shared pairs would otherwise never be labelled); and the μ rule applies the plan's 3.0″ second-model radius
  explicitly (`mu_second_radius_arcsec`). On MACS1149 each system has one shared pair and both second-model images are
  within 0.9″, so the labels are the same.

**Limits.**
- This is the rules' own defining field, so the pass shows the functions reproduce the hand reading, not that the
  rules are validated; adopting either label as a class needs a second-field gate (plan, "What changes the plan").
- 1192 passes the μ rule narrowly (ratio 2.07 against the 2× ASSUMPTION) with the pinned `macs1149-sharon` entry;
  D-037 quoted μ 2.2 from scratch code (z_lens 0.542, image 0.56″ away), which was further from the threshold.
- The fitted-redshift and BCG-annulus rules of D-037 are not coded.

No flag, no candidate; the MACS1149 tally is unchanged.
