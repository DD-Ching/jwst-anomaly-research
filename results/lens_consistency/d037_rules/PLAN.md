# D-037 forced-photometry rules as functions (plan, frozen before the run)

**Question.** Can two of D-037's hand rules for `forced_check` flags be computed, so that `absent` / `confused` rows
are labelled automatically (TASKS Supporting 1, `forced_check` item: "the D-037 rules as code")?

**Rules (thresholds are D-037's ASSUMPTIONs).**
- **Repeated pair** (`model_topology`): a non-`observed` predicted row of system S whose catalogued images share one
  predicted image (`validate`'s `shared_match`) and which lies within 8″ of **every** image of that shared pair.
- **Model-dependent μ** (`mu_model_dependent`): the second model's nearest image of the same source (D-041 cross-model
  run, within its 3.0″ loose radius) has |μ| differing from the first model's by more than 2×.
- Only rows with `forced_class` in {`absent`, `confused`} receive a label; other rows are reported as a diagnostic.
- The fitted-redshift rule (needs a z scan per row) and the BCG-annulus rule (photometry) are not coded here.

**Inputs.** Stored CATS table `results/lens_consistency/macs1149-cats/images_predicted.ecsv`; a re-run of
`lens_consistency.py --model macs1149-cats validate` for `shared_match` (CATS maps from archive.stsci.edu, not
stored); stored Sharon μ in `results/lens_consistency/macs1149_crossmodel/summary.json`.

**Known-case gate (R).** First check the route: the re-run `validate` reproduces rms 0.673″ over 145 images and the
shared pairs recorded in docs/fields/macs1149.md (2.2/2.3, 16.1/16.2, 21.1/21.2, 1015.3/1015.4, 1050.2/1050.3,
1052.4/1052.5). Then R passes when the two functions reproduce D-037's hand labels on its 3 flags exactly:
system 16 `absent` → `model_topology`; system 2 `confused` → `model_topology`; knot 1192 `confused` →
`mu_model_dependent` and not `model_topology`. Otherwise the rules stay by-eye and the mismatches are listed.

**Ordinary explanations targeted.** Model topology / position error near critical curves; model-dependent μ.
**What changes the plan.** A pass makes both labels a reported column for future forced runs (gated next on a second
field before any use as a class); a fail keeps D-037's per-row reading by eye and closes this item.
