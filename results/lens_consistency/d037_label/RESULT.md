# D-037 `d037_label` diagnostic column in `images --forced-image` (result)

Code: `lens_consistency.py images --forced-image ... --shared-imageplane <validate imageplane.ecsv>` adds
`d037_pair` (shared-match ids, `/`-joined) and `d037_label` (`model_topology` or empty) to `images_predicted.ecsv`
and to the forced rows of `images.json`, through `d037_columns` (D-037 repeated-pair rule, `absent` / `confused`
rows only, radius 8″ = `PAIR_RADIUS_ARCSEC`, ASSUMPTION). `forced_class` is not changed: the label is a diagnostic,
not a class (D-037 addendum b). The model-dependent-μ rule needs a second model and stays in the gate scripts.

**Check** (`scripts/d037_label_check.py`, `summary.json`): on the stored MACS1149 CATS forced table
(`results/lens_consistency/macs1149-cats/images_predicted.ecsv`) with the stored `validate` imageplane
(`results/lens_consistency/d037_rules/imageplane_macs1149-cats.ecsv`), the column reproduces every gate-R
repeated-pair reading (`results/lens_consistency/d037_rules/summary.json`): 2 `confused` → `model_topology`
(2.2/2.3), 16 `absent` → `model_topology` (16.1/16.2), 1192 `confused` → none; no other `absent` / `confused`
row; `passes: true`. Model files from archive.stsci.edu (proxy CA + certifi bundle).

**Limits.** Same-route reproduction on the defining field only: it checks the column's wiring, not the rule (the
rule's second-field gate is D-037 addendum b). The column needs `validate`'s `shared_match` for the same model;
a model whose image ids differ from its published table (Abell 2744 Bergamini `JD1` vs `46.1`) needs that alias
fixed first. Provenance: forced classes derived (stored 2026-10-14 run); labels derived. No flag, no candidate.
