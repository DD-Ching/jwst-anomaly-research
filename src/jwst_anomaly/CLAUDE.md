# src/jwst_anomaly — module rules

Root CLAUDE.md (owner decisions, merge policy, budget) still applies; these rules govern the package.

## Stage contracts
- One module per pipeline stage; docs/architecture.md has the stage table.
- Stages exchange `astropy.table.Table`; column contracts live in `schema.py`.
- Every table sets `meta["provenance"]` (a `schema.Provenance` value: observed / derived / simulated /
  model_prediction / assumption / hypothesis) and `meta["source"]` (what it was derived from), including tables
  returned by helper functions.
- Changing a public stage signature requires updating docs/architecture.md and `pipeline.py` in the same PR.

## Exotic signatures (D-047, D-054)
- `signatures.py` is the registry: one `Signature` per D-047 code (W1, W2, W3, W5) with predict / inject (bound
  `exotic_sim` functions), screens, ordinary mimics (cheapest first), limits doc and decisions. A new survey enters
  through a thin adapter satisfying `LightCurveSurvey` or `CatalogueSurvey`; update the registry entry when a screen
  or limit is added.
- `exotic_sim.py` is the only source of exotic predictions and injections (`simulated` / `model_prediction`).
  Exotic physics is a hypothesis; a screen flag is an anomaly, not evidence.
- Survey adapters (`ogle.py`, `gaia_mulens.py`, `moa.py`, `lenscats.py`) pin every input by sha256 and record it in
  a manifest (`data/manifests/CLAUDE.md`).

## Code
- Thresholds are ASSUMPTIONs: name them in a `Params` dataclass or a module constant with a comment.
- Default tests are offline (synthetic fixtures); live-service tests are marked for `--run-network`.
