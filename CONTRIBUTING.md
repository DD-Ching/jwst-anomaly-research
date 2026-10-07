# Contributing

Thanks for your interest. The project is maintained by its owner ([@DD-Ching](https://github.com/DD-Ching))
together with AI agents, and the same rules apply to humans and agents. Read these first:

- [docs/methodology.md](docs/methodology.md): provenance labels and how anomaly candidates are vetted.
- [docs/architecture.md](docs/architecture.md): pipeline stages and table contracts.
- [CLAUDE.md](CLAUDE.md) and [docs/agent-charter.md](docs/agent-charter.md): the operating rules that
  agents follow. They are a useful summary of project norms for human contributors too.
- [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) (Contributor Covenant 2.1).

## Development setup

Use [uv](https://docs.astral.sh/uv/) and Python ≥ 3.11 (CI tests 3.11 and 3.12):

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/Scripts/python.exe -e ".[dev,cloud]"   # Linux/macOS: .venv/bin/python
```

The commands below use the Windows path `.venv/Scripts/python.exe`. On Linux or macOS, use `.venv/bin/python`.

## Tests

```bash
.venv/Scripts/python.exe -m pytest -q                 # offline tests (what CI runs on every PR)
.venv/Scripts/python.exe -m pytest -q --run-network   # also live MAST/CDS/... tests
```

- Mark any test that needs the internet with `@pytest.mark.network`. These tests are skipped unless
  you pass `--run-network` or set `JWST_ANOMALY_NETWORK=1`, so PR CI stays offline and deterministic.
- `.github/workflows/network-tests.yml` runs the network tests weekly and can be started manually
  from the Actions tab. A failed weekly run opens or updates the issue "Weekly network tests failing".
  If you change archive or catalog access code, run the network tests locally too.
- Offline tests use small synthetic or recorded fixtures (well under 1 MB), never downloaded survey data.

## Lint, format and pre-commit

CI runs ruff on `src`, `tests` and `scripts`:

```bash
.venv/Scripts/python.exe -m ruff check src tests scripts
.venv/Scripts/python.exe -m ruff format --check src tests scripts
```

[pre-commit](https://pre-commit.com/) runs ruff on the same paths, plus file-hygiene hooks: files over
1 MB, private keys, YAML/TOML syntax, whitespace and merge markers. Its ruff version is pinned and can lag
behind the latest ruff that CI uses, so if they disagree, CI wins. Install it once:

```bash
uv tool install pre-commit
pre-commit install            # run the hooks on every commit
pre-commit run --all-files    # run them on the whole repository
```

## Provenance and scientific integrity

- Every table a stage produces sets `meta["provenance"]` to one of observed, derived, simulated,
  model_prediction, assumption or hypothesis, and sets `meta["source"]` to the inputs it came from
  (see `schema.py`). Written claims in PRs, issues and docs use the same labels.
- **An anomaly score is not evidence of new physics.** It means "unusual relative to this sample
  under this feature set and model". Follow the explanation order in
  [docs/methodology.md](docs/methodology.md). Exotic interpretations are hypotheses of last resort.
- Never fabricate citations, measurements, metadata or conclusions. Cite only sources you have opened,
  with enough detail to find them again (URL plus DOI or version). If you are unsure, say so.
- Keep a simple baseline. A more sophisticated model has to beat it on injection-recovery before it
  replaces the baseline.

## Reuse decisions and sources

Before building a subsystem, check whether a mature official or maintained tool already does the job,
and read [DECISIONS.md](DECISIONS.md) first, because the question may already be settled.

- Record important reuse decisions as `D-NNN` entries in DECISIONS.md, with **Decision** /
  **Alternatives rejected** (and why) / **Evidence** / **Revisit if**.
- Record external data, software, papers and models in [SOURCES.md](SOURCES.md) with the URL, the
  version or DOI, and the date you checked them.

## Data and secrets

- Never commit astronomical data. Downloads go under `$JWST_ANOMALY_DATA` (default `<repo>/data`) and run
  outputs under `$JWST_ANOMALY_OUTPUTS` (default `<repo>/outputs`). `.gitignore` covers `data/cache/`,
  `data/raw/`, `data/derived/`, `outputs/` and common binary formats, but not every file you might put
  in `data/` (for example `.ecsv` or `.parquet`), so check `git status` before committing. Tracked
  manifests in `data/manifests/` (URI, sha256, size, pipeline version) and the scripts in `scripts/`
  are the reproducibility record.
- No file over 1 MB, no secrets, no personal data. A MAST token, if needed, comes only from the
  `MAST_API_TOKEN` environment variable. See [SECURITY.md](SECURITY.md).

## Pull requests

- Every change goes through a pull request, and the owner reviews and merges it. Nobody pushes to
  `main`, force-pushes or rewrites published history. External contributors work from a fork.
- Keep a PR to one coherent unit and fill in the template (Summary / Decision / Evidence /
  Limitations / Follow-ups, plus the checklist). CI (ruff and offline tests on Linux and Windows) must pass.
- Labels: PRs opened by AI agents carry `agent`. Add `needs-human` when a PR or issue needs an owner
  decision, meaning one that is scientific, irreversible, costly, legal or credential-related. The
  other labels are `candidate`, `science`, `reuse-decision` and `infra`.
- State files: TASKS.md is the work queue and CHANGELOG.md holds results and handoff notes. Update
  them when your change affects them. In parallel batch work, list the updates under "Follow-ups"
  instead and the coordinator folds them in.
- If you change a public stage signature, update [docs/architecture.md](docs/architecture.md) and
  `pipeline.py` in the same PR.

## Reporting an anomaly candidate

Open an issue with the **Anomaly candidate** form. It asks for:

- the `source_uid`, RA/Dec, the run ID and the exact data products (archive IDs, pipeline version);
- the evidence, with each statement labeled (observed, derived, model_prediction, ...);
- an outcome for each ordinary explanation, following steps 1–3 of docs/methodology.md: artifacts,
  catalog effects, and known populations (including SIMBAD/NED/Gaia matches). Say how you tested each one;
- links to cutouts, notebooks or runs (link to data, don't attach it).

A candidate is reported as "unexplained under tests X, Y, Z", never as a discovery. Maintainers vet
it further before any interpretation, following docs/methodology.md. A `/vet-candidate` agent skill
for this is planned in CLAUDE.md.

## Bugs, questions, security and conduct

- Bugs: use the **Bug report** form. Proposals and questions: use the **Task** form.
- Security issues: report them privately as described in [SECURITY.md](SECURITY.md), never in a public issue.
- Code of Conduct concerns: see the contact in [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).

By contributing, you agree that your contributions are licensed under the project's
[BSD-3-Clause license](LICENSE).
