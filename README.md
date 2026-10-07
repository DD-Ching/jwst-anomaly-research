# jwst-anomaly-research

Open, reproducible search for **anomalous astronomical observations in public JWST data** (via MAST).
The goal is to rank unusual sources efficiently and transparently, then vet them against ordinary
instrumental and astrophysical explanations. Later milestones target unusual gravitational-lensing
signatures.

> **Scientific stance.** An anomaly score only ranks sources. It is not evidence of new physics. Exotic
> interpretations are hypotheses of last resort, considered only after artifacts and known populations
> have been ruled out ([methodology](docs/methodology.md)).

**Status:** M0 bootstrap is complete. The M1 catalog-level slice runs end to end on real data; first
results and limitations are in [CHANGELOG.md](CHANGELOG.md). See [TASKS.md](TASKS.md) and [ROADMAP.md](ROADMAP.md).

## Pipeline (v0)

MAST query → cached acquisition with checksummed manifests → JWST level-3 pipeline catalogs →
cross-band source table → derived features → baseline anomaly ranking → image cutouts for top-k →
SIMBAD/NED/Gaia cross-check → candidate store and report. Details are in [docs/architecture.md](docs/architecture.md).

## Setup

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/Scripts/python.exe -e ".[dev,cloud]"   # Linux/macOS: .venv/bin/python
.venv/Scripts/python.exe -m pytest -q                                 # offline tests
.venv/Scripts/python.exe -m pytest -q --run-network                  # + live MAST/CDS tests
```

## Run the reference sample

```bash
export JWST_ANOMALY_DATA=~/jwst-anomaly-data             # downloaded catalogs (never committed)
export JWST_ANOMALY_OUTPUTS=~/jwst-anomaly-data/outputs  # run outputs and the candidate DB
python scripts/fetch_reference_sample.py --config configs/reference_sample.yaml --catalogs-only  # 17 catalogs, checksummed
jwst-anomaly run --config configs/reference_sample.yaml   # ~3 min: rank, S3 cutouts, cross-match, report
jwst-anomaly candidates list --run <run_id>
```

Each run writes `outputs/runs/<run_id>/report.md`, with provenance, a stage log and the top-k table,
plus a `contact_sheet.png` per sample. `data/manifests/` records exactly what was fetched: sha256,
size and pipeline version. Nothing in a report is vetted until a `/vet-candidate` record says so.

## Agent operations

AI agents and the owner maintain this project together. Agents follow [CLAUDE.md](CLAUDE.md) and
the [charter](docs/agent-charter.md), working through `/research-cycle`. How to run them is in
[docs/operations.md](docs/operations.md): a single cycle, the local `/loop`, cloud routines, or
`@claude` on issues.

## Project memory

The repository doubles as the long-term memory of its maintainers, both human and AI agents:
[TASKS.md](TASKS.md) (queue), [CHANGELOG.md](CHANGELOG.md) (progress and handoff),
[DECISIONS.md](DECISIONS.md) (why we chose what we chose), [SOURCES.md](SOURCES.md) (provenance),
[CLAUDE.md](CLAUDE.md) + [docs/agent-charter.md](docs/agent-charter.md) (agent operating rules).

## License

BSD-3-Clause. See [LICENSE](LICENSE).
