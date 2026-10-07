# jwst-anomaly-research

Open, reproducible search for **anomalous astronomical observations in public JWST data** (via MAST).
The goal is to rank unusual sources efficiently and transparently, then vet them against ordinary
instrumental and astrophysical explanations. Later milestones target unusual gravitational-lensing
signatures.

> **Scientific stance.** An anomaly score only ranks sources. It is not evidence of new physics. Exotic
> interpretations are hypotheses of last resort, considered only after artifacts and known populations
> have been ruled out ([methodology](docs/methodology.md)).

**Status:** bootstrap (M0). The first vertical slice is under construction. See [TASKS.md](TASKS.md) and [ROADMAP.md](ROADMAP.md).

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

Downloaded data goes to `$JWST_ANOMALY_DATA` (default `./data`) and is never committed.
`data/manifests/` records exactly what was fetched.

## Project memory

The repository doubles as the long-term memory of its maintainers, both human and AI agents:
[TASKS.md](TASKS.md) (queue), [CHANGELOG.md](CHANGELOG.md) (progress and handoff),
[DECISIONS.md](DECISIONS.md) (why we chose what we chose), [SOURCES.md](SOURCES.md) (provenance),
[CLAUDE.md](CLAUDE.md) + [docs/agent-charter.md](docs/agent-charter.md) (agent operating rules).

## License

BSD-3-Clause. See [LICENSE](LICENSE).
