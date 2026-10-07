"""Run provenance capture (code version, config, environment). Owner: bootstrap unit 6."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def capture_run_context(config_path: str | Path | None = None) -> dict[str, Any]:
    """Return git commit/dirty state, config hash, package versions and UTC time for a run."""
    raise NotImplementedError("bootstrap unit 6: candidate store + runner")
