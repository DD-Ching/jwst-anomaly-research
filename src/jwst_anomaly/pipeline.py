"""Config-driven end-to-end runner chaining the stage modules. Owner: bootstrap unit 6."""

from __future__ import annotations

from pathlib import Path


def run(config_path: str | Path) -> str:
    """Run the pipeline described by a YAML config and return the run id."""
    raise NotImplementedError("bootstrap unit 6: candidate store + runner")
