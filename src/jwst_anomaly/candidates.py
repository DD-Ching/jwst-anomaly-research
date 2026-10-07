"""Persistent candidate store with vetting status. Owner: bootstrap unit 6."""

from __future__ import annotations

from pathlib import Path


class CandidateStore:
    """SQLite-backed store of runs, ranked candidates and vetting notes."""

    def __init__(self, path: str | Path) -> None:
        raise NotImplementedError("bootstrap unit 6: candidate store + runner")
